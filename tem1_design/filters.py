"""Property, alert and novelty filters applied to every generated molecule.

Generative models emit plenty of nonsense, so the triage is explicit and every
rejection reason is recorded rather than silently dropped -- the rejection
table is one of the more informative outputs of the whole pipeline.
"""
from __future__ import annotations

import os
import sys
from functools import lru_cache

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize

from .config import REFERENCE_SMILES, Config

RDLogger.DisableLog("rdApp.*")

# SA score ships in RDKit's contrib tree, not the main API.
sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

_lfc = rdMolStandardize.LargestFragmentChooser()
_unchg = rdMolStandardize.Uncharger()


@lru_cache(maxsize=1)
def _pains_catalog() -> FilterCatalog.FilterCatalog:
    params = FilterCatalog.FilterCatalogParams()
    for cat in (
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK,
    ):
        params.AddCatalog(cat)
    return FilterCatalog.FilterCatalog(params)


def standardize(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = _unchg.uncharge(_lfc.choose(rdMolStandardize.Cleanup(mol)))
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return Chem.MolToSmiles(mol) if mol.GetNumHeavyAtoms() else None


def describe(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    ring_info = mol.GetRingInfo()
    return {
        "mw": Descriptors.MolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "qed": Descriptors.qed(mol),
        "sa_score": sascorer.calculateScore(mol),
        "max_ring_size": max((len(r) for r in ring_info.AtomRings()), default=0),
        "formal_charge": Chem.GetFormalCharge(mol),
    }


def _fpgen(cfg: Config):
    return rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg.ecfp_radius, fpSize=cfg.ecfp_bits
    )


def reference_fingerprints(cfg: Config):
    gen = _fpgen(cfg)
    return [
        gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in REFERENCE_SMILES.values()
    ]


def max_similarity_to_references(smiles: list[str], cfg: Config) -> np.ndarray:
    gen = _fpgen(cfg)
    refs = reference_fingerprints(cfg)
    out = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            out.append(np.nan)
            continue
        fp = gen.GetFingerprint(mol)
        out.append(max(DataStructs.BulkTanimotoSimilarity(fp, refs)))
    return np.asarray(out, dtype=float)


def reasons(row: pd.Series, cfg: Config) -> list[str]:
    """Every filter a molecule fails, so rejections can be audited in bulk."""
    bad = []
    if not (cfg.mw_range[0] <= row["mw"] <= cfg.mw_range[1]):
        bad.append("mw")
    if row["clogp"] > cfg.clogp_max:
        bad.append("clogp")
    if row["hbd"] > cfg.hbd_max:
        bad.append("hbd")
    if row["hba"] > cfg.hba_max:
        bad.append("hba")
    if row["rotb"] > cfg.rotb_max:
        bad.append("rotb")
    if row["tpsa"] > cfg.tpsa_max:
        bad.append("tpsa")
    if row["qed"] < cfg.qed_min:
        bad.append("qed")
    if row["sa_score"] > cfg.sa_max:
        bad.append("sa_score")
    if row["max_ring_size"] > cfg.max_ring_size:
        bad.append("ring_size")
    if cfg.require_pains_clean and row.get("alerts"):
        bad.append("structural_alert")
    if row.get("max_sim_reference", 0) > cfg.novelty_max_tanimoto:
        bad.append("too_similar_to_reference")
    return bad


def triage(smiles: list[str], cfg: Config, source: str = "") -> pd.DataFrame:
    """Standardise, describe, flag. Returns one row per unique molecule."""
    std = [standardize(s) for s in smiles]
    df = pd.DataFrame({"raw_smiles": smiles, "smiles": std, "source": source})
    df = df.dropna(subset=["smiles"]).drop_duplicates("smiles").reset_index(drop=True)
    if df.empty:
        return df

    props = pd.DataFrame([describe(s) for s in df["smiles"]])
    df = pd.concat([df, props], axis=1)

    catalog = _pains_catalog()
    df["alerts"] = [
        ";".join(sorted({e.GetDescription() for e in catalog.GetMatches(Chem.MolFromSmiles(s))}))
        for s in df["smiles"]
    ]
    df["max_sim_reference"] = max_similarity_to_references(df["smiles"].tolist(), cfg)
    df["fail_reasons"] = [";".join(reasons(r, cfg)) for _, r in df.iterrows()]
    df["passes"] = df["fail_reasons"] == ""
    return df


def rejection_summary(df: pd.DataFrame) -> pd.DataFrame:
    exploded = (
        df.loc[~df["passes"], "fail_reasons"].str.split(";").explode().value_counts()
    )
    return exploded.rename("n_molecules").to_frame().reset_index(names="filter")
