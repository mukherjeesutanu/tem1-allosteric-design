"""Protein-ligand interaction fingerprints, used as a pose-fidelity filter.

Validation showed Vina's score alone is not trustworthy in this pocket: on the
naked receptor it ranks a wrong pose above the crystallographic one. So a
generated molecule is not judged by score alone -- its pose must also
*reproduce the interaction pattern* of the reference ligand. That is what this
module measures.

LIMITS OF THIS METRIC, measured on the committed run — read before quoting it:

* ``_load_protein`` selects ``"protein"``, which DROPS the 227 crystallographic
  waters. The validation argument is that those waters define this site, so the
  fidelity score is blind to every water-mediated contact.
* The prepared receptor PDB carries no hydrogens, so ProLIF's angle-dependent
  HBond/XBond criteria cannot be evaluated properly.
* Consequently the FTA reference row sets only **4 bits, all Hydrophobic**
  (Leu221, Leu225, Val261, Ile287) of 24 columns, although FTA makes 16 residue
  contacts in the crystal. No PiStacking or ionic bit is set anywhere.

A Tanimoto against a 4-bit reference is coarse: 2 shared bits plus 2 extras
gives exactly 0.333, just over the 0.3 threshold, and 77.9% of the top-800
pass. Treat "mode fidelity" as a weak hydrophobic-contact overlap, not as
evidence that a molecule reproduces the reference binding mode.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from .config import ALLOSTERIC_POCKET_RESIDUES, Config

RDLogger.DisableLog("rdApp.*")


def _load_protein(receptor_pdb: Path):
    import MDAnalysis as mda
    import prolif as plf

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(receptor_pdb))
        sel = u.select_atoms("protein")
        return plf.Molecule.from_mda(sel, NoImplicit=False)


def fingerprint(
    poses: list[Chem.Mol],
    receptor_pdb: Path,
    interactions: list[str] | None = None,
) -> tuple[pd.DataFrame, object]:
    """ProLIF fingerprint for a list of docked poses against one receptor."""
    import prolif as plf

    protein = _load_protein(receptor_pdb)
    fp = plf.Fingerprint(
        interactions
        or [
            "Hydrophobic",
            "HBDonor",
            "HBAcceptor",
            "PiStacking",
            "Anionic",
            "Cationic",
            "CationPi",
            "PiCation",
            "XBDonor",
        ]
    )
    ligands = []
    for mol in poses:
        if mol is None:
            continue
        try:
            ligands.append(plf.Molecule.from_rdkit(mol))
        except Exception:
            continue
    if not ligands:
        return pd.DataFrame(), fp
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fp.run_from_iterable(ligands, protein, progress=False)
    return fp.to_dataframe(), fp


def tanimoto_to_reference(df: pd.DataFrame, reference_row: int = 0) -> np.ndarray:
    """Per-pose Tanimoto over interaction bits, against one reference row."""
    if df.empty:
        return np.array([])
    bits = df.to_numpy().astype(bool)
    ref = bits[reference_row]
    out = []
    for row in bits:
        inter = np.logical_and(row, ref).sum()
        union = np.logical_or(row, ref).sum()
        out.append(inter / union if union else 0.0)
    return np.asarray(out)


def pocket_residue_contacts(df: pd.DataFrame) -> pd.Series:
    """Fraction of the reference pocket's core residues each pose touches.

    A blunt but robust check that the pose is in the right place, independent
    of interaction typing.
    """
    if df.empty:
        return pd.Series(dtype=float)
    core = set(ALLOSTERIC_POCKET_RESIDUES)
    hit_counts = []
    for _, row in df.iterrows():
        touched = set()
        for (_, residue, _), value in row.items():
            if not value:
                continue
            digits = "".join(ch for ch in str(residue) if ch.isdigit())
            if digits:
                touched.add(int(digits))
        hit_counts.append(len(touched & core) / len(core))
    return pd.Series(hit_counts, index=df.index)


def summarize_interactions(df: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    """Most frequently formed interactions across a pose set."""
    if df.empty:
        return pd.DataFrame()
    freq = df.mean(axis=0).sort_values(ascending=False).head(top_n)
    return (
        freq.rename("frequency")
        .to_frame()
        .reset_index()
        .rename(columns={"ligand": "ligand", "protein": "residue", "interaction": "type"})
    )
