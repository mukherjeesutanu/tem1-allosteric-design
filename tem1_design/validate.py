"""Self- and cross-docking validation of the docking protocol.

This runs *before* any design work and decides whether docking scores into
this pocket mean anything at all. The answer turned out to be conditional:

  * 1PZP / FTA, crystallographic waters retained ....... 0.95 A  (success)
  * 1PZP / FTA, waters stripped ........................ 3.42 A  (failure)
  * 1PZO / CBT300, with the second CBT copy present .... 1.75 A  (success)
  * 1PZO / CBT301, naked protein ....................... 7.54 A  (failure)

The TEM-1 allosteric cleft is shallow, hydrophobic and solvent-mediated: its
binding-competent boundary is defined by ordered waters and, in 1PZO, by a
second copy of the ligand. Docking the naked protein reproduces neither
crystal pose, and -- worse -- scores the crystallographic pose (-6.2 kcal/mol)
below the wrong pose it prefers (-7.8). Any design campaign that skipped this
check would have optimised against an artefact.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from .config import PREPARED, REFERENCE_SMILES, Box, Config
from .dock import (
    embed_3d,
    ligand_pdbqt_string,
    pdbqt_to_mol,
    pose_rmsd,
    prepare_receptor_pdbqt,
)
from .structure import write_ligand_pdb

RDLogger.DisableLog("rdApp.*")

SUCCESS_RMSD = 2.0   # the conventional redocking success criterion


def crystal_ligand(pdb: str, resname: str, resid: int, smiles_key: str) -> Chem.Mol:
    """Crystal coordinates with bond orders restored from the SMILES template."""
    path = write_ligand_pdb(
        Path(f"data/raw/{pdb}.pdb"), resname, resid, PREPARED / f"{pdb}_{resname}{resid}.pdb"
    )
    raw = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=False)
    template = Chem.MolFromSmiles(REFERENCE_SMILES[smiles_key])
    return AllChem.AssignBondOrdersFromTemplate(template, raw)


def build_receptor(
    pdb: str,
    keep_waters: bool = True,
    keep_ligand: tuple[str, int] | None = None,
    chain: str = "A",
) -> Path:
    """Receptor PDBQT, optionally retaining waters and one ligand copy."""
    tag = f"{pdb}{'_wat' if keep_waters else ''}"
    if keep_ligand:
        tag += f"_{keep_ligand[0]}{keep_ligand[1]}"
    out_pdb = PREPARED / f"{tag}_receptor.pdb"
    lines = []
    for line in Path(f"data/raw/{pdb}.pdb").read_text().splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            lines.append(line)
        elif line.startswith("HETATM"):
            res = line[17:20].strip()
            if keep_waters and res == "HOH":
                lines.append(line)
            elif keep_ligand and res == keep_ligand[0] and int(line[22:26]) == keep_ligand[1]:
                lines.append(line)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(lines) + "\nEND\n")
    return prepare_receptor_pdbqt(out_pdb, PREPARED / f"{tag}_receptor.pdbqt")


def box_around(mol: Chem.Mol, padding: float = 5.0) -> Box:
    P = mol.GetConformer().GetPositions()
    lo, hi = P.min(axis=0) - padding, P.max(axis=0) + padding
    return Box(tuple(((lo + hi) / 2).round(2)), tuple((hi - lo).round(1)))


def redock(
    smiles: str,
    receptor_pdbqt: Path,
    box: Box,
    reference: Chem.Mol,
    sf_name: str = "vina",
    exhaustiveness: int = 32,
    seed: int = 42,
    n_cpu: int = 0,
) -> dict:
    from vina import Vina

    mol = embed_3d(smiles, seed=seed)
    text = ligand_pdbqt_string(mol) if mol is not None else None
    if text is None:
        return {"status": "prep_failed", "rmsd": np.nan, "score": np.nan}
    v = Vina(sf_name=sf_name, cpu=n_cpu, seed=seed, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_string(text)
    v.compute_vina_maps(center=list(box.center), box_size=list(box.size))
    v.dock(exhaustiveness=exhaustiveness, n_poses=5)
    pose = pdbqt_to_mol(v.poses(n_poses=1))
    ref_c = reference.GetConformer().GetPositions().mean(axis=0)
    pose_c = Chem.RemoveHs(pose).GetConformer().GetPositions().mean(axis=0)
    return {
        "status": "ok",
        "score": float(v.energies(n_poses=1)[0][0]),
        "rmsd": pose_rmsd(pose, reference),
        "centroid_shift": float(np.linalg.norm(pose_c - ref_c)),
        "pose_pdbqt": v.poses(n_poses=1),
    }


def score_crystal_pose(
    reference: Chem.Mol, receptor_pdbqt: Path, sf_name: str = "vina"
) -> dict:
    """Score the crystal pose where it sits - separates scoring-function error
    from sampling error."""
    from vina import Vina

    mol = Chem.AddHs(reference, addCoords=True)
    Chem.SanitizeMol(mol)
    text = ligand_pdbqt_string(mol)
    if text is None:
        return {"crystal_score": np.nan, "crystal_score_minimised": np.nan}
    v = Vina(sf_name=sf_name, cpu=0, seed=42, verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_string(text)
    v.compute_vina_maps(
        center=list(reference.GetConformer().GetPositions().mean(axis=0)),
        box_size=[24, 24, 24],
    )
    return {
        "crystal_score": float(v.score()[0]),
        "crystal_score_minimised": float(v.optimize()[0]),
    }


def run_validation(cfg: Config, n_cpu: int = 0) -> pd.DataFrame:
    """The full matrix reported in outputs/VALIDATION.md."""
    rows = []
    cases = [
        ("FTA", "1PZP", "FTA", 300, None),
        ("CBT", "1PZO", "CBT", 301, None),
        ("CBT", "1PZO", "CBT", 301, ("CBT", 300)),
        ("CBT", "1PZO", "CBT", 300, ("CBT", 301)),
    ]
    for smiles_key, pdb, resname, resid, co_ligand in cases:
        ref = crystal_ligand(pdb, resname, resid, smiles_key)
        box = box_around(ref, padding=5.0)
        for sf, waters in itertools.product(("vina", "vinardo"), (False, True)):
            receptor = build_receptor(pdb, keep_waters=waters, keep_ligand=co_ligand)
            r = redock(
                REFERENCE_SMILES[smiles_key], receptor, box, ref,
                sf_name=sf, exhaustiveness=32, n_cpu=n_cpu,
            )
            c = score_crystal_pose(ref, receptor, sf_name=sf)
            rows.append(
                {
                    "ligand": f"{resname}{resid}",
                    "receptor": pdb,
                    "scoring_function": sf,
                    "waters": waters,
                    "co_ligand": f"{co_ligand[0]}{co_ligand[1]}" if co_ligand else "-",
                    "dock_score": r.get("score"),
                    "rmsd": r.get("rmsd"),
                    "centroid_shift": r.get("centroid_shift"),
                    **c,
                    "success": bool(r.get("rmsd", np.inf) < SUCCESS_RMSD),
                }
            )
    return pd.DataFrame(rows)
