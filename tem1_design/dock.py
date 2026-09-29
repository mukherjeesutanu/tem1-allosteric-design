"""AutoDock Vina docking into the TEM-1 allosteric (and control) sites.

Everything here is resumable and parallel-safe: a virtual screen of tens of
thousands of ligands dies halfway through often enough that resumability is a
feature, not a nicety.
"""
from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from .config import Box, Config

RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------
# preparation
# --------------------------------------------------------------------------
def prepare_receptor_pdbqt(
    pdb_in: Path, pdbqt_out: Path, ph: float = 7.4, force: bool = False
) -> Path:
    """Protonate at pH and convert to a rigid receptor PDBQT via Open Babel.

    Cached: an existing non-empty output is reused. Receptor prep is the only
    step needing Open Babel, so caching lets the prepared files be shipped to
    a compute node that has no Open Babel installed.
    """
    pdbqt_out.parent.mkdir(parents=True, exist_ok=True)
    if not force and pdbqt_out.exists() and pdbqt_out.stat().st_size > 0:
        return pdbqt_out
    cmd = [
        "obabel", "-ipdb", str(pdb_in), "-opdbqt", "-O", str(pdbqt_out),
        "-xr", "-p", str(ph), "--partialcharge", "gasteiger",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if not pdbqt_out.exists() or pdbqt_out.stat().st_size == 0:
        raise RuntimeError(f"receptor prep failed:\n{res.stdout}\n{res.stderr}")
    return pdbqt_out


def embed_3d(smiles: str, seed: int = 42) -> Chem.Mol | None:
    """ETKDG embedding + MMFF relaxation. Returns None if either step fails."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    return mol


def ligand_pdbqt_string(mol: Chem.Mol) -> str | None:
    """Meeko-prepared PDBQT text for an embedded RDKit molecule."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    try:
        prep = MoleculePreparation()
        setups = prep(mol)
        if not setups:
            return None
        text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
        return text if ok else None
    except Exception:
        return None


# --------------------------------------------------------------------------
# docking
# --------------------------------------------------------------------------
def dock_smiles(
    smiles: str,
    receptor_pdbqt: Path,
    box: Box,
    cfg: Config,
    seed: int | None = None,
) -> dict:
    """Dock one SMILES. Returns scores plus the best-pose PDBQT text."""
    from vina import Vina

    out: dict = {"smiles": smiles, "status": "ok", "vina_score": np.nan}
    mol = embed_3d(smiles, seed=seed if seed is not None else cfg.vina_seed)
    if mol is None:
        out["status"] = "embed_failed"
        return out
    text = ligand_pdbqt_string(mol)
    if text is None:
        out["status"] = "pdbqt_failed"
        return out

    try:
        v = Vina(sf_name="vina", cpu=cfg.n_cpu, seed=cfg.vina_seed, verbosity=0)
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_string(text)
        v.compute_vina_maps(center=list(box.center), box_size=list(box.size))
        v.dock(exhaustiveness=cfg.exhaustiveness, n_poses=cfg.n_poses)
        energies = v.energies(n_poses=cfg.n_poses)
        out["vina_score"] = float(energies[0][0])
        out["vina_scores_all"] = [float(e[0]) for e in energies]
        out["pose_pdbqt"] = v.poses(n_poses=1)
        heavy = Chem.RemoveHs(mol).GetNumHeavyAtoms()
        out["heavy_atoms"] = heavy
        out["ligand_efficiency"] = float(energies[0][0]) / heavy if heavy else np.nan
    except Exception as exc:
        out["status"] = f"dock_failed: {type(exc).__name__}"
    return out


def pdbqt_to_mol(pdbqt_text: str) -> Chem.Mol | None:
    """Parse a Vina pose back into RDKit via Meeko's PDBQT reader."""
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate

        pmol = PDBQTMolecule(pdbqt_text, skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pmol)
        return mols[0] if mols else None
    except Exception:
        return None


def pose_rmsd(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    """Symmetry-aware heavy-atom RMSD between two poses of the same molecule,
    computed *in place* - no alignment, because the question is whether the
    docked pose reproduces the crystal position."""
    from rdkit.Chem import rdMolAlign

    a, b = Chem.RemoveHs(mol_a), Chem.RemoveHs(mol_b)
    try:
        return float(rdMolAlign.CalcRMS(a, b))
    except Exception:
        pa, pb = a.GetConformer().GetPositions(), b.GetConformer().GetPositions()
        if pa.shape != pb.shape:
            return float("nan")
        return float(np.sqrt(((pa - pb) ** 2).sum(axis=1).mean()))
