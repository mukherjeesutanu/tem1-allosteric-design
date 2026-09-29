"""Receptor and reference-ligand preparation, plus the pocket analysis that
decides which crystallographic ligand copy is the real allosteric one.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import ALLOSTERIC_LIGANDS, CATALYTIC_RESIDUES, Config


@dataclass
class Atom:
    record: str
    name: str
    resname: str
    chain: str
    resid: int
    xyz: np.ndarray
    element: str
    line: str


def read_pdb(path: Path) -> list[Atom]:
    atoms = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atoms.append(
            Atom(
                record=line[:6].strip(),
                name=line[12:16].strip(),
                resname=line[17:20].strip(),
                chain=line[21],
                resid=int(line[22:26]),
                xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                element=(line[76:78].strip() or line[12:16].strip()[:1]),
                line=line,
            )
        )
    return atoms


def coords(atoms: list[Atom]) -> np.ndarray:
    return np.array([a.xyz for a in atoms]) if atoms else np.empty((0, 3))


def ligand_atoms(atoms: list[Atom], resname: str, resid: int | None = None) -> list[Atom]:
    return [
        a
        for a in atoms
        if a.resname == resname and a.record == "HETATM" and (resid is None or a.resid == resid)
    ]


def contact_residues(
    atoms: list[Atom], ligand: list[Atom], cutoff: float = 4.5
) -> list[tuple[int, str]]:
    prot = [a for a in atoms if a.record == "ATOM"]
    P, L = coords(prot), coords(ligand)
    if len(P) == 0 or len(L) == 0:
        return []
    dmin = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=2).min(axis=1)
    return sorted({(prot[i].resid, prot[i].resname) for i in np.flatnonzero(dmin < cutoff)})


def distance_to_catalytic_serine(atoms: list[Atom], ligand: list[Atom]) -> float:
    ser = [a for a in atoms if a.resid == 70 and a.name == "OG"]
    if not ser or not ligand:
        return float("nan")
    return float(np.linalg.norm(coords(ligand).mean(axis=0) - ser[0].xyz))


def identify_allosteric_copy(pdb_path: Path, resname: str) -> dict:
    """Rank crystallographic copies of a ligand by how well they match the
    published allosteric site: buried, >12 A from Ser70, hydrophobic lining."""
    atoms = read_pdb(pdb_path)
    out = []
    for resid in sorted({a.resid for a in atoms if a.resname == resname and a.record == "HETATM"}):
        lig = ligand_atoms(atoms, resname, resid)
        contacts = contact_residues(atoms, lig)
        hydrophobic = sum(
            1 for _, rn in contacts if rn in {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}
        )
        out.append(
            {
                "resid": resid,
                "n_contacts": len(contacts),
                "frac_hydrophobic": hydrophobic / max(len(contacts), 1),
                "dist_ser70": distance_to_catalytic_serine(atoms, lig),
                "center": coords(lig).mean(axis=0).round(2).tolist(),
                "contacts": [f"{rn}{ri}" for ri, rn in contacts],
            }
        )
    return {"pdb": Path(pdb_path).stem, "resname": resname, "copies": out}


def write_receptor(
    pdb_path: Path, out_path: Path, keep_waters: bool = False, chain: str = "A"
) -> Path:
    """Strip heteroatoms (and by default waters) to leave a clean receptor."""
    lines = []
    for line in Path(pdb_path).read_text().splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            lines.append(line)
        elif keep_waters and line.startswith("HETATM") and line[17:20].strip() == "HOH":
            lines.append(line)
    lines.append("END")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def write_ligand_pdb(pdb_path: Path, resname: str, resid: int, out_path: Path) -> Path:
    atoms = read_pdb(pdb_path)
    lig = ligand_atoms(atoms, resname, resid)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(a.line for a in lig) + "\nEND\n")
    return out_path


def pocket_box(
    pdb_paths: list[Path], padding: float = 5.0
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Box enclosing every allosteric ligand copy across the given structures."""
    pts = []
    for p in pdb_paths:
        key = Path(p).stem
        if key not in ALLOSTERIC_LIGANDS:
            continue
        resname, resid = ALLOSTERIC_LIGANDS[key]
        pts.append(coords(ligand_atoms(read_pdb(p), resname, resid)))
    P = np.vstack(pts)
    lo, hi = P.min(axis=0) - padding, P.max(axis=0) + padding
    center = tuple(((lo + hi) / 2).round(2))
    size = tuple((hi - lo).round(1))
    return center, size


def ca_rmsd(pdb_a: Path, pdb_b: Path) -> float:
    """CA RMSD without superposition - tests whether two PDBs share a frame."""
    def cas(p):
        return {a.resid: a.xyz for a in read_pdb(p) if a.name == "CA" and a.record == "ATOM"}

    A, B = cas(pdb_a), cas(pdb_b)
    keys = sorted(set(A) & set(B))
    X = np.array([A[k] for k in keys])
    Y = np.array([B[k] for k in keys])
    return float(np.sqrt(((X - Y) ** 2).sum(axis=1).mean()))
