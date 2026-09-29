"""Parallel docking of a molecule library, resumable across restarts.

Each worker runs a single-CPU Vina instance; N workers saturate the node far
better than one Vina using N threads, because Vina's internal parallelism
scales poorly past a handful of cores on a small box.
"""
from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from .config import Box, Config

_WORKER: dict = {}


THREAD_ENV_VARS = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
)


def limit_threads() -> None:
    """Pin every worker to a single thread.

    Without this each of N workers also spawns BLAS threads, the node ends up
    oversubscribed several-fold, and throughput drops by roughly half.
    """
    for var in THREAD_ENV_VARS:
        os.environ[var] = "1"


def _init_worker(receptor: str, box_center, box_size, sf_name: str, seed: int):
    """One Vina object and one precomputed map set per worker process.

    Recomputing the affinity maps for every ligand dominates runtime; building
    them once per worker is the single biggest speedup in the screen.
    """
    limit_threads()
    from vina import Vina

    v = Vina(sf_name=sf_name, cpu=1, seed=seed, verbosity=0)
    v.set_receptor(receptor)
    v.compute_vina_maps(center=list(box_center), box_size=list(box_size))
    _WORKER["vina"] = v


def _dock_one(args) -> dict:
    smiles, name, exhaustiveness, n_poses, seed = args
    from rdkit import Chem

    from .dock import embed_3d, ligand_pdbqt_string

    rec = {"name": name, "smiles": smiles, "status": "ok", "vina_score": None}
    mol = embed_3d(smiles, seed=seed)
    if mol is None:
        rec["status"] = "embed_failed"
        return rec
    text = ligand_pdbqt_string(mol)
    if text is None:
        rec["status"] = "pdbqt_failed"
        return rec
    try:
        v = _WORKER["vina"]
        v.set_ligand_from_string(text)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = v.energies(n_poses=n_poses)
        heavy = Chem.RemoveHs(mol).GetNumHeavyAtoms()
        rec["vina_score"] = float(energies[0][0])
        rec["ligand_efficiency"] = float(energies[0][0]) / heavy if heavy else None
        rec["heavy_atoms"] = heavy
        rec["pose_pdbqt"] = v.poses(n_poses=1)
    except Exception as exc:
        rec["status"] = f"dock_failed: {type(exc).__name__}"
    return rec


def screen(
    library: pd.DataFrame,
    receptor_pdbqt: Path,
    box: Box,
    cfg: Config,
    out_jsonl: Path,
    n_workers: int | None = None,
    exhaustiveness: int | None = None,
    smiles_col: str = "smiles",
    name_col: str | None = None,
) -> pd.DataFrame:
    """Dock every row of `library`, appending results to `out_jsonl`.

    Already-docked molecules (by SMILES) are skipped, so an interrupted screen
    resumes where it stopped.
    """
    out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    done: set[str] = set()
    if out_jsonl.exists():
        for line in out_jsonl.read_text().splitlines():
            try:
                done.add(json.loads(line)["smiles"])
            except Exception:
                continue

    todo = library[~library[smiles_col].isin(done)]
    n_workers = n_workers or max(1, (os.cpu_count() or 4) - 1)
    exhaustiveness = exhaustiveness or cfg.exhaustiveness

    tasks = [
        (
            row[smiles_col],
            str(row[name_col]) if name_col else f"mol_{i}",
            exhaustiveness,
            cfg.n_poses,
            cfg.vina_seed,
        )
        for i, row in todo.iterrows()
    ]
    if not tasks:
        return read_results(out_jsonl)

    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker,
        initargs=(
            str(receptor_pdbqt),
            box.center,
            box.size,
            cfg.scoring_function,
            cfg.vina_seed,
        ),
    ) as pool, out_jsonl.open("a") as fh:
        futures = [pool.submit(_dock_one, t) for t in tasks]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="docking"):
            rec = fut.result()
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
    return read_results(out_jsonl)


def read_results(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    rows = []
    for line in path.read_text().splitlines():
        try:
            rows.append(json.loads(line))
        except Exception:
            continue
    return pd.DataFrame(rows)
