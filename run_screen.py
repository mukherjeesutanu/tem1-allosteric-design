#!/usr/bin/env python
"""Dock the filtered library into the allosteric site and two control sites.

    python run_screen.py [--workers 24] [--exhaustiveness 16] [--limit N]

Three sites are used, all on the same receptor and with the same protocol:

  allosteric   the validated design site (1PZP FTA300 pocket)
  orthosteric  the catalytic Ser70 site -> selectivity control
  surface      a crystallographic lattice site 29 A from Ser70, i.e. a real
               but non-functional surface depression -> site-specificity control

A molecule that scores well everywhere is a sticky compound, not a hit. The
controls are what separate the two.
"""
from __future__ import annotations

import os

# Must precede numpy/RDKit import: worker processes each get one thread.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tem1_design.config import CFG, GENERATED, OUTPUTS, REFERENCE_SMILES, Box
from tem1_design.screen import read_results, screen

# The allosteric and orthosteric boxes come from the config so there is a
# single source of truth: an oversized box here silently breaks the validated
# protocol (Vina drifts into the neighbouring surface groove) while still
# producing plausible-looking numbers.
SITES = {
    "allosteric": CFG.allosteric_box,
    "orthosteric": CFG.orthosteric_box,
    # crystallographic lattice site from 1PZP (FTA copy 301), 29 A from Ser70:
    # a real surface depression that has no functional role.
    "surface": Box((13.85, -1.87, 10.18), (18.0, 18.0, 18.0)),
}


def diverse_subset(lib: pd.DataFrame, n: int) -> pd.DataFrame:
    """MaxMin pick over ECFP4, proportionally per generator.

    References are always retained: they are the yardstick every generated
    molecule is compared against.
    """
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.SimDivFilters import MaxMinPicker

    refs = lib[lib["source"] == "reference"]
    pool = lib[lib["source"] != "reference"]
    if len(pool) <= n:
        return lib

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    picked = []
    for source, group in pool.groupby("source"):
        quota = max(1, round(n * len(group) / len(pool)))
        quota = min(quota, len(group))
        fps, idx = [], []
        for i, smi in zip(group.index, group["smiles"]):
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                fps.append(gen.GetFingerprint(mol))
                idx.append(i)
        if not fps:
            continue
        chosen = MaxMinPicker().LazyBitVectorPick(fps, len(fps), min(quota, len(fps)), seed=42)
        picked.extend(idx[c] for c in chosen)
    return pd.concat([lib.loc[picked], refs]).drop_duplicates("smiles")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0, help="0 = cpu_count-1")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument(
        "--subset", type=int, default=0,
        help="dock a MaxMin diversity-selected subset of this size, stratified "
             "by generator. Docking is the pipeline bottleneck; a diverse "
             "subset covers the generated chemical space far better than an "
             "arbitrary head() of it.",
    )
    ap.add_argument(
        "--top-from-allosteric", type=int, default=0,
        help="dock only the N best-scoring molecules from the allosteric run. "
             "Selectivity and specificity only matter for molecules that look "
             "good at the target site, so the controls run on the funnel's "
             "survivors rather than the whole library.",
    )
    ap.add_argument(
        "--random-from-allosteric", type=int, default=0,
        help="dock a RANDOM sample of N already-docked molecules, chosen "
             "independently of their allosteric score. --top-from-allosteric "
             "selects on the very quantity the control is meant to test, so "
             "its apparent selectivity is partly regression to the mean; this "
             "is the unbiased comparison.",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--tag", default="",
        help="suffix for the output files, e.g. _random",
    )
    ap.add_argument("--sites", default=",".join(SITES))
    ap.add_argument(
        "--receptor",
        default="data/prepared/1PZP_wat_receptor.pdbqt",
        help="validated receptor: 1PZP with crystallographic waters",
    )
    args = ap.parse_args()

    lib_path = GENERATED / "library_filtered.parquet"
    if not lib_path.exists():
        raise SystemExit(f"{lib_path} missing - run run_generate.py first")
    lib = pd.read_parquet(lib_path)
    if args.limit:
        lib = lib.head(args.limit)

    # reference ligands go through the identical protocol as the yardstick
    refs = pd.DataFrame(
        {"smiles": list(REFERENCE_SMILES.values()), "name": list(REFERENCE_SMILES)}
    )
    lib = pd.concat(
        [lib.assign(name=lib["smiles"]), refs.assign(source="reference")],
        ignore_index=True,
    )

    if args.subset:
        lib = diverse_subset(lib, args.subset)
        print(f"diversity-selected subset: {len(lib):,} molecules")

    if args.random_from_allosteric:
        prev = read_results(GENERATED / "docked_allosteric.jsonl")
        if prev.empty:
            raise SystemExit("--random-from-allosteric needs a finished allosteric run")
        pool = prev[prev["status"] == "ok"]["smiles"].tolist()
        rng = np.random.default_rng(args.seed)
        n = min(args.random_from_allosteric, len(pool))
        chosen = set(rng.choice(pool, n, replace=False).tolist())
        chosen |= set(REFERENCE_SMILES.values())
        lib = lib[lib["smiles"].isin(chosen)]
        print(f"controls on a RANDOM sample of {len(lib):,} molecules "
              f"(seed {args.seed}), independent of allosteric score")

    if args.top_from_allosteric:
        prev = read_results(GENERATED / "docked_allosteric.jsonl")
        if prev.empty:
            raise SystemExit("--top-from-allosteric needs a finished allosteric run")
        best = (
            prev[prev["status"] == "ok"]
            .nsmallest(args.top_from_allosteric, "vina_score")["smiles"]
            .tolist()
        )
        keep = set(best) | set(REFERENCE_SMILES.values())
        lib = lib[lib["smiles"].isin(keep)]
        print(f"controls restricted to {len(lib):,} molecules "
              f"(top {args.top_from_allosteric} + references)")

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    receptor = Path(args.receptor)
    for site in args.sites.split(","):
        site = site.strip()
        if site not in SITES:
            continue
        t0 = time.time()
        out = GENERATED / f"docked_{site}{args.tag}.jsonl"
        print(f"\n=== {site} site ({len(lib):,} molecules) -> {out.name}")
        screen(
            lib, receptor, SITES[site], CFG, out,
            n_workers=args.workers or None,
            exhaustiveness=args.exhaustiveness,
            smiles_col="smiles", name_col="name",
        )
        df = read_results(out)
        ok = df[df["status"] == "ok"]
        print(
            f"    {len(ok):,}/{len(df):,} docked in {time.time() - t0:.0f}s | "
            f"median {ok['vina_score'].median():.2f} | best {ok['vina_score'].min():.2f}"
        )


if __name__ == "__main__":
    main()
