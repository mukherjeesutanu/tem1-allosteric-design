#!/usr/bin/env python
"""Generate a candidate library with every available generator, then triage.

    python run_generate.py [--n 3000] [--generators crem,ga,safe]

Writes data/generated/library.parquet plus a per-generator comparison, which
is one of the project's actual results: the three paradigms differ enormously
in how much of their raw output survives basic medicinal-chemistry filters.
"""
from __future__ import annotations

import argparse
import time

import pandas as pd

from tem1_design.config import CFG, GENERATED, OUTPUTS, REFERENCE_SMILES
from tem1_design.filters import rejection_summary, triage

GENERATORS = ("crem", "ga", "safe")


def run_generator(name: str, n: int) -> list[str]:
    t0 = time.time()
    if name == "crem":
        from tem1_design.generate.crem_gen import grow_and_mutate

        out = grow_and_mutate(CFG, n_target=n)
    elif name == "ga":
        from tem1_design.generate.selfies_ga import run_ga

        out = run_ga(CFG, n_target=n)
    elif name == "safe":
        from tem1_design.generate.safe_gen import generate, is_available

        if not is_available():
            print("  safe: package/checkpoint unavailable, skipping")
            return []
        out = generate(CFG, n_target=n)
    else:
        raise ValueError(name)
    print(f"  {name}: {len(out)} raw molecules in {time.time() - t0:.0f}s")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=CFG.n_per_generator)
    ap.add_argument("--generators", default=",".join(GENERATORS))
    args = ap.parse_args()

    GENERATED.mkdir(parents=True, exist_ok=True)
    OUTPUTS.mkdir(parents=True, exist_ok=True)

    frames, stats = [], []
    for name in args.generators.split(","):
        name = name.strip()
        if not name:
            continue
        raw = run_generator(name, args.n)
        if not raw:
            stats.append({"generator": name, "raw": 0, "unique_valid": 0, "passed": 0})
            continue
        df = triage(raw, CFG, source=name)
        frames.append(df)
        stats.append(
            {
                "generator": name,
                "raw": len(raw),
                "unique_valid": len(df),
                "passed": int(df["passes"].sum()),
                "pass_rate": float(df["passes"].mean()),
                "median_qed": float(df["qed"].median()),
                "median_sa": float(df["sa_score"].median()),
                "median_mw": float(df["mw"].median()),
                "median_sim_to_ref": float(df["max_sim_reference"].median()),
            }
        )

    if not frames:
        raise SystemExit("no generator produced anything")

    lib = pd.concat(frames, ignore_index=True)
    # a molecule found by two generators is kept once, crediting both
    dup = lib.groupby("smiles")["source"].apply(lambda s: ",".join(sorted(set(s))))
    lib = lib.drop_duplicates("smiles").drop(columns=["source"]).merge(
        dup.rename("source"), on="smiles"
    )
    lib.to_parquet(GENERATED / "library.parquet", index=False)

    passed = lib[lib["passes"]].copy()
    passed.to_parquet(GENERATED / "library_filtered.parquet", index=False)

    stat_df = pd.DataFrame(stats)
    stat_df.to_csv(OUTPUTS / "generator_comparison.csv", index=False)
    rejection_summary(lib).to_csv(OUTPUTS / "rejection_summary.csv", index=False)

    print("\n=== generator comparison ===")
    print(stat_df.to_string(index=False))
    print(f"\nlibrary: {len(lib):,} unique molecules, {len(passed):,} pass all filters")
    print("\nrejection reasons:")
    print(rejection_summary(lib).to_string(index=False))

    ref = triage(list(REFERENCE_SMILES.values()), CFG, source="reference")
    ref.to_parquet(GENERATED / "references.parquet", index=False)


if __name__ == "__main__":
    main()
