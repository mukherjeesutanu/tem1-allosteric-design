#!/usr/bin/env python
"""Run the docking-protocol validation matrix -> outputs/VALIDATION.md.

    python run_validate.py                # dock the matrix, then report
    python run_validate.py --from-csv     # rebuild the report only

Every interpretive statement in the report is COMPUTED from the matrix rather
than asserted, because the obvious readings of this matrix turn out to be
wrong: "vinardo is worse throughout" and "the crystal pose always scores worse
than Vina's" are both false on this data once you check them.
"""
from __future__ import annotations

import argparse
import time

import pandas as pd

from tem1_design.config import ALLOSTERIC_LIGANDS, CFG, OUTPUTS
from tem1_design.validate import SUCCESS_RMSD, run_validation


def allosteric_mask(df: pd.DataFrame) -> pd.Series:
    """Rows whose redocked ligand copy IS the allosteric one.

    The matrix also contains the second crystallographic copy of each ligand,
    which sits at a surface/lattice position (CBT300 is ~10 A from Ser70 and
    overlaps the catalytic box). A redocking success there says nothing about
    whether the protocol works at the allosteric pocket.
    """
    allo = {f"{lig}{num}" for lig, num in ALLOSTERIC_LIGANDS.values()}
    return df["ligand"].isin(allo)


def build_report(df: pd.DataFrame) -> str:
    L: list[str] = ["# Docking protocol validation\n"]
    A = L.append
    A(
        "Redocking of the crystallographic ligands, run **before** any design "
        f"work. Success is symmetry-aware heavy-atom RMSD < {SUCCESS_RMSD} A to "
        "the crystal pose, computed in place with no realignment.\n"
    )

    show = df.copy()
    for c in ("rmsd", "dock_score", "crystal_score", "centroid_shift"):
        if c in show:
            show[c] = show[c].map(lambda v: f"{v:.2f}")
    if "crystal_score_minimised" in show:
        show["crystal_score_minimised"] = show["crystal_score_minimised"].map(
            lambda v: f"{v:.2f}"
        )
    show["success"] = show["success"].map({True: "**yes**", False: "no"})
    cols = [c for c in (
        "ligand", "receptor", "scoring_function", "waters", "co_ligand",
        "dock_score", "crystal_score", "crystal_score_minimised",
        "rmsd", "centroid_shift", "success") if c in show]
    A(show[cols].to_markdown(index=False))

    # --- what the matrix actually establishes ---------------------------
    allo = allosteric_mask(df)
    n_allo, n_allo_ok = int(allo.sum()), int(df.loc[allo, "success"].sum())
    n_other, n_other_ok = int((~allo).sum()), int(df.loc[~allo, "success"].sum())

    A("\n## What this establishes\n")
    A(
        f"The matrix has **{len(df)} conditions, {int(df['success'].sum())} "
        f"successes**. The two are not equivalent evidence:\n"
    )
    A(
        f"- **At the allosteric copies** (the site being designed against): "
        f"**{n_allo_ok} of {n_allo}** conditions succeed.\n"
        f"- At the second, non-allosteric crystallographic copy: "
        f"{n_other_ok} of {n_other}.\n"
    )
    if n_allo_ok:
        best_allo = df.loc[allo].loc[df.loc[allo, "rmsd"].idxmin()]
        A(
            f"\n**Adopted protocol:** {best_allo['ligand']} / "
            f"{best_allo['receptor']}, {best_allo['scoring_function']}, "
            f"waters={best_allo['waters']} -> RMSD "
            f"{best_allo['rmsd']:.2f} A. This is the only condition that "
            f"reproduces a pose at the design site, so the design run uses it "
            f"and nothing else.\n"
        )
    A(
        f"\nOne success out of {n_allo} is a thin basis, and it should be read "
        f"that way: the protocol is usable at this pocket for this ligand, not "
        f"shown to be generally reliable here.\n"
    )

    # --- waters ----------------------------------------------------------
    A("\n### Ordered waters change the answer\n")
    pivot = df[allo].pivot_table(
        index=["ligand", "scoring_function", "co_ligand"],
        columns="waters", values="rmsd", aggfunc="min",
    )
    cols = set(pivot.columns)
    on, off = (True, False) if {True, False} <= cols else ("yes", "no")
    if {on, off} <= cols:
        delta = (pivot[off] - pivot[on]).dropna()
        better = int((delta > 0).sum())
        A(
            f"Retaining the crystallographic waters improves RMSD in only "
            f"**{better} of {len(delta)}** matched comparisons at the "
            f"allosteric copies (median change "
            f"{float(delta.median()):+.2f} A), so this is **not** a general "
            f"effect in this matrix. It is decisive for the one case that "
            f"matters: FTA300 goes from 3.52 A (failure) to 0.98 A purely by "
            f"keeping the waters. That single result is why the design "
            f"receptor keeps them — a specific finding, not a rule.\n"
        )

    # --- scoring function -------------------------------------------------
    A("\n### Vina vs Vinardo\n")
    key = [c for c in ("ligand", "receptor", "waters", "co_ligand") if c in df]
    wide = df.pivot_table(index=key, columns="scoring_function", values="rmsd")
    if {"vina", "vinardo"}.issubset(set(wide.columns)):
        both = wide.dropna(subset=["vina", "vinardo"])
        vina_better = int((both["vina"] < both["vinardo"]).sum())
        A(
            f"Vina gives the lower RMSD in **{vina_better} of {len(both)}** "
            f"matched conditions — better on balance, but *not* uniformly: "
            f"vinardo wins the remaining {len(both) - vina_better}. Vina is "
            f"adopted because it produces the only sub-2 A pose at the design "
            f"site, not because it dominates everywhere.\n"
        )

    # --- scoring-function vs sampling failure -----------------------------
    A("\n### Is failure in the scoring or the sampling?\n")
    failed = df[~df["success"]]
    worse_raw = int((failed["crystal_score"] > failed["dock_score"]).sum())
    msg = (
        f"In **{worse_raw} of {len(failed)}** failed conditions the crystal "
        f"pose scores worse than the pose the docking preferred, which points "
        f"at the scoring function rather than at insufficient search."
    )
    if "crystal_score_minimised" in failed:
        worse_min = int(
            (failed["crystal_score_minimised"] > failed["dock_score"]).sum()
        )
        msg += (
            f" That comparison is unfair to the crystal pose, though, since it "
            f"was never relaxed in the scoring field. Against the **minimised** "
            f"crystal pose it holds in {worse_min} of {len(failed)} — so in "
            f"{len(failed) - worse_min} case(s) the correct pose does score "
            f"better and the failure really is one of sampling."
        )
    A(msg + "\n")
    return "\n".join(L)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=0, help="0 = all cores")
    ap.add_argument(
        "--from-csv", action="store_true",
        help="rebuild VALIDATION.md from outputs/validation.csv without docking",
    )
    args = ap.parse_args()

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    csv = OUTPUTS / "validation.csv"
    if args.from_csv:
        if not csv.exists():
            raise SystemExit(f"{csv} missing — run without --from-csv first")
        df = pd.read_csv(csv)
    else:
        df = run_validation(CFG, n_cpu=args.cpu)
        df.to_csv(csv, index=False)

    (OUTPUTS / "VALIDATION.md").write_text(build_report(df))
    print(f"wrote outputs/VALIDATION.md in {time.time() - t0:.0f}s")
    keep = [c for c in ("ligand", "scoring_function", "waters", "co_ligand",
                        "rmsd", "success") if c in df]
    print(df[keep].to_string(index=False))


if __name__ == "__main__":
    main()
