#!/usr/bin/env python
"""Analyse the screen: PLIF fidelity, consensus ranking, figures, RESULTS.md.

    python run_analyze.py [--top-plif 800]
"""
from __future__ import annotations

import argparse
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tem1_design import plots
from tem1_design.config import CFG, FIGURES, GENERATED, OUTPUTS, PREPARED, REFERENCE_SMILES
from tem1_design.dock import pdbqt_to_mol
from tem1_design.plif import fingerprint, pocket_residue_contacts, tanimoto_to_reference
from tem1_design.score import apply_hurdles, consensus_rank, hurdle_funnel, merge_sites
from tem1_design.screen import read_results


def compute_plif(docked: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Interaction-fingerprint similarity to the reference ligand's pose.

    Restricted to the best-scoring molecules: PLIF is the expensive step and
    only matters for molecules already past the affinity hurdle.
    """
    ok = docked[(docked["status"] == "ok") & docked["pose_pdbqt"].notna()]
    ref_smiles = set(REFERENCE_SMILES.values())
    refs = ok[ok["smiles"].isin(ref_smiles)]
    best = ok.nsmallest(top_n, "vina_score")
    subset = pd.concat([refs, best]).drop_duplicates("smiles")
    if subset.empty:
        return pd.DataFrame()

    poses, smiles = [], []
    for _, row in subset.iterrows():
        mol = pdbqt_to_mol(row["pose_pdbqt"])
        if mol is not None:
            poses.append(mol)
            smiles.append(row["smiles"])
    if not poses:
        return pd.DataFrame()

    receptor = PREPARED / "1PZP_wat_receptor.pdb"
    df, _ = fingerprint(poses, receptor)
    if df.empty:
        return pd.DataFrame()

    # reference row = FTA if present, else the best-scoring molecule
    ref_idx = 0
    for i, smi in enumerate(smiles):
        if smi == REFERENCE_SMILES["FTA"]:
            ref_idx = i
            break
    sim = tanimoto_to_reference(df, reference_row=ref_idx)
    cov = pocket_residue_contacts(df)
    return pd.DataFrame(
        {
            "smiles": smiles[: len(sim)],
            "plif_similarity": sim,
            "pocket_coverage": cov.to_numpy()[: len(sim)],
            "n_interactions": df.to_numpy().astype(bool).sum(axis=1)[: len(sim)],
        }
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-plif", type=int, default=800)
    args = ap.parse_args()

    FIGURES.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    allo = read_results(GENERATED / "docked_allosteric.jsonl")
    if allo.empty:
        raise SystemExit("no allosteric docking results - run run_screen.py first")
    ortho = read_results(GENERATED / "docked_orthosteric.jsonl")
    surf = read_results(GENERATED / "docked_surface.jsonl")

    print(f"docked: allosteric {len(allo):,} | orthosteric {len(ortho):,} | surface {len(surf):,}")

    merged = merge_sites(allo, ortho, surf)
    lib = pd.read_parquet(GENERATED / "library.parquet")
    merged = merged.merge(
        lib[["smiles", "source", "mw", "clogp", "qed", "sa_score", "max_sim_reference"]],
        on="smiles", how="left",
    )
    merged["source"] = merged["source"].fillna("reference")

    print("computing interaction fingerprints…")
    plif = compute_plif(allo, args.top_plif)
    if not plif.empty:
        merged = merged.merge(plif, on="smiles", how="left")

    ref_rows = merged[merged["smiles"].isin(REFERENCE_SMILES.values())]
    ref_score = float(ref_rows["score_allosteric"].min()) if not ref_rows.empty else None

    scored = apply_hurdles(merged, CFG, reference_score=ref_score)
    ranked = consensus_rank(scored, CFG)
    ranked.to_parquet(OUTPUTS / "ranked_candidates.parquet", index=False)
    ranked.head(200).drop(columns=[c for c in ("pose_pdbqt",) if c in ranked]).to_csv(
        OUTPUTS / "top_candidates.csv", index=False
    )
    funnel = hurdle_funnel(scored)
    funnel.to_csv(OUTPUTS / "funnel.csv", index=False)

    print("plotting…")
    plots.fig_site_comparison(merged, FIGURES / "fig_site_comparison.png")
    plots.fig_generator_quality(merged, FIGURES / "fig_generator_quality.png")
    plots.fig_score_vs_size(merged, FIGURES / "fig_score_vs_size.png")
    plots.fig_funnel(funnel, FIGURES / "fig_funnel.png")
    if "plif_similarity" in merged:
        plots.fig_plif(merged, FIGURES / "fig_plif.png")

    write_results(merged, ranked, funnel, ref_rows, time.time() - t0)
    CFG.to_json(OUTPUTS / "run_config.json")
    print(f"wrote outputs/RESULTS.md in {time.time() - t0:.0f}s")


def write_results(merged, ranked, funnel, ref_rows, runtime) -> None:
    L: list[str] = []
    A = L.append
    A("# Results — de novo design against the TEM-1 allosteric pocket\n")
    A(f"*Generated in {runtime:.0f}s on {platform.node()}. "
      "Docking protocol validated in [VALIDATION.md](VALIDATION.md).*\n")

    A("## What was screened\n")
    A(f"- **{len(merged):,} molecules** docked at the allosteric site")
    for col, name in (("score_orthosteric", "orthosteric (catalytic) site"),
                      ("score_surface", "non-functional surface site")):
        if col in merged:
            A(f"- **{int(merged[col].notna().sum()):,}** also docked at the {name}")
    A("")

    if not ref_rows.empty:
        A("### Reference ligands under the same protocol\n")
        cols = [c for c in ("smiles", "score_allosteric", "le_allosteric",
                            "score_orthosteric", "score_surface", "plif_similarity")
                if c in ref_rows]
        A(ref_rows[cols].to_markdown(index=False))
        A("")

    A("## Per-generator outcome\n")
    grp = merged.groupby("source").agg(
        n=("smiles", "size"),
        median_score=("score_allosteric", "median"),
        best_score=("score_allosteric", "min"),
        median_le=("le_allosteric", "median"),
    ).reset_index()
    A(grp.to_markdown(index=False))

    A("\n## Triage funnel\n")
    A("Each hurdle is independent; a molecule that docks well everywhere is a "
      "sticky compound, not a designed binder.\n")
    A(funnel.to_markdown(index=False))

    n_cand = int(merged.get("is_candidate", pd.Series(dtype=bool)).sum()) if "is_candidate" in merged else 0
    A(f"\n## Top candidates\n")
    show_cols = [c for c in (
        "smiles", "source", "score_allosteric", "le_allosteric", "selectivity",
        "site_specificity", "plif_similarity", "qed", "sa_score",
        "max_sim_reference", "consensus_z") if c in ranked]
    A(ranked.head(CFG.top_n_report)[show_cols].to_markdown(index=False))

    A("\n## Figures\n")
    for f, cap in (
        ("fig_site_comparison.png", "Score distributions at the three sites"),
        ("fig_funnel.png", "Triage funnel"),
        ("fig_generator_quality.png", "Docking outcome by generator"),
        ("fig_score_vs_size.png", "Score vs ligand size"),
        ("fig_plif.png", "Interaction-fingerprint fidelity"),
    ):
        if (FIGURES / f).exists():
            A(f"![{cap}](figures/{f})\n")

    A("\n## Honest limitations\n")
    A("- **Nothing here is experimentally validated.** These are docking-ranked "
      "proposals, not inhibitors.\n"
      "- Vina's scoring function **fails to reproduce the CBT crystal pose** on "
      "the naked receptor and scores that pose below a wrong one; the protocol "
      "only passes validation with crystallographic waters retained. Scores "
      "outside that setup should not be trusted.\n"
      "- The receptor is **rigid**. The TEM-1 allosteric site is a "
      "conformationally responsive, core-disrupting pocket, so rigid docking "
      "cannot capture the induced fit that makes it druggable.\n"
      "- Selectivity and site-specificity are **docking-score differences**, "
      "not measured selectivity.\n"
      "- Generated molecules are novel by fingerprint distance from the two "
      "reference ligands only; no patent or full-ChEMBL novelty check was run.\n")
    (OUTPUTS / "RESULTS.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
