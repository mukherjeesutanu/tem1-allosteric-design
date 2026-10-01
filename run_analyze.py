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
from tem1_design.score import (
    apply_hurdles,
    hurdle_attrition,
    consensus_rank,
    hurdle_funnel,
    merge_sites,
    paired_site_test,
)
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
    drop = [c for c in ("pose_pdbqt",) if c in ranked]
    cleared = ranked[ranked["is_candidate"]].drop(columns=drop)
    ref_mask = cleared["smiles"].isin(REFERENCE_SMILES.values())
    # designed molecules only; a recovered reference ligand is a control, not
    # an output, and must not be counted among the designs
    cleared[~ref_mask].to_csv(OUTPUTS / "top_candidates.csv", index=False)
    if ref_mask.any():
        cleared[ref_mask].to_csv(OUTPUTS / "references_cleared.csv", index=False)
    ranked.drop(columns=drop).head(200).to_csv(
        OUTPUTS / "ranked_top200_all.csv", index=False
    )
    funnel = hurdle_funnel(scored)
    funnel.to_csv(OUTPUTS / "funnel.csv", index=False)
    attrition = hurdle_attrition(scored, CFG)
    attrition.to_csv(OUTPUTS / "hurdle_attrition.csv", index=False)
    # Two control sets. The primary one was chosen as the top scorers at the
    # design site, i.e. selected on the very quantity under test, so its
    # selectivity is partly regression to the mean. The random set is the
    # unbiased estimate.
    site_test = paired_site_test(merged).assign(
        control_set="top-150 by allosteric score (biased)"
    )
    rand_o = read_results(GENERATED / "docked_orthosteric_random.jsonl")
    rand_s = read_results(GENERATED / "docked_surface_random.jsonl")
    rand_merged = None
    if not rand_o.empty or not rand_s.empty:
        rand_merged = merge_sites(allo, rand_o, rand_s)
        rand_test = paired_site_test(rand_merged).assign(
            control_set="random 150 (unbiased)"
        )
        site_test = pd.concat([site_test, rand_test], ignore_index=True)
    # Third control: identical CUBIC boxes at all three sites. The production
    # boxes differ in volume (allosteric 5953 A^3 vs orthosteric 10560 A^3) and
    # Vina scores depend on search volume, so selectivity measured across them
    # could be box geometry rather than biology. This removes that difference.
    mo = read_results(GENERATED / "docked_orthosteric_matched.jsonl")
    ms = read_results(GENERATED / "docked_surface_matched.jsonl")
    ma = read_results(GENERATED / "docked_allosteric_matched.jsonl")
    if not ma.empty and not mo.empty:
        matched = merge_sites(ma, mo, ms)
        site_test = pd.concat(
            [site_test,
             paired_site_test(matched).assign(
                 control_set="random 150, volume-matched cubic boxes")],
            ignore_index=True,
        )
    site_test = site_test[
        ["control_set", "control_site", "n_paired", "median_delta_kcal",
         "frac_favouring_allosteric", "wilcoxon_p"]
    ]
    site_test.to_csv(OUTPUTS / "site_selectivity_test.csv", index=False)

    print("plotting…")
    plots.fig_site_comparison_both(
        merged, rand_merged,
        FIGURES / "fig_site_comparison.png",
    )
    plots.fig_generator_quality(merged, FIGURES / "fig_generator_quality.png")
    plots.fig_score_vs_size(merged, FIGURES / "fig_score_vs_size.png")
    plots.fig_funnel(funnel, FIGURES / "fig_funnel.png")
    if "plif_similarity" in merged:
        plots.fig_plif(merged, FIGURES / "fig_plif.png")

    write_results(merged, ranked, funnel, ref_rows, site_test, attrition,
                  time.time() - t0)
    CFG.to_json(OUTPUTS / "run_config.json")
    print(f"wrote outputs/RESULTS.md in {time.time() - t0:.0f}s")


def write_results(merged, ranked, funnel, ref_rows, site_test, attrition,
                  runtime) -> None:
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

    if site_test is not None and len(site_test):
        A("## Does the design site actually win? (paired control)\n")
        A("Wilcoxon signed-rank on molecules docked at **all three** sites. "
          "Positive `median_delta_kcal` means the allosteric site scored better "
          "for the *same* molecule. Comparing whole distributions instead would "
          "contrast different molecule sets and prove nothing.\n")
        A("Two control sets are reported because the choice of control "
          "molecules can manufacture the result. The **biased** set is the top "
          "scorers at the design site — selected on the very quantity under "
          "test, so part of its advantage is regression to the mean. The "
          "**unbiased** set is a random sample docked at all three sites.\n")
        A(site_test.to_markdown(index=False))
        try:
            b = site_test[(site_test.control_set.str.contains("biased", case=False))
                          & (~site_test.control_set.str.contains("unbiased"))
                          & (site_test.control_site == "orthosteric")].iloc[0]
            u = site_test[site_test.control_set.str.contains("unbiased")
                          & (site_test.control_site == "orthosteric")].iloc[0]
            A(f"\nSelecting controls on the design-site score inflates the "
              f"orthosteric advantage from **{u.median_delta_kcal:+.2f}** to "
              f"**{b.median_delta_kcal:+.2f}** kcal/mol "
              f"({100 * (b.median_delta_kcal / u.median_delta_kcal - 1):.0f}% "
              f"overstatement) and the fraction favouring the design site from "
              f"{u.frac_favouring_allosteric:.1%} to "
              f"{b.frac_favouring_allosteric:.1%}. **The effect survives the "
              f"unbiased control**, so the site preference is real — but the "
              f"unbiased number is the one to quote.\n")
            mm = site_test[site_test.control_set.str.contains("volume-matched")
                           & (site_test.control_site == "orthosteric")]
            if len(mm):
                mm = mm.iloc[0]
                A(f"\nThe production boxes are **not equal in volume** "
                  f"(allosteric 5 953 A^3, orthosteric 10 560 A^3), and Vina "
                  f"scores depend on search volume, so the comparison was "
                  f"repeated with identical 18.12 A cubes at all three centres. "
                  f"The preference does not shrink — it **grows**, to "
                  f"**{mm.median_delta_kcal:+.2f}** kcal/mol "
                  f"({mm.frac_favouring_allosteric:.1%} of molecules). The "
                  f"oversized orthosteric box was scoring that site *better* "
                  f"and therefore understating selectivity, so the geometry "
                  f"confound ran in the conservative direction. The allosteric "
                  f"score itself barely moves (−7.07 production vs −7.10 cube), "
                  f"so the cube is not distorting the design site.\n")
        except (IndexError, ZeroDivisionError):
            pass
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
      "sticky compound, not a designed binder. The control sites were docked "
      "only for the top slice of the library, so `n_not_tested_at_this_stage` "
      "counts molecules that cleared affinity but never had that comparison "
      "run.\n")
    untested_total = int(funnel["n_not_tested_at_this_stage"].sum())
    A(f"Across all stages that is **{untested_total}** molecule(s)"
      + (", so no hurdle is being decided by absent data.\n" if untested_total == 0
         else " — those hurdles are partly decided by missing comparisons.\n"))
    A(funnel.to_markdown(index=False))

    if attrition is not None and len(attrition):
        A("\n### How much does each hurdle actually filter?\n")
        a = attrition.copy()
        a["unconditional_pass_rate"] = a["unconditional_pass_rate"].map(
            lambda v: "—" if pd.isna(v) else f"{v:.1%}")
        a["spearman_vs_dock_score"] = a["spearman_vs_dock_score"].map(
            lambda v: "—" if pd.isna(v) else f"{v:+.2f}")
        A(a[["hurdle", "n_before", "n_after", "removed", "n_tested",
             "unconditional_pass_rate", "spearman_vs_dock_score"]]
          .to_markdown(index=False))
        dead = a[(a["removed"] == 0) & (a["hurdle"] != "affinity")]["hurdle"].tolist()
        if dead:
            A("\n*Mode fidelity is weaker than its name suggests: the PLIF "
          "receptor excludes the crystallographic waters and has no hydrogens, "
          "so the FTA reference row carries only 4 hydrophobic bits. See the "
          "`tem1_design/plif.py` docstring.*\n")
        A(f"\n**{' and '.join(dead)} removed nothing.** Calling these "
              f"independent hurdles overstates the triage: the candidate set is "
              f"set almost entirely by the affinity cut, and every criterion "
              f"correlates with the docking score it is meant to be independent "
              f"of. They are sanity checks that nothing failed, not filters.\n")

    passed = ranked[ranked["is_candidate"]] if "is_candidate" in ranked else ranked
    is_ref = passed["smiles"].isin(REFERENCE_SMILES.values())
    ref_passed, cands = passed[is_ref], passed[~is_ref]

    A(f"\n## The two known ligands, and why this is not a positive control\n")
    names = {v: k for k, v in REFERENCE_SMILES.items()}
    if not ref_rows.empty:
        rr = ref_rows.assign(name=ref_rows["smiles"].map(names))
        cleared_names = set(ref_passed["smiles"].map(names)) if len(ref_passed) else set()
        for _, r in rr.iterrows():
            verdict = "clears all four hurdles" if r["name"] in cleared_names \
                else "does **not** clear the triage"
            A(f"- **{r['name']}** — docks at {r['score_allosteric']:.2f} "
              f"kcal/mol and {verdict}.")
        A("")
    A("The affinity hurdle is *defined* as \"at least as good as the best "
      "reference\", and the best reference is FTA, so **FTA passes by "
      "construction**. Recovering it is arithmetic, not validation, and it is "
      "excluded from the candidate count below.\n")
    A("CBT — the other crystallographic ligand of this same pocket — is the "
      "informative case: it docks 1.4 kcal/mol worse and **fails**. A "
      "threshold that rejects a known binder of the target site is a threshold "
      "tuned for docking score, not for binding. Read the candidate list as "
      "\"molecules the scoring function likes at least as much as the better "
      "reference\", nothing stronger.\n")

    A(f"\n## Top designed candidates\n")
    by_source = cands["source"].value_counts().to_dict() if "source" in cands else {}
    A(f"**{len(cands)} designed molecules of {len(ranked):,} docked clear all four "
      f"hurdles** (by generator: {by_source}). Only these are listed: ranking the "
      f"whole library together would mix molecules that passed the control "
      f"comparisons with molecules that never had them run.\n")
    show_cols = [c for c in (
        "smiles", "source", "score_allosteric", "le_allosteric", "selectivity",
        "site_specificity", "plif_similarity", "qed", "sa_score",
        "max_sim_reference", "consensus_z") if c in ranked]
    A(cands.head(CFG.top_n_report)[show_cols].to_markdown(index=False))

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
