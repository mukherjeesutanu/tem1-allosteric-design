"""Consensus triage of docked candidates.

No single number decides a hit here, because validation showed the Vina score
alone is unreliable in this pocket. A candidate has to clear four independent
hurdles:

  1. affinity        - Vina score at the allosteric site, and ligand efficiency
  2. selectivity     - better at the allosteric site than at the catalytic site
  3. site specificity- better there than at a non-functional surface site
  4. mode fidelity   - its pose reproduces the reference interaction pattern

Hurdles 2-4 are the ones that separate a designed binder from a sticky
molecule that docks well everywhere.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Config


def merge_sites(
    allosteric: pd.DataFrame,
    orthosteric: pd.DataFrame | None = None,
    surface: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """One row per molecule with the score at each site."""
    keep = ["smiles", "vina_score", "ligand_efficiency", "heavy_atoms", "status"]
    df = allosteric[allosteric["status"] == "ok"][keep].rename(
        columns={"vina_score": "score_allosteric", "ligand_efficiency": "le_allosteric"}
    )
    for other, tag in ((orthosteric, "orthosteric"), (surface, "surface")):
        if other is None or other.empty:
            continue
        sub = other[other["status"] == "ok"][["smiles", "vina_score"]].rename(
            columns={"vina_score": f"score_{tag}"}
        )
        df = df.merge(sub, on="smiles", how="left")

    if "score_orthosteric" in df:
        # positive = prefers the allosteric site (Vina scores are negative)
        df["selectivity"] = df["score_orthosteric"] - df["score_allosteric"]
    if "score_surface" in df:
        df["site_specificity"] = df["score_surface"] - df["score_allosteric"]
    return df.drop(columns=["status"])


def zscore(series: pd.Series) -> pd.Series:
    sd = series.std(ddof=0)
    return (series - series.mean()) / sd if sd else series * 0.0


def consensus_rank(
    df: pd.DataFrame,
    cfg: Config,
    weights: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Rank by a transparent weighted z-score of the available criteria.

    Weights are explicit rather than fitted: there is no experimental data to
    fit them to, and pretending otherwise would be the usual generative-design
    overclaim.
    """
    weights = weights or {
        "score_allosteric": 0.35,   # lower is better -> negated below
        "le_allosteric": 0.20,
        "selectivity": 0.20,
        "site_specificity": 0.10,
        "plif_similarity": 0.15,
    }
    out = df.copy()
    total = pd.Series(0.0, index=out.index)
    weight_present = pd.Series(0.0, index=out.index)
    n_used = pd.Series(0, index=out.index)
    used = {}
    for col, w in weights.items():
        if col not in out or out[col].isna().all():
            continue
        # Missing criteria are NOT imputed. The control sites were only docked
        # for the top slice of the library, so filling a missing selectivity
        # with the median would hand un-tested molecules a fabricated score and
        # let them outrank molecules that were actually measured. Instead each
        # molecule is scored on the criteria it has, and the weights are
        # renormalised by how much weight was actually available to it.
        v = out[col]
        present = v.notna()
        z = -zscore(v[present]) if col in ("score_allosteric", "le_allosteric") else zscore(v[present])
        total = total.add(w * z.reindex(out.index).fillna(0.0))
        weight_present = weight_present + w * present.astype(float)
        n_used = n_used + present.astype(int)
        used[col] = w

    total_weight = sum(used.values())
    out["consensus_z"] = np.where(
        weight_present > 0, total / weight_present.replace(0, np.nan) * total_weight, np.nan
    )
    out["n_criteria"] = n_used
    out["criteria_weight_fraction"] = weight_present / total_weight if total_weight else 0.0
    out.attrs["weights_used"] = used
    return out.sort_values("consensus_z", ascending=False).reset_index(drop=True)


def apply_hurdles(
    df: pd.DataFrame,
    cfg: Config,
    reference_score: float | None = None,
) -> pd.DataFrame:
    """Flag which hurdles each molecule clears; nothing is silently dropped."""
    out = df.copy()
    out["passes_affinity"] = (
        out["score_allosteric"] <= reference_score
        if reference_score is not None
        else out["score_allosteric"] <= out["score_allosteric"].median()
    )
    out["passes_selectivity"] = (
        out["selectivity"] > 0 if "selectivity" in out else True
    )
    out["passes_specificity"] = (
        out["site_specificity"] > 0 if "site_specificity" in out else True
    )
    out["passes_mode_fidelity"] = (
        out["plif_similarity"] >= cfg.min_plif_similarity
        if "plif_similarity" in out
        else True
    )
    hurdles = [
        "passes_affinity",
        "passes_selectivity",
        "passes_specificity",
        "passes_mode_fidelity",
    ]
    out["n_hurdles_cleared"] = out[hurdles].sum(axis=1)
    out["is_candidate"] = out[hurdles].all(axis=1)
    return out


def hurdle_funnel(df: pd.DataFrame) -> pd.DataFrame:
    """How many molecules survive each successive hurdle.

    `n_not_tested` matters as much as `n`: the control sites were docked only
    for the top slice of the library, so a molecule can fail a hurdle either
    because it lost the comparison or because the comparison was never run.
    Collapsing those two into one number would overstate how discriminating
    the hurdle is.
    """
    a = df["passes_affinity"]
    stages = [
        ("docked", np.ones(len(df), dtype=bool), None),
        ("affinity", a.to_numpy(), None),
        ("+ selectivity", (a & df["passes_selectivity"]).to_numpy(),
         "selectivity"),
        ("+ site specificity",
         (a & df["passes_selectivity"] & df["passes_specificity"]).to_numpy(),
         "site_specificity"),
        ("+ mode fidelity", df["is_candidate"].to_numpy(), "plif_similarity"),
    ]
    rows = []
    for name, mask, col in stages:
        untested = (
            int((a & df[col].isna()).sum()) if col and col in df else 0
        )
        rows.append({"stage": name, "n": int(mask.sum()),
                     "n_not_tested_at_this_stage": untested})
    return pd.DataFrame(rows)


def paired_site_test(df: pd.DataFrame) -> pd.DataFrame:
    """Wilcoxon signed-rank test of design site vs each control, paired.

    Only molecules docked at *all* three sites are used. An unpaired
    comparison here would contrast the whole library at the design site with a
    different, score-selected handful at the controls, which is not a control
    at all. The effect size is the median per-molecule difference, in kcal/mol.
    """
    from scipy import stats

    cols = ["score_allosteric", "score_orthosteric", "score_surface"]
    have = [c for c in cols if c in df]
    paired = df.dropna(subset=have)
    rows = []
    for control in [c for c in have if c != "score_allosteric"]:
        a = paired["score_allosteric"].to_numpy(dtype=float)
        b = paired[control].to_numpy(dtype=float)
        if len(a) < 5:
            continue
        diff = b - a                      # positive = allosteric site is better
        try:
            stat, p = stats.wilcoxon(a, b)
        except ValueError:
            stat, p = float("nan"), float("nan")
        rows.append(
            {
                "control_site": control.replace("score_", ""),
                "n_paired": int(len(a)),
                "median_delta_kcal": float(np.median(diff)),
                "frac_favouring_allosteric": float((diff > 0).mean()),
                "wilcoxon_p": float(p),
            }
        )
    return pd.DataFrame(rows)
