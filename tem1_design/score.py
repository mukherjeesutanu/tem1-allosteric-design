"""Consensus triage of docked candidates.

No single number decides a hit here, because validation showed the Vina score
alone is unreliable in this pocket. A candidate has to clear four hurdles:

  1. affinity        - Vina score at the allosteric site, and ligand efficiency
  2. selectivity     - better at the allosteric site than at the catalytic site
  3. site specificity- better there than at a non-functional surface site
  4. mode fidelity   - its pose reproduces the reference interaction pattern

Hurdles 2-4 are intended to separate a designed binder from a sticky molecule
that docks well everywhere. **On the committed run they largely do not.**
Measured unconditional pass rates, over every molecule that has the datum:

    selectivity        263/267 = 98.5%   -> removed 0 molecules
    site specificity   151/151 = 100.0%  -> removed 0 molecules
    mode fidelity      623/800 = 77.9%   -> removed 2 molecules

and all three correlate with the raw docking score they are meant to be
independent of (Spearman rho = -0.28, -0.43, -0.22 against score_allosteric,
which itself tracks ligand size at rho = -0.48). So the candidate set is set
almost entirely by hurdle 1, and hurdles 2-3 are better described as sanity
checks that nothing failed than as filters. `hurdle_attrition()` recomputes
these numbers for any run rather than trusting this comment.
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

    Molecules missing any weighted criterion are **not ranked at all**. Two
    earlier attempts were both wrong: filling gaps with the column median
    invents data, and renormalising the weights by the fraction present is
    algebraically the same thing - it imputes each gap with the molecule's own
    mean z and then inflates the surviving terms by 1/fraction, so a molecule
    measured on two easy criteria can outrank one measured on five.

    The z-scores are computed over the complete-data subset only: standardising
    each column over a different population (1222 molecules for the docking
    score, 151 for site specificity) makes the columns incommensurable and
    their weighted sum meaningless.
    """
    weights = weights or {
        "score_allosteric": 0.35,   # lower is better -> negated below
        "le_allosteric": 0.20,
        "selectivity": 0.20,
        "site_specificity": 0.10,
        "plif_similarity": 0.15,
    }
    out = df.copy()
    cols = [c for c in weights if c in out and not out[c].isna().all()]
    complete = out[cols].notna().all(axis=1) if cols else pd.Series(False, index=out.index)
    out["n_criteria"] = out[cols].notna().sum(axis=1) if cols else 0
    out["fully_scored"] = complete

    total = pd.Series(np.nan, index=out.index)
    sub = out.loc[complete]
    if len(sub):
        acc = pd.Series(0.0, index=sub.index)
        for col in cols:
            z = zscore(sub[col])
            if col in ("score_allosteric", "le_allosteric"):
                z = -z                      # more negative is better
            acc = acc + weights[col] * z
        total.loc[sub.index] = acc
    out["consensus_z"] = total
    out.attrs["weights_used"] = {c: weights[c] for c in cols}
    out.attrs["n_fully_scored"] = int(complete.sum())
    # unranked molecules sort to the bottom and keep a NaN score, never a
    # fabricated one
    return out.sort_values(
        "consensus_z", ascending=False, na_position="last"
    ).reset_index(drop=True)



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
    # An ABSENT column means the control was never run, which is not evidence
    # that the molecule passed. Defaulting these to True would report every
    # affinity passer as a four-hurdle candidate on a run with no controls.
    out["passes_selectivity"] = (
        out["selectivity"] > 0 if "selectivity" in out else False
    )
    out["passes_specificity"] = (
        out["site_specificity"] > 0 if "site_specificity" in out else False
    )
    out["passes_mode_fidelity"] = (
        out["plif_similarity"] >= cfg.min_plif_similarity
        if "plif_similarity" in out
        else False
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


def hurdle_attrition(df: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """How discriminating each hurdle actually is, measured not assumed.

    `removed` is what the hurdle costs *given the ones before it*; `pass_rate`
    is unconditional over every molecule that has the datum. A hurdle with a
    100% pass rate is not filtering anything, however sensible it sounds.
    """
    from scipy import stats

    spec = [
        ("affinity", "passes_affinity", "score_allosteric", None),
        ("selectivity", "passes_selectivity", "selectivity", 0.0),
        ("site specificity", "passes_specificity", "site_specificity", 0.0),
        ("mode fidelity", "passes_mode_fidelity", "plif_similarity",
         cfg.min_plif_similarity),
    ]
    rows, running = [], pd.Series(True, index=df.index)
    for name, flag, col, thr in spec:
        before = int(running.sum())
        running = running & df[flag]
        tested = df[col].notna().sum() if col in df else 0
        passed = int(df.loc[df[col].notna(), flag].sum()) if col in df else 0
        rho = p = float("nan")
        if col in df and col != "score_allosteric":
            m = df[["score_allosteric", col]].dropna()
            if len(m) > 2:
                rho, p = stats.spearmanr(m["score_allosteric"], m[col])
        rows.append({
            "hurdle": name,
            "n_before": before,
            "n_after": int(running.sum()),
            "removed": before - int(running.sum()),
            "n_tested": int(tested),
            "unconditional_pass_rate": passed / tested if tested else float("nan"),
            "spearman_vs_dock_score": rho,
            "p": p,
        })
    return pd.DataFrame(rows)


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
