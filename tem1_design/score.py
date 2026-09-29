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
    used = {}
    for col, w in weights.items():
        if col not in out or out[col].isna().all():
            continue
        v = out[col].fillna(out[col].median())
        # Vina scores and LE are "more negative is better"
        z = -zscore(v) if col in ("score_allosteric", "le_allosteric") else zscore(v)
        total = total + w * z
        used[col] = w
    out["consensus_z"] = total
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
    """How many molecules survive each successive hurdle."""
    stages = [
        ("docked", np.ones(len(df), dtype=bool)),
        ("affinity", df["passes_affinity"].to_numpy()),
        ("+ selectivity", (df["passes_affinity"] & df["passes_selectivity"]).to_numpy()),
        (
            "+ site specificity",
            (df["passes_affinity"] & df["passes_selectivity"] & df["passes_specificity"]).to_numpy(),
        ),
        ("+ mode fidelity", df["is_candidate"].to_numpy()),
    ]
    return pd.DataFrame(
        [{"stage": name, "n": int(mask.sum())} for name, mask in stages]
    )
