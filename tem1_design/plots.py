"""Figures. One idea per panel, legible in greyscale."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "figure.facecolor": "white",
    }
)

TEAL, CORAL, SLATE, AMBER, PURPLE = "#2a9d8f", "#e76f51", "#4a5568", "#e9c46a", "#8367c7"
GEN_COLORS = {"crem": TEAL, "ga": CORAL, "safe": PURPLE, "reference": SLATE}


def fig_site_comparison_both(biased: pd.DataFrame, random_set, path: Path) -> None:
    """Both control sets side by side.

    Left: controls chosen as the best scorers at the design site - selected on
    the quantity under test, so part of the gap is regression to the mean.
    Right: a random sample, the unbiased estimate. Showing only the left panel
    would overstate the site preference.
    """
    cols = ["score_allosteric", "score_orthosteric", "score_surface"]
    labels = ["allosteric\n(design)", "orthosteric\n(catalytic)", "surface\n(non-func.)"]
    panels = [(biased, "controls = top 150 by design-site score\n(biased)")]
    if random_set is not None and not random_set.empty:
        panels.append((random_set, "controls = random 150\n(unbiased)"))

    fig, axes = plt.subplots(1, len(panels), figsize=(5.0 * len(panels), 3.9),
                             sharey=True)
    axes = np.atleast_1d(axes)
    for ax, (frame, title) in zip(axes, panels):
        have = [c for c in cols if c in frame]
        paired = frame.dropna(subset=have)
        if paired.empty:
            continue
        x = np.arange(1, len(have) + 1)
        for row in paired[have].to_numpy():
            ax.plot(x, row, color=SLATE, alpha=0.06, lw=0.6, zorder=1)
        parts = ax.violinplot([paired[c] for c in have], positions=x,
                              showmedians=True, widths=0.7)
        for pc, color in zip(parts["bodies"], (TEAL, AMBER, CORAL)):
            pc.set_facecolor(color); pc.set_alpha(0.6); pc.set_zorder(2)
        for key in ("cmedians", "cbars", "cmins", "cmaxes"):
            if key in parts:
                parts[key].set_color(SLATE); parts[key].set_zorder(3)
        ax.set_xticks(x)
        ax.set_xticklabels(labels[: len(have)], fontsize=8)
        for i, c in zip(x, have):
            ax.text(i + 0.07, paired[c].median(), f" {paired[c].median():.2f}",
                    va="center", fontsize=8.5)
        ax.set_title(f"{title}   n = {len(paired)}", loc="left", fontsize=8.5)
    axes[0].set_ylabel("Vina score (kcal/mol)")
    axes[0].invert_yaxis()
    fig.suptitle("Does the design site actually win? (paired, same molecules)",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_site_comparison(df: pd.DataFrame, path: Path) -> None:
    """The control figure: is the allosteric site actually preferred?

    The comparison is PAIRED. The control sites were docked only for a subset
    of the library, so plotting each site's full distribution would compare
    1,222 molecules at the design site against ~150 different ones at the
    controls - three different populations, not a control. Only molecules
    docked at all three sites are shown, and each grey line is one molecule.
    """
    cols = [c for c in ("score_allosteric", "score_orthosteric", "score_surface") if c in df]
    labels = {"score_allosteric": "allosteric\n(design site)",
              "score_orthosteric": "orthosteric\n(catalytic)",
              "score_surface": "surface\n(non-functional)"}
    paired = df.dropna(subset=cols)
    if paired.empty:
        paired = df

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    data = [paired[c] for c in cols]
    x = np.arange(1, len(cols) + 1)

    # per-molecule trajectories
    vals = paired[cols].to_numpy()
    for row in vals:
        ax.plot(x, row, color=SLATE, alpha=0.06, lw=0.6, zorder=1)

    parts = ax.violinplot(data, positions=x, showmedians=True, widths=0.7)
    for pc, color in zip(parts["bodies"], (TEAL, AMBER, CORAL)):
        pc.set_facecolor(color)
        pc.set_alpha(0.6)
        pc.set_zorder(2)
    for key in ("cmedians", "cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_color(SLATE)
            parts[key].set_zorder(3)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[c] for c in cols], fontsize=8)
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.invert_yaxis()
    for i, d in zip(x, data):
        ax.text(i + 0.06, d.median(), f" {d.median():.2f}", va="center", fontsize=8.5)
    ax.set_title(f"Docking score by site — paired, n = {len(paired)} molecules "
                 f"docked at all three", loc="left", fontsize=9.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_generator_quality(df: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2))
    sources = [s for s in ("crem", "ga", "safe") if (df["source"].fillna("") == s).any()]
    for ax, (col, label) in zip(
        axes, (("score_allosteric", "Vina score (kcal/mol)"),
               ("le_allosteric", "ligand efficiency"))
    ):
        data, keep = [], []
        for s in sources:
            v = df.loc[df["source"] == s, col].dropna()
            if len(v):
                data.append(v)
                keep.append(s)
        if not data:
            continue
        bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False)
        for patch, s in zip(bp["boxes"], keep):
            patch.set_facecolor(GEN_COLORS.get(s, SLATE))
            patch.set_alpha(0.7)
        for med in bp["medians"]:
            med.set_color("black")
        ax.set_xticklabels(keep)
        ax.set_ylabel(label)
        ax.invert_yaxis()
    axes[0].set_title("Docking outcome by generator", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_score_vs_size(df: pd.DataFrame, path: Path) -> None:
    """Vina scores scale with ligand size; ligand efficiency is the honest axis."""
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    sub = df.dropna(subset=["score_allosteric", "heavy_atoms"])
    for s in ("crem", "ga", "safe"):
        d = sub[sub["source"] == s]
        if d.empty:
            continue
        ax.scatter(d["heavy_atoms"], d["score_allosteric"], s=7, alpha=0.4,
                   color=GEN_COLORS[s], label=s, edgecolors="none")
    refs = sub[sub["source"] == "reference"]
    if not refs.empty:
        ax.scatter(refs["heavy_atoms"], refs["score_allosteric"], s=90, marker="*",
                   color="black", label="reference ligands", zorder=5)
    if len(sub) > 2:
        r = np.corrcoef(sub["heavy_atoms"], sub["score_allosteric"])[0, 1]
        ax.text(0.98, 0.04, f"Pearson r = {r:.2f}", transform=ax.transAxes,
                ha="right", fontsize=8, color=SLATE)
    ax.set_xlabel("heavy atoms")
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Score tracks ligand size", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_funnel(funnel: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.0))
    y = np.arange(len(funnel))[::-1]
    ax.barh(y, funnel["n"], color=[TEAL, AMBER, AMBER, CORAL, CORAL][: len(funnel)])
    ax.set_yticks(y)
    ax.set_yticklabels(funnel["stage"], fontsize=8.5)
    for yi, n in zip(y, funnel["n"]):
        ax.text(n, yi, f" {n:,}", va="center", fontsize=8)
    ax.set_xlabel("molecules surviving")
    ax.set_title("Triage funnel", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_plif(df: pd.DataFrame, path: Path) -> None:
    sub = df.dropna(subset=["plif_similarity"])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.scatter(sub["plif_similarity"], sub["score_allosteric"], s=9, alpha=0.45,
               color=TEAL, edgecolors="none")
    refs = sub[sub["source"] == "reference"]
    if not refs.empty:
        ax.scatter(refs["plif_similarity"], refs["score_allosteric"], s=95,
                   marker="*", color="black", zorder=5, label="reference")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel("interaction-fingerprint similarity to reference pose")
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.invert_yaxis()
    ax.set_title("Affinity is not binding-mode fidelity", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
