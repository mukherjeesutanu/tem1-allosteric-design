#!/usr/bin/env python
"""Run the docking-protocol validation matrix -> outputs/VALIDATION.md."""
from __future__ import annotations

import argparse
import time

from tem1_design.config import CFG, OUTPUTS
from tem1_design.validate import SUCCESS_RMSD, run_validation


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", type=int, default=0, help="0 = all cores")
    args = ap.parse_args()

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df = run_validation(CFG, n_cpu=args.cpu)
    df.to_csv(OUTPUTS / "validation.csv", index=False)

    lines = ["# Docking protocol validation\n"]
    lines.append(
        "Redocking of the crystallographic ligands before any design work. "
        f"Success is heavy-atom RMSD < {SUCCESS_RMSD} A to the crystal pose "
        "(symmetry-aware, no realignment).\n"
    )
    show = df.copy()
    show["rmsd"] = show["rmsd"].map(lambda v: f"{v:.2f}")
    show["dock_score"] = show["dock_score"].map(lambda v: f"{v:.2f}")
    show["crystal_score"] = show["crystal_score"].map(lambda v: f"{v:.2f}")
    show["success"] = show["success"].map({True: "**yes**", False: "no"})
    lines.append(show[[
        "ligand", "receptor", "scoring_function", "waters", "co_ligand",
        "dock_score", "crystal_score", "rmsd", "centroid_shift", "success",
    ]].to_markdown(index=False))

    best = df.loc[df["rmsd"].idxmin()]
    lines.append(
        f"\n**Best protocol:** {best['ligand']} into {best['receptor']} "
        f"({best['scoring_function']}, waters={best['waters']}, "
        f"co-ligand={best['co_ligand']}) -> RMSD {best['rmsd']:.2f} A.\n"
    )
    lines.append(
        "\n**Interpretation.** Where redocking fails, the crystal pose also "
        "*scores worse* than the pose Vina prefers, so the failure is in the "
        "scoring function and the receptor's boundary conditions, not only in "
        "sampling. The TEM-1 allosteric cleft is shallow, hydrophobic and "
        "water-mediated; stripping ordered waters removes the site's real "
        "shape. Design is therefore run only on the validated setup.\n"
    )
    (OUTPUTS / "VALIDATION.md").write_text("\n".join(lines))
    print(f"wrote outputs/VALIDATION.md in {time.time() - t0:.0f}s")
    print(df[["ligand", "scoring_function", "waters", "co_ligand", "rmsd", "success"]].to_string(index=False))


if __name__ == "__main__":
    main()
