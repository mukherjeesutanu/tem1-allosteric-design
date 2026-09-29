#!/usr/bin/env python
"""Decide which crystallographic ligand copy occupies the allosteric pocket.

Both 1PZO and 1PZP contain two copies of their ligand. Only one per structure
is the allosteric site; the others are surface/lattice positions. This script
makes that call from the structures themselves rather than by assertion.
"""
from pathlib import Path
import json

from tem1_design.config import ALLOSTERIC_LIGANDS, OUTPUTS
from tem1_design.structure import identify_allosteric_copy

if __name__ == "__main__":
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    report = [
        identify_allosteric_copy(Path(f"data/raw/{pdb}.pdb"), resname)
        for pdb, (resname, _) in ALLOSTERIC_LIGANDS.items()
    ]
    (OUTPUTS / "pocket_identification.json").write_text(json.dumps(report, indent=2))
    for r in report:
        print(f"\n{r['pdb']} {r['resname']}")
        for c in r["copies"]:
            print(f"  copy {c['resid']}: {c['n_contacts']} contacts, "
                  f"{c['frac_hydrophobic']:.0%} hydrophobic, "
                  f"{c['dist_ser70']:.1f} A from Ser70-OG")
