"""Tests that do not need the 5.8 GB CReM database or a docking run.

The important ones are the pocket-identification and scoring-direction tests:
both encode facts that, if silently wrong, would invalidate every downstream
number while still producing a plausible-looking report.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from tem1_design.config import (
    ALLOSTERIC_LIGANDS,
    ALLOSTERIC_POCKET_RESIDUES,
    CFG,
    REFERENCE_SMILES,
)
from tem1_design.filters import describe, standardize, triage
from tem1_design.score import consensus_rank, merge_sites
from tem1_design.structure import (
    ca_rmsd,
    contact_residues,
    distance_to_catalytic_serine,
    identify_allosteric_copy,
    ligand_atoms,
    read_pdb,
)

RAW = Path("data/raw")
HAVE_PDB = (RAW / "1PZO.pdb").exists() and (RAW / "1PZP.pdb").exists()
needs_pdb = pytest.mark.skipif(not HAVE_PDB, reason="PDB files not downloaded")


# --------------------------------------------------------------------------
# structure / pocket
# --------------------------------------------------------------------------
@needs_pdb
def test_reference_ligands_present():
    for pdb, (resname, resid) in ALLOSTERIC_LIGANDS.items():
        atoms = read_pdb(RAW / f"{pdb}.pdb")
        assert ligand_atoms(atoms, resname, resid), f"{resname}{resid} missing from {pdb}"


@needs_pdb
def test_allosteric_copy_is_distal_and_hydrophobic():
    """The designated copy must match the published site: far from the
    catalytic serine and in a hydrophobic lining. The other copy must not."""
    report = identify_allosteric_copy(RAW / "1PZO.pdb", "CBT")
    copies = {c["resid"]: c for c in report["copies"]}
    chosen = copies[ALLOSTERIC_LIGANDS["1PZO"][1]]
    assert chosen["dist_ser70"] > 12.0
    assert chosen["frac_hydrophobic"] > 0.6


@needs_pdb
def test_pocket_residues_actually_contact_the_ligand():
    atoms = read_pdb(RAW / "1PZP.pdb")
    resname, resid = ALLOSTERIC_LIGANDS["1PZP"]
    contacts = {ri for ri, _ in contact_residues(atoms, ligand_atoms(atoms, resname, resid))}
    core = set(ALLOSTERIC_POCKET_RESIDUES)
    # the documented core must be mostly recovered from the structure itself
    assert len(core & contacts) / len(core) >= 0.7


@needs_pdb
def test_structures_share_a_frame():
    """1PZO and 1PZP are used together without superposition, so they must
    already be in the same crystal frame."""
    assert ca_rmsd(RAW / "1PZO.pdb", RAW / "1PZP.pdb") < 1.0


@needs_pdb
def test_sites_are_well_separated():
    """The boxes must be distinct searches, not two views of one pocket.

    Centre distance alone does not establish that - two large boxes 20 A apart
    still overlap - so the overlap VOLUME is measured. A small overlap is
    tolerable (the shared slab here is a fraction of a percent); a large one
    would make the selectivity control meaningless.
    """
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from run_screen import SITES

    def overlap_volume(b1, b2) -> float:
        v = 1.0
        for c1, s1, c2, s2 in zip(b1.center, b1.size, b2.center, b2.size):
            lo = max(c1 - s1 / 2, c2 - s2 / 2)
            hi = min(c1 + s1 / 2, c2 + s2 / 2)
            v *= max(0.0, hi - lo)
        return v

    allo = SITES["allosteric"]
    allo_vol = float(np.prod(allo.size))
    for name in ("orthosteric", "surface"):
        other = SITES[name]
        assert np.linalg.norm(
            np.array(allo.center) - np.array(other.center)
        ) > 15.0, name
        frac = overlap_volume(allo, other) / allo_vol
        assert frac < 0.05, f"{name} box overlaps the allosteric box by {frac:.1%}"


@needs_pdb
def test_allosteric_site_is_distal_from_catalytic_serine():
    atoms = read_pdb(RAW / "1PZP.pdb")
    lig = ligand_atoms(atoms, *ALLOSTERIC_LIGANDS["1PZP"])
    assert distance_to_catalytic_serine(atoms, lig) > 15.0


# --------------------------------------------------------------------------
# filters
# --------------------------------------------------------------------------
def test_reference_ligands_are_chemically_sane():
    for name, smi in REFERENCE_SMILES.items():
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, name
        d = describe(smi)
        assert 250 < d["mw"] < 400
        assert d["sa_score"] < 4.0


def test_standardize_strips_salt():
    out = standardize("CC(=O)Oc1ccccc1C(=O)[O-].[Na+]")
    assert out and "Na" not in out


def test_triage_flags_reference_acrylonitrile():
    """FTA carries an acrylonitrile Michael acceptor; the alert catalogue must
    catch it. If this stops firing, the filters have been silently weakened."""
    df = triage([REFERENCE_SMILES["FTA"]], CFG, source="test")
    assert df.loc[0, "alerts"]


def test_triage_rejects_garbage_without_crashing():
    df = triage(["not_a_smiles", "C", REFERENCE_SMILES["CBT"]], CFG, source="test")
    assert len(df) >= 1
    assert df["smiles"].map(lambda s: Chem.MolFromSmiles(s) is not None).all()


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def _fake_docking(scores, tag):
    return pd.DataFrame(
        {
            "smiles": [f"C{i}" for i in range(len(scores))],
            "vina_score": scores,
            "ligand_efficiency": [s / 20 for s in scores],
            "heavy_atoms": [20] * len(scores),
            "status": ["ok"] * len(scores),
        }
    )


def test_selectivity_sign_convention():
    """Vina scores are negative; a molecule better at the allosteric site must
    get POSITIVE selectivity. A sign flip here would invert every conclusion."""
    allo = _fake_docking([-9.0, -6.0], "a")
    ortho = _fake_docking([-6.0, -9.0], "o")
    m = merge_sites(allo, ortho)
    assert m.loc[0, "selectivity"] > 0   # better allosterically
    assert m.loc[1, "selectivity"] < 0   # better orthosterically


def test_consensus_prefers_stronger_binder():
    df = pd.DataFrame(
        {
            "smiles": ["a", "b"],
            "score_allosteric": [-9.5, -5.0],
            "le_allosteric": [-0.45, -0.20],
            "selectivity": [2.0, -1.0],
            "site_specificity": [1.5, -0.5],
            "plif_similarity": [0.8, 0.1],
        }
    )
    ranked = consensus_rank(df, CFG)
    assert ranked.iloc[0]["smiles"] == "a"


# --- consensus ranking must not invent data -------------------------------
def test_consensus_rank_does_not_rank_incomplete_molecules():
    """A molecule with no control docking must not be ranked at all.

    The fixture is deliberately adversarial: the untested molecule has the
    BEST docking score, so any scheme that imputes its missing criteria - with
    the column median, or equivalently by renormalising the weights - lets it
    outrank molecules that were actually measured. An earlier version of this
    test gave all three the same score, which could not detect that.
    """
    import numpy as np
    import pandas as pd
    from tem1_design.config import CFG
    from tem1_design.score import consensus_rank

    df = pd.DataFrame(
        {
            "smiles": ["measured_bad", "untested_but_best_score", "measured_good"],
            "score_allosteric": [-8.0, -9.5, -8.2],
            "le_allosteric": [-0.35, -0.45, -0.36],
            "selectivity": [-2.0, None, 2.0],
            "site_specificity": [-2.0, None, 2.0],
            "plif_similarity": [0.5, 0.9, 0.6],
        }
    )
    out = consensus_rank(df, CFG)
    untested = out[out["smiles"] == "untested_but_best_score"].iloc[0]
    assert np.isnan(untested["consensus_z"]), "incomplete molecule was scored"
    assert not untested["fully_scored"]
    # it must not occupy a ranked position ahead of a measured molecule
    ranked = out[out["consensus_z"].notna()]["smiles"].tolist()
    assert ranked == ["measured_good", "measured_bad"], ranked


def test_absent_control_columns_do_not_pass_hurdles():
    """With no control docking at all, nothing may be called a candidate."""
    import pandas as pd
    from tem1_design.config import CFG
    from tem1_design.score import apply_hurdles

    df = pd.DataFrame({"smiles": ["x"], "score_allosteric": [-9.0],
                       "le_allosteric": [-0.4]})
    out = apply_hurdles(df, CFG, reference_score=-8.0)
    assert bool(out["passes_affinity"].iloc[0]) is True
    assert bool(out["is_candidate"].iloc[0]) is False


def test_untested_molecule_cannot_become_a_candidate():
    """Missing control data must fail the hurdle, never pass it silently."""
    import pandas as pd
    from tem1_design.config import CFG
    from tem1_design.score import apply_hurdles

    df = pd.DataFrame(
        {
            "smiles": ["untested"],
            "score_allosteric": [-9.0],
            "le_allosteric": [-0.4],
            "selectivity": [None],
            "site_specificity": [None],
            "plif_similarity": [0.9],
        }
    )
    out = apply_hurdles(df, CFG, reference_score=-8.0)
    assert bool(out["passes_affinity"].iloc[0]) is True
    assert bool(out["passes_selectivity"].iloc[0]) is False
    assert bool(out["is_candidate"].iloc[0]) is False


def test_funnel_separates_failed_from_untested():
    import pandas as pd
    from tem1_design.config import CFG
    from tem1_design.score import apply_hurdles, hurdle_funnel

    df = pd.DataFrame(
        {
            "smiles": ["a", "b"],
            "score_allosteric": [-9.0, -9.0],
            "le_allosteric": [-0.4, -0.4],
            "selectivity": [-1.0, None],   # a failed, b never tested
            "site_specificity": [1.0, None],
            "plif_similarity": [0.9, 0.9],
        }
    )
    f = hurdle_funnel(apply_hurdles(df, CFG, reference_score=-8.0))
    row = f[f["stage"] == "+ selectivity"].iloc[0]
    assert row["n"] == 0
    assert row["n_not_tested_at_this_stage"] == 1
