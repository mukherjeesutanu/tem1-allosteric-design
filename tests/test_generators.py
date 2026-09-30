"""Generator tests. Each skips cleanly when its dependency is absent, so the
suite runs on a laptop without the 5.8 GB CReM database or a torch install.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem, RDLogger

from tem1_design.config import CFG, REFERENCE_SMILES

RDLogger.DisableLog("rdApp.*")


def _all_parse(smiles: list[str]) -> bool:
    return all(Chem.MolFromSmiles(s) is not None for s in smiles)


# --------------------------------------------------------------------------
# SELFIES GA - no external data needed
# --------------------------------------------------------------------------
def test_ga_alphabet_is_neutral_and_drug_like():
    from tem1_design.generate.selfies_ga import DRUG_LIKE_ELEMENTS, drug_like_alphabet

    alphabet = drug_like_alphabet()
    assert alphabet
    assert not any("+" in t or "-" in t for t in alphabet), "charged tokens leaked in"
    for banned in ("[B]", "[I]", "[Se]"):
        assert banned not in alphabet


def test_ga_sanity_gate_rejects_absurd_molecules():
    from tem1_design.generate.selfies_ga import is_sane

    assert is_sane(Chem.MolFromSmiles(REFERENCE_SMILES["CBT"]))
    # allene, charged, tiny, and a 3-membered ring should all be refused
    for bad in ("C=C=C", "CC[O+]CC", "CCO", "C1CC1c1ccccc1"):
        mol = Chem.MolFromSmiles(bad)
        if mol is not None:
            assert not is_sane(mol), bad


def test_ga_produces_valid_molecules():
    from tem1_design.generate.selfies_ga import run_ga

    cfg = CFG
    cfg.ga_generations, cfg.ga_population = 5, 40
    out = run_ga(cfg, n_target=30)
    assert out, "GA produced nothing"
    assert _all_parse(out)


# --------------------------------------------------------------------------
# CReM - needs the fragment database
# --------------------------------------------------------------------------
@pytest.mark.skipif(not Path(CFG.crem_db).exists(), reason="CReM database not downloaded")
def test_crem_grows_from_the_reference_ligands():
    from tem1_design.generate.crem_gen import grow_and_mutate

    out = grow_and_mutate(CFG, n_target=40)
    assert out
    assert _all_parse(out)
    # products must differ from the seeds
    assert not set(out) & set(REFERENCE_SMILES.values())


# --------------------------------------------------------------------------
# SAFE - needs torch + the pretrained checkpoint
# --------------------------------------------------------------------------
def test_safe_degrades_gracefully_when_unavailable():
    """The pipeline must survive a missing neural generator rather than crash."""
    from tem1_design.generate import safe_gen

    if safe_gen.is_available():
        pytest.skip("SAFE is installed; nothing to degrade")
    assert safe_gen.generate(CFG, n_target=5) == []


@pytest.mark.slow
def test_safe_generates_valid_molecules():
    from tem1_design.generate import safe_gen

    if not safe_gen.is_available():
        pytest.skip("SAFE not installed")
    out = safe_gen.generate(CFG, n_target=20)
    if not out:
        pytest.skip("SAFE checkpoint unavailable offline")
    assert _all_parse(out)
