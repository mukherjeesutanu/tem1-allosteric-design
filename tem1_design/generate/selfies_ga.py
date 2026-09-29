"""Generator 2 - SELFIES genetic algorithm: evolutionary string design.

Paradigm: *stochastic evolutionary search in a robust string space*. SELFIES
are constructed so that essentially any token sequence decodes to a valid
molecule, which lets mutation and crossover be brutally simple without
producing chemical garbage at the syntax level (chemical *sense* is still the
filters' job).

The fitness function is pluggable. By default it is a cheap ligand-efficiency
proxy so the GA can run standalone; `run_generate.py` can pass a
docking-derived fitness instead, at much higher cost.
"""
from __future__ import annotations

import random
import re
from typing import Callable

import numpy as np
import selfies as sf
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator

from ..config import REFERENCE_SMILES, Config

RDLogger.DisableLog("rdApp.*")

Fitness = Callable[[list[str]], np.ndarray]


def default_fitness(cfg: Config) -> Fitness:
    """Reward drug-likeness, penalise hard synthesis, and reward *some*
    resemblance to the reference chemotype without copying it.

    This is a deliberately weak prior: the real selection pressure in this
    project is docking, applied downstream.
    """
    import os, sys
    sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, "SA_Score"))
    import sascorer

    gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=cfg.ecfp_radius, fpSize=cfg.ecfp_bits
    )
    refs = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in REFERENCE_SMILES.values()]

    def score(smiles: list[str]) -> np.ndarray:
        out = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if not is_sane(mol):
                out.append(-10.0)
                continue
            qed = Descriptors.qed(mol)
            sa = sascorer.calculateScore(mol)
            sim = max(DataStructs.BulkTanimotoSimilarity(gen.GetFingerprint(mol), refs))
            logp = Crippen.MolLogP(mol)
            # a hydrophobic pocket wants logP ~3-4.5, not maximal logP
            logp_term = -abs(logp - 3.8) / 3.0
            # peak reward at moderate similarity: novel but chemotype-aware
            novelty_term = 1.0 - abs(sim - 0.35) / 0.35
            out.append(
                qed + logp_term + 0.5 * novelty_term - 0.25 * sa
                - 1.5 * strain_penalty(mol)
            )
        return np.asarray(out, dtype=float)

    return score


# SELFIES' fully robust alphabet includes boron, iodine and exotic valences,
# which the GA happily exploits into molecules no chemist would draw. Restrict
# it to the elements that actually populate lead-like space.
DRUG_LIKE_ELEMENTS = ("C", "N", "O", "S", "F", "Cl", "Br", "P")
_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")


def drug_like_alphabet() -> list[str]:
    """Neutral, common-element tokens only.

    Charged tokens are dropped too: the GA otherwise discovers that
    [O+]/[C-] satisfy SELFIES' valence rules and litters the population with
    zwitterions that no medicinal chemist would propose.
    """
    allowed = set(DRUG_LIKE_ELEMENTS)
    out = []
    for token in sorted(sf.get_semantic_robust_alphabet()):
        if "+" in token or "-" in token:
            continue
        elements = set(_ELEMENT_RE.findall(token)) - {"Ring", "Branch", "Expl"}
        if elements <= allowed:
            out.append(token)
    return out


# Motifs that are valence-legal but synthetically or metabolically absurd.
STRAIN_SMARTS = [
    "[*]=[*]=[*]",            # allene / cumulene
    "[#6]=[#16]=[#6]",        # thioketene-like
    "[NX2]=[NX2]=[NX2]",      # azo-cumulene
    "[OX2][OX2]",             # peroxide
    "[NX3][NX3][NX3]",        # triazane
    "[#6]#[#6]#[#6]",         # poly-yne
]
_STRAIN = [Chem.MolFromSmarts(s) for s in STRAIN_SMARTS]


def strain_penalty(mol) -> float:
    return sum(1.0 for patt in _STRAIN if patt is not None and mol.HasSubstructMatch(patt))


def is_sane(mol) -> bool:
    """Hard gate applied during selection.

    Without it the GA drifts into valence-legal but absurd space (exotic ring
    systems, hypervalent sulfur, radical centres). Gating in the fitness keeps
    the population inside chemistry a synthetic chemist would recognise.
    """
    if mol is None or mol.GetNumHeavyAtoms() < 12:
        return False
    if strain_penalty(mol) > 0:
        return False
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        return False
    if Chem.GetFormalCharge(mol) != 0:
        return False
    for ring in mol.GetRingInfo().AtomRings():
        if not 4 <= len(ring) <= 7:
            return False
    # phosphorus and sulfur only in common oxidation environments
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "S" and atom.GetTotalValence() not in (2, 4, 6):
            return False
        if atom.GetSymbol() == "P" and atom.GetTotalValence() not in (3, 5):
            return False
    return True


def _mutate(tokens: list[str], alphabet: list[str], rng: random.Random) -> list[str]:
    t = list(tokens)
    op = rng.random()
    if op < 0.4 and t:                      # substitute
        t[rng.randrange(len(t))] = rng.choice(alphabet)
    elif op < 0.75:                          # insert
        t.insert(rng.randrange(len(t) + 1), rng.choice(alphabet))
    elif t and len(t) > 6:                   # delete
        del t[rng.randrange(len(t))]
    else:
        t.append(rng.choice(alphabet))
    return t


def _crossover(a: list[str], b: list[str], rng: random.Random) -> list[str]:
    if not a or not b:
        return a or b
    i, j = rng.randrange(len(a)), rng.randrange(len(b))
    return a[:i] + b[j:]


def run_ga(
    cfg: Config,
    fitness: Fitness | None = None,
    seeds: list[str] | None = None,
    n_target: int | None = None,
) -> list[str]:
    rng = random.Random(cfg.random_seed)
    fitness = fitness or default_fitness(cfg)
    seeds = seeds or list(REFERENCE_SMILES.values())
    n_target = n_target or cfg.n_per_generator

    alphabet = drug_like_alphabet()
    population = []
    for smi in seeds:
        try:
            population.append(list(sf.split_selfies(sf.encoder(smi))))
        except Exception:
            continue
    if not population:
        raise RuntimeError("no seed could be encoded as SELFIES")
    while len(population) < cfg.ga_population:
        population.append(_mutate(rng.choice(population), alphabet, rng))

    seen: set[str] = set()
    for _ in range(cfg.ga_generations):
        smiles = []
        for toks in population:
            try:
                smiles.append(sf.decoder("".join(toks)))
            except Exception:
                smiles.append("")
        scores = fitness(smiles)
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi) if smi else None
            if is_sane(mol):
                seen.add(Chem.MolToSmiles(mol))

        # tournament selection
        elite_idx = np.argsort(scores)[::-1][: max(2, cfg.ga_population // 10)]
        elites = [population[i] for i in elite_idx]
        new_pop = list(elites)
        while len(new_pop) < cfg.ga_population:
            a, b = rng.choice(elites), rng.choice(population)
            child = _crossover(a, b, rng)
            if rng.random() < 0.8:
                child = _mutate(child, alphabet, rng)
            new_pop.append(child[:120])
        population = new_pop
        if len(seen) >= n_target * 3:
            break

    seen.difference_update(Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in seeds)
    return sorted(seen)
