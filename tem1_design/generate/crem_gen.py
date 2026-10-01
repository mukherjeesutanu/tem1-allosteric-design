"""Generator 1 - CReM: knowledge-based fragment replacement.

Paradigm: *context-constrained enumeration*. Every new fragment comes from a
database of fragments observed in real ChEMBL molecules, matched on the
chemical context that surrounds the attachment point. Nothing is invented, so
output is synthetically plausible by construction -- the opposite failure mode
to a neural sampler.

Seeded from the crystallographic ligands, so the search starts inside the
pocket's known chemotype and walks outward.
"""
from __future__ import annotations

import random
from pathlib import Path

from rdkit import Chem, RDLogger

from ..config import REFERENCE_SMILES, Config

RDLogger.DisableLog("rdApp.*")


def _require_db(cfg: Config) -> Path:
    db = Path(cfg.crem_db)
    if not db.exists():
        raise FileNotFoundError(
            f"CReM fragment database missing at {db}. "
            "Run scripts/download_crem_db.sh first."
        )
    return db


def grow_and_mutate(
    cfg: Config,
    seeds: list[str] | None = None,
    n_target: int | None = None,
) -> list[str]:
    """Mutate, grow and link-decorate the seed ligands.

    Three CReM operations are used because they explore different directions:
      mutate_mol  - swap a fragment for a context-compatible one
      grow_mol    - append a fragment at an implicit hydrogen
      (grow is applied to mutation products too, giving a 2-step walk)
    """
    from crem.crem import grow_mol, mutate_mol

    db = str(_require_db(cfg))
    rng = random.Random(cfg.random_seed)
    seeds = seeds or list(REFERENCE_SMILES.values())
    n_target = n_target or cfg.n_per_generator

    out: set[str] = set()
    frontier = list(seeds)
    generation = 0

    while len(out) < n_target and frontier and generation < 4:
        generation += 1
        next_frontier: list[str] = []
        for smi in frontier:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                muts = list(
                    mutate_mol(
                        mol,
                        db_name=db,
                        radius=cfg.crem_radius,
                        min_freq=cfg.crem_min_freq,
                        max_replacements=200,
                        return_mol=False,
                    )
                )
                grows = list(
                    grow_mol(
                        mol,
                        db_name=db,
                        radius=cfg.crem_radius,
                        min_freq=cfg.crem_min_freq,
                        max_replacements=100,
                        return_mol=False,
                    )
                )
            except Exception:
                continue
            products = [p for p in muts + grows if p]
            rng.shuffle(products)
            out.update(products)
            # carry a sample forward so later generations drift further from
            # the seed without the frontier exploding
            next_frontier.extend(products[:12])
            if len(out) >= n_target:
                break
        frontier = next_frontier[:40]

    out.difference_update(seeds)
    # Truncate by RANDOM sample, not by sorting. `sorted(out)[:n]` keeps the
    # alphabetically-first SMILES, which silently drops every molecule whose
    # canonical string sorts late - sulfur-first, aromatic-first ("c", "n"),
    # bracket-atom-first - and makes the resulting pass-rate a statement about
    # a biased slice rather than about the generator.
    ordered = sorted(out)                      # deterministic starting order
    if len(ordered) <= n_target:
        return ordered
    # fresh RNG so the sample does not depend on how much the search consumed
    return sorted(random.Random(cfg.random_seed).sample(ordered, n_target))
