"""Generator 3 - SAFE: pretrained deep generative model.

Paradigm: *learned sequence model*. SAFE (Sequential Attachment-based Fragment
Embedding, Noutahi et al. 2024) reorders SMILES into fragment blocks so that a
GPT-style model can do scaffold decoration, linker design and superstructure
generation as plain autoregressive completion. Unlike CReM it can propose
fragments never seen in any lookup table, and unlike the GA it has a learned
prior over what molecules look like.

Used here in *scaffold-decoration* mode on the shared pharmacophore of the two
crystallographic ligands - an N-linked tetrazole - so the generated series
keeps the anchor that the pocket recognises while the rest is redesigned.
"""
from __future__ import annotations

import warnings

from rdkit import Chem, RDLogger

from ..config import Config

RDLogger.DisableLog("rdApp.*")

# Both reference ligands present a tetrazole attached to an sp2/sp3 linker.
# Decorating from this anchor keeps the recognition element fixed.
DEFAULT_SCAFFOLDS = [
    "[*]c1nnn[nH]1",                 # bare tetrazole anchor
    "[*]N([*])c1nnn[nH]1",           # CBT-like amino-tetrazole
    "[*]Cc1ccc(Cl)cc1",              # CBT chlorobenzyl arm
]


def _load_designer_class():
    """SAFEDesign moved between `safe.sample` and `safe.generative` across
    releases; accept either rather than pinning to one layout."""
    for module in ("safe.sample", "safe.generative"):
        try:
            return __import__(module, fromlist=["SAFEDesign"]).SAFEDesign
        except Exception:
            continue
    return None


def is_available() -> bool:
    return _load_designer_class() is not None


def generate(
    cfg: Config,
    scaffolds: list[str] | None = None,
    n_target: int | None = None,
    model_name: str = "datamol-io/safe-gpt",
    device: str = "cpu",
) -> list[str]:
    """Scaffold decoration with the pretrained SAFE-GPT model.

    Returns [] if the package or checkpoint is unavailable, so the pipeline
    degrades to the other two generators rather than failing.
    """
    if not is_available():
        return []
    SAFEDesign = _load_designer_class()
    if SAFEDesign is None:
        return []

    scaffolds = scaffolds or DEFAULT_SCAFFOLDS
    n_target = n_target or cfg.n_per_generator
    per_scaffold = max(50, n_target // max(len(scaffolds), 1))

    try:
        designer = SAFEDesign.load_default(model_dir=model_name, device=device, verbose=False)
    except Exception as exc:  # network / checkpoint problems
        warnings.warn(f"SAFE model unavailable: {exc}")
        return []

    out: set[str] = set()
    for scaffold in scaffolds:
        try:
            gen = designer.scaffold_decoration(
                scaffold=scaffold,
                n_samples_per_trial=per_scaffold,
                n_trials=1,
                sanitize=True,
                do_not_fragment_further=True,
combine_strategy=None,
            )
        except TypeError:
            # older/newer signature: fall back to the minimal call
            gen = designer.scaffold_decoration(
                scaffold=scaffold, n_samples_per_trial=per_scaffold, sanitize=True
            )
        except Exception as exc:
            warnings.warn(f"SAFE decoration failed for {scaffold}: {exc}")
            continue
        for smi in gen:
            mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else smi
            if mol is not None:
                out.add(Chem.MolToSmiles(mol))
    return sorted(out)
