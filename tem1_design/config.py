"""Configuration for de novo design against the TEM-1 allosteric pocket.

Target
------
TEM-1 beta-lactamase carries a well-validated allosteric site ~19 A from the
catalytic Ser70, between helix H11 and the Omega-loop. Horn & Shoichet
(J. Mol. Biol. 2004, 336:1283) solved two complexes with non-beta-lactam
inhibitors that bind there and block the conformational change catalysis
needs -- "core disruption":

    1PZO   CBT   N,N-bis(4-chlorobenzyl)-1H-tetrazol-5-amine
    1PZP   FTA   3-(4-phenylamino-phenylamino)-2-(1H-tetrazol-5-yl)acrylonitrile

Both crystals contain two copies of the ligand. Only one copy per structure
sits in the allosteric pocket; the others are surface/lattice sites. The copies
used here were selected by contact analysis (see `scripts/identify_pocket.py`)
and agree with the published site definition.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RAW = DATA / "raw"
PREPARED = DATA / "prepared"
GENERATED = DATA / "generated"
OUTPUTS = ROOT / "outputs"
FIGURES = OUTPUTS / "figures"

# Ligand copy that occupies the allosteric pocket, per structure.
ALLOSTERIC_LIGANDS = {"1PZO": ("CBT", 301), "1PZP": ("FTA", 300)}
# Copies that do NOT: kept documented so the choice is auditable.
NON_POCKET_LIGANDS = {"1PZO": ("CBT", 300), "1PZP": ("FTA", 301)}

REFERENCE_SMILES = {
    # from the PDB chemical component dictionary
    "CBT": "Clc1ccc(CN(Cc2ccc(Cl)cc2)c2nnn[nH]2)cc1",
    "FTA": "N#C/C(=C\\Nc1ccc(Nc2ccccc2)cc1)c1nnn[nH]1",
}

# Pocket lining shared by both complexes (contact analysis, 4.5 A).
# Strikingly hydrophobic, which is why the site tolerates the greasy
# bis-chlorobenzyl chemotype.
ALLOSTERIC_POCKET_RESIDUES = (
    221, 224, 225, 261, 263, 279, 280, 283, 286,
)
CATALYTIC_RESIDUES = (70, 73, 130, 166, 170, 234, 237)


@dataclass
class Box:
    center: tuple[float, float, float]
    size: tuple[float, float, float]


@dataclass
class Config:
    # --- validated docking protocol -------------------------------------
    # Established by tem1_design/validate.py, not assumed. Redocking succeeds
    # (0.95 A) only on 1PZP with crystallographic waters retained; the naked
    # protein fails on both ligands. Design therefore runs against 1PZP+waters.
    design_pdb: str = "1PZP"
    design_ligand: tuple = ("FTA", 300)
    keep_crystal_waters: bool = True
    scoring_function: str = "vina"
    validated_redock_rmsd: float = 0.95

    # --- structures -----------------------------------------------------
    holo_pdb: str = "1PZO"           # design template (CBT complex)
    cross_pdb: str = "1PZP"          # second holo, used for cross-docking control
    apo_pdb: str = "1BTL"            # apo TEM-1, pocket-collapse control
    receptor_pdb: Path = RAW / "1PZO.pdb"

    # Boxes derived from the union of both allosteric ligand copies (+5 A).
    # 1PZO and 1PZP share a crystal frame (0.57 A CA-RMSD as deposited, no
    # superposition applied), so boxes transfer between them unchanged.
    # Box around the FTA(300) crystal pose + 5 A - the pose the protocol was
    # validated to reproduce. Deliberately NOT the union of both ligand copies:
    # the wider box lets Vina drift into an adjacent surface groove.
    allosteric_box: Box = field(
        default_factory=lambda: Box((3.51, 13.78, 17.65), (14.4, 21.2, 19.5))
    )
    # Catalytic site, for the selectivity control. Centres are 22.3 A apart.
    # NOTE: this box is 10560 A^3 vs the allosteric box's 5953 A^3. Vina scores
    # depend on search volume, so run_screen.py --cubic-box repeats the
    # comparison with identical cubes; see outputs/RESULTS.md.
    orthosteric_box: Box = field(
        default_factory=lambda: Box((9.80, 7.31, 38.07), (24.0, 20.0, 22.0))
    )

    # --- docking --------------------------------------------------------
    exhaustiveness: int = 16
    n_poses: int = 5
    vina_seed: int = 42
    n_cpu: int = 0                   # 0 = all available

    # --- generation -----------------------------------------------------
    crem_db: Path = RAW / "chembl22_sa25_hac12.db"
    crem_radius: int = 3
    crem_min_freq: int = 10
    n_per_generator: int = 2000
    ga_generations: int = 30
    ga_population: int = 200
    random_seed: int = 42

    # --- property filters ----------------------------------------------
    # Deliberately lead-like, not drug-like: the reference ligands are small
    # and greasy, and the pocket is buried and hydrophobic.
    mw_range: tuple[float, float] = (200.0, 480.0)
    clogp_max: float = 5.5
    hbd_max: int = 3
    hba_max: int = 8
    rotb_max: int = 8
    tpsa_max: float = 120.0
    qed_min: float = 0.30
    sa_max: float = 4.5              # SA score, 1 (easy) - 10 (hard)
    max_ring_size: int = 7
    require_pains_clean: bool = True

    # --- novelty --------------------------------------------------------
    novelty_max_tanimoto: float = 0.6   # vs the two reference ligands
    ecfp_radius: int = 2
    ecfp_bits: int = 2048

    # --- selection ------------------------------------------------------
    top_n_report: int = 25
    min_plif_similarity: float = 0.3

    outputs: Path = OUTPUTS
    figures: Path = FIGURES

    def to_json(self, path: Path) -> None:
        d = {}
        for k, v in asdict(self).items():
            d[k] = str(v) if isinstance(v, Path) else v
        path.write_text(json.dumps(d, indent=2, default=str))


CFG = Config()
