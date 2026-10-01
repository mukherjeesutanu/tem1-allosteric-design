# tem1-allosteric-design — de novo design against the TEM-1 allosteric pocket

Three different generative paradigms, one experimentally validated allosteric
site, and a docking protocol that has to **earn the right to be used** before
any molecule is designed.

> **Attribution.** The scientific design of this project is my own. I used
> Claude Code (Anthropic) as a coding and documentation assistant. Every
> number below was produced by the code in this repository and can be
> reproduced with the commands given.

---

## The target

TEM-1 β-lactamase — the archetypal antibiotic-resistance enzyme — has an
allosteric pocket about 19 Å from the catalytic Ser70, between helix H11 and
the Ω-loop. Horn & Shoichet
([*J. Mol. Biol.* **336**:1283, 2004](https://doi.org/10.1016/j.jmb.2003.12.068))
showed that non-β-lactam ligands binding there inhibit by **core disruption**:
they block the conformational change catalysis requires, rather than competing
with substrate.

| PDB | Ligand | Chemistry |
|---|---|---|
| [1PZO](https://www.rcsb.org/structure/1PZO) | **CBT** | N,N-bis(4-chlorobenzyl)-1H-tetrazol-5-amine |
| [1PZP](https://www.rcsb.org/structure/1PZP) | **FTA** | 3-(4-phenylamino-phenylamino)-2-(1H-tetrazol-5-yl)acrylonitrile |

Both crystals contain **two copies** of their ligand, and only one per
structure is the allosteric site. `scripts/identify_pocket.py` makes that call
from the structures themselves — contact count, hydrophobicity of the lining
and distance from Ser70 — rather than by assertion:

```
1PZO CBT  copy 300: 12 contacts, 58% hydrophobic, 10.2 Å from Ser70-OG
         copy 301: 12 contacts, 75% hydrophobic, 17.6 Å from Ser70-OG   <- allosteric
1PZP FTA  copy 300: 16 contacts, 88% hydrophobic, 21.6 Å from Ser70-OG  <- allosteric
         copy 301:  6 contacts, 50% hydrophobic, 29.3 Å from Ser70-OG   (lattice site)
```

The shared lining is strikingly greasy — Leu221, Ala224, Leu225, Val261,
Ile263, Ile279, Ala280, Gly283, Leu286 — which is why the site tolerates the
bis-chlorobenzyl chemotype.

---

## The part most pipelines skip: validating the docking

Before designing anything, the protocol was tested by redocking the
crystallographic ligands. It **failed** in the obvious setup, and the failure
mode is instructive ([full matrix](outputs/VALIDATION.md)):

| Ligand | Receptor | Waters | Co-ligand | Vina score | Crystal-pose score | RMSD | Success |
|---|---|---|---|---|---|---|---|
| FTA300 | 1PZP | ✗ | – | −8.08 | −5.71 | 3.52 Å | no |
| **FTA300** | **1PZP** | **✓** | – | **−8.44** | −6.52 | **0.98 Å** | **yes** |
| CBT301 | 1PZO | ✓ | – | −8.04 | −6.66 | 7.80 Å | no |
| CBT301 | 1PZO | ✓ | CBT300 | −8.42 | −7.60 | 3.11 Å | no |
| CBT300 | 1PZO | ✓ | CBT301 | −7.80 | −7.25 | 1.73 Å | yes¹ |

**16 conditions, 2 successes — but only 1 of them is at the design site.**
¹CBT300 is the *non*-allosteric copy: it sits 10.2 Å from Ser70 and 15 of its
22 atoms fall inside the catalytic control box. Counting it as validation of
the allosteric protocol would be wrong. Restricted to the allosteric copies
(FTA300, CBT301), the matrix is **1 success in 12 conditions**. That single
condition is the one the design run uses, and one in twelve is a thin basis —
stated plainly rather than averaged away.

Three further readings, each computed from the matrix rather than asserted
(`run_validate.py` recomputes them, after an earlier hand-written version of
this section got two of them wrong):

1. **Waters are decisive for FTA, but not a general effect.** FTA300 goes from
   3.52 Å to 0.98 Å purely by keeping them — which is why the design receptor
   keeps them. Across all matched comparisons at the allosteric copies, though,
   waters improve RMSD in only 2 of 6 (median change −0.00 Å).
2. **Vina is better on balance, not everywhere.** It wins 5 of 8 matched
   conditions; Vinardo wins 3. Vina is adopted because it produces the only
   sub-2 Å pose at the design site.
3. **Failure is mostly scoring, partly sampling.** Against the un-minimised
   crystal pose, the crystal scores worse than Vina's choice in 14 of 14
   failures. That comparison is unfair — the crystal pose was never relaxed in
   the Vina field — and against the *minimised* crystal pose it holds in 12 of
   14. So in 2 cases the correct pose does score better and more search could
   have found it.

Design runs **only** on the validated setup: 1PZP with crystallographic waters,
Vina, box around the FTA site.

---

## Three generators, deliberately different

| | Paradigm | What it can and cannot invent |
|---|---|---|
| **CReM** | knowledge-based fragment replacement | every fragment comes from a ChEMBL-derived database, matched on chemical context — synthetically plausible by construction, but cannot leave known fragment space |
| **SELFIES GA** | evolutionary search in a robust string space | mutation/crossover on a representation where almost any token string decodes to a valid molecule — maximal freedom, minimal chemical taste |
| **SAFE-GPT** | pretrained autoregressive sequence model | a learned prior over drug-like chemistry; used in scaffold-decoration mode from the shared tetrazole anchor |

All three are seeded from, or anchored on, the crystallographic ligands, then
put through the identical filter cascade.

### Result: the paradigms differ enormously in usable yield

| Generator | Raw | Unique valid | Passed filters | **Pass rate** | Median QED | Median SA | Median MW |
|---|---|---|---|---|---|---|---|
| CReM | 3000 | 3000 | 1826 | **60.9%** | 0.71 | 2.09 | 341 |
| SAFE-GPT | 2206 | 2147 | 500 | **23.3%** | 0.61 | 3.51 | 247 |
| SELFIES GA | 1206 | 1206 | 93 | **7.7%** | 0.50 | 4.20 | 319 |

**6,352 unique molecules, 2,419 pass every filter.** The ordering —
knowledge-based ≫ learned prior ≫ naive evolutionary — is a yield statement,
not a quality-of-best-hit statement.

⚠️ **Caveat on CReM's row.** The committed run truncated CReM's output with
`sorted(out)[:3000]`, which keeps the alphabetically-first SMILES and silently
drops everything sorting after `O=…` — all sulfur-first, aromatic-first and
bracket-atom-first canonical strings. CReM's 60.9% is therefore measured on a
biased slice of its output. The code now samples randomly instead
(`crem_gen.py`), but the committed numbers predate that fix and the ordering
should be treated as indicative until the generation is re-run.

A sanity check on the filters themselves: they flag **FTA**, one of the two
literature ligands, for its acrylonitrile Michael acceptor. That is chemically
correct, and a useful reminder that "passes the filters" and "is a good
molecule" are different claims.

---

## Triage: four independent hurdles

Because the validation showed the Vina score alone cannot be trusted in this
pocket, no candidate is selected on score:

1. **Affinity** — Vina score and ligand efficiency at the allosteric site
2. **Selectivity** — better there than at the catalytic (orthosteric) site
3. **Site specificity** — better there than at the second crystallographic FTA
   site, a surface position 29 Å from Ser70. This is a deliberately hard
   control: it is a *real* ligand-binding depression (FTA occupies it in the
   1PZP lattice), just one with no catalytic or allosteric role — so beating it
   means more than beating a random patch of surface.
4. **Binding-mode fidelity** — its ProLIF interaction fingerprint reproduces
   the reference ligand's contact pattern

Hurdles 2–4 are *intended* to separate a designed binder from a greasy
molecule that docks well everywhere. **Measured, they largely do not:**

| hurdle | removed | pass rate where tested | Spearman vs dock score |
|---|---|---|---|
| 1 affinity | 1,176 | 3.8% | — |
| 2 selectivity | **0** | 98.5% | −0.28 |
| 3 site specificity | **0** | 100% | −0.43 |
| 4 mode fidelity | 2 | 77.9% | −0.22 |

Hurdles 2 and 3 eliminated nothing, and every criterion correlates with the
docking score it is supposed to be independent of (which in turn tracks ligand
size, ρ = −0.48). The candidate set is set almost entirely by the affinity cut;
the rest are sanity checks that nothing failed. Reporting them as four
independent filters would overstate the triage, so the numbers are printed in
`outputs/hurdle_attrition.csv` on every run.

### Does the design site actually win?

This is the claim the whole project rests on, so it is tested **paired** — the
same molecule docked at all three sites — and with two different control sets,
because the choice of control molecules can manufacture the answer.

![Paired site comparison](outputs/figures/fig_site_comparison.png)

| control set | vs orthosteric | vs surface |
|---|---|---|
| top 150 by design-site score (**biased**) | +1.42 kcal/mol, 99.3% favour | +2.82, 100% |
| random 150 (**unbiased**) | **+0.84 kcal/mol, 86.8% favour** | **+2.17, 100%** |
| random 150, **volume-matched cubic boxes** | **+1.32 kcal/mol, 89.5% favour** | **+2.19, 100%** |

Two confounds, tested rather than argued:

**Selection.** Choosing the control molecules by their design-site score — the
obvious thing to do, since those are the molecules you care about — selects on
the very quantity under test and **overstates the orthosteric advantage by
68%**. The preference survives the unbiased control (p = 1×10⁻²⁰).

**Box geometry.** The production boxes are not equal in volume (allosteric
5,953 Å³ vs orthosteric 10,560 Å³), and Vina scores depend on search volume, so
the comparison was repeated with identical 18.12 Å cubes at all three centres.
The preference does not shrink — it **grows to +1.32 kcal/mol**. The oversized
orthosteric box had been scoring that site *better*, so the geometry confound
was running in the conservative direction. The allosteric score itself barely
moves (−7.07 production vs −7.10 cube).

A paired exhaustiveness test (8 vs 16) at all three sites changes the median
score by 0.000 kcal/mol, so search effort is not driving this either.

**The honest caveat is the surface control, not the orthosteric one.** The
surface box is the least enclosed of the three (56% free volume vs 28% at the
design site), so beating it is close to automatic — 100% of molecules do, under
every control set. The orthosteric comparison is the one carrying weight.

### Outcome

| | |
|---|---|
| Docked at the allosteric site | 1,222 |
| Clear the affinity hurdle (≥ FTA, −8.42) | 46 |
| Clear all four hurdles | 44 |
| — of which a **known ligand recovered** (FTA) | 1 |
| — **designed** candidates | **43**, all from CReM |

**This is not a positive control.** The affinity hurdle is defined as "at least
as good as the best reference", and the best reference is FTA — so FTA passes
by construction. The informative case is **CBT, the other crystallographic
ligand of this same pocket, which docks 1.4 kcal/mol worse and fails the
triage.** A bar that rejects a known binder of the target site is a bar tuned
to the docking score, not to binding. Read the candidate list as "molecules the
scoring function likes at least as much as the better reference".

Neither the SELFIES-GA nor the SAFE-GPT molecules produced a single candidate —
no GA or SAFE molecule cleared even hurdle 1.

The 43 designed candidates score −9.24 to −8.44 kcal/mol (43 of 43 better than
the best reference), with median QED 0.67, median SA 1.91, and median Tanimoto
to the references of 0.20 — chemically distinct from the literature ligands.

Full funnel, per-generator breakdown and the ranked table are in
[outputs/RESULTS.md](outputs/RESULTS.md).

---

## Reproducing

```bash
conda create -n tem1 -c conda-forge python=3.11 rdkit numpy scipy pandas \
    matplotlib openbabel -y
conda activate tem1
pip install -r requirements.txt

bash scripts/fetch_structures.sh      # 1PZO, 1PZP, 1BTL from RCSB
bash scripts/download_crem_db.sh      # CReM fragment DB (~846 MB -> 5.8 GB)

python scripts/identify_pocket.py     # which ligand copy is the real site
python run_validate.py                # docking protocol validation  -> VALIDATION.md
python run_generate.py --n 3000       # three generators + filters
python run_screen.py  --workers 22 --subset 1100 --sites allosteric
python run_screen.py  --workers 22 --sites orthosteric,surface --top-from-allosteric 150
# unbiased control: same two sites on a RANDOM sample, independent of score
python run_screen.py  --workers 22 --sites orthosteric,surface \
                      --random-from-allosteric 150 --seed 0 --tag _random
python run_analyze.py                 # ranking, figures -> RESULTS.md
pytest -q
```

Receptor preparation is the only step needing the Open Babel CLI, and its
output is cached, so a compute node without Open Babel can run everything else.

### Layout

```
tem1_design/
  config.py      pocket definition, validated protocol, filter thresholds
  structure.py   PDB parsing, pocket identification, receptor prep
  validate.py    redocking / cross-docking benchmark
  filters.py     properties, PAINS+Brenk alerts, novelty
  generate/      crem_gen.py · selfies_ga.py · safe_gen.py
  dock.py        Vina wrapper, 3D embedding, pose RMSD
  screen.py      resumable parallel screening
  plif.py        ProLIF fingerprints, binding-mode fidelity
  score.py       four-hurdle triage and consensus ranking
  plots.py       figures
```

Docking was run on a 24-core node (`solaris/node5`); everything else runs on a
laptop.

---

## What this is not

- **No compound here is experimentally validated.** These are docking-ranked
  proposals. β-lactamase inhibition is straightforward to assay, and until
  someone does, "hit" is the wrong word.
- The receptor is **rigid**, while the whole point of this site is that it is
  conformationally responsive. Rigid docking cannot model the induced fit that
  makes core disruption work.
- Selectivity and site specificity are **differences between docking scores**,
  not measured selectivity. They are reported paired and against an unbiased
  random control, but a Vina score difference is still not an affinity.
- The three production boxes are **not identical volumes**; this was tested
  with volume-matched cubes and runs in the conservative direction, but the
  surface control remains an easy one (least enclosed of the three).
- **Mode fidelity is weaker than it sounds.** The PLIF receptor drops the
  crystallographic waters and has no hydrogens, so the FTA reference row sets
  only 4 hydrophobic bits of 24 and 78% of the top-800 clear the threshold.
- **Provenance.** The committed `docked_allosteric.jsonl` has 1,225 entries
  where the documented `--subset 1100` yields 1,101; the extra molecules and
  116 of the orthosteric rows are carry-over from an earlier run reused by the
  resume logic. They do not enter the paired tests, which use only molecules
  docked at all three sites.
- Novelty is measured against the two reference ligands only — no patent or
  full-ChEMBL search was run.
- The screen docks a MaxMin diversity-selected subset of the filtered library,
  not all of it; docking is the bottleneck.

## License

MIT — see [LICENSE](LICENSE). Cite via [CITATION.cff](CITATION.cff).
