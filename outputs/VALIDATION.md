# Docking protocol validation

Redocking of the crystallographic ligands, run **before** any design work. Success is symmetry-aware heavy-atom RMSD < 2.0 A to the crystal pose, computed in place with no realignment.

| ligand   | receptor   | scoring_function   | waters   | co_ligand   |   dock_score |   crystal_score |   crystal_score_minimised |   rmsd |   centroid_shift | success   |
|:---------|:-----------|:-------------------|:---------|:------------|-------------:|----------------:|--------------------------:|-------:|-----------------:|:----------|
| FTA300   | 1PZP       | vina               | False    | -           |        -8.08 |           -5.71 |                     -6.61 |   3.52 |             2.33 | no        |
| FTA300   | 1PZP       | vina               | True     | -           |        -8.44 |           -6.52 |                     -7.6  |   0.98 |             0.57 | **yes**   |
| FTA300   | 1PZP       | vinardo            | False    | -           |        -6.45 |           -3.18 |                     -4.77 |   4.03 |             2.53 | no        |
| FTA300   | 1PZP       | vinardo            | True     | -           |        -6.63 |           -3.67 |                     -5.73 |   2.98 |             1.33 | no        |
| CBT301   | 1PZO       | vina               | False    | -           |        -7.77 |           -6.15 |                     -6.66 |   7.79 |             6.14 | no        |
| CBT301   | 1PZO       | vina               | True     | -           |        -8.04 |           -6.66 |                     -7.71 |   7.8  |             6.15 | no        |
| CBT301   | 1PZO       | vinardo            | False    | -           |        -5.8  |           -3.66 |                     -4.3  |   7.86 |             6.23 | no        |
| CBT301   | 1PZO       | vinardo            | True     | -           |        -5.97 |           -4.41 |                     -5.61 |   7.86 |             6.23 | no        |
| CBT301   | 1PZO       | vina               | False    | CBT300      |        -8.92 |           -7.09 |                     -7.39 |   2.4  |             1.41 | no        |
| CBT301   | 1PZO       | vina               | True     | CBT300      |        -8.42 |           -7.6  |                     -8.52 |   3.11 |             2.14 | no        |
| CBT301   | 1PZO       | vinardo            | False    | CBT300      |        -6.7  |           -4.63 |                     -5.23 |   2.36 |             1.37 | no        |
| CBT301   | 1PZO       | vinardo            | True     | CBT300      |        -5.54 |           -5.37 |                     -6.06 |   2.8  |             1.55 | no        |
| CBT300   | 1PZO       | vina               | False    | CBT301      |        -7.73 |           -6.88 |                     -6.97 |   4.42 |             0.86 | no        |
| CBT300   | 1PZO       | vina               | True     | CBT301      |        -7.8  |           -7.25 |                     -7.55 |   1.73 |             0.31 | **yes**   |
| CBT300   | 1PZO       | vinardo            | False    | CBT301      |        -5.62 |           -5.36 |                     -5.58 |   2.12 |             1.02 | no        |
| CBT300   | 1PZO       | vinardo            | True     | CBT301      |        -6.43 |           -5.81 |                     -6.15 |   2.05 |             0.89 | no        |

## What this establishes

The matrix has **16 conditions, 2 successes**. The two are not equivalent evidence:

- **At the allosteric copies** (the site being designed against): **1 of 12** conditions succeed.
- At the second, non-allosteric crystallographic copy: 1 of 4.


**Adopted protocol:** FTA300 / 1PZP, vina, waters=True -> RMSD 0.98 A. This is the only condition that reproduces a pose at the design site, so the design run uses it and nothing else.


One success out of 12 is a thin basis, and it should be read that way: the protocol is usable at this pocket for this ligand, not shown to be generally reliable here.


### Ordered waters change the answer

Retaining the crystallographic waters improves RMSD in only **2 of 6** matched comparisons at the allosteric copies (median change -0.00 A), so this is **not** a general effect in this matrix. It is decisive for the one case that matters: FTA300 goes from 3.52 A (failure) to 0.98 A purely by keeping the waters. That single result is why the design receptor keeps them — a specific finding, not a rule.


### Vina vs Vinardo

Vina gives the lower RMSD in **5 of 8** matched conditions — better on balance, but *not* uniformly: vinardo wins the remaining 3. Vina is adopted because it produces the only sub-2 A pose at the design site, not because it dominates everywhere.


### Is failure in the scoring or the sampling?

In **14 of 14** failed conditions the crystal pose scores worse than the pose the docking preferred, which points at the scoring function rather than at insufficient search. That comparison is unfair to the crystal pose, though, since it was never relaxed in the scoring field. Against the **minimised** crystal pose it holds in 12 of 14 — so in 2 case(s) the correct pose does score better and the failure really is one of sampling.
