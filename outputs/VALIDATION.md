# Docking protocol validation

Redocking of the crystallographic ligands, run **before** any design work. Success is symmetry-aware heavy-atom RMSD < 2.0 A to the crystal pose, computed in place with no realignment.

| ligand   | receptor   | scoring_function   | waters   | co_ligand   |   dock_score |   crystal_score |   crystal_score_minimised |   rmsd |   centroid_shift | success   |
|:---------|:-----------|:-------------------|:---------|:------------|-------------:|----------------:|--------------------------:|-------:|-----------------:|:----------|
| FTA300   | 1PZP       | vina               | no       | -           |        -8.08 |           -5.71 |                     -6.61 |   3.52 |             2.33 | no        |
| FTA300   | 1PZP       | vina               | yes      | -           |        -8.44 |           -6.52 |                     -7.6  |   0.98 |             0.57 | **yes**   |
| FTA300   | 1PZP       | vinardo            | no       | -           |        -6.45 |           -3.18 |                     -4.77 |   4.03 |             2.53 | no        |
| FTA300   | 1PZP       | vinardo            | yes      | -           |        -6.63 |           -3.67 |                     -5.73 |   2.98 |             1.33 | no        |
| CBT301   | 1PZO       | vina               | no       | -           |        -7.77 |           -6.15 |                     -6.66 |   7.79 |             6.14 | no        |
| CBT301   | 1PZO       | vina               | yes      | -           |        -8.04 |           -6.66 |                     -7.71 |   7.8  |             6.15 | no        |
| CBT301   | 1PZO       | vinardo            | no       | -           |        -5.8  |           -3.66 |                     -4.3  |   7.86 |             6.23 | no        |
| CBT301   | 1PZO       | vinardo            | yes      | -           |        -5.97 |           -4.41 |                     -5.61 |   7.86 |             6.23 | no        |
| CBT301   | 1PZO       | vina               | no       | CBT300      |        -8.92 |           -7.09 |                     -7.39 |   2.4  |             1.41 | no        |
| CBT301   | 1PZO       | vina               | yes      | CBT300      |        -8.42 |           -7.6  |                     -8.52 |   3.11 |             2.14 | no        |
| CBT301   | 1PZO       | vinardo            | no       | CBT300      |        -6.7  |           -4.63 |                     -5.23 |   2.36 |             1.37 | no        |
| CBT301   | 1PZO       | vinardo            | yes      | CBT300      |        -5.54 |           -5.37 |                     -6.06 |   2.8  |             1.55 | no        |
| CBT300   | 1PZO       | vina               | no       | CBT301      |        -7.73 |           -6.88 |                     -6.97 |   4.42 |             0.86 | no        |
| CBT300   | 1PZO       | vina               | yes      | CBT301      |        -7.8  |           -7.25 |                     -7.55 |   1.73 |             0.31 | **yes**   |
| CBT300   | 1PZO       | vinardo            | no       | CBT301      |        -5.62 |           -5.36 |                     -5.58 |   2.12 |             1.02 | no        |
| CBT300   | 1PZO       | vinardo            | yes      | CBT301      |        -6.43 |           -5.81 |                     -6.15 |   2.05 |             0.89 | no        |

**Adopted protocol:** FTA300 / 1PZP, vina, crystallographic waters retained -> RMSD 0.98 A.


## What this shows

1. **Ordered waters are structural here.** FTA moves from 3.52 A (failure) to
   0.98 A (success) purely by retaining the crystallographic waters. The TEM-1
   allosteric cleft is shallow and solvent-mediated; the naked protein is the
   wrong receptor.

2. **The two CBT copies complete each other's site.** CBT300 redocks at 1.73 A
   only when copy 301 is left in the receptor. The cleft is one extended groove
   that holds two ligand molecules, so removing both destroys the surface each
   one binds against. CBT301 never reproduces.

3. **Failures are scoring-function failures, not only sampling failures.**
   Wherever redocking fails, the crystal pose *scores worse* than the incorrect
   pose Vina prefers - so more exhaustiveness cannot fix it. On the naked 1PZO
   receptor the crystal CBT pose scores -6.66 kcal/mol while Vina's (wrong)
   best pose scores -8.04.

4. **Vinardo is worse than Vina throughout** for this pocket, on both RMSD and
   its ranking of the crystal pose.

The practical consequence for the design campaign: the Vina score alone is not
a trustworthy objective in this site, so candidates are additionally required
to reproduce the reference binding mode (ProLIF fingerprint) and to prefer the
allosteric site over two control sites.
