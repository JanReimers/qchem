# TeF₆: is ε a PP tool? (lane B, 2026-11-05)

**Answer: no, not on the XC mesh.**  With the PP terms' own quadrature taken out, all-electron needs 1.9-2.7× the mesh points of the PP run to meet the same ε, and nothing about that cost is ridiculous.  The term that fails on both arms is the Coulomb fit, not the mesh: stock A1_coul misses ε·N by 142× all-electron and 18× on the PP density.  ε is not consumed on molecules at all, so every number here says which settings WOULD meet ε, not that the machinery gets there.

Run: TeF₆, Te-F 1.82 Å (stretch 1.87), LDA, no symmetry.  ε = 1e-5, so |ΔE| ≤ ε·N = 7e-5 Ha, max|ΔH_ij| ≤ 1e-5.
AE: DZVP, 129 Cartesian functions.  PP: GTH-LDA Te q6 / F q7, DZVP-MOLOPT-SR-GTH from `~/Code/cp2k/data/BASIS_MOLOPT` (PBE-optimised), 98 Cartesian functions.

## 1. XC only: PP against all-electron

Fixed density on OUR mesh points, PySCF `numint`, reference 400_77_gl.  The PP matrices do not enter: this is the "hold PP at the converged mesh, vary only XC" arm, to first order in the density, so no patch was needed.  Checks: at 150_47_gl the audit reproduces our `Eloc`, `EenNL`, `Exc` to 1e-8, and its summed dE reproduces the SCF ladder (default −6.105e-4 against −6.11e-4).  The AE reference is converged (300_65_gl against it: dExc 2.3e-9, max|dH| 3.1e-7).

| mesh | pts | ΔE AE | ΔE PP | max\|ΔH\| AE | max\|ΔH\| PP | stretch AE | stretch PP |
|---|---|---|---|---|---|---|---|
| 50_29 | 105,436 | −4.24e-3 | −5.54e-4 | 6.5e-4 | 1.72e-4 | −2.30e-4 | −2.99e-4 |
| **80_29** (default) | 168,688 | −4.48e-5 | −2.19e-5 | 1.12e-5 | 7.36e-6 | +1.36e-4 | −1.97e-5 |
| 100_29 | 210,860 | +7.34e-5 | −3.76e-6 | 1.13e-5 | 7.54e-6 | +1.75e-5 | +7.5e-6 |
| 150_29 | 316,284 | +3.26e-5 | −3.65e-6 | 1.12e-5 | 7.53e-6 | −8.0e-6 | +4.2e-6 |
| 80_35 | 240,592 | −6.4e-5 | −2.14e-5 | 1.73e-5 | 5.19e-6 | +1.60e-4 | −2.29e-5 |
| 150_35 | 451,092 | +2.2e-6 | −2.1e-6 | 4.3e-6 | 1.27e-6 | −1.4e-6 | −1.3e-7 |
| 150_47_gl | 1,194,096 | +1.0e-6 | −2.1e-7 | 4.1e-6 | 1.1e-7 | +2.3e-7 | −1.3e-7 |
| 200_53_gl | 2,017,052 | −5e-8 | −9e-8 | 1.6e-6 | 2.0e-8 | +5e-8 | +8e-9 |

- L = 29 has an angular floor on both arms: 1.12e-5 AE, 7.5e-6 PP.  The AE floor sits just above ε and the PP floor just below it, and that alone sets the max|ΔH| row below.  At ε = 5e-6 both would need L = 35.
- Worst elements.  AE at the default: Te 4s×4s (D 3.4), mid-core; past nR 80 the L floor sits on F polarization d with D ~ 3e-4.  The absolute core diagonals are not the hard case.  PP: Te d (D 0.05) and F d (D 1.4e-3).
- All-electron does not hide its error in the core: free atoms are converged by nR 80 (e_Te 9e-7), so no core error exists to cancel.  The error depends on R, so it comes from the molecular Becke region.  That the cause is neighbour-nucleus sampling on Te's radial mesh is a hypothesis, unmeasured.

## 2. The PP terms' own quadrature error (the structure gap, not ε's)

| mesh | ΔE_loc | ΔE_NL | max\|ΔV_loc\| | max\|ΔV_NL\| |
|---|---|---|---|---|
| 80_29 | −5.65e-4 | −2.3e-5 | 4.65e-4 | 5.5e-6 |
| 150_29 | −2.04e-4 | +6.7e-6 | 4.75e-4 | 1.8e-6 |
| 150_35 | −8.5e-5 | +4.7e-6 | 7.9e-5 | 4.5e-7 |
| 150_47_gl | −7.8e-6 | +2.5e-6 | 1.4e-6 | 5.1e-7 |
| 200_53_gl | −1.9e-6 | +2.1e-7 | 5.8e-7 | 1.4e-7 |

V_loc is the term: 26× the XC error at the default, 8.7× over ε·N in total.  With it counted, the PP run would look HARDER than all-electron, which is the false reading the split removes.

## 3. Cost per target

| target | all-electron | PP (XC only) |
|---|---|---|
| \|ΔE\| ≤ 7e-5 | 150_29: 316,284 pts, 1:59 (80_29 passes by luck, 100_29 fails) | 80_29: 168,688 pts, 1:09 |
| stretch difference ≤ 7e-5 | 100_29: 210,860 pts, 1:54 | 80_29 |
| occupied valence block ≤ 1e-5 | 80_29: 168,688 pts, 1:04 | 80_29 |
| AO max\|ΔH\| ≤ 1e-5 | 150_35: 451,092 pts, 3:12 | 80_29: 168,688 pts |
| full MO max ≤ 1e-5 | 150_47_gl: 1,194,096 pts, 8:28 | not measured (no virtuals dumped) |
| Coulomb fit \|ΔE\| ≤ 7e-5 | ~2400 functions (l≤3, ratio 2); 2835 comfortably (l≤4, ratio 3) | 680 (l≤3, ratio 3) |
| Coulomb fit max\|ΔJ\| ≤ 1e-5 | ~4200 (l≤4, ratio 2) | 1855 (l≤4, ratio 2) |

Walls are ours, 2 SCFs, 4 threads, on a loaded box.  Fit sizes are measured in PySCF only; our fit cost at those sizes is unmeasured.

## 4. The Coulomb fit: the term that fails first, on both arms

Exact Coulomb is the reference.  Density: exact-J for AE; for PP our A1-fit SCF density, so its numbers are fit error at that density.

| | AE | PP |
|---|---|---|
| stock A1 ΔE | −9.92e-3 (142×) | −1.24e-3 (18×) |
| stock A1 stretch difference | 4.3e-5, passes | −9.4e-5, fails 1.3× |
| stock A1 AO max\|ΔJ\| | 2.8e-3 | 7.9e-3 |
| stock A1 occupied block | 6.7e-4 | 2.8e-3 |

- Ours minus exact: AE −9.94e-3 (exact-J −7194.67283302, ours −7194.6828193).  PP −1.216e-3 (exact-J −153.4071868, ours −153.4084030).  Both match the fixed-density residual, so the fit explains the gap.
- A1 is an all-electron-shaped set, but it is no better on valence density in the H matrix: max|ΔJ| is 2.8× WORSE on PP.
- Even-tempered sweep: l_max matters more than ratio.  On AE l≤2 never passes, with a floor of 1.6e-4; plain Cholesky stops near ratio 2 at l≤4 (cond ~1e13) and an eigen-cut still works.
- Bases with l_max 2 need fit l ≥ 3-4, so U (orbital l_max 4) needs well above A1's l ≤ 2.

## 5. (b) against (a) on the mesh, and on the seed

Sent earlier; unchanged.  (b) does not track (a) in value: at nR 80 the sign is opposite and |a| is 2.4-3× |b|.  At 7e-5, converged (b) disagrees with (a) on 5 of 16 rungs, seed (b) on 1, and both pass the default where (a) fails 2×.  With a 2× margin both first pass at 150_29 and pass nothing (a) fails; that margin was chosen after seeing the data, on one molecule.  The seed loses little: 1-2e-5 on the L ≤ 29 rungs.  The seed here is PySCF's spherical neutral atoms, not our SAD.  (b) does NOT size the fit: A1 gives (a) 4.3e-5 and (b) 2.24e-3.  The free-atom UKS densities were reported unconverged.  No PP free-atom runs, so (b) has no PP number; PP (a) is in §1.

## 6. TeF₆ in a box (solid GPW path, analytic PP)

20 bohr cube, Γ, spherical, ε 1e-5 through the pilot.  It is cheap: 23:41 wall, 2.0 GB, converged, `[basis trim]` clean (min eig S 1.48e-2, cond 268).

- **Rasters: the machinery works.**  The pilot settled after 2 re-convergences, density ball cutoffFactor 8.  Audit: V_long −5.2e-5, Hartree −2.8e-6, max|ΔH| 3.07e-6.  Every raster term meets ε.
- **The XC mesh is not sized by the pilot** (`XC L0` = 0).  Against (200,35): default (40,29) ΔE −1.13e-4 (1.6×), max|ΔH| 1.23e-4 (12×); stated (80,29) ΔE −6.7e-5 (passes), max|ΔH| 6.7e-5 (6.7×), worst on Te d.
- **Unresolved:** at nominally (80,29) the box's XC max|ΔH| is 9× the molecule's 7.4e-6.  They are not the same mesh or basis: the box drops 24% of its points as tails (beckeEps 1e-6, 127,980 against 168,688) and is spherical.  The cause is unmeasured.
- **Agreement with the molecule:** box (80,29) −153.4040315.  Exact-J spherical molecule −153.4038870 (PP at the converged mesh).  Raw gap 1.4e-4; removing the audited XC and V_long errors leaves 2.5e-5, inside ε·N, IF the audit's ΔE is run minus reference (sign not verified).  Cartesian against spherical differs by 3.3e-3 (−153.4071868), so compare only like with like.
- **Two things needed to run it at all:** `gth-lda-valence` has no Te, so a deck cannot name it (I added `"Te": 6` to the data file for the launch and restored it; nothing committed); `atomic_valence_densities.json` has no Te q6, so `IonicSAD` throws and the run used `seed: Uniform`, which converged without trouble.

## 7. Other findings

- **ε on the molecule facade:** none consumed.  Fixed mesh (MHL 80, Lebedev 29), `screenEps` 1e-10, no H audit, no term audit; `AccuracyTarget accuracy` is only carried.  All H audits here are external (PySCF on our mesh).
- **Seed:** AE and PP ran from the core guess, not the default.  Default SAD stops at Ne (row filed); it did not make the SCF hard: AE 18 + 20 iterations, no cap hit.
- **One cap hit:** PP at (20,11), 100 iterations, not converged, E −150.34.  That mesh is far below any default; every rung from (30,17) up converged.
- **Spherical molecular PP:** NaN in H from nR ≥ 50 (row filed), which is why the PP ladder is Cartesian.
- **jorge-DZP TeF₆:** reader takes a normalised copy only (BSE's raw header does not parse).  (80,29) −7198.2872610510 and (200,35) −7198.2872664890, 5.4e-6 apart; 35 iterations; min eig S 5.5e-5, cond 1.1e5, 0 drop; 2:54 and 9:58, 2.6 GB.  No exact-J number on jorge yet.
- **UF₆ blocker (from reading the code, not run):** A1_coul has no U, and `Gaussian94Reader::FindAtom` throws `element U (Z=92) has no entry…` in the Hamiltonian build.  Ham_DFT and Ham_PP always use `FittedVee`; only HF has exact `Vee`, so a DFT exact-Coulomb route would be new code.

Scripts and logs: this directory (`pp_audit.py`, `pp_fit.py`, `pp_exact*.py`, `xc_audit.py`, `diff_audit.py`, `seed_audit.py`, `fit_arms.py`, `fit_diff.py`, `gpw/`).  `pp_audit182_wrongTe.log` is void in its V_loc and KB columns (O's entry matched the q6 tag); its XC column was right.
