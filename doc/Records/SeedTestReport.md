# Seed test — size the XC mesh from the seed, before iteration 1 (lane B, 2026-11-07)

Rows: "Size every mesh from the seed density, before iteration 1" and "Two-step runs" (doc/OpenWork.md).
Data, scripts and logs: `~/Code/qchem6-runs/SeedSizing/` (`seed_mol.py`, `analyse_mol.py`, `solid_seed.py`, `analyse_solid.py`,
`systems.json`; one directory per system and step).  Solids ran on the radial batch's snapshot `snap/laneb-rr-20261106-085824`
(its `RR_SCAN` patch: one frozen D, a list of candidate meshes, each rung's E_xc and V_xc matrix dumped).

## The question, and how it is measured

A **ladder** of Becke meshes (nR × angular rung) is audited at a FIXED density D against a fine reference at the same D
((400, 77 GL) on TeF₆, (300, 77 GL) on solids).  A rung **passes** at ε when |ΔE_xc| ≤ ε·N_atoms **and** max|ΔH_ij| ≤ ε
(H = the V_xc matrix).  The **pick** is persistent: at each angular rung, the smallest nR that passes with every larger nR on the
ladder passing too; the pick is the cheapest of those over angular rungs.  The **truth** is the pick made at the converged D (Df).

Three densities per system:
- **seed** — the density before iteration 1.  TeF₆: **seed = PySCF spherical atoms** (`get_init_guess(key='atom')`; our SAD stops at
  Ne, so it cannot seed Te).  Solids: the D of the FIRST FILL — one Fock build at our SAD (MnO: IonicSAD, and SAD beside it),
  diagonalized.  LiF at that D: Etot −23.40 vs −24.32 converged, Δρ 0.90.
- **D4** — converged at ε = 1e-4 on the rung the seed picks at 1e-4 (the two-step run's step 1).
- **Df** — converged on a fine mesh, restarted from D5 (D5 = step 2: restart D4 on the rung D4 picks at 1e-5).

Rules: today's defaults throughout (the radial rule is not settled): molecule MHL m3 α2 Becke order 2; solid MHL m2 α2 order 3.
Poisson energies: none used (the molecular ρ client is PySCF's own integrals plus the s-fit V; the solids audit XC only), so no
pre-fix Poisson number enters.

**Clients.**  XC (E and H, both codes).  ρ (the Hartree client) on TeF₆ only: lane S's δ Coulomb route is molecular and REFUSED on a
cell, so the solids audit XC only.  On TeF₆ the ρ client is measured as the first-order Hartree error dJ1 = Σ w ρ V − (same on the
reference), V from an s-only Coulomb-metric fit of ρ at that D (its error enters at second order; it agrees with the full A1 fit to a
median 1e-4 relative), and judged against ε·N.  Beside it the cheap prediction **2 E_H δq / N** (δq = Σwρ − N).

## Verdicts

### (a) Does the seed pick the right rung at ε = 1e-4?

| system | client | truth (Df) | seed pick | gap (nR, L) | points | at Df: \|ΔE\|/εN, max\|ΔH\|/ε |
|---|---|---|---|---|---|---|
| TeF₆ PP (seed = PySCF spherical atoms) | XC | (70,23) 93,992 | (70,23) | 0, 0 | ×1.00 | 0.12, 0.48 |
| TeF₆ PP (seed = PySCF spherical atoms) | ρ | (70,29) 147,608 | (80,35) | +1, +1 | ×1.63 | 0.74 |
| TeF₆ AE (seed = PySCF spherical atoms) | XC | (70,29) 147,608 | (70,29) | 0, 0 | ×1.00 | 0.87, 0.47 |
| TeF₆ AE (seed = PySCF spherical atoms) | ρ | (150,29) 316,284 | (150,35) | 0, +1 | ×1.43 | 0.32 |
| LiF (seed = first fill at SAD) | XC | (60,GL41) 67,996 | (60,GL41) | 0, 0 | ×1.00 | 0.01, 0.26 |
| NaF (first fill at SAD) | XC | (40,35) 22,936 | (60,35) | +1, 0 | ×1.50 | 0.03, 0.17 |
| MnO (first fill at IonicSAD) | XC | (60,GL41) 138,384 | (60,GL41) | 0, 0 | ×1.00 | 0.03, 0.36 |
| MnO (first fill at SAD) | XC | (60,GL41) 138,384 | (60,GL41) | 0, 0 | ×1.00 | 0.03, 0.36 |

### (b1–b2) Does the rung picked at the 1e-4 D deliver ε = 1e-5 at the final D, and how far off is it?

| system | client | truth at 1e-5 (Df) | D4 pick | gap | points | at Df: \|ΔE\|/εN, max\|ΔH\|/ε |
|---|---|---|---|---|---|---|
| TeF₆ PP (seed = PySCF spherical atoms) | XC | (80,35) 240,832 | (80,35) | 0, 0 | ×1.00 | 0.31, 0.57 |
| TeF₆ PP (seed = PySCF spherical atoms) | ρ | (100,29) 210,860 | (100,29) | 0, 0 | ×1.00 | 0.41 |
| TeF₆ AE (seed = PySCF spherical atoms) | XC | (100,35) 301,212 | (100,35) | 0, 0 | ×1.00 | 0.61, 0.48 |
| TeF₆ AE (seed = PySCF spherical atoms) | ρ | (150,47) 803,588 | (150,47) | 0, 0 | ×1.00 | 0.02 |
| LiF | XC | (100,GL47) 149,276 | (100,GL47) | 0, 0 | ×1.00 | 0.02, 0.58 |
| NaF | XC | (60,GL41) 69,728 | (60,GL41) | 0, 0 | ×1.00 | 0.05, 0.22 |
| MnO | XC | (80,L65) 307,288 | (80,L65) | 0, 0 | ×1.00 | 0.05, 0.27 |

### (b3) The same from the seed (the comparison: the 1e-4 D against the seed)

| system | client | truth at 1e-5 | seed pick | gap | points | at Df |
|---|---|---|---|---|---|---|
| TeF₆ PP (seed = PySCF spherical atoms) | XC | (80,35) | (80,35) | 0, 0 | ×1.00 | 0.31, 0.57 |
| TeF₆ PP (seed = PySCF spherical atoms) | ρ | (100,29) | (100,41) | 0, +2 | ×1.94 | 0.16 |
| TeF₆ AE (seed = PySCF spherical atoms) | XC | (100,35) | (100,35) | 0, 0 | ×1.00 | 0.61, 0.48 |
| TeF₆ AE (seed = PySCF spherical atoms) | ρ | (150,47) | (150,47) | 0, 0 | ×1.00 | 0.02 |
| LiF (first fill at SAD) | XC | (100,GL47) | **none** | — | — | — |
| NaF (first fill at SAD) | XC | (60,GL41) | (60,GL41) | 0, 0 | ×1.00 | 0.05, 0.22 |
| MnO (first fill at IonicSAD) | XC | (80,L65) | **(80,GL47)** | 0, **−2** | ×0.79 | 0.01, **2.33 — FAILS** |
| MnO (first fill at SAD) | XC | (80,L65) | **none** | — | — | — |

LiF's seed picks NOTHING at 1e-5: at the first-fill D no rung's H gets below ε persistently — the floor at the finest rungs is
1.5–5 ε (L65 nR 100–150: 1.53–1.59 ε; GL47: 3.3 ε), against 0.56–0.61 ε at Df.  (L65 nR 80 passes at 0.37 ε, but nR 100 then fails,
so the persistent rule refuses it.)  The cause is not measured; what is measured is that the seed's H error at the finest rungs is
3–6× the converged D's, so a 1e-5 audit at the seed asks for a mesh finer than the ladder holds.

## Per-rung tables: ΔE against ε·N and max|ΔH| against ε (controller ask (1))

Units: dE/(1e-5·N), max|dH|/1e-5 — a cell ≤ 1 meets ε = 1e-5, ≤ 10 meets 1e-4.  **(40,29) Lebedev is today's solid default.**

### LiF (N = 2)

| rung | points | seed: ΔE, ΔH | D4: ΔE, ΔH | Df: ΔE, ΔH |
|---|---|---|---|---|
| (30,L23) | 7472 | -218.1, 1566.2 | 17.0, 341.2 | 17.0, 341.2 |
| (40,L23) | 10008 | -342.4, 1471.0 | -10.5, 494.0 | -10.5, 494.0 |
| (60,L23) | 15180 | -281.7, 1299.7 | -5.9, 385.8 | -5.9, 385.8 |
| (80,L23) | 20272 | -287.1, 1322.6 | -5.9, 396.0 | -5.9, 396.0 |
| (100,L23) | 25350 | -291.0, 1338.3 | -6.0, 402.3 | -6.0, 402.3 |
| (120,L23) | 30490 | -291.1, 1338.5 | -6.0, 402.3 | -6.0, 402.3 |
| (150,L23) | 38150 | -290.9, 1338.0 | -6.0, 402.1 | -6.0, 402.1 |
| (30,L29) | 11680 | -128.8, 1049.7 | 18.8, 223.5 | 18.8, 223.5 |
| **(40,L29) default** | 15560 | -115.5, 459.1 | -5.6, 177.4 | -5.6, 177.4 |
| (60,L29) | 23564 | -91.6, 417.5 | -2.0, 120.4 | -2.0, 120.4 |
| (80,L29) | 31424 | -96.1, 435.6 | -2.0, 129.6 | -2.0, 129.6 |
| (100,L29) | 39294 | -98.8, 446.5 | -2.1, 134.2 | -2.1, 134.2 |
| (120,L29) | 47394 | -98.8, 446.8 | -2.1, 134.2 | -2.1, 134.2 |
| (150,L29) | 59318 | -98.8, 446.5 | -2.1, 134.1 | -2.1, 134.1 |
| (30,L35) | 16664 | -101.1, 898.1 | 19.5, 192.7 | 19.5, 192.7 |
| (40,L35) | 22288 | -45.1, 148.4 | -4.1, 82.2 | -4.1, 82.2 |
| (60,L35) | 33604 | -31.9, 144.5 | -0.8, 40.1 | -0.8, 40.1 |
| (80,L35) | 45192 | -33.8, 150.9 | -0.7, 45.1 | -0.7, 45.1 |
| (100,L35) | 56478 | -35.5, 157.7 | -0.7, 48.2 | -0.7, 48.2 |
| (120,L35) | 67906 | -35.5, 157.8 | -0.7, 48.2 | -0.7, 48.2 |
| (150,L35) | 85078 | -35.4, 157.6 | -0.7, 48.1 | -0.7, 48.1 |
| (30,GL41) | 33552 | -87.2, 830.8 | 19.9, 180.7 | 19.9, 180.7 |
| (40,GL41) | 45052 | -7.8, 18.8 | -3.3, 29.5 | -3.3, 29.5 |
| (60,GL41) | 67996 | -0.2, 1.5 | -0.1, 2.6 | -0.1, 2.6 |
| (80,GL41) | 91116 | 1.7, 6.7 | 0.0, 3.0 | 0.0, 3.0 |
| (100,GL41) | 114020 | 1.2, 4.9 | 0.0, 1.8 | 0.0, 1.8 |
| (120,GL41) | 136996 | 1.2, 4.8 | 0.0, 1.8 | 0.0, 1.8 |
| (150,GL41) | 171536 | 1.2, 4.8 | 0.0, 1.8 | 0.0, 1.8 |
| (30,GL47) | 43956 | -85.5, 821.8 | 19.9, 177.9 | 19.9, 177.9 |
| (40,GL47) | 58996 | -10.9, 9.1 | -3.4, 36.2 | -3.4, 36.2 |
| (60,GL47) | 89116 | -1.5, 7.7 | -0.1, 0.4 | -0.1, 0.4 |
| (80,GL47) | 119120 | 1.2, 4.3 | 0.0, 1.5 | 0.0, 1.5 |
| (100,GL47) | 149276 | 0.9, 3.3 | 0.0, 0.6 | 0.0, 0.6 |
| (120,GL47) | 179272 | 0.9, 3.3 | 0.0, 0.6 | 0.0, 0.6 |
| (150,GL47) | 224472 | 0.9, 3.3 | 0.0, 0.6 | 0.0, 0.6 |
| (30,L53) | 37320 | -86.3, 825.2 | 19.9, 177.9 | 19.9, 177.9 |
| (40,L53) | 50136 | -14.0, 9.6 | -3.4, 40.2 | -3.4, 40.2 |
| (60,L53) | 75644 | -4.0, 17.9 | -0.2, 3.3 | -0.2, 3.3 |
| (80,L53) | 101264 | -1.7, 7.4 | -0.0, 1.9 | -0.0, 1.9 |
| (100,L53) | 126966 | -2.1, 9.0 | -0.0, 3.1 | -0.0, 3.1 |
| (120,L53) | 152394 | -2.1, 8.9 | -0.0, 3.0 | -0.0, 3.0 |
| (150,L53) | 190638 | -2.1, 8.8 | -0.0, 3.0 | -0.0, 3.0 |
| (30,L65) | 55712 | -85.2, 820.8 | 19.9, 176.6 | 19.9, 176.6 |
| (40,L65) | 74648 | -13.1, 6.5 | -3.4, 38.9 | -3.4, 38.9 |
| (60,L65) | 112668 | -2.9, 13.0 | -0.1, 1.8 | -0.1, 1.8 |
| (80,L65) | 150888 | -0.1, 0.4 | -0.0, 0.4 | -0.0, 0.4 |
| (100,L65) | 188998 | -0.4, 1.5 | -0.0, 0.6 | -0.0, 0.6 |
| (120,L65) | 227090 | -0.4, 1.6 | -0.0, 0.6 | -0.0, 0.6 |
| (150,L65) | 284046 | -0.4, 1.6 | -0.0, 0.6 | -0.0, 0.6 |

On LiF the default (40,29) misses at every density: at Df ΔE = −5.6 εN (1.1e-4 Ha) and max|ΔH| = 177 ε — the ladder REFUSES it
at the seed already (−115 εN, 459 ε), so a seed-sized run never lands on it.  H, not E, sets the rung: at (40, L35) ΔE is −4.1 εN but
H is 82 ε; the L23–L35 ladders never meet ε on H at any nR (an angular floor: 400 ε, 134 ε, 48 ε at Df).

### NaF (N = 2)

| rung | points | seed: ΔE, ΔH | D4: ΔE, ΔH | Df: ΔE, ΔH |
|---|---|---|---|---|
| (30,L23) | 7708 | -25.7, 34.6 | -15.5, 26.7 | -15.5, 26.7 |
| (40,L23) | 10312 | -52.7, 33.2 | -2.7, 22.6 | -2.7, 22.6 |
| (60,L23) | 15536 | -39.0, 24.3 | -2.0, 15.3 | -2.0, 15.3 |
| (80,L23) | 20772 | -40.8, 25.3 | -1.9, 15.4 | -1.9, 15.4 |
| (100,L23) | 26068 | -40.5, 25.1 | -1.9, 15.4 | -1.9, 15.4 |
| (120,L23) | 31216 | -40.5, 25.2 | -1.9, 15.4 | -1.9, 15.4 |
| (150,L23) | 39048 | -40.5, 25.2 | -1.9, 15.4 | -1.9, 15.4 |
| (30,L29) | 11884 | 6.6, 38.9 | -13.9, 24.5 | -13.9, 24.5 |
| **(40,L29) default** | 15928 | -22.0, 13.7 | -1.6, 11.6 | -1.6, 11.6 |
| (60,L29) | 24056 | -13.3, 8.0 | -0.8, 4.5 | -0.8, 4.5 |
| (80,L29) | 32172 | -14.6, 8.7 | -0.7, 4.6 | -0.7, 4.6 |
| (100,L29) | 40444 | -14.3, 8.5 | -0.7, 4.6 | -0.7, 4.6 |
| (120,L29) | 48352 | -14.3, 8.6 | -0.7, 4.6 | -0.7, 4.6 |
| (150,L29) | 60696 | -14.3, 8.5 | -0.7, 4.6 | -0.7, 4.6 |
| (30,L35) | 16924 | 17.6, 40.1 | -13.4, 27.5 | -13.4, 27.5 |
| (40,L35) | 22936 | -10.8, 8.1 | -1.1, 8.5 | -1.1, 8.5 |
| (60,L35) | 34448 | -4.9, 2.9 | -0.3, 1.7 | -0.3, 1.7 |
| (80,L35) | 46260 | -5.5, 3.2 | -0.3, 1.5 | -0.3, 1.5 |
| (100,L35) | 57844 | -5.2, 3.1 | -0.3, 1.5 | -0.3, 1.5 |
| (120,L35) | 69568 | -5.3, 3.1 | -0.3, 1.5 | -0.3, 1.5 |
| (150,L35) | 86856 | -5.3, 3.1 | -0.3, 1.5 | -0.3, 1.5 |
| (30,GL41) | 34540 | 25.0, 41.0 | -13.1, 29.4 | -13.1, 29.4 |
| (40,GL41) | 46240 | -3.5, 7.2 | -0.8, 7.2 | -0.8, 7.2 |
| (60,GL41) | 69728 | 0.1, 0.4 | -0.1, 0.2 | -0.1, 0.2 |
| (80,GL41) | 93244 | 0.2, 0.4 | 0.0, 0.2 | 0.0, 0.2 |
| (100,GL41) | 116692 | 0.2, 0.4 | 0.0, 0.2 | 0.0, 0.2 |
| (120,GL41) | 140408 | 0.2, 0.4 | 0.0, 0.2 | 0.0, 0.2 |
| (150,GL41) | 175504 | 0.2, 0.4 | 0.0, 0.2 | 0.0, 0.2 |
| (30,GL47) | 45136 | 24.0, 40.8 | -13.1, 29.4 | -13.1, 29.3 |
| (40,GL47) | 60352 | -3.3, 7.1 | -0.8, 7.1 | -0.8, 7.1 |
| (60,GL47) | 91136 | -0.2, 0.2 | -0.1, 0.1 | -0.1, 0.1 |
| (80,GL47) | 121840 | 0.0, 0.2 | 0.0, 0.2 | 0.0, 0.2 |
| (100,GL47) | 152688 | 0.1, 0.2 | 0.0, 0.1 | 0.0, 0.1 |
| (120,GL47) | 183424 | 0.1, 0.2 | 0.0, 0.1 | 0.0, 0.1 |
| (150,GL47) | 229504 | 0.1, 0.2 | 0.0, 0.1 | 0.0, 0.1 |
| (30,L53) | 38332 | 23.5, 40.8 | -13.1, 29.1 | -13.1, 29.1 |
| (40,L53) | 51448 | -4.1, 7.2 | -0.8, 7.1 | -0.8, 7.1 |
| (60,L53) | 77480 | -0.5, 0.3 | -0.1, 0.3 | -0.1, 0.3 |
| (80,L53) | 103644 | -0.4, 0.2 | -0.0, 0.1 | -0.0, 0.1 |
| (100,L53) | 129724 | -0.3, 0.2 | -0.0, 0.1 | -0.0, 0.1 |
| (120,L53) | 156016 | -0.3, 0.2 | -0.0, 0.1 | -0.0, 0.1 |
| (150,L53) | 195000 | -0.3, 0.2 | -0.0, 0.1 | -0.0, 0.1 |
| (30,L65) | 56764 | 23.8, 40.8 | -13.1, 29.1 | -13.1, 29.1 |
| (40,L65) | 76552 | -3.6, 7.2 | -0.8, 7.0 | -0.8, 7.0 |
| (60,L65) | 115304 | -0.3, 0.2 | -0.1, 0.2 | -0.1, 0.2 |
| (80,L65) | 154236 | -0.1, 0.0 | -0.0, 0.0 | -0.0, 0.0 |
| (100,L65) | 193036 | -0.1, 0.0 | -0.0, 0.0 | -0.0, 0.0 |
| (120,L65) | 232336 | -0.1, 0.0 | -0.0, 0.0 | -0.0, 0.0 |
| (150,L65) | 290664 | -0.1, 0.0 | -0.0, 0.0 | -0.0, 0.0 |

On NaF the default (40,29) misses ε = 1e-5 (ΔE −1.6 εN, max|ΔH| 11.6 ε at Df) and sits on the edge of 1e-4 (H 1.16×); the
seed ladder refuses it too (−22 εN, 13.7 ε).

### MnO (N = 4)

| rung | points | seed: ΔE, ΔH | seed_SAD: ΔE, ΔH | D4: ΔE, ΔH | Df: ΔE, ΔH |
|---|---|---|---|---|---|
| **(40,L29) default** | 31664 | -46.9, 131.1 | -68.4, 145.6 | -11.9, 105.4 | -11.9, 105.4 |
| (60,L29) | 47920 | -48.7, 146.6 | -82.5, 167.2 | -18.3, 117.5 | -18.3, 117.5 |
| (80,L29) | 63928 | -43.7, 130.4 | -73.4, 148.5 | -16.4, 104.6 | -16.4, 104.6 |
| (100,L29) | 79992 | -44.6, 133.6 | -75.2, 152.1 | -16.8, 107.1 | -16.8, 107.1 |
| (120,L29) | 96320 | -44.4, 133.2 | -74.9, 151.6 | -16.7, 106.8 | -16.7, 106.8 |
| (40,L35) | 45680 | -19.1, 54.5 | -25.0, 58.7 | -1.4, 44.1 | -1.4, 44.1 |
| (60,L35) | 68944 | -18.3, 56.5 | -31.7, 64.4 | -6.8, 45.3 | -6.8, 45.3 |
| (80,L35) | 91864 | -15.5, 47.8 | -26.8, 54.4 | -5.8, 38.3 | -5.8, 38.3 |
| (100,L35) | 115080 | -16.0, 49.7 | -27.9, 56.7 | -6.0, 39.9 | -6.0, 39.9 |
| (120,L35) | 138224 | -15.9, 49.4 | -27.7, 56.2 | -6.0, 39.6 | -6.0, 39.6 |
| (40,GL41) | 91824 | -4.4, 11.3 | -0.1, 11.2 | 4.1, 11.2 | 4.1, 11.2 |
| (60,GL41) | 138384 | -0.7, 0.7 | -0.5, 4.8 | -0.2, 3.6 | -0.2, 3.6 |
| (80,GL41) | 185400 | 0.6, 2.5 | 1.4, 4.8 | 0.2, 3.6 | 0.2, 3.6 |
| (100,GL41) | 231928 | 0.5, 2.1 | 1.1, 4.8 | 0.2, 3.6 | 0.2, 3.6 |
| (120,GL41) | 278528 | 0.6, 2.2 | 1.2, 4.8 | 0.2, 3.6 | 0.2, 3.6 |
| (40,GL47) | 119808 | -3.5, 11.3 | 1.3, 11.2 | 4.5, 11.2 | 4.5, 11.2 |
| (60,GL47) | 181056 | -1.1, 2.3 | -1.4, 4.2 | -0.4, 2.3 | -0.4, 2.3 |
| (80,GL47) | 242208 | 0.1, 0.8 | 0.2, 4.2 | 0.0, 2.3 | 0.0, 2.3 |
| (100,GL47) | 303488 | -0.0, 0.7 | 0.1, 4.2 | -0.0, 2.3 | -0.0, 2.3 |
| (120,GL47) | 364640 | 0.0, 0.7 | 0.1, 4.2 | -0.0, 2.3 | -0.0, 2.3 |
| (40,L53) | 102320 | -4.2, 11.3 | -0.1, 11.2 | 4.2, 11.2 | 4.2, 11.2 |
| (60,L53) | 154096 | -2.4, 6.7 | -3.8, 7.6 | -0.9, 5.3 | -0.9, 5.3 |
| (80,L53) | 206104 | -0.9, 2.9 | -1.6, 3.3 | -0.3, 2.3 | -0.3, 2.3 |
| (100,L53) | 258168 | -1.0, 3.3 | -1.8, 3.8 | -0.4, 2.6 | -0.4, 2.6 |
| (120,L53) | 309920 | -0.9, 3.2 | -1.8, 3.6 | -0.4, 2.5 | -0.4, 2.5 |
| (40,L65) | 152528 | -3.4, 11.3 | 1.5, 11.2 | 4.5, 11.2 | 4.5, 11.2 |
| (60,L65) | 229552 | -1.5, 4.0 | -2.3, 4.5 | -0.6, 3.2 | -0.6, 3.2 |
| (80,L65) | 307288 | -0.1, 0.3 | -0.2, 1.1 | -0.1, 0.3 | -0.1, 0.3 |
| (100,L65) | 384888 | -0.2, 0.6 | -0.3, 1.1 | -0.1, 0.5 | -0.1, 0.5 |
| (120,L65) | 461984 | -0.2, 0.6 | -0.3, 1.1 | -0.1, 0.5 | -0.1, 0.5 |

On MnO the default (40,29) misses by 105 ε on H at Df (ΔE −11.9 εN), and every density's ladder refuses it (131–146 ε).  The
MnO seeds are far from the answer: the IonicSAD first fill has Etot −48.85 Ha (converged −61.30) and site moment 1.03 e (seed 4.79,
converged 4.46); the SAD first fill −59.73 Ha, 3.93 e.  At the IonicSAD fill the GL47 ladder's H floor is 0.7 ε, below the converged
D's 2.3 ε, so its 1e-5 audit accepts (80,GL47), which delivers 2.33 ε at Df: **the one unsafe pick in the test.**  At the SAD fill the
floor is 4.2 ε (GL47) and 1.1 ε (L65), so nothing passes.

## The ρ client's cheap prediction 2 E_H δq / N — a calibration, on TWO systems (TeF₆ PP and AE)

| density | PP: dJ1/pred median, 10–90% | AE: median, 10–90% |
|---|---|---|
| seed (PySCF spherical atoms) | 0.89, [0.77, 1.11] | 0.62, [0.26, 1.06] |
| D4 | 1.17, [0.63, 3.40] | 0.59, [0.15, 1.23] |
| Df | 1.17, [0.63, 3.40] | 0.59, [0.15, 1.23] |

Picking on the prediction instead of dJ1: PP at Df picks (60,23) at 1e-4 (measured 0.29 εN, passes) and (100,35) at 1e-5 (0.74, passes);
AE at Df picks (150,35) at 1e-4 (0.32) and (150,47) at 1e-5 (0.02).  So on these two the prediction picks a passing rung, but its ratio
to the measurement spans 0.15–3.4 per rung: a factor of 3 under-prediction sits inside the 10–90% band on PP.  Any margin built on it
is a named calibration resting on **2 systems** (one element pair, two Hamiltonians) — not enough to ship without NaF/MnO-class
evidence, which the solids cannot give until the δ route exists on a cell.

## Step-2 restart cost (the two-step run)

| system | step 1 (seed → D4 on its 1e-4 pick) | step 2 (D4 → D5 on its 1e-5 pick) | Df (D5 → fine mesh) |
|---|---|---|---|
| TeF₆ PP (seed = PySCF spherical atoms; PySCF DIIS) | 12 cycles | 5 | 4 |
| TeF₆ AE (seed = PySCF spherical atoms; PySCF DIIS) | 9 | 4 | 2 |
| LiF (ours) | 13 iterations | 2 | 2 |
| NaF (ours) | 15 iterations | 2 | 3 |
| MnO (ours) | 15 iterations | 3 | 1 |

## Audit cost against the SCF it replaces

Largest case MnO (4 atoms, 26 e, 4 OpenMP threads; D4 ran with no other run of mine beside it, other lanes' load ~8 on the box).
Wall from the bin files' mtimes and the run report:

| system, density | the SCF (setup + iterations) | reference (300, 77 GL) | ladder | per rung, median / max |
|---|---|---|---|---|
| **MnO, D4** | **220 s** (62 s setup + 15 iterations) | **752 s** | 2444 s (30 rungs) | 79 s / 201 s |
| MnO, Df | 224 s (201 s setup + 1 iteration) | 830 s | 1787 s | 49 s / 136 s |
| LiF, D4 | 329 s | 202 s | 1027 s (49 rungs) | 18 s / 44 s |
| NaF, D4 | 24 s | 201 s | 716 s (49 rungs) | 10 s / 41 s |

On MnO the audit as run (reference + 30 rungs) costs **53 min against a 3.7 min SCF, ×15**.  The REFERENCE alone is ×3.4 the SCF
it would size, and a single rung is a third to nearly the whole of it: each rung rebuilds the periodic Becke weights (order 3, the
ε-converged image series), which is the build cost the radial batch already measured.  A production audit would not walk 30 rungs
(walk L at one nR, then nR at the chosen L: about 8), but even then it is reference + 8 × 79 s ≈ 23 min, ×6 the SCF.  So with
today's Becke build an audit-sized run is not cheaper than running the SCF on a generous mesh; the audit pays only if the
reference is avoided (judge each rung against the next rung up, not a fixed fine reference) or the weight build gets cheaper
(the 'Becke order 2 on solids' row).

## What this does and does not show

1. **The two-step run works on all five systems.**  The rung the 1e-4 D picks at ε = 1e-5 IS the truth (the pick at the converged
   D) on every system and client: TeF₆ PP and AE (XC and ρ), LiF, NaF, MnO.  Gap 0 rungs, ×1.00 points, and at the converged D
   it delivers |ΔE| ≤ 0.61 εN and max|ΔH| ≤ 0.58 ε (ρ: ≤ 0.41 εN).  Step 2 restarts in 2–5 iterations (ours: 2, 2, 3).
2. **The seed alone is not reliable at ε = 1e-5.**  Of the five XC cases it picks the truth on three (TeF₆ PP, AE, NaF), NOTHING on
   LiF and on MnO-from-SAD (no rung passes at the seed), and UNDER-sizes MnO-from-IonicSAD by two angular rungs, a pick that
   delivers 2.33 ε on H at the converged D.  The seed's own H floor sits on either side of the converged D's (LiF 3–6× above, MnO-
   IonicSAD 3× below), so no fixed margin on the seed's numbers corrects it: one case needs a looser test, another a tighter one.
3. **At ε = 1e-4 the seed is good enough.**  XC: right on TeF₆ PP and AE, LiF and MnO (both seeds); one nR rung over on NaF (×1.50
   points, safe side).  ρ (TeF₆): one rung over on each (×1.63, ×1.43, safe side).  No pick at 1e-4 under-sizes.
4. **So the answer to 'does a 1e-4 D beat SAD': yes, on every system.**  Seed-size the 1e-4 step, converge, re-audit, restart on the
   1e-5 pick.  No calibration margin is needed for the XC client on this evidence (5 systems, 0 misses from D4); none is proposed.
5. **The ladder catches today's solid default (40,29) everywhere**, at every density, mostly on H: at the converged D it misses
   ε = 1e-5 by 177 ε (LiF), 11.6 ε (NaF), 105 ε (MnO) on max|ΔH|, against 5.6, 1.6 and 11.9 εN on E.  H, not E, sets the rung on
   all three solids.
6. **The ρ client's cheap prediction 2 E_H δq / N is not yet a sizing rule**: measured/predicted spans 0.15–3.4 per rung on TeF₆ and
   the calibration rests on two systems (one element pair); the solids cannot add to it until the δ Coulomb route exists on a cell.
7. **The audit is the expensive part** (×15 the SCF as run on MnO, ×6 for a lean walk), its reference mesh above all.

Not measured: the ρ client on solids (the δ route is refused on a cell); a seed of ours for Te (TeF₆ used PySCF spherical atoms);
a geometry or +U change between the two steps; the radial rule under test (all of this is today's MHL m2/m3 α2 default); kT > 0.
