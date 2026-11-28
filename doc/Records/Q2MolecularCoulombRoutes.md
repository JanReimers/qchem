# Q2 RECORD: molecular Coulomb route × accuracy (lane S, 2026-11-27)

What the user decides from: which Coulomb route a molecular DFT run gets by default.  **No default has changed.**

- **Runs.** LDA, free (no imposed symmetry), converged SCF, `accuracy` ε = 1e-5 (e5) or 1e-6 (e6) Ha, 4 threads on a shared box (so timings are loose: re-time on a quiet box before quoting).
- **Columns.**
  - `E−E_exact`: the converged energy minus the exact-J converged energy at the same ε (same mesh).
  - `fixed-D ΔE_H`: the arm's Hartree energy minus exact J **at one fixed density** (the arm's converged D), so the route's own error with no SCF relaxation in it.  The one to judge against ε.
- **Arms.**
  - `exact`: analytic 4-index J (`coulomb.route: Exact`), the reference.  Cross-checked against PySCF J at our D: 46.652772136820 vs 46.6527721368.
  - `poisson`: the δ fit, Poisson-solved Hartree on the Becke mesh (**today's default**).
  - `analytic`: the δ fit with analytic (ij|δ_g) integrals.
    - Cost: 90 % generation, 2 % contraction; unscreened, matrix-free, 2 passes per iteration.
    - Not a production route yet.
  - `auxQ` / `auxNone`: a Gaussian auxiliary-basis fit (RI-J), with and without the charge constraint.
- **Binaries.** Fitted arms from the snapshot built at 6d87dc765; the exact column and the audits at 76536bd2f (the step-2 landing HEAD).  The fitted energies did not move between the two (the audit reruns reproduce them).

## Findings

1. **The Gaussian fits (auxQ/auxNone) miss ε everywhere, and do not improve at e6.**  0.1 mHa on water, 0.3 on CH₄, 2.9 on SF₆, 1.2 on TeF₆ PP, 10 on TeF₆ AE.  Fit-basis incompleteness: ε does not drive them.  They are fast and small.
2. **The δ Poisson default meets ε on the light molecules and on the PP cases, but misses on heavy all-electron TeF₆.**
   - TeF₆ AE: −1.4e-4 at e5 (14× ε) and −1.8e-4 at e6 (180× ε).  It does not improve with ε, so the accuracy knob is not reaching the term that limits it: an open defect, not a cost choice.
   - SF₆ AE: −3e-5 at e5 (3× ε), −4.6e-6 at e6 (inside ε).
   - CH₄ AE: −1.6e-5 at e5, −1.5e-6 at e6 (both ~1.5× ε).
3. **Analytic δ is the most accurate fitted route on the small cases (1e-7 on H₂O PP), but too slow.**  It is 20-500× Poisson's SCF time, and SF₆/TeF₆ were killed at the 3 h timeout.
4. **Exact J is cheap at these sizes.**  It has the fastest SCF in most rows.  Its RAM is mixed: higher than Poisson on SF₆ (0.9-1.1 GB vs 0.6-0.8 GB), lower on TeF₆ PP (1.4-1.6 GB vs 2.7-3.9 GB) and on TeF₆ AE (3.9-4.2 GB vs 4.8-5.6 GB).  Its 4-index cost grows as N⁴ in the basis, so this does not carry to large molecules.

| case | ε | arm | E−E_exact (Ha) | fixed-D ΔE_H (Ha) | iters | SCF s | RSS MB |
|---|---|---|---|---|---|---|---|
| H2O_PP | e5 | exact | 0 (E=-17.1822824961) | — | 11 | 3.56 | 163.5 |
| H2O_PP | e5 | poisson | +1.73e-05 | +1.73e-05 | 7 | 3.237 | 164.9 |
| H2O_PP | e5 | analytic | +1.26e-07 | +1.26e-07 | 11 | 1087 | 134.9 |
| H2O_PP | e5 | auxQ | -1.21e-04 | -1.22e-04 | 7 | 2.646 | 180.4 |
| H2O_PP | e5 | auxNone | -1.19e-04 | -1.20e-04 | 7 | 2.672 | 179.4 |
| H2O_PP | e6 | exact | 0 (E=-17.1822824961) | — | 11 | 4.093 | 163 |
| H2O_PP | e6 | poisson | +2.33e-06 | +2.31e-06 | 10 | 22.95 | 470 |
| H2O_PP | e6 | analytic | +1.26e-07 | +1.26e-07 | 11 | 1233 | 143 |
| H2O_PP | e6 | auxQ | -1.21e-04 | -1.22e-04 | 7 | 3.024 | 177.9 |
| H2O_PP | e6 | auxNone | -1.19e-04 | -1.20e-04 | 7 | 3.37 | 177.3 |
| H2O_AE | e5 | exact | 0 (E=-75.8772997308) | — | 7 | 0.8725 | 135.2 |
| H2O_AE | e5 | poisson | -8.92e-06 | -9.53e-06 | 11 | 14.26 | 173.7 |
| H2O_AE | e5 | analytic | +8.76e-06 | +8.76e-06 | 7 | 76.08 | 140 |
| H2O_AE | e5 | auxQ | -9.44e-05 | -9.48e-05 | 8 | 3.218 | 152.9 |
| H2O_AE | e5 | auxNone | -9.44e-05 | -9.48e-05 | 8 | 2.894 | 152.3 |
| H2O_AE | e6 | exact | 0 (E=-75.8772997308) | — | 7 | 0.9335 | 136.1 |
| H2O_AE | e6 | poisson | +3.89e-06 | +3.10e-06 | 9 | 23.5 | 363 |
| H2O_AE | e6 | analytic | +1.08e-06 | +8.51e-07 | 8 | 213.2 | 284 |
| H2O_AE | e6 | auxQ | -9.44e-05 | -9.48e-05 | 8 | 2.291 | 153.4 |
| H2O_AE | e6 | auxNone | -9.44e-05 | -9.48e-05 | 8 | 2.256 | 152.5 |
| CH4_PP | e5 | exact | 0 (E=-8.03523969481) | — | 10 | 9.751 | 326.7 |
| CH4_PP | e5 | poisson | +2.61e-06 | +2.61e-06 | 7 | 15.06 | 309.7 |
| CH4_PP | e5 | analytic | -8.25e-06 | — | 10 | 4353 | 304.3 |
| CH4_PP | e5 | auxQ | -3.61e-04 | -3.64e-04 | 7 | 8.408 | 344.3 |
| CH4_PP | e5 | auxNone | -3.60e-04 | -3.63e-04 | 7 | 8.722 | 343.2 |
| CH4_PP | e6 | exact | 0 (E=-8.03523969481) | — | 10 | 9.211 | 327 |
| CH4_PP | e6 | poisson | +3.26e-06 | +4.63e-06 | 8 | 72.12 | 769.2 |
| CH4_PP | e6 | analytic | -6.97e-07 | — | 11 | 7333 | 521.3 |
| CH4_PP | e6 | auxQ | -3.61e-04 | -3.64e-04 | 7 | 9.644 | 343.3 |
| CH4_PP | e6 | auxNone | -3.60e-04 | -3.63e-04 | 7 | 9.175 | 343.8 |
| CH4_AE | e5 | exact | 0 (E=-40.1052899352) | — | 6 | 1.9 | 231 |
| CH4_AE | e5 | poisson | -1.57e-05 | -1.57e-05 | 6 | 13.07 | 235.4 |
| CH4_AE | e5 | analytic | -1.17e-05 | — | 6 | 191.5 | 209 |
| CH4_AE | e5 | auxQ | -2.98e-04 | -3.00e-04 | 6 | 5.514 | 254.4 |
| CH4_AE | e5 | auxNone | -2.97e-04 | -2.99e-04 | 6 | 5.443 | 253.9 |
| CH4_AE | e6 | exact | 0 (E=-40.1052899352) | — | 6 | 2.219 | 231.3 |
| CH4_AE | e6 | poisson | -5.67e-08 | -1.51e-06 | 7 | 79.4 | 654.8 |
| CH4_AE | e6 | analytic | -1.52e-06 | — | 9 | 466.1 | 414.5 |
| CH4_AE | e6 | auxQ | -2.98e-04 | -3.00e-04 | 6 | 5.713 | 254.3 |
| CH4_AE | e6 | auxNone | -2.97e-04 | -2.99e-04 | 6 | 5.027 | 254.1 |
| SF6_AE | e5 | exact | 0 (E=-983.452628971) | — | 10 | 14.15 | 896.4 |
| SF6_AE | e5 | poisson | +1.49e-05 | -2.96e-05 | 13 | 165.6 | 627.7 |
| SF6_AE | e5 | analytic | killed at the 3 h timeout, no energy | | | | |
| SF6_AE | e5 | auxQ | -2.82e-03 | -2.92e-03 | 9 | 10.11 | 214.5 |
| SF6_AE | e5 | auxNone | -2.82e-03 | -2.92e-03 | 9 | 11.1 | 214.6 |
| SF6_AE | e6 | exact | 0 (E=-983.452584488) | — | 11 | 90.91 | 1097 |
| SF6_AE | e6 | poisson | -4.66e-06 | -4.64e-06 | 13 | 183.2 | 833.7 |
| SF6_AE | e6 | analytic | killed at the 3 h timeout, no energy | | | | |
| SF6_AE | e6 | auxQ | -2.82e-03 | -2.92e-03 | 12 | 110 | 599.7 |
| SF6_AE | e6 | auxNone | -2.82e-03 | -2.92e-03 | 12 | 77.8 | 597.6 |
| TeF6_PP | e5 | exact | 0 (E=-153.407208931) | — | 11 | 122 | 1362 |
| TeF6_PP | e5 | poisson | +1.74e-06 | -2.10e-05 | 15 | 113.9 | 2748 |
| TeF6_PP | e5 | analytic | killed at the 3 h timeout, no energy | | | | |
| TeF6_PP | e5 | auxQ | -1.22e-03 | -1.24e-03 | 9 | 14.46 | 937.8 |
| TeF6_PP | e5 | auxNone | -1.19e-03 | -1.21e-03 | 9 | 14.78 | 935 |
| TeF6_PP | e6 | exact | 0 (E=-153.407190752) | — | 14 | 120.1 | 1554 |
| TeF6_PP | e6 | poisson | +5.92e-06 | +1.99e-06 | 13 | 222.5 | 3932 |
| TeF6_PP | e6 | analytic | killed at the 3 h timeout, no energy | | | | |
| TeF6_PP | e6 | auxQ | -1.22e-03 | -1.24e-03 | 10 | 45.2 | 1102 |
| TeF6_PP | e6 | auxNone | -1.20e-03 | -1.21e-03 | 10 | 49.63 | 1102 |
| TeF6_AE | e5 | exact | 0 (E=-7194.67283823) | — | 13 | 83.98 | 4179 |
| TeF6_AE | e5 | poisson | -1.40e-04 | -1.39e-04 | 15 | 220.8 | 4827 |
| TeF6_AE | e5 | analytic | killed at the 3 h timeout, no energy | | | | |
| TeF6_AE | e5 | auxQ | -9.94e-03 | -9.97e-03 | 13 | 73.87 | 2806 |
| TeF6_AE | e5 | auxNone | -9.94e-03 | -9.97e-03 | 13 | 77.19 | 2810 |
| TeF6_AE | e6 | exact | 0 (E=-7194.67283823) | — | 13 | 216.4 | 3910 |
| TeF6_AE | e6 | poisson | -1.75e-04 | -1.78e-04 | 15 | 469.2 | 5592 |
| TeF6_AE | e6 | analytic | killed at the 3 h timeout, no energy | | | | |
| TeF6_AE | e6 | auxQ | -9.94e-03 | -9.97e-03 | 13 | 109.5 | 2772 |
| TeF6_AE | e6 | auxNone | -9.94e-03 | -9.97e-03 | 13 | 145.3 | 2734 |

## Provenance

- **Decks and logs** (run data, not in the repo): `~/Code/qchem6-runs/Q2/` (`decks/`, `logs/`, `rows.txt` fitted arms, `rows_exact.txt` exact column + audits, `run_fitted.sh`, `run_exact.sh`).
- **Driver:** `CLIapps/q2probe` (`--audit` prints `auditExact_dEH` / `auditExact_maxdJ`: the arm's J and E_H against the four-centre J at the arm's converged D).
- **PySCF cross-checks:** `~/Code/qchem6-runs/Q2/xcheck/` (`jcheck.py`: J at our D; `hybrid.py`: the LDA-hybrid water references that `M_DFT.WaterLDAHybrid_*` gates).
- **Open work** cut from this record: OpenWork §3 "The δ Poisson Hartree misses ε on heavy all-electron molecules"; OpenWork §4 "The Analytic (ij|δ_g) Hartree is ~300× Poisson".
