# CP2K reference inputs

A small library of CP2K (Quickstep GPW) input/control files used as **independent
oracles** for the qchem GPW work — CP2K is the reference Gaussian-and-Plane-Waves
implementation (Lippert–Hutter), so its total energy + per-term breakdown validate
our own GPW (see `doc/GPWPlan.md`). These are checked in as reference data, not run
by the C++ test suite.

Collected results (energies, breakdowns, runtimes) are tabulated in **`doc/CP2Kresults.md`**.

## Files
- `si_fcc_gpw.inp` — FCC-Si at Γ, GPW/LDA, matching our `GPW_SCF` setup: primitive cell
  a=10.26 bohr, 2 Si at `(0,0,0)+(¼,¼,¼)`, `GTH-PADE-q4` (== our GTH-LDA q4),
  `LDA_X + LDA_C_VWN` (Slater/Dirac exchange + VWN5). Corner atom kept at `(0,0,0)`
  deliberately (the raster-bug trigger).
- `si_fcc_gpw_222.inp` — same, with a `&KPOINTS MONKHORST-PACK 2 2 2` mesh (multi-k reference).
- `SIPP-SR-BASIS` — our SIPP_SR Si valence basis (`BasisSetData/sipp_sr.bsd`,
  3s3p uncontracted) transcribed into CP2K `BASIS_SET` format.
- (NaF/CsI: blocked on q-mismatched bases — CP2K ships no q1 Na/Cs basis, no I basis. See `doc/CP2Kresults.md`.)

## Reference result (CP2K 2026.1, serial ssmp)
FCC-Si Γ, SIPP_SR, GTH-PADE-q4, LDA_X+VWN5:
- **Total energy = −7.11506 Ha, charge = 8** (converged; −7.115058 by `CUTOFF` 80 Ry ≈ 40 Ha).
- Breakdown: Core-H (kin+PP) +5.565, Hartree +10.380, XC −2.544; PP total −7.548
  (local −8.489, nonlocal +0.941); core self-energy −20.516.
- `CUTOFF` convergence: −7.115107 (40 Ry) → −7.115058 (80 Ry) → flat.
- **This is the number our GPW+SIPP_SR should hit at Γ** once the collocation is off the
  FFT raster (`densityEcut` ≈ 30–40 Ha). NOT our PW −7.2273 (a different, plane-wave basis).

## How to run (built at ~/Code/cp2k, serial)
```sh
source ~/Code/cp2k/tools/toolchain/install/setup
export LD_LIBRARY_PATH=~/Code/cp2k/install/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
cd UnitTests/CP2K   # BASIS_SET_FILE_NAME is relative (./SIPP-SR-BASIS)
~/Code/cp2k/install/bin/cp2k.ssmp -i si_fcc_gpw.inp -o si_fcc_gpw.out
grep "ENERGY| Total FORCE_EVAL" si_fcc_gpw.out
```
Note: `si_fcc_gpw.inp` points `POTENTIAL_FILE_NAME` at CP2K's shipped
`~/Code/cp2k/data/GTH_POTENTIALS` — adjust that path on another machine. Parameter-matching
table (qchem ↔ CP2K) is in `doc/GPWPlan.md` (TODO 2).
