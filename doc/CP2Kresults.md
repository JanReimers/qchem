# CP2K reference results

Independent GPW oracles from **CP2K 2026.1** (serial `ssmp`, built at `~/Code/cp2k`) for the qchem GPW/PW
work. Input decks live in `UnitTests/CP2K/`; run recipe + the qchem↔CP2K parameter map are in that folder's
README and in `doc/GPWPlan.md` TODO 2. All runs: `METHOD GPW`, `LDA_X + LDA_C_VWN` (Slater/Dirac exchange +
VWN5), GTH-PADE PP (== our GTH-LDA), FCC/rocksalt/CsCl cells matching the `GPW_SCF`/`PlaneWaveDFT` tests.
`CUTOFF` is CP2K's density-grid cutoff in **Ry** (= 2× our `densityEcut` in Ha); converged values (Si: flat by
~80 Ry). CP2K has no orbital `Ecut` (Gaussians) and no `Rcut` knob (neighbour lists / `EPS_PGF_ORB`).

## Results

Decks use **`CUTOFF 80` Ry** — the minimal CONVERGED cutoff (the total is flat from 80 Ry: −7.115058 at
80/150/300/600 Ry; 40 Ry gives −7.115107, still ~5e-5 under). 80 Ry is ~4–7× faster than 300 for the same
number, so it's the right choice for a reference oracle.

| system | k-mesh | basis (CP2K) | CUTOFF (Ry) | **Etot (Ha)** | Core-H | Hartree | XC | PP loc | PP nonloc | time | inp |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | 80 | **−7.11506** | +5.565 | +10.380 | −2.544 | −8.489 | +0.941 | 3.5 s | `si_fcc_gpw.inp` |
| Si (FCC) | 2×2×2 | SIPP_SR | 80 | **−7.86744** | +4.384 | +10.671 | −2.407 | −8.489 | +0.941 | 4.8 s | `si_fcc_gpw_222.inp` |

Charge = 8 for both; both SCF-converged. (Self-energy of the core charge −20.516 and the PP local/nonlocal
totals are k-independent; Core-H/Hartree/XC carry the k-dispersion.) CP2K's GPW electrostatic split differs
from ours (it uses a compensating-core-charge scheme) — compare the **total** and the cleaner sub-terms
(nonlocal-PP, XC).

## qchem comparison
- **Si Γ, SIPP_SR — the tight BASIS-MATCHED gate:** our GPW **−7.11505** vs CP2K **−7.11506** (1e-5), Exc
  −2.544 = −2.544, nonlocal-PP → +0.9406. This validated the bulk fix (see `doc/GPWPlan.md`, "Bulk
  over-binding FIXED").
- **Si 2×2×2:** qchem GPW has no multi-k yet; our PW 2×2×2 is −7.7613 (Ecut=4, under-converged plane waves) —
  a loose check vs CP2K's converged −7.86744 (different basis). A future multi-k GPW at `Rcut ≥ 2a` + SIPP_SR
  should approach −7.86744.

## Blocked: NaF, CsI (basis / PP-q mismatch)
Both were requested but **cannot be run with CP2K's shipped data** because our qchem PPs use LOW valence q
that CP2K doesn't ship matching bases for:
- **Na q1, Cs q1:** CP2K ships only the semicore basis (`Na/Cs DZVP-MOLOPT-SR-GTH`, optimised for q9) — CP2K
  aborts: *"Basis-set and pseudo-potential were optimized for different valence electron numbers."* No q1
  Gaussian basis is shipped.
- **I (iodine):** no GTH-optimised Gaussian basis shipped at all.
- (F q7 is fine — `F DZVP-GTH-PADE`.)

Two ways forward (either lets NaF/CsI join this table):
1. **Hand-roll SIPP-style low-q valence bases** for Na(q1)/F(q7)/Cs(q1)/I(q7) — a few uncontracted s/p (and d
   for Cs/I) Gaussians, transcribed to CP2K `BASIS_SET` format (as `SIPP-SR-BASIS` was for Si). These would
   also be the qchem GPW bases when the multi-species GPW path is built — so it's shared work.
2. **Re-anchor the qchem NaF/CsI tests to CP2K's standard q** (Na q9, etc.) — but that changes the physics
   (semicore) and the PW anchors.

Recommended: (1), together with extending qchem GPW to multi-species (needs those bases anyway).

## How to run
See `UnitTests/CP2K/README.md`. In short: `source ~/Code/cp2k/tools/toolchain/install/setup`,
`export LD_LIBRARY_PATH=~/Code/cp2k/install/lib:$LD_LIBRARY_PATH`, then
`cp2k.ssmp -i UnitTests/CP2K/si_fcc_gpw.inp -o si.out` from a dir where `./SIPP-SR-BASIS` is visible.
