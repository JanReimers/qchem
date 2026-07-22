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
| NaF (rocksalt) | Γ | VALENCE-LOWQ-SR | 320 | **−27.93128**¹ | +7.247 | +32.135 | −3.740 | — | — | 9m32s | `naf_gpw_sr_diag.inp` |

**Oracle re-validated on the new machine (2026-07-21, conda-forge CP2K 2026.1 — see `UnitTests/CP2K/README.md`
for the new run recipe):** Si Γ **−7.11505788** and NaF **−27.9312751** reproduce the rows above exactly
(NaF grid-charge leak 1.95e-4 e, same class as recorded).  **The NaF deck's actual multigrid (from the
`PW_GRID|` log blocks — the numbers the qchem grid-matched validation must force):**

| level | cutoff (Ha) | N (points/axis) | pair count (`gaussian_gridlevel`) |
|---|---|---|---|
| 1 | 160.0 | 36 | 3973 |
| 2 | 53.33 | 24 | 3149 |
| 3 | 17.78 | 12 | 3532 |
| 4 | 5.926 | 8 | 2342 |

(`REL_CUTOFF 60` Ry = 30 Ha: pair→level rule `needed = (αᵢ+αⱼ)·30` Ha, coarsest level ≥ needed, finest as
fallback — `gaussian_gridlevels.F`; progression factor 3.  Si deck at CUTOFF 80 Ry: levels 40/13.3/4.4/1.5 Ha,
N=24/12/8/4.)  Note CP2K's N=36 at 160 Ha vs our radix-2 FFT's N=128 at the same cutoff — same {G} ball,
45× the raster points on our side (their N is mixed-radix 2²·3²).  CP2K grid-side match knobs in qchem:
`GPW_MGRID_ECUTS` (explicit sub-level cutoffs) + `GPW_RELCUTOFF` (absolute Ha-per-exponent assignment rule).

**GRID-MATCHED qchem runs (2026-07-21; full record `doc/GPWPlan.md` §0e★+§0f, logs `~/Code/naf_gridmatched.log`
et al.): FINAL VERDICT = GPW VALIDATED — sub-mHa agreement on the same basis.**  The first matched run gave
−23.6739 vs the −27.9313 oracle (Δ=4.26 Ha), initially mis-attributed to the collocation method; two probes
decomposed it exactly: **0.76 Ha = a MOM-pinned EXCITED state** (the `AdoptMOMReference` coarse→fine transfer;
pure-aufbau rerun → −24.4317, clean) **+ 3.50 Ha = a BASIS MISMATCH** (qchem ran `VALENCE_LOWQ_SR2`; the
oracle deck ran `VALENCE-LOWQ-SR` — it predates the SR2 trim).  A 480-Ha-ball probe also falsified the
ball/Gibbs hypothesis (0.15 mHa effect).  **Apples-to-apples row (SR2, `naf_gpw_sr2_diag.inp`):**

| system | basis | qchem GPW (matched grids, aufbau) | CP2K | Δ |
|---|---|---|---|---|
| NaF Γ | VALENCE-LOWQ-SR2 | **−24.4316608** (22 iters, Δρ 2e-6) | **−24.4312134** (converged, 16 steps) | **0.45 mHa** |
| NaF Γ (2026-07-22: analytic short + κ rule + 5-smooth rasters) | VALENCE-LOWQ-SR2 | **−24.4314027** (22 iters, fine raster 72³, whole 2-stage run **190 s**) | −24.4312134 (5.8 s) | **0.19 mHa** |

Exc −4.95597 vs −4.95531 (0.7 mHa).  Note CP2K's SCF **fully converges** on SR2 (the eternal density
limit-cycle was the near-singular SR basis, cond≈8e3 — both codes struggled on SR, both are clean on SR2).
The SR↔SR2 3.5 Ha is REAL variational physics (both codes agree on both bases): Na's diffuse s 0.0857 +
p 0.05 are not physically null despite λ_min~1e-6 — see GPWPlan §1 (rank-reduction) for the
conditioning-vs-completeness resolution.  CP2K detailed components of the SR rerun: Core-H +7.2466, Hartree
+32.1349, XC −3.7398, core-charge self −63.5730, core-charge overlap +0.0000172.

¹ **NaF (2026-07-15, q1/q7 on OUR transcribed low-q SR basis): the ENERGY is settled, the DENSITY is not.**
Broyden(α=0.2, Kerker β=1.5) + traditional diagonalization: Etot flat at −27.9312754 to ~1e-6 from iteration
~130 on, but the density RMS gradient LIMIT-CYCLES at 0.03–0.12 forever (never reaches EPS_SCF 1e-6;
`IGNORE_CONVERGENCE_FAILURE` used for the clean exit).  **The SAME charge-transfer limit cycle our GPW shows**
— it is the system+basis (overlap cond ≈ 8e3), not either implementation.  OT is WORSE: the earlier OT run
(`naf_gpw_sr.out`) bounced between −25.7 and +253 and never settled E at all.  CP2K core-self −63.573
(its compensating-core split; compare TOTALS).  Grid charge: CP2K loses 2.0e-4 e at 320 Ry (our readout
loses ~1e-3-e-scale on the same system — same class).  **Our GPW must converge to −27.9313 to pass;
our 60-iter capped runs still swing (last E −24.03 is NOT settled) — mixing work needed (CP2K's working
recipe: plain damped Broyden, α=0.2, no DIIS-from-start).**

Charge = 8 for both; both SCF-converged. (Self-energy of the core charge −20.516 and the PP local/nonlocal
totals are k-independent; Core-H/Hartree/XC carry the k-dispersion.) CP2K's GPW electrostatic split differs
from ours (it uses a compensating-core-charge scheme) — compare the **total** and the cleaner sub-terms
(nonlocal-PP, XC).

## qchem comparison
- **Si Γ, SIPP_SR — the tight BASIS-MATCHED gate:** our GPW **−7.11505** vs CP2K **−7.11506** (1e-5), Exc
  −2.544 = −2.544, nonlocal-PP → +0.9406. This validated the bulk fix (see `doc/GPWPlan.md`, "Bulk
  over-binding FIXED").
- **Grid convergence, matched cutoff** (`CUTOFF` Ry = 2× our `densityEcut` Ha; our FFT `N` is `NextPow2`, so
  it jumps 32→64):

  | our grid | nominal | qchem GPW | CP2K (same cutoff) | vs CP2K-converged |
  |---|---|---|---|---|
  | N=32 (dE=20 Ha) | 40 Ry | −7.11467 | −7.115107 | +0.39 mHa |
  | N=64 (dE=30 Ha) | 60 Ry | −7.11505 | (flat, ≈−7.11506) | +0.01 mHa |
  | — | ≥80 Ry | — | **−7.11506** | reference |

  Our point-collocation total converges to CP2K-converged **from above**: N=64 matches to ~1e-5. At the matched
  40 Ry, N=32 sits ~0.4 mHa above CP2K's own 40 Ry (−7.115107) — CP2K itself only moves 5e-5 from 40→80 Ry, so
  that 0.4 mHa is OUR grid, not CP2K convergence: `NextPow2` point-collocation is a touch coarser than CP2K's
  analytic Gaussian-to-grid mapping at equal cutoff, closed by N=64. The fast gate
  (`GPW_SCF.DISABLED_SR_GammaRcut2a_CP2KReference`, N=32, ~45 s) pins −7.11467 with a 2e-3 tolerance that
  absorbs this; the tight match is N=64 (~5× slower).
- **Si 2×2×2 — multi-k GPW VALIDATED (2026-07-10).** Dispersive multi-k GPW now runs (KB Bloch-orbital fix +
  `Rcut = 2a`, SIPP_SR): charge 8, the total drops with k-sampling (Γ −7.11467 → 2×1×1 −7.451 → 2×2×2
  −7.7778 — real dispersion). **Grid-for-grid at the SAME Γ-centred mesh: our −7.7778 vs CP2K −7.77846
  (~0.7 mHa, the N=32 grid gap).** The 90 mHa vs CP2K's DEFAULT −7.86744 is purely the **k-mesh convention**:
  our GPW path is Γ-CENTRED only (`GPW_BasisSet` lrounds `k·N` → integer `ik/N`), while CP2K's
  `MONKHORST-PACK 2 2 2` is the classic SHIFTED grid (k at ±¼, confirmed from its k-point list). Matching the
  shifted grid needs the fractional k threaded through `BlochFactory` (a follow-up); the general-k PHYSICS is
  validated here. Decks: `si_fcc_gpw_222.inp` (shifted, −7.86744) + `si_fcc_gpw_222_gamma.inp` (Γ-centred,
  −7.77846). Test: `GPW_SCF.DISABLED_SR_2x2x2GammaCentred_vs_CP2K`.

  | mesh | convention | qchem GPW | CP2K (80 Ry) |
  |---|---|---|---|
  | 2×2×2 | Γ-centred (0, ½) | −7.7778 | −7.77846 |
  | 2×2×2 | shifted (±¼, CP2K default) | *(needs shifted-k support)* | −7.86744 |

## NaF: UNBLOCKED (2026-07-15) — own-basis transcription; CsI still blocked
NaF now runs (row above): the fix was route (1) below — `VALENCE-LOWQ-BASIS` (our
`valence_lowq{,_sr}.bsd` transcribed) + `BASIS_SET_FILE_NAME`, which carries no q tag so the
q1/q7 kinds pair cleanly.  The historical blocker, kept for CsI (iodine still has no basis anywhere):

Both were requested but **could not be run with CP2K's shipped data** because our qchem PPs use LOW valence q
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
