# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**Step 1 of `doc/OpenWork.md`.  This is an INSTRUMENT, not a report.**  Steps 2–4 (folds, Φ-screening, RAM)
have no finish line without it: everything optimised in the 2026-08 runtime campaign was measured against a
single hand-run disabled test at `GPW_MNO_NMAX=3`, which is how two charter premises survived for months
while being wrong.  A row here is only meaningful with its PROVENANCE columns — same span held full-rank by
both codes, and the fold/threading/BLAS state that produced the number.

## How to produce a row — the SAME wrapper for both codes

```bash
scripts/bench "Si Gamma qchem" -- build/Release/IntegrationTests/ITMain --gtest_filter=GPW_SCF.SiliconGammaConverges
scripts/bench "Si Gamma cp2k"  -- cp2k -i IntegrationTests/CP2K/si_fcc_gpw.inp
```

**Peak RAM is measured from OUTSIDE, identically for both codes** (`/usr/bin/time -v` → "Maximum resident
set size" = the kernel's `VmHWM` high-water mark).  That is the whole reason the RAM column is now
obtainable without touching CP2K: no in-program hook was ever needed, and measuring both sides through one
wrapper is what makes the two columns comparable.  On a qchem row the wrapper also prints qchem's OWN
internal `VmHWM` as a cross-check — they should agree (measured: 266 vs 266.7 MB on Si Γ).

For the qchem detail (per-bucket ledger, folds, stream coverage) add `GPW_REPORT=1`:

```bash
GPW_REPORT=1 build/Release/IntegrationTests/ITMain --gtest_filter=<test> --gtest_also_run_disabled_tests
```

Every GPW run now prints, without extra flags:

| what | where it comes from |
|---|---|
| `Etot` at **10 s.f.**, iterations, convergence verdict | the run fingerprint + summary lines |
| wall time, and the per-bucket ledger | the `timing` section |
| **PEAK RSS (MB, process high-water)** | `timing`, from Linux `VmHWM` (added 2026-08-19) |
| `[fold] <site>: … = F×` for every fold site | `EmitFold`, Step 0b |
| `[stream cache] … pts64/pts32/runs/meanRun/dropped` | the stream-cache readout |
| `[site moments] … [e]` (polarized) | `QCHEM_SITE_MOMENTS=1`, Step 0a |

…**including the ANNEALED driver**, which had none of the summary half until 2026-08-19 — and every MnO row
runs through it, so the "every GPW run reports PEAK RSS" claim above was false for precisely the runs whose
RAM this table most needs.  Two precision leaks were repaired at the same time: the `[fold]` line left
`cout` at 2 s.f. for the rest of the run (`std::defaultfloat` restores the format flag, not the precision),
and the fingerprint's `Efinal` inherited whatever a verbose table had left behind.  **A number's precision
must not depend on which other diagnostics were switched on** — every line above now states and restores
its own.

**Read the provenance, not just the number.**  `PEAK RSS` is the process high-water mark, so it is only a
clean per-config figure when the process runs ONE config — a `--gtest_filter` naming several tests reports
the watermark of the whole process.  `GPW_OMP_THREADS` governs the GPW **pair loops only** and says nothing
about the BLAS: with it unset, blaze still ran these rows at 115–239% CPU.  BLAS routing is
`QCHEM_BLAZE_BLAS` (default ON).  **qchem has no equivalent of CP2K's `GLOBAL| Number of threads` banner
line — a run cannot state how parallel it actually was**, which is why `scripts/bench` reports measured CPU%
and CPU seconds rather than trusting a knob.  Worth closing: one line at run start stating the effective
BLAS/OpenMP thread counts would make every future row self-describing.

## The table

Energies in Ha.  **Both columns measured on this box (14 GB, 16 cores) through `scripts/bench`, 2026-08-19** —
the CP2K side is no longer banked prose: `apt`'s CP2K 2025.2 reproduces every banked 2026.1 deck value to the
printed digits (`doc/CP2KBuild.md`), so both codes are measured under one wrapper.  Provenance per row is in
*How each row was produced*.

> **⚠ THE TWO SIDES DO NOT USE THE SAME NUMBER OF CORES, SO READ THE CPU COLUMN.**  CP2K is genuinely serial
> here (`OMP_NUM_THREADS=1`, measured 97–99% CPU).  qchem is NOT: `GPW_OMP_THREADS` is unset, but blaze runs
> the BLAS multi-threaded regardless — measured 115–239% CPU.  **A thread-count knob is not a measurement**,
> and an earlier cut of this table printed "1 thr" on the qchem rows on the strength of the unset knob.
> Wall time therefore FLATTERS qchem by the number of cores it took; **CPU time (user+sys) is the honest
> ratio** and it is roughly 2× the wall ratio.

| system | k-mesh | span | qchem Etot | CP2K Etot | Δ (qchem−CP2K) | wall q / c | **CPU q / c** | **CPU ×** | peak RSS q / c |
|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | −7.115068508 | −7.115057882 | **−10.6 µHa** | 6.1 s / 5.2 s | 7.4 / 5.0 s | 1.5× | 267 / 148 MB |
| Si (FCC) | 2×2×2 Γ-centred | SIPP_SR | −7.778472814 | −7.778457865 | **−14.9 µHa** | 8.3 s / 5.8 s | 9.6 / 5.6 s | 1.7× | 269 / 153 MB |
| Si (FCC) | 2×2×2 shifted MP | SIPP_SR | ⛔ −3.735074784 ¹ | −7.867436530 | ¹ | 29.1 s / 6.1 s | — ¹ | — | 269 / 153 MB |
| NaF (rocksalt) | Γ | LOWQ_SR2 (both) | −24.430336482 | −24.431213375 | **+0.877 mHa** | 39.5 s / 7.4 s | 94.5 / 7.2 s | **13.2×** | 577 / 173 MB |
| NaF (rocksalt) | 2×2×2 Γ-centred | LOWQ_SR2 | −24.546883793 | — ² | | 57.4 s / — | 112 s / — | — | 590 / — MB |
| NaF (rocksalt) | Γ | LOWQ_SR (full) | −24.430944039 | −24.432293467 | **+1.349 mHa** | 2m44s / 1m42s | 219 / 102 s | 2.1× | 3090 / 186 MB |
| MnO AFM-II | Γ | **VA (N=118)** | −61.402976200 | −61.303325178 | **−99.65 mHa** | 20m05s / 6m14s | 2240 / 373 s | **6.0×** | 4947 / 217 MB |
| MnO FM | Γ | **VA (N=118)** | −61.441583060 | −61.304782531 | **−136.80 mHa** | 21m45s / 3m13s | 2321 / 192 s | **12.1×** | 4947 / 217 MB |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | | ❓ |

Δ(AFM−FM) on the VA span: **qchem +38.61 mHa, CP2K +1.46 mHa** — both order FM first, and the
**configuration-SELECTIVE part of the offset is −37.15 mHa** (`OpenWork` Step 5).  Every one of these four
energies reproduces its banked value (runs 61/62 and the CP2K VA pair) to the digits those were recorded at.

¹ ⛔ **BLOCKED BY A DEFECT AT EXACTLY k=¼, diagnosed 2026-08-19.**  `DISABLED_SR_2x2x2ShiftedMP_vs_CP2K`
banked −7.86673 against CP2K's −7.86744 when the complex-k path landed (`doc/GPWHistory.md`, `745d03ff`);
today it returns −3.7351.  **The shifted 2×2×2 mesh's k are (±¼,±¼,±¼) — and ¼ is the one place this fails.**
`GPW_SCF.DISABLED_SingleKSweepProbe` runs ONE k-point per run (mesh 1×1×1, so k = `GPW_KSHIFT` exactly), no
weights, no IBZ, no symmetry — E(k) must be smooth, and it is, except for a single-point blow-up:

| k | 0 | ⅛ | 0.249 | **¼** | 0.2500001 | 0.251 | ⅓ | ⅜ | ⅝ | ¾ | ⅞ | ½ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Etot | −7.1151 | −7.3121 | −7.5638 | **−5.0278** | −7.5660 | −7.5672 | −7.6830 | −7.7189 | −7.7189 | **−5.0174** | −7.3121 | −7.7763 |

**The failing set is exactly {¼, ¾}** — where the Bloch phase \f$e^{2\pi ik}\f$ is purely imaginary and its
REAL part underflows to \f$\cos(\pi/2)\approx6\times10^{-17}\f$ instead of being O(1).  ⅛, ⅜, ⅝, ⅞ and ⅓ are
all equally complex and all clean (and ⅞/⅝ reproduce ⅛/⅜ to every digit, so k↔−k is sound).  Eliminated, each
by measurement: **symmetry imposition** (the free run is wrong the same way), **the collocation streams**
(budgets off → analytic path → same energy to 8 digits, converged in 18 iterations), **frontier smearing**
(kT up to 0.02 moves the answer only −5.028 → −5.116, nowhere near the −7.57 the curve wants), **the
TRIM/real typing** (`BlochQN::IsReal()` is exact and correctly answers *false* at ¼; the blocks are built
complex).  **And it is not physics: k = ¼ + 1e-9 lands back on the
smooth curve** (−7.5224, still descending) — a 1e-9 change in k cannot move a converged energy by 2.5 Ha, so
this is an exact-value branch.  Signature at the bad k: Ekin 5.67 where the curve wants 3.63, Een −0.02
where it wants −0.71.

**★ THE STATIC OPERATORS ARE EXONERATED** — `GPW.GeneralK_OneElectronMatricesAreContinuousAtQuarterK` (new,
ENABLED, 2 s) builds the matrices at k = 0.249 / 0.250 / 0.251 with no SCF, no density and no occupations,
and every one is smooth and monotone through the bad k:

| | 0.249 | **0.250** | 0.251 |
|---|---|---|---|
| ‖S‖ | 6.72656 | 6.72899 | 6.73143 |
| ‖T‖ | 21.9720 | 21.9719 | 21.9718 |
| ‖V‖ | 121.742 | 121.754 | 121.766 |
| ‖V_loc-long‖ | 11.9330 | 11.9333 | 11.9336 |
| ‖V_loc-short‖ | 8.14542 | 8.14551 | 8.14560 |
| ‖V_KB‖ | 17.8310 | 17.8315 | 17.8319 |

So S, T, V, both local-PP pieces and the KB nonlocal are all CORRECT at k=¼ — including the KB projector,
the historic home of complex-only phase bugs.  Two more gates (also new, also enabled) close the rest of the
static path:

- **the collocation/integrate-back**: \f$\langle\chi_i^k|V_0|\chi_j^k\rangle=V_0S^{Bloch}_{ij}\f$ holds at any
  k, and its residual is 7.54e-8 / 7.67e-8 / 7.80e-8 at 0.249 / 0.250 / 0.251 — smooth and tiny.  (The
  pre-existing version of this invariant runs at Γ only and compares only REAL parts, so at complex k, where
  the imaginary part IS the physics, it had never been checked.)
- **the 1E spectrum**: solving \f$HC=\varepsilon SC\f$ for the static H, every level moves by ~1e-6 over
  dk=1e-5 and the frontier gap is a smooth **0.129648 / 0.129649 / 0.129650**.

**★ SO THE OPERATORS AND THEIR SPECTRUM ARE CLEAN, AND THE SCF STILL DIVERGES AT ITERATION 1.**  From the
same k-independent uniform seed: E₁ = −5.62164 at k=0.24999, **−4.86010 at k=0.25**, −5.62148 at 0.25001 —
a 0.76 Ha jump for a 1e-5 change in k.  All three gates stay enabled as the first coverage quarter-integer
k has ever had.

**★★ THE SYMPTOM, NAMED (2026-08-19).**  With the frontier window now available on SOLIDS and printing
DEGENERACY (`ε(occ/degen)`, `GPW_BANDGAP=1`), iteration 1 reads:

| k | occupied levels ε(occ/degen) | pattern |
|---|---|---|
| 0.25001 | −0.2197(2/2)  0.0592(2/2)  **0.1604(4.0/4)** | `1,1,2` ✓ |
| **0.25** | −0.7197(2/2) −0.1013(2/2) −0.0237(2/2) **0.1449(2.0/4)** | `1,1,1` + half-filled doublet ✗ |
| 0.5 (L) | −0.2056(2/2) −0.1339(2/2)  0.1023(4/4) | `1,1,2` ✓ |
| 0 (Γ) | −0.1315(2/2)  0.3080(6/6) | `1,3` ✓ (cubic) |

**Si's four valence bands along Λ are Λ₁, Λ₁, Λ₃(×2) — the `1,1,2` its neighbours and the L point show.
At exactly ¼ an EXTRA SINGLET appears below the frontier and the doublet is left straddling E_F with 2 of
its 4 states filled** — that is the `m` flag, and a half-filled shell spreads 2 electrons over 4 states,
which is a symmetry-broken density by construction.  The whole spectrum shifts with it (the lowest level
moves −0.2197 → −0.7197, half a hartree, for dk=1e-5).

Both spectra carry 24 orbitals / 48 states, so no function is lost.  **Also ruled out: the
orthogonalisation** — `GPW_ORTHO=cholesky|eigen|svd` all give ≈−5.02 at the bad k.

**★★★ AND THE MATRICES ARE SMOOTH ELEMENT-WISE, not merely in norm.**  The continuity gate originally
compared Frobenius norms, which cannot see the bug class it exists for — \f$\|\overline{M}\|=\|M\|\f$, so a
conjugated Bloch phase (the historic complex-k defect) passes a norm check at every k.  It now compares
M(¼) against the midpoint of its neighbours **entry by entry**, and covers one more path: the potential
matrix at **G≠0**, which the constant-field gate does not reach (its closure is nonzero only at
`dm=(0,0,0)`, i.e. G=0 alone).  At dk=1e-4 the worst entry of each is

| S | T | V_loc-long | V_loc-short | V_KB | V(G≠0) collocation |
|---|---|---|---|---|---|
| 5.03e-8 | 4.20e-8 | 2.41e-8 | 1.53e-8 | 4.27e-8 | 8.29e-8 |

and **these are curvature, not defect**: shrinking dk ten-fold twice drops every entry by exactly 100×
(5.03e-6 → 5.03e-8 → 5.03e-10), the signature of a smooth function sampled at three points.

**So every matrix that enters H is smooth at k=¼, the ortho is irrelevant, and the assembled H's spectrum
still jumps.**  That is a sharp contradiction, and it puts the defect in the ASSEMBLY step itself — how the
Hamiltonian combines these (correct) matrices into the block it hands the solver — or in a k-dependent
ingredient not yet enumerated.  **Caveat that bounds the inference**: the static gate's H
(T/2 + V_loc + V_KB) is a stand-in; its occupied-level PATTERN matches the SCF's at the good k (`1,1,2`)
and not at ¼, which is suggestive, but it is not the SCF's own H.
The CP2K half of this row is measured and sound.  The only test covering shifted k is DISABLED, so nothing
caught it — re-enable it, or add a cheap quarter-integer gate, once fixed.
**Why the two Si 2×2×2 rows differ by 89 mHa** (asked 2026-08-19; settled by running CP2K to the k-limit,
seconds per point, so it never has to be argued again).  They are two samplings of ONE integral and converge
to the same answer — the gap is sampling error at a very coarse mesh, not a setup discrepancy:

| CP2K mesh | E | error vs limit |
|---|---|---|
| 2×2×2 Γ-centred (8 pts) | −7.778458 | 95.8 mHa |
| 3×3×3 MP (odd n ⇒ Γ-centred) | −7.851866 | 22.3 mHa |
| **2×2×2 shifted (8 pts)** | **−7.867437** | **6.8 mHa** |
| 4×4×4 shifted | −7.874138 | 0.08 mHa |
| 6×6×6 / 8×8×8 | −7.874213 / −7.874215 | converged |

At identical cost (8 points) the shifted grid is **14× closer to the limit** — the textbook Monkhorst–Pack
result, since a Γ-centred grid on FCC spends its points on Γ and the zone boundary where the bands are
extremal.  So qchem's Γ-centred −7.778473 is being compared against the right reference for ITS sampling,
and the shifted row (once footnote ¹ is fixed) is the one that would say something about k-convergence.

² CP2K's NaF decks carry no `&KPOINTS` section, i.e. they are Γ.  The qchem test's own default is a **2×2×2**
mesh (8 k → 3 irreducible) — 116 mHa of band dispersion below its Γ value — so the row that compares to CP2K
is the `NAF_KMESH=1` one, and the 2×2×2 row is a qchem-only cost/size datapoint until a k-point CP2K deck
exists.  This mismatch was live in the test until 2026-08-19: it carried a Γ-era anchor (−24.4304) against a
2×2×2 configuration and simply failed.
³ the full-SR span runs FULL RANK at Γ ("kept 32 of 32", λ_min 4.35e-4) — the historical near-null trouble
was the multi-k cell.  The retracted −27.93128 "oracle" (a 3.5 Ha `EPS_PGF_ORB` screening artifact) is gone
from this table for good.

### How each row was produced

```bash
# qchem  (GPW_REPORT=1 for the ledger; the energy line now prints 10 s.f.)
GPW_REPORT=1 scripts/bench "Si Gamma qchem"  -- build/Release/IntegrationTests/ITMain --gtest_filter=GPW_SCF.SiliconGammaConverges
# NaF: NAF_KMESH picks the mesh (1 = the CP2K-comparable Γ), NAF_SPAN the basis (sr2 default, sr = full)
GPW_REPORT=1 NAF_KMESH=1 NAF_SPAN=sr2 scripts/bench "NaF SR2 Gamma qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.DISABLED_NaFRocksaltGamma --gtest_also_run_disabled_tests
# MnO, VA span, the runs-61/62 recipe -- now selected BY NAME instead of by overwriting a committed file.
# ONE ARM PER INVOCATION (MNO_SKIP_FM / MNO_SKIP_AFM): peak RSS is a PROCESS watermark, so a run that does
# both arms reports one number for the pair.
GPW_SPHERICAL=1 GPW_BASIS_SPAN=va MNO_ANNEAL="5e-3,0" MNO_ACC="Ladder,GDM" MNO_MOM=0 \
MNO_ORTHO_TOL=1e-3 MNO_SHARED_MU=1 MNO_IMPOSE=1 MNO_SKIP_FM=1 GPW_REPORT=1 \
    scripts/memsafe scripts/bench "MnO AFM2 VA qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.DISABLED_MnO_AFM2_RhombohedralGamma --gtest_also_run_disabled_tests

# CP2K -- from the DECK'S directory (relative BASIS_SET_FILE_NAME) and at one thread
cd IntegrationTests/CP2K
OMP_NUM_THREADS=1 ../../scripts/bench "Si Gamma cp2k"    -- mpirun -np 1 cp2k.psmp -i si_fcc_gpw.inp
OMP_NUM_THREADS=1 ../../scripts/bench "MnO AFM2 VA cp2k" -- mpirun -np 1 cp2k.psmp -i mno_afm2_gpw_va.inp
```

Verify each side reproduces its own history before reading a Δ: the CP2K decks against `doc/CP2Kresults.md`
(all five re-validated 2026-08-19, `doc/CP2KBuild.md`), and the qchem runs against the tests' own anchors.

### What is still missing

- **★ REPEAT THE WHOLE TABLE AT 12 THREADS** (user, 2026-08-19).  This cut is the SERIAL baseline — CP2K
  genuinely serial, qchem serial-except-BLAS — which is the right reference for "how much work does each
  code do", and the wrong one for "how fast is it on this box".  Both sides scale from the environment:
  CP2K's `psmp` takes `OMP_NUM_THREADS`, qchem takes `GPW_OMP_THREADS` for the pair loops plus whatever
  blaze does with the BLAS.  Run every row at `OMP_NUM_THREADS=12` / `GPW_OMP_THREADS=12` and keep BOTH
  cuts: the serial one is the algorithmic comparison, the threaded one is the user-facing time, and the
  ratio between them is qchem's parallel efficiency — which nothing currently measures.
  **Deliberately deferred until the qchem times come down** (user): re-measuring an 20-minute row now, when
  Steps 2–3 are about to change it by an order of magnitude, buys a number with a short shelf life.
- **The multi-k MnO row** (`MNO_KMESH=2`, ALL-TRIM so the whole k-sum runs real) — cost unmeasured, and CP2K
  needs a matching `&KPOINTS` deck.  At 20 min for the Γ row, budget for it before starting.
- **A k-point CP2K deck for NaF**, so the 2×2×2 row becomes a comparison instead of a qchem-only datapoint.
- **The Si shifted-MP row**, which needs the fractional-k regression fixed first (footnote ¹).
- **An A/B on NaF's XC route.**  Today's direct runs agree with CP2K to +0.88 mHa (SR2) and +1.35 mHa (SR),
  whereas July's numbers were 0.19 / 0.10 mHa — but those came from the GRID-CONTINUATION test
  (`DISABLED_NaFGridContinuation`, coarse→fine seeding, uniform-XC era) and not from this one, so it is not
  a like-for-like drift.  One A/B would say whether ~1 mHa is NaF's honest agreement at today's defaults or
  the price of a recipe change; until then do not quote NaF as "0.1 mHa class".

### Like-for-like: what "same span" costs, and how it is now held

CP2K can be forced down to the spherical spans qchem holds at FULL RANK — `IntegrationTests/CP2K/`
`mno_{afm2,fm}_gpw_v{a,b}.inp` with the `VALENCE-LOWQ-V{A,B}` entries in `VALENCE-LOWQ-BASIS`; **VA = N 118**
(full rank in both codes), **VB = N 128**.  So the MnO comparison does NOT wait on the 136-span question
(`OpenWork` Step 6).

On the qchem side that span used to be produced by `doc/scripts/bisect_valence_sph.py` **overwriting the
committed `BasisSetData/valence_lowq_sph.bsd` in the working tree** — so a run could not state which span it
had used, and no row over it was reproducible after the file was restored.  VA and VB are now committed
basis sets (`valence_lowq_v{a,b}.bsd`, `BasisSetData::VALENCE_LOWQ_V{A,B}`) selected by
**`GPW_BASIS_SPAN=va|vb|sph|sr`**, and the run prints its own `nFunctions` — reproducing run 61's basis block
exactly (118 functions, λ_min 1.29e-3, cond 4.41e3).  Si and NaF already share one span per material.

### Fold state and cost profile of the MnO rows (their provenance, and Steps 2–3's target)

The MnO rows above run `MNO_IMPOSE=1`, so unlike the FREE production run they DO fold — printed by the run
itself (Step 0b):

| site | this row | free production run |
|---|---|---|
| XC mesh (Becke star-average) | **23.03×** (24 ops, magnetic/Shubnikov) | NONE |
| `V_loc`-long {G}-star | **10.37×** (12 ops) | NONE |
| collocation streams (T3 pairs) | **NONE** (8778 items in full) `[opt-in: GPW_STREAM_FOLD=1]` | NONE |

So the imposed row already claims two of the three folds, and the one still unclaimed is the real-space pair
stream — `OpenWork` Step 2's biggest item, measured at 71× on the cache.

Where the 20 minutes go (AFM arm, `GPW_REPORT=1` ledger, seconds of 1205):

| bucket | s | share |
|---|---|---|
| setup: XC-mesh **Φ tables** | 370.0 | **31%** |
| scf: collocate density (pair scatter) | 263.4 | 22% |
| scf: integrate-back (pair gather) | 167.6 | 14% |
| setup: collocation stream build | 110.1 | 9% |
| setup: Becke mesh build | 70.6 | 6% |
| scf: XC-mesh ρ sampling (matrix-free) | 57.1 | 5% |
| scf: XC-mesh H_xc quadrature | 33.8 | 3% |

**The single largest bucket is the Φ-table build**, not the scatter/gather that rounds 3–4 were spent on —
and it is precisely `OpenWork` Step 3's named "Φ-table SCREENING" lever (Φ is stored dense, npts × n, while
the true object is sparse).  On this cell it is bigger than either pair loop.

## Standing observations

- **Si Γ agrees to 1.1e-5 Ha** at 1.5× CP2K's CPU time and 1.8× its RAM.  On the small cells the gap is
  modest; the diffuse and magnetic cells are where it opens up.
- **NaF Γ agrees to 0.9 mHa (SR2) / 1.3 mHa (full SR)** at 13.2× / 2.1× the CPU time and 3.3× / 16.6× the
  RAM.  The full span costs us 3090 MB where CP2K needs 186 MB — the same RAM signature as MnO, and the
  reason the diffuse cells are where this table bites.  (The SR2 row's 13.2× against the full span's 2.1× is
  not a typo and is not yet explained: the SMALLER basis is the worse ratio.  Both are Γ, same cell, same
  code path — worth one look, because whatever it is scales the wrong way.)
- **MnO is where both columns hurt: 6.0× (AFM) / 12.1× (FM) the CPU time, and 23× the RAM** — 4947 MB
  against CP2K's 217 MB, on an identical 118-function span.  CP2K caches nothing and its kernels are just
  fast; the RAM gap is our streams and Φ tables, which is why Steps 2–3 (fold, then screen Φ) are the
  answer and Step 4 is their readout, not a campaign of its own.
- **MnO carries a ~100 mHa configuration-BLIND offset and a −37.15 mHa configuration-SELECTIVE one** on that
  identical span (`OpenWork` Step 5) — and qchem sits BELOW a variational reference there, which convicts an
  operator or convention rather than a basis.  The run now prints its own term breakdown
  (`Ekin/Een/Eee/Exc/Enn/E_alphaZ`), which is where Step 5's term-by-term comparison starts.
- **The FREE production MnO run still folds NOTHING** — `[fold]` prints `NONE` at all three sites on a cell
  whose magnetic group has 12–24 ops (`OpenWork` Step 2).  Any MnO runtime row taken from a free run is
  measuring an unfolded run, and must say so.
