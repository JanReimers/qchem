# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**Step 1 of `doc/OpenWork.md`.  This is an INSTRUMENT, not a report.**  Steps 2–4 (folds, Φ-screening, RAM)
have no finish line without it: everything optimised in the 2026-08 runtime campaign was measured against a
single hand-run disabled test at `GPW_MNO_NMAX=3`, which is how two charter premises survived for months
while being wrong.  A row here is only meaningful with its PROVENANCE columns — same span held full-rank by
both codes, and the fold/threading/BLAS state that produced the number.

## ★★★ BENCHMARK PROTOCOL — the three rules a comparable row must satisfy (user, 2026-08-25)

A timing number is a claim about two codes doing the same work on the same hardware.  Three things break
that silently, and all three have already bitten this table, so they are rules rather than advice.

### 1. EVERY TIMING TABLE STATES ITS THREAD STATE — per table, for BOTH codes

Not in prose somewhere above it: in the caption or a column, so a row cannot be quoted out of context.
Today's tables carry a prose warning and it is not enough — the warning is right (CP2K is genuinely serial
here at 97–99% CPU; qchem is not, at 115–239%) but a reader lifting one row will not carry it.

⚠ **qchem has no `GLOBAL| Number of threads` banner** — this document has asked for one twice.  Until it
exists the thread state is asserted by the person taking the row, which is exactly the discipline that
fails.  **Build the banner** (index item 2): thread counts AND the feature flags below, printed by the run
itself, so every future row is self-describing.

### 2. SINGLE-THREAD PARITY FIRST, THEN THREADS — in that order (user)

The systematic procedure, and the reason for it: a threaded comparison against a serial code measures two
different things at once — the algorithm and the parallel efficiency — and if the single-thread gap is
unknown you cannot tell which one you are looking at.  So:

1. **Get qchem's SINGLE-THREAD time roughly in line with CP2K.**  `OMP_NUM_THREADS=1 GPW_OMP_THREADS=1`,
   which pins both our own OpenMP regions and the BLAS underneath blaze.  This is the honest
   algorithm-to-algorithm number and the one that should drive optimisation.
2. **THEN turn on N=8 or 16 and look for OMP-SHAPED GAPS** — poor scaling, barrier waits, regions that do
   not thread at all.  Those are a different class of defect with different fixes, and they only read
   cleanly once step 1 is settled.

⚠ Related and already measured: our OpenMP threads **BUSY-WAIT at the barrier**, so a threaded qchem run
bills far more CPU than it uses (a 294 s serial build billed ~590 s CPU at 16 threads).  That is why the
CPU column overstates qchem wherever it threads — and another reason the single-thread row is the clean one.

### 3. qchem-ONLY ACCELERATIONS ARE **OFF** FOR A HEAD-TO-HEAD ROW (user)

If qchem runs an algorithm CP2K does not, the row stops being a comparison and becomes an advertisement.
Turn them off, take the comparable number, and report the acceleration SEPARATELY as a qchem-vs-qchem
delta — which is the more useful statement anyway.

| flag | default | for a head-to-head row |
|---|---|---|
| `QCHEM_DM_LOWRANK` — the factored/low-rank \f$\rho\f$ (\f$D=LL^\dagger\f$, so \f$\rho_g=\lVert L^\dagger\Phi_g\rVert^2\f$)  | **ON** | **`=0`.**  CP2K collocates PAIRS; this is a singles/low-rank route it does not run.  ⚠ Every row taken since `07d13bf6` has it ON. |
| `GPW_XC_DM_SOURCE` — XC fed the retained DM ρ instead of the mixed ρ̃ | off | leave off; and see the pin — it is a TRAJECTORY change, so it must be pinned one way for Step 5's term-by-term breakdown either way |
| `GPW_OMP_THREADS` / `OMP_NUM_THREADS` | 1 / inherited | `=1` for the parity row (rule 2) |
| the T3 pair-stream orbit fold | **ARMED** | **DECLARE it, decide per row.**  Do not assume CP2K has no equivalent — state the fold state in the row and let the reader judge. |

**The general rule, so this table does not have to be exhaustive:** anything that makes qchem faster and is
not in CP2K's algorithm gets declared in the row and defaults to OFF in the comparison.  A new accelerator
is not finished until it is on this list.

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
| Si (FCC) | 2×2×2 shifted MP | SIPP_SR | −7.868473428 ¹ | −7.867436530 | **−1.04 mHa** | 14.3 s / 6.1 s | 15.6 / 6.0 s | 2.6× | 269 / 153 MB |
| NaF (rocksalt) | Γ | LOWQ_SR2 (both) | −24.430336482 | −24.431213375 | **+0.877 mHa** | 39.5 s / 7.4 s | 94.5 / 7.2 s | **13.2×** | 577 / 173 MB |
| NaF (rocksalt) | 2×2×2 Γ-centred | LOWQ_SR2 | −24.546883793 | — ² | | 57.4 s / — | 112 s / — | — | 590 / — MB |
| NaF (rocksalt) | Γ | LOWQ_SR (full) | −24.430944039 | −24.432293467 | **+1.349 mHa** | 2m44s / 1m42s | 219 / 102 s | 2.1× | 3090 / 186 MB |
| MnO AFM-II | Γ | **VA (N=118)** | −61.40297621 ⁴ | −61.303325178 | **−99.65 mHa** | 6m56s / 6m14s | 663 / 373 s | **1.8×** | **1323** / 217 MB |
| MnO FM | Γ | **VA (N=118)** | −61.441583060 ⁵ | −61.304782531 | **−136.80 mHa** | 21m45s / 3m13s | 2321 / 192 s | **12.1×** | 4947 / 217 MB |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | | ❓ |

Δ(AFM−FM) on the VA span: **qchem +38.61 mHa, CP2K +1.46 mHa** — both order FM first, and the
**configuration-SELECTIVE part of the offset is −37.15 mHa** (`OpenWork` Step 5).  Every one of these four
energies reproduces its banked value (runs 61/62 and the CP2K VA pair) to the digits those were recorded at.

¹ **SOLVED 2026-08-19 — and it was a SCREEN, not a phase.**  This test had rotted to −3.7351 while DISABLED
(it is the suite's ONLY fractional-k SCF coverage; every other k is TRIM, where the defect is structurally
invisible).  It is now ENABLED at ~14 s and reads −7.868473428, converged in 16 iterations at Δρ 1.0e-9.

**The bug** (`PG_Cart_MnD/Evaluator.C`, the D-aware integrate-back screen): the term was dropped when
`|Re(D_ij · conj(phase))| · maxv < eps` — a REAL PART used as if it were a magnitude.  `Re[D e^{-ikR}]` is
the right coefficient on the COLLOCATION side, where it multiplies a real pair product and a zero means a
genuinely zero contribution to ρ; the integrate-back's term is `phase·b`, whose size is `|b|` however the
phase is oriented.  **At a quarter-integer k, \f$e^{2\pi ikn}=i^n\f$ is purely imaginary for every ODD
offset**, so for real-ish D the screen discarded every odd-offset term and the Hartree/XC matrix came out
EXACTLY REAL — measured `maxIm(dV) = 0` at k=¼ against 0.067 at k=0.25001.  An H missing its imaginary part
has the wrong spectrum, so the SCF converged 2.5 Ha high.  Fix: screen on the true magnitude `|D_ij|`
(= `|D_ij conj(phase)|`, since `|phase|=1`) — strictly more conservative, and what the project's own
"the magnitude screen is the only truncation" rule always meant.

**How it was found**, since the sequence is the reusable part: a single-k sweep showed E(k) smooth except at
exactly ¼ and ¾; k=¼+1e-9 was fine, so it was an exact-value branch and not physics; every operator and the
1E spectrum were proven element-wise continuous (three gates, still enabled); the symptom was an extra
singlet with the Λ₃ doublet straddling E_F; a Fock-matrix fingerprint then showed the density-dependent
potential losing its imaginary part **only** at ¼; and `GPW_DENSITY_EPS=1e-30` recovered the right answer,
naming the screen.  Verification: E(k) is now smooth (0.249 → −7.563844, ¼ → −7.565529, 0.251 → −7.567208)
and **k=¾ equals k=¼ exactly**, as time reversal requires.  Si Γ, Si 2×2×2 Γ-centred and NaF Γ are unchanged
to every digit — at TRIM k the phase is real and the old test agreed with the new one.

² CP2K's NaF decks carry no `&KPOINTS` section, i.e. they are Γ.  The qchem test's own default is a **2×2×2**
mesh (8 k → 3 irreducible) — 116 mHa of band dispersion below its Γ value — so the row that compares to CP2K
is the `NAF_KMESH=1` one, and the 2×2×2 row is a qchem-only cost/size datapoint until a k-point CP2K deck
exists.  This mismatch was live in the test until 2026-08-19: it carried a Γ-era anchor (−24.4304) against a
2×2×2 configuration and simply failed.
³ the full-SR span runs FULL RANK at Γ ("kept 32 of 32", λ_min 4.35e-4) — the historical near-null trouble
was the multi-k cell.  The retracted −27.93128 "oracle" (a 3.5 Ha `EPS_PGF_ORB` screening artifact) is gone
from this table for good.

⁴ **RE-MEASURED TWICE, 2026-08-19/20** — same command, same box.  All three cuts agree on the ENERGY to
nine significant figures (−61.402976200 → −61.40297623 → −61.40297622, a spread of ~1e-8 Ha) with `m_stag`
±0.6667 over 17 iterations: both changes below are cost changes, not physics.

| cut | wall | CPU | peak RSS |
|---|---|---|---|
| banked (no folds, dense Φ) | 20m05s | 2240 s | 4947 MB |
| + **Step 2**: T3 stream fold ARMED (plan T3.5) | 13m25s | 1809 s | **1349 MB** |
| + **Step 3**: Φ-table build (screen + sparse spherical transform) | 8m36s | 1554 s | 1350 MB |
| + **Step 3**: Becke partition ε 1e-8 → 1e-6 | **6m56s** | **663 s** | 1323 MB |
| | **2.90×** | **3.38×** | **3.7×** |

**Step 2** (seconds, banked → armed): pair **scatter 263.4 → 41.0** (6.4×), pair **gather 167.6 → 23.1**
(7.3×), **stream build 110.1 → 28.2** (3.9×) — a larger factor than the 4.60× rep-pair reduction, because
the pairs the fold drops are the expensive ones.  **THE RAM CAME FROM HERE**, and it was the streams.
**Step 3**: the **Φ-table build 379.2 → 83.3 s** (4.6×; 190 → 36.8 per anneal stage).  Its dominant cause
was NOT Φ's density but a dense 122×118 cart→spherical transform applied INSIDE the lattice-image loop —
~150 mat-vecs per mesh point where one would do, 10.6% of all cycles — plus a missing magnitude screen on
the pointwise sweep.  Details and the rejected third hypothesis in `doc/OpenWork.md` Step 3.
**Becke ε** (2026-08-20): the partition's ε-converged competitor series ran at ε=1e-8, which had only ever
been probed TIGHTER.  ε fixes |im| (3183 competitor images per live point on this cell) and the partition
costs O(|P-set|·|im|) per point, so ε scales the dominant loop rather than shaving it: **1e-8 → 1e-6 takes
the becke build 36.95 → 8.33 s threaded (294 → ~66 s SERIAL), 4.44×**, at an energy identical to nine
significant figures (−61.40297622 → −61.40297621).  Margin: the binding gate is
`BeckeEquivalentSitesOwnEqualShares` (site shares equal to 1e-8 RELATIVE), which survives 1e-6 and fails at
1e-5, while the quadrature error itself is ~1e-4.  ⚠ This is a **TOLERANCE trade, not a bit-identical
restructuring** like the Φ work — the weights move at ~1e-6 relative.  `GPW_BECKE_EPS` overrides for A/B.
Against CP2K the CPU gap narrows 6.0× → 4.2× → **1.8×** and RAM 23× → **6.1×**; on WALL this row is 1.11×.

> **⚠ AND IT EXPOSED A FLAW IN THIS TABLE'S OWN CPU COLUMN.**  The becke build is 294 s SERIAL but bills
> ~590 s of CPU when threaded 16-way (36.95 s wall × 16) — because the OpenMP threads BUSY-WAIT at the
> barrier, so CPU time counts spinning as work.  Two anneal stages of that was ~1180 s of the banked row's
> 1554 s CPU (76%), which is why removing 4.44× of a "38%" bucket cut the row by 2.34×.  The user's pin
> "compare CPU, not wall" assumes CPU tracks work done; against a serial CP2K it does not, wherever qchem
> threads.  **The serial column is the honest algorithmic comparison** — which is exactly why the whole
> table was cut at one thread.  Anyone reading a threaded CPU number here should divide by the parallel
> efficiency, which still nothing measures.

**What is hot now:** the top bucket is `scf: XC-mesh ρ sampling (matrix-free)` at 84.1 s.
⚠ **THIS ROW PREVIOUSLY MIS-NAMED THAT BUCKET "the Φ-shaped GEMM, i.e. the Φ-SPARSITY item" — IT IS
NEITHER** (corrected 2026-08-20).  *Matrix-free* means exactly "carries no density matrix", so it cannot be
the DM GEMM.  Per `PWTerms.C:698`, it is the **ρ̃-MIXED density sampled on the XC mesh — a batched inverse
FT over the whole {G}, on EVERY Kerker/Pulay iteration** (plus the iteration-0 seed).  The real DM GEMM is
the separate `scf: XC-mesh ρ sampling (all iterations)` line, and on this recipe the mixer hands XC a
ρ̃-backed density from iteration 1 on, so the GEMM is nearly BYPASSED: measured 1.70 s against the
matrix-free bucket's 35.0 s on a 6-iteration cut.  The code had already split these two buckets for this
exact reason — *"lumping it into the GEMM hid the fact that the mixed-density sampling, not the GEMM, was
the iteration's largest XC cost"* — and this table re-lumped them in prose.  **The per-iteration lever is
the ρ̃-mixed sampling, not anything Φ- or D-shaped.**
⁵ the FM row is still the PRE-Step-2 measurement (its energy is unaffected, its cost is not) — re-run it
with the AFM command when the threaded cut of this whole table is taken.  **The Si and NaF rows likewise
predate Step 3**: their ENERGIES are unchanged (verified — Si Γ still reads −7.115067665 to every digit),
but any row whose XC runs on the Becke mesh may now be faster than its cost columns say.  Only the MnO AFM
row has been re-measured end to end.  Nothing in the table is stale in the Δ column, which is the column it
exists for.

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
| collocation streams (T3 pairs) | **4.60×** (12 ops, 8778 → 1909 rep pairs) ⁴ | NONE |

So the imposed row now claims all three folds.  The stream fold is **12 ops, not 24** — the S3 pin: a
magnetic imposition may fold the PER-CHANNEL streams only under the σ=None (sublattice-preserving) subgroup,
since a flip op relates D↑ to D↓.  And 4.60× on a 4-atom cell is not the 71× the diamond gate cell showed:
the orbit factor is a property of the cell's symmetry, so the headline number belongs to a cell, never to
the feature.  **The FREE production run still folds NOTHING at any of the three sites** — that is a
`MNO_IMPOSE` decision, not a plumbing gap.

Where the time goes (AFM arm, `GPW_REPORT=1` ledger) — BEFORE the fold was armed (of 1205 s) and AFTER
(of 805 s):

| bucket | s (banked) | s (+Step 2 fold) | s (+Step 3 Φ) |
|---|---|---|---|
| setup: XC-mesh **Φ tables** | 370.0 (31%) | 379.2 (47%) | **83.3 (16%)** |
| scf: collocate density (pair scatter) | 263.4 (22%) | 41.0 (5%) | 42.1 (8%) |
| scf: integrate-back (pair gather) | 167.6 (14%) | 23.1 (3%) | 23.3 (5%) |
| setup: collocation stream build | 110.1 (9%) | 28.2 (4%) | 29.0 (6%) |
| setup: Becke mesh build | 70.6 (6%) | 69.4 (9%) | 71.5 (14%) |
| **scf: XC-mesh ρ sampling (matrix-free)** | 57.1 (5%) | 76.8 (10%) | **84.1 (16%)** |
| scf: XC-mesh H_xc quadrature | 33.8 (3%) | 44.1 (5%) | 28.3 (5%) |

**The profile has flattened and the head has moved.**  The pair loops (Step 2) and the Φ BUILD (Step 3) are
each down to single-digit-to-16% shares, and the largest bucket is now the **Φ-shaped ρ GEMM** — the
Φ-SPARSITY item, which is what `OpenWork` Step 3's "Φ-table screening" was always aimed at and which the
build's cost used to hide.  Second is the **Becke mesh build**, whose 71 s of wall is ~50% of the run's
CYCLES (it is the one loop threaded by default) — so on the CPU column, which is the honest one, the Becke
partition is now the single biggest item in the code.  ⁶ the XC-mesh buckets are not comparable one-to-one
across cuts (iteration counts differ per cut); per-iteration cost is what matters there.
⁶ the XC-mesh buckets are not comparable one-to-one across the two runs (the folded run converged its second
stage in 17 iterations; per-iteration cost is what matters there, and it did not change).

## Standing observations

- **Si Γ agrees to 1.1e-5 Ha** at 1.5× CP2K's CPU time and 1.8× its RAM.  On the small cells the gap is
  modest; the diffuse and magnetic cells are where it opens up.
- **NaF Γ agrees to 0.9 mHa (SR2) / 1.3 mHa (full SR)** at 13.2× / 2.1× the CPU time and 3.3× / 16.6× the
  RAM.  The full span costs us 3090 MB where CP2K needs 186 MB — the same RAM signature as MnO, and the
  reason the diffuse cells are where this table bites.  (The SR2 row's 13.2× against the full span's 2.1× is
  not a typo and is not yet explained: the SMALLER basis is the worse ratio.  Both are Γ, same cell, same
  code path — worth one look, because whatever it is scales the wrong way.)
- **MnO is where both columns hurt — but Steps 2 and 3 took a large bite out of both (2026-08-19/20).**  The
  AFM row went 6.0× → **4.2×** CPU, 23× → **6.2×** RAM (4947 → 1350 MB against CP2K's 217 MB) and 3.2× →
  **1.38×** wall on an identical 118-function span; the FM row still carries the pre-Step-2 numbers.  **The
  RAM half was the streams** (Step 2), so Step 4 was answered by Step 2 rather than by a campaign of its
  own, exactly as predicted.  CP2K caches nothing and its kernels are just fast — what is left on our side
  is the Φ-shaped ρ GEMM and the Becke partition, in that order for wall and the reverse for CPU.
- **MnO carries a ~100 mHa configuration-BLIND offset and a −37.15 mHa configuration-SELECTIVE one** on that
  identical span (`OpenWork` Step 5) — and qchem sits BELOW a variational reference there, which convicts an
  operator or convention rather than a basis.  The run now prints its own term breakdown
  (`Ekin/Een/Eee/Exc/Enn/E_alphaZ`), which is where Step 5's term-by-term comparison starts.
- **The FREE production MnO run still folds NOTHING** — `[fold]` prints `NONE` at all three sites on a cell
  whose magnetic group has 12–24 ops (`OpenWork` Step 2).  Any MnO runtime row taken from a free run is
  measuring an unfolded run, and must say so.  Arming T3 does not change this: it is a `MNO_IMPOSE`
  decision, i.e. physics, and the 1.5× wall / 3.7× RAM measured above is what a free run is paying for the
  freedom.
