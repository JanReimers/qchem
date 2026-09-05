# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**This is an INSTRUMENT, not a report.**  A row here is a claim that two codes did the same work on the
same hardware; the sections below are what makes that claim checkable.
📖 **The reasoning, the wrong turns and every superseded number are in `doc/BenchmarkHistory.md`** (split
out 2026-09-04).  Come back here for what is TRUE NOW; go there for why.

**Read in this order:** §1 the process → §2 what parity actually means → §3 the rules → §4 the commands →
§5 the rows.

---

## 1. THE SYSTEMATIC PROCESS — single thread first, then threads; and the four bins

### 1a. Single-thread parity FIRST, then threads — in that order (user)

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


### 1b. The four gap-close bins (user, 2026-08-28) — priority order

> *"1) Per iteration CPU time, 2) Init (pre iteration) time, 3) top RAM usage (mostly solved I think),
> 4) very roughly match the # of iterations.  I think items that fall into bin 4 can just be documented."*

| bin | the axis | where it stands (MnO AFM-II VA, 2026-09-04) |
|---|---|---|
| **1** | **per-iteration CPU** | MnO: ALL DEFAULTS **1.47×** CP2K, `QCHEM_BECKE_XC=0` **0.47×** (we WIN), `CP2K_COMPAT=1` **1.99×**.  Si: **0.15×** at Γ but **1.01×** at 8 k — the k-scaling gap (§5a) is where the small-cell story is.  ⇒ **NOT closed** |
| **2** | **init / pre-iteration** | ✅ **~1% of wall** on the parity routes.  On the DEFAULT route it is still the Becke mesh + Φ tables (~57 s of 328 s), so bin 2 is a Becke-mesh question and only there |
| **3** | **peak RAM** | ✅ **solved, and we win**: ~470 MB defaults, **~105 MB on the parity routes against CP2K's 217 MB** |
| **4** | **iteration count** | 31 (defaults) / capped (parity) against CP2K's 44 — ⇒ DOCUMENT, do not chase.  The two codes do not run the same ρ-mixing algorithm (doc/OpenWork.md) |

The live tracker for these is `doc/OpenWork.md`; this file holds the measurements behind them.

---

## 2. WHAT `CP2K_COMPAT=1` ENCOMPASSES — ⚠ AN EMERGING LIST, NOT A FINISHED ONE

**Every deviation found so far is an ACCELERATION, not physics**: turning them all off moves the MnO total
by **3e-8 Ha** (agreeing to 10 s.f.).  That is the property the switch most needed to demonstrate about
itself, and it is measured rather than asserted (history §2).

⚠ **THIS LIST HAS GROWN EVERY TIME SOMEONE LOOKED.**  It started at four items; it is seven.  Assume it is
still incomplete — a parity row is only as honest as the last thing we noticed we were doing differently.

| # | knob | what qchem does that CP2K does not | found |
|---|---|---|---|
| 1 | `QCHEM_DM_LOWRANK` | factored/low-rank ρ (\f$D=LL^\dagger\f$) — a SINGLES route; CP2K collocates PAIRS | 08-25 |
| 2 | `GPW_STREAM_FOLD` | orbit fold on the collocation pair streams (5.2× on MnO's pair count) | 08-25 |
| 3 | `QCHEM_MIX_RHO_M` | (ρ,m) mixing channels instead of (up,dn) | 08-25 |
| 4 | `GPW_XC_DM_SOURCE` | \f$V_{xc}\f$ fed ρ[D] wholesale instead of ρ_mix | 08-25 |
| 5 | `QCHEM_IMPOSE_SYMMETRY` | space-group imposition: BZ fold + ρ star-average + site-adapted XC mesh.  ⚠ **OVERRULES the caller** | 08-26 |
| 6 | `QCHEM_BECKE_XC` | atom-centred (Becke) XC quadrature instead of the uniform grid.  ⚠ **OVERRULES the caller**; it was **43% of the MnO row** | 08-28 |
| 7 | `GPW_DAWARE_SCREEN` | D-aware collocation box tolerance \f$\varepsilon/|c_{ij}|\f$ instead of flat \f$\varepsilon\f$ | 09-04 |

**Verified locally** that CP2K does none of the symmetry work: the 1129-line `bench_MnO_AFM2_VA_cp2k.log`
contains **zero** occurrences of "irrep", "symmetry" or "point group"; QuickStep keeps K and P as DBCSR
sparse ATOM-BLOCK matrices and blocks by atom-pair SPARSITY, not by irrep.  The one symmetry knob that
exists is BZ-side and defaults OFF (`BRILLOUIN| K-Point point group symmetrization  OFF`).

⇒ **These rows compare a SYMMETRY-exploiting code against a SPARSITY-exploiting one on a small,
high-symmetry cell — the regime that most favours us.**  A 100-water box would invert it; CP2K's design
centre is large disordered systems where every group has order 1.

**THE OVERRIDE RULE**: an explicitly-set knob WINS over `CP2K_COMPAT=1`, and the banner marks it `(stated)`
so a run that thinks it has parity and does not says so out loud.  ⇒ **A new accelerator is NOT FINISHED
until it appears on that deviation line** (`qchem::RunPolicy`).

---

## 3. THE RULES a comparable row must satisfy

### 3a. COPY the command from §4 — never reconstruct it

⛔ The MnO rows need **eight** env vars.  Drop them and you silently get a different basis
(`VALENCE_LOWQ_SR`, Cartesian, not VA/spherical) and a single SCF stage instead of the two-stage anneal —
a different system, on which CP2K's 8.5 s/step means nothing.  This cost a full retraction on 2026-09-04
(history §1).
★ **THE TELL**: that probe is documented as *"~6 minutes"*; the bad runs took **90 seconds**.  **A large
discrepancy against the DOCUMENTED cost of the SAME probe is a CONFIGURATION difference until proven
otherwise — it is not your speedup.**
✅ **The check that catches it in one run**: a correct recipe reproduces the banked \f$E_{tot}\f$ AND
iteration count to all digits.  If it does not, you are not running the banked system.

### 3b. EVERY timing table states its thread state — per table, for BOTH codes

Not in prose somewhere above it: in the caption or a column, so a row cannot be quoted out of context.
Today's tables carry a prose warning and it is not enough — the warning is right (CP2K is genuinely serial
here at 97–99% CPU; qchem is not, at 115–239%) but a reader lifting one row will not carry it.

✅ **THE BANNER EXISTS AS OF 2026-08-26** (doc/OpenWork.md T5/N5).  `qchem::SolidCalculation` prints it
UNCONDITIONALLY — no verbose flag, no opt-in — four lines at construction plus one per `Converge`:

```
[<label> run] system: 4 atoms, 26 valence e, multiplicity 1 (POLARIZED), seed=IonicSAD
[<label> run] grids: densityEcut=auto C=2 raster=BallOnly xcMesh=Becke (nR=40 L=29)
[<label> run] symmetry: IMPOSED (Shubnikov from the decoration);  threads: OMP_NUM_THREADS=1 GPW_OMP_THREADS=1 (BLAS pinned to 1)
[<label> run] CP2K_COMPAT=0 -> DEVIATING;  QCHEM_DM_LOWRANK=on*  GPW_STREAM_FOLD=on*  QCHEM_MIX_RHO_M=off  GPW_XC_DM_SOURCE=off   [* = differs from CP2K]
[<label> scf] mixer: Kerker(G0=1.000000) alpha=0.45;  XC rho source: rho_mix;  accel: Ladder;  kT=0.005 MOM=on NMaxIter=80
```

**So a row is now taken by COPYING those lines beside it, not by remembering what was set.**  The
deviation line is generated from the SAME policy object the factories consult (`qchem.RunPolicy`), so it
cannot drift from what was actually built, and `CP2K_COMPAT=1` turns the whole set off in one place.
⚠ The standing rule is unchanged and now has somewhere to point: **a new accelerator is not finished until
it appears on that deviation line.**


### 3c. qchem-ONLY accelerations are OFF for a head-to-head row (user)

If qchem runs an algorithm CP2K does not, the row stops being a comparison and becomes an advertisement.
Turn them off, take the comparable number, and report the acceleration SEPARATELY as a qchem-vs-qchem
delta — which is the more useful statement anyway.

| flag | default | for a head-to-head row |

The full list is §2; the mechanism is `CP2K_COMPAT=1`.  ⇒ Report the acceleration SEPARATELY as a
qchem-vs-qchem delta — which is the more useful statement anyway (§6).

### 3d. The `CPU` column is a WHOLE-RUN TOTAL — divide by the iterations

Two runs in the table differ by 3× in iteration count, so a run that needs twice the steps looks twice as
slow on a per-ITERATION basis it may actually be winning.  ⇒ Quote **CPU/iteration** for bin 1.
⚠ And per-iteration is not comparable across different ITERATION CAPS either (the stage mix changes the
GDM line-search); the cap-independent measure is **per CALL**, straight off the ledger.

### 3e. Read gather MISSES, not closure calls

The ledger prints both and **only misses are work**.  A cost estimate taken off call counts is wrong
whenever a memo sits underneath — that is how a 2026-09-04 fix was over-estimated **70×** (history §1b).

---

## 4. HOW TO PRODUCE A ROW — the same wrapper for both codes

```bash
scripts/bench "Si Gamma qchem" -- build/Release/IntegrationTests/ITMain --gtest_filter=GPW_SCF.SiliconGammaConverges
scripts/bench "Si Gamma cp2k"  -- cp2k -i IntegrationTests/CP2K/si_fcc_gpw.inp
```

**Peak RAM is measured from OUTSIDE, identically for both codes** (`/usr/bin/time -v` → "Maximum resident
set size" = the kernel's `VmHWM` high-water mark).  That is the whole reason the RAM column is now
obtainable without touching CP2K: no in-program hook was ever needed, and measuring both sides through one
wrapper is what makes the two columns comparable.  On a qchem row the wrapper also prints qchem's OWN
internal `VmHWM` as a cross-check — they should agree (measured: 266 vs 266.7 MB on Si Γ).

For the qchem detail (per-bucket ledger, folds, task-list geometry) add `GPW_REPORT=1`:

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
| `[collocation] kernel=… ;  task list: … tasks, … MB` | the kernel + task-list readout (unconditional) |
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

# Si k-MESH ROWS -- ⚠ THESE WERE MISSING FROM THIS BLOCK UNTIL 2026-09-04, i.e. two table rows had no
# printed recipe at all (rule 3a's own failure mode).  Recovered by matching Etot to the banked value.
GPW_REPORT=1 scripts/bench "Si 222g qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.DISABLED_SR_2x2x2GammaCentred_vs_CP2K --gtest_also_run_disabled_tests
GPW_REPORT=1 scripts/bench "Si 222shift qchem" -- build/Release/IntegrationTests/ITMain \
    --gtest_filter=GPW_SCF.SR_2x2x2ShiftedMP_vs_CP2K

# CP2K -- from the DECK'S directory (relative BASIS_SET_FILE_NAME) and at one thread
cd IntegrationTests/CP2K
OMP_NUM_THREADS=1 ../../scripts/bench "Si Gamma cp2k"    -- mpirun -np 1 cp2k.psmp -i si_fcc_gpw.inp
OMP_NUM_THREADS=1 ../../scripts/bench "MnO AFM2 VA cp2k" -- mpirun -np 1 cp2k.psmp -i mno_afm2_gpw_va.inp
```

Verify each side reproduces its own history before reading a Δ: the CP2K decks against `doc/CP2Kresults.md`
(all five re-validated 2026-08-19, `doc/CP2KBuild.md`), and the qchem runs against the tests' own anchors.


---

## 5. THE ROWS — single thread

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

**★ RE-TAKEN 2026-08-27/28**, most recently on the spin-native XC pair route.  ⚠ **THE ROWS ARE NOT ALL
FROM THE SAME BINARY**, and this session moved fast enough that saying so matters more than usual — five
changes landed in two days (cache deletion, `template<int LP>`, the two memo fixes, the exp recurrence, the
XC decoupling), and a row is only as current as the last time it was run.  The `taken` column says when.

| row | last taken | why it may have moved since |
|---|---|---|
| Si Γ, NaF SR2 Γ, NaF full-SR Γ, MnO AFM-II, **MnO parity** | **08-28, current** | — |
| Si 2×2×2 (both), NaF SR2 2×2×2 | 08-28, mid-session | predate the exp recurrence (~1.4× on the box walk); ⚠ CHEAP to re-take (8–68 s) |
| **MnO `CP2K_COMPAT=1`** | **08-28, current** | the old row's *"does not converge"* verdict is RETRACTED — see ⁷ |
| **MnO FM** | **08-19** | predates EVERYTHING; ~20 min to re-take |

⇒ **The `CP2K_COMPAT=1` row is now real rather than aspirational**, because the XC pair route made the
parity ROUTE affordable — see footnote ⁷ (§5d).  CP2K column untouched throughout.

| system | k-mesh | span | qchem Etot | CP2K Etot | Δ (qchem−CP2K) | wall q / c | **CPU q / c** | **CPU ×** | peak RSS q / c |
|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | −7.115067844 | −7.115057882 | **−10.0 µHa** | **0.8 s** / 5.2 s | **2.2** / 5.0 s | **0.44×** | **29** / 148 MB |
| Si (FCC) | 2×2×2 Γ-centred | SIPP_SR | −7.778472833 | −7.778457865 | **−15.0 µHa** | 7.5 s / 5.8 s | 8.9 / 5.6 s | 1.6× | **30** / 153 MB |
| Si (FCC) | 2×2×2 shifted MP | SIPP_SR | −7.868473428 ¹ | −7.867436530 | **−1.04 mHa** | 16.2 s / 6.1 s | 17.5 / 6.0 s | 2.9× | **31** / 153 MB |
| NaF (rocksalt) | Γ | LOWQ_SR2 (both) | −24.4303364755 | −24.431213375 | **+0.877 mHa** | 21.1 s / 7.4 s | 37.9 / 7.2 s | **5.3×** | **61** / 173 MB |
| NaF (rocksalt) | 2×2×2 Γ-centred | LOWQ_SR2 | −24.5468834873 | — ² | | 1m07.8s / — | 84.5 s / — | — | **68** / — MB |
| NaF (rocksalt) | Γ | LOWQ_SR (full) | −24.4309472653 | −24.432293467 | **+1.346 mHa** | **25.0 s** / 1m42s | **41.4** / 102 s | **0.41×** | **65** / 186 MB |
| **MnO AFM-II — ALL DEFAULTS** ⁴ | Γ | **VA (N=118)** | −61.40297551 ⁴ | −61.303325178 | **−99.65 mHa** | **5m16.6s** / 6m14s | **563** / 373 s | **1.51×** | **490** / 217 MB |
| **MnO AFM-II, `QCHEM_BECKE_XC=0`** ⁶ | Γ | **VA (N=118)** | −61.40358773 | −61.303325178 | −100.26 mHa | **4m05.4s** / 6m14s | **246** / 373 s | **0.66×** | **105** / 217 MB |
| ⚠ STALE MnO FM | Γ | **VA (N=118)** | −61.441583060 ⁵ | −61.304782531 | **−136.80 mHa** | 21m45s / 3m13s | 2321 / 192 s | **12.1×** | 4947 / 217 MB |
| **MnO AFM-II, `CP2K_COMPAT=1`** ⁷ | Γ | **VA (N=118)** | −61.39789688 ⁷ | −61.303325178 | −94.57 mHa | 45m40s / 6m14s | **2736** / 373 s | **7.3×** | **112** / 217 MB |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | | ❓ |


### 5a. Per-ITERATION CPU — the bin-1 ladder

★★★ **THE SMALL ROWS, DECOMPOSED (2026-09-04).**  Splitting `CPU×` into (CPU/iteration) × (iteration count)
dissolves most of the column's outliers.  All four re-taken with NATIVE codegen; **every one reproduces its
banked \f$E_{tot}\f$ exactly**, which is what says they are the same systems:

| row | CPU banked | **CPU now** | qchem iters | CP2K iters | **qchem s/iter** | CP2K s/iter | **per-iter ×** | old whole-run × |
|---|---|---|---|---|---|---|---|---|
| Si Γ | 2.2 s | **0.69 s** | 11 | 12 | **0.063** | 0.417 | **0.15×** | 0.44× |
| Si 2×2×2 Γ-centred | 8.9 s | **3.04 s** | 7 | 13 | **0.434** | 0.431 | **1.01×** | 1.6× |
| Si 2×2×2 shifted MP | 17.5 s | **6.57 s** | 16 | 14 | **0.411** | 0.429 | **0.96×** | 2.9× |
| NaF SR2 Γ | 37.9 s | **26.3 s** | 29 | — | 0.906 | — | — | 5.3× |

⇒ **PER ITERATION WE ARE AT OR BETTER THAN CP2K ON ALL THREE Si ROWS.**  The banked 1.6× and 2.9× were two
artifacts stacked: rows stale by 2.7–3.2×, plus iteration-count differences (bin 4).

⛔ **BUT THE k-SCALING IS THE REAL BIN-1 STORY HERE.**  Same system, same span, only the k-mesh changing:

| | qchem s/iter | CP2K s/iter |
|---|---|---|
| Γ | 0.063 | 0.417 |
| 8 k-points | 0.434 | 0.431 |
| **rise** | **6.9×** | **1.03×** |

We start **6.7× ahead** of CP2K at Γ and spend the entire lead on k-points.  Cause and the two refuted
fixes: doc/OpenWork.md's k-scaling entry (the gather memo is bypassed whenever a density screen is passed,
and 32 of 51 gathers on the 8-k run are the SAME FIELD).

★ **AND NaF's 5.3× IS A BECKE COST, NOT A GPW ONE.**  Its ledger: Becke mesh build **6.96 s (27%)** +
XC-mesh ρ sampling **8.16 s (31%)** + hamiltonian ctor 2.0 s (8%) — against a box walk of **~1 s (4%)**.
⇒ On the system where Becke should be at its best, Becke IS the runtime.  That makes §6's open
\f$\lVert V_{xc}-V_{xc}^{fit}\rVert\f$ study the deciding measurement for this row too.



★ **RE-TAKEN 2026-09-04 WITH THE PRINTED COMMAND** (serial; `-O3` and `-O3 -march=native` arms).  Every
row reproduces its banked \f$E_{tot}\f$ and iteration count EXACTLY -- which is the correctness statement
for the whole 09-04 change set (screener seam + one-gather XC term + CacheSpin + the codegen flip):

| the row | iterations | banked 08-28 | **O3 now** | **NATIVE now** | \f$E_{tot}\f$ |
|---|---|---|---|---|---|
| **ALL DEFAULTS** | 14+17 = **31** | 18.2 s (2.14×) | **13.78 s (1.62×)** | **12.45 s (1.47×)** | −61.40297551 = banked |
| **`QCHEM_BECKE_XC=0`** | 14+25 = **39** | 6.31 s (0.74×) ⁸ᵇ | **4.36 s (0.51×)** | **3.99 s (0.47×)** | −61.40358773 = banked |
| **`CP2K_COMPAT=1`**, `GPW_MNO_NMAX=10` probe | 10+10 = **20** | 19.3 s (2.3×) | **17.89 s (2.11×)** | **16.94 s (1.99×)** | (capped, both arms equal) |
| ⚠ `CP2K_COMPAT=1`, the full table row | 13+80 = **93**, both CAPPED | 29.4 s ⁸ | not re-taken (45 min) | — | |
| CP2K (its own log) | 44 | 8.5 s | — | — | |

⁸ᵇ the banked `QCHEM_BECKE_XC=0` row recorded 246 s CPU but only stage 2's iteration count; stage 1 is
assumed 14 (this re-take's value, and its \f$E_{tot}\f$ matches to all digits), so 6.31 s/iter is derived,
not banked.

⇒ **WHAT MOVED, honestly**: ALL DEFAULTS **−24.3%** CPU against 08-28, `BECKE_XC=0` **−30.9%**, the parity
probe only **−7.3%**.  ⚠ THIS IS "CURRENT vs BANKED", NOT AN ATTRIBUTION: other work landed between 08-28
and now, and no A/B against the parent commit was run on this recipe, so it is NOT established how much of
it today's three changes account for.
⇒ **BIN 1 IS NOT CLOSED.**  Best case 1.47× (defaults, NATIVE); parity still 1.99×.  The `BECKE_XC=0` row
is the standout -- **0.47×, i.e. we are ~2× FASTER than CP2K per step on the XC-parity route.**


★ **AND THE BOX WALK IS ~95% OF THIS RUN** (109 s of 114.6 s), of which the GATHER is 74% (81.5 s over 65
calls, against 27.9 s over 22).  ⇒ Amdahl leaves nothing outside the walk worth touching, and any per-call
win is worth ~3× more on the gather side than on the collocate side.

★ **STANDING PROBE for bin 1** (user, 2026-08-28: *"we just cut off at ~10 or so iterations, just to get a
decent average"*): `CP2K_COMPAT=1 GPW_MNO_NMAX=10` — ~6 minutes, and quote per-call beside per-iteration.

⇒ Two separate facts the whole-run column blurs together: we are **1.47× per step on the default route and
1.99× at parity** (09-04, VA recipe; it was 2.14×/3.5× when this line was written), and at parity we take **93 steps to CP2K's 44** (bin 4, and both of
ours hit the cap rather than converging).  The parity row's per-step figure is the WORSE of the two because
parity also removes the stream fold — 5.2× on MnO's pair count.


⚠ **AND THE PER-STEP COMPARISON MAY NOT BE APPLES-TO-APPLES AT ALL — OPEN, 2026-09-04.**  The 09-04 ledger
shows **~9 distinct KS-field integrations per SCF iteration** on a 2-channel system where the physics needs
~3 (one \f$V_H\f$ gather + one \f$V_{xc}\f$ per channel); `GPW_INTEGRATE_CENSUS=1` says all of them are
genuinely NEW fields, so it is NOT redundancy a memo could remove.  Two readings, and they call for
opposite responses:
- **the accelerator is buying bin 4 with bin 1 work** — if the Ladder line-search builds \f$H\f$ several
  times per "iteration", then our iteration is not CP2K's step and per-iteration is the WRONG metric; the
  honest one is total \f$H\f$-builds (or wall) to convergence;
- **or the term assembly asks more often than it needs to**, in which case it is a ~3× on the dominant
  bucket and would put the run BELOW CP2K.

⇒ **Distinguish them before optimising either way**: count \f$H\f$-builds per SCF step directly, and read
CP2K's own per-step \f$H\f$ count out of its log.  Cheap, and it decides whether bin 1 is finished.


### 5d. Footnotes to the table

Compact here; the full stories are in `doc/BenchmarkHistory.md` at the section named after each.

- **¹** Si 2×2×2 shifted MP (−1.04 mHa) — the residual after a **D-aware integrate-back SCREEN defect** was
  fixed 2026-08-19.  It is the suite's ONLY fractional-k SCF coverage (every other k is TRIM, where the
  defect is structurally invisible), which is why it had rotted to −3.7351 while DISABLED.  *(history: the
  screen was reading \f$\mathrm{Re}[D\overline{e^{ikR}}]\f$ where it needed \f$|D|\f$.)*
- **²** No CP2K counterpart: CP2K's NaF decks carry no `&KPOINTS` section, i.e. they are Γ, while the qchem
  test defaults to 2×2×2.
- **⁴** MnO ALL DEFAULTS — **re-measured twice** (08-19/20), same command, same box; all cuts agree on the
  energy to the printed digits, and re-taken again 09-04 with the same \f$E_{tot}\f$.
- **⁵** ⚠ The MnO **FM** row is the only PRE-Step-2 measurement left in the table: its ENERGY is unaffected
  but its COST is stale by everything since 08-19.  ~20 min to re-take.
- **⁶** The XC-parity row — **the first MnO row on which qchem beat CP2K on BOTH axes**, and still the
  standout (0.47× per step, ~105 MB against 217 MB).
- **⁷** The parity row: an earlier *"at true parity our recipe does not converge"* verdict was **RETRACTED**
  (08-28).  It still hits the iteration cap, but for a far more benign reason — 5.1 mHa short, not 3.8 Ha,
  with the AFM order surviving both stages.
- **⁸** The 93-iteration parity figure is pre-fix and its stage mix differs — **do not** read 29.4 → 19.3 as
  a 1.5×; per-iteration is not comparable across iteration caps (rule 3d).

### 5b. The three MnO rows are ONE system and ONE recipe

★ **THE THREE `MnO AFM-II` ROWS ARE ONE SYSTEM AND ONE RECIPE, with progressively more of OUR deviations
switched off.**  Read them as a ladder, not as three experiments — the only thing changing is which of the
six declared deviations are on:

| the row | what is off | CPU | RSS |
|---|---|---|---|
| **ALL DEFAULTS** | nothing — every deviation at its qchem default (Becke XC, imposition, low-rank ρ, stream fold) | 584 s | 491 MB |
| **`QCHEM_BECKE_XC=0`** | the Becke XC mesh only — so XC runs CP2K's way, everything else still ours | **246 s** | **105 MB** |
| **`CP2K_COMPAT=1`** | ALL SIX — the honest algorithm-to-algorithm row | 2736 s | 112 MB |

⚠ The middle row is the fastest because Becke's atom-centred quadrature is expensive machinery, and the
bottom row is the slowest because parity also removes the STREAM FOLD (5.2× on MnO's pair count) and the
low-rank ρ.  ⇒ **Always say WHICH of the three** — "the MnO row" has meant all three of these in the space
of one day, and they differ by 11× in CPU.


### 5c. Like-for-like: what "same span" costs


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


---

## 6. qchem ACCELERATIONS NOT IN CP2K — the qchem-vs-qchem deltas

These are what §2's knobs BUY.  They are excluded from a head-to-head row by rule 3c and reported here
instead, which is the more useful statement.

| acceleration | measured worth | provenance |
|---|---|---|
| Becke XC mesh (`QCHEM_BECKE_XC`) | ⚠ **NEGATIVE on MnO**: the default row is 2.8× the CPU of the `BECKE_XC=0` row (12.45 vs 3.99 s/iter, NATIVE).  It buys atom-centred accuracy, not speed | §5a |
| symmetry imposition | buys CONVERGENCE, not accuracy and not the magnetic basin — the AFM order survives a free run | history §2 |
| stream fold (`GPW_STREAM_FOLD`) | 5.2× on MnO's pair COUNT | history §7 |
| D-aware screen (`GPW_DAWARE_SCREEN`) | +14.5% wall on the collocation route | doc/ScreeningPlan.md §6 |
| collocation memo depth 5 | 3.74× on its bucket (60 → 16 collocations) | history §4 |
| gather memo on (V, screen) | catches the exact duplicates; the remainder are genuinely distinct fields | history §5 |
| `-march=native` (now the Release DEFAULT) | −9.6% / −8.4% / −5.3% CPU on the three MnO rows, \f$E_{tot}\f$ bit-identical | CMakeLists.txt |

⚠ **NOT an acceleration and NOT a deviation**: `GPW_CONTRACT_CUBE` (the separable-contraction collocation
kernel).  CP2K collocates exactly this way (`grid_cpu_collint.h`, Mathieu's three 2-D tables).  It is on
this page only so a row can state which kernel produced it — every run prints `[collocation] kernel=…`.

---

## 6b. ACCELERATIONS **CP2K** HAS THAT **WE** DO NOT — ⚠ ALSO AN EMERGING LIST

§2 and §6 are one direction; this is the other, and it had no home until 2026-09-04.  A row is only fair
if BOTH lists are known, and this one started at zero because nobody had looked.

| # | what CP2K does | what it costs us | found |
|---|---|---|---|
| 1 | **TIME-REVERSAL k FOLDING** — a Monkhorst-Pack mesh is folded \f$8\to4\f$ (`BRILLOUIN\| List of Kpoints ... 4`, weights 0.25, with `K-Point point group symmetrization OFF`).  We run all 8. | up to **2×** on any non-TRIM mesh.  ⚠ It folds NOTHING on a Γ-centred 2×2×2, where every k is its own inverse — which is why the two Si k-rows behave differently | 09-04 |
| 2 | **k-INDEPENDENT per-step cost** — their per-iteration time is flat from Γ to 8 k (0.417 → 0.431 s) where ours rises 6.9×.  Not a "feature" so much as a consequence of collocating the k-summed density once per step | the whole Γ lead (§5a) | 09-04 |

⚠ Item 1 is a REAL algorithmic advantage they hold, not a deviation to switch off; item 2 is a gap of ours
with a diagnosed cause.  ⇒ **Neither belongs in `CP2K_COMPAT`** — that switch turns OUR accelerations off,
and turning theirs off is not available to us.

⚠ CP2K also runs `Wavefunction type COMPLEX` on BOTH Si k-rows, including the all-TRIM Γ-centred one where
we use REAL blocks.  So on that row we hold an advantage they decline to take — and are still only at
parity per iteration.

---

## 7. THE ROWS — threaded (OMP)

⛔ **EMPTY BY DESIGN, AND THAT IS THE POINT.**  Rule 1a says single-thread parity comes FIRST: a threaded
comparison against a serial code measures the algorithm and the parallel efficiency at once, and if the
single-thread gap is unknown you cannot tell which you are looking at.  Bin 1 is not closed at parity
(1.99×), so this table stays empty on purpose.

What it will need when bin 1 is settled:

| system | threads q / c | wall q / c | CPU q / c | **speedup vs own serial** | parallel efficiency |
|---|---|---|---|---|---|
| — | — | — | — | — | — |

⚠ **KNOWN BEFORE WE START**: our OpenMP threads **BUSY-WAIT at the barrier**, so a threaded qchem run bills
far more CPU than it uses (a 294 s serial build billed ~590 s CPU at 16 threads).  ⇒ On this table the
honest column is **speedup against our OWN serial time**, not CPU; and CP2K must be given the same core
count, not left at `OMP_NUM_THREADS=1`.

---

## 8. WHAT IS STILL MISSING

★ **AS OF 2026-09-04**, in rough priority order:
- **the \f$\lVert V_{xc}-V_{xc}^{fit}\rVert\f$ study** — every Becke-vs-uniform cost comparison in this
  file is taken at UNKNOWN-EQUAL accuracy (nR=40 × L=29 against a 20³ raster), which is not a controlled
  comparison.  Becke may well WIN at equal fit quality, especially on NaF.  Until that is measured, §6's
  "Becke is a negative acceleration" row is a cost statement only, NOT a verdict;
- **the k-scaling gap** (§5a) — diagnosed, two fixes refuted, next step recorded in doc/OpenWork.md;
- **the MnO FM row** (footnote ⁵) — the last pre-Step-2 measurement in the table, ~20 min to re-take;
- **an attribution A/B** — the 09-04 deltas are "current vs banked" across everything that landed since
  08-28, with no parent-commit A/B on the VA recipe;
- **the threaded table** (§7), which needs the busy-wait barrier understood first;
- **NaF's CP2K iteration count** — its log was not to hand, so that row has no per-iteration column.




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


---

## 9. STANDING OBSERVATIONS


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