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
| `QCHEM_BECKE_XC` — the atom-centred (Becke) XC quadrature | **ON** (default since 2026-08-02) | **`=0`, and `CP2K_COMPAT=1` now does it for you.**  CP2K evaluates \f$V_{xc}\f$ on the uniform realspace grid.  ⚠ **This is 43% of the MnO row** (158 s of 367 s wall, the single largest block) — by far the biggest item ever added to this table, and it sat outside it until 2026-08-28 because it is a TYPED option rather than an env flag. |
| `GPW_CONTRACT_CUBE` — the separable-contraction collocation kernel | **ON** (default since 2026-08-27) | **LEAVE IT ON.**  It is NOT a deviation: CP2K collocates exactly this way (`grid_cpu_collint.h`, Mathieu's three 2-D tables).  It is on this list only so a row states which kernel produced it — every run prints `[collocation] kernel=…` unconditionally.  `=0` reverts to the reference box walk, which is an investigation opt-out, not a comparison setting. |

**The general rule, so this table does not have to be exhaustive:** anything that makes qchem faster and is
not in CP2K's algorithm gets declared in the row and defaults to OFF in the comparison.  A new accelerator
is not finished until it is on this list.

★ **AND THE RULE IS A PERMISSION, NOT A BAN** (user, 2026-08-28: *"it's ok to code up non-CP2K
accelerations as long as they are guarded with the `CP2K_COMPAT` flag"*).  Nothing here says do not build
it — it says build it THROUGH `qchem.RunPolicy` so it appears on the deviation line and one switch turns it
off.  An acceleration that cannot be switched off is the only kind this table cannot live with.

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
parity ROUTE affordable — see the parity row and the note below it.  CP2K column untouched throughout.

| system | k-mesh | span | qchem Etot | CP2K Etot | Δ (qchem−CP2K) | wall q / c | **CPU q / c** | **CPU ×** | peak RSS q / c |
|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | −7.115067844 | −7.115057882 | **−10.0 µHa** | **0.8 s** / 5.2 s | **2.2** / 5.0 s | **0.44×** | **29** / 148 MB |
| Si (FCC) | 2×2×2 Γ-centred | SIPP_SR | −7.778472833 | −7.778457865 | **−15.0 µHa** | 7.5 s / 5.8 s | 8.9 / 5.6 s | 1.6× | **30** / 153 MB |
| Si (FCC) | 2×2×2 shifted MP | SIPP_SR | −7.868473428 ¹ | −7.867436530 | **−1.04 mHa** | 16.2 s / 6.1 s | 17.5 / 6.0 s | 2.9× | **31** / 153 MB |
| NaF (rocksalt) | Γ | LOWQ_SR2 (both) | −24.4303364755 | −24.431213375 | **+0.877 mHa** | 21.1 s / 7.4 s | 37.9 / 7.2 s | **5.3×** | **61** / 173 MB |
| NaF (rocksalt) | 2×2×2 Γ-centred | LOWQ_SR2 | −24.5468834873 | — ² | | 1m07.8s / — | 84.5 s / — | — | **68** / — MB |
| NaF (rocksalt) | Γ | LOWQ_SR (full) | −24.4309472653 | −24.432293467 | **+1.346 mHa** | **25.0 s** / 1m42s | **41.4** / 102 s | **0.41×** | **65** / 186 MB |
| **MnO AFM-II — ALL DEFAULTS** ⁴ | Γ | **VA (N=118)** | −61.40297551 ⁴ | −61.303325178 | **−99.65 mHa** | **5m28.1s** / 6m14s | **584** / 373 s | **1.57×** | **491** / 217 MB |
| **MnO AFM-II, `QCHEM_BECKE_XC=0`** ⁶ | Γ | **VA (N=118)** | −61.40358773 | −61.303325178 | −100.26 mHa | **4m05.4s** / 6m14s | **246** / 373 s | **0.66×** | **105** / 217 MB |
| ⚠ STALE MnO FM | Γ | **VA (N=118)** | −61.441583060 ⁵ | −61.304782531 | **−136.80 mHa** | 21m45s / 3m13s | 2321 / 192 s | **12.1×** | 4947 / 217 MB |
| **MnO AFM-II, `CP2K_COMPAT=1`** ⁷ | Γ | **VA (N=118)** | −61.39789688 ⁷ | −61.303325178 | −94.57 mHa | 45m40s / 6m14s | **2736** / 373 s | **7.3×** | **112** / 217 MB |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | | ❓ |

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

⁷ **THE PARITY ROW — IT EXISTS AGAIN, AND THE OLD VERDICT IS RETRACTED (re-measured 2026-08-28).**
`doc/OpenWork.md` and this file have carried *"AT TRUE PARITY OUR MnO RECIPE DOES NOT CONVERGE ... stage 2
caps at −57.620, 3.8 Ha short"* since 2026-08-26.  Re-run on today's tree, same recipe, `CP2K_COMPAT=1`
(so the imposition, the low-rank ρ, the stream fold AND the Becke mesh are all off):

| | 2026-08-26 | **2026-08-28** |
|---|---|---|
| stage 1 | caps at 80, −60.431 | UNSETTLED at 13, −61.41070717 |
| stage 2 | caps at 80, **−57.620** | **FIT-FLOOR STALL at 80, −61.39789688** |
| against the imposed answer (−61.40297551) | **3.8 Ha short** | **5.1 mHa short** |
| the AFM order | — | **SURVIVED both stages** (m_stag 0.636, 0.610) |
| peak RSS | 5034 MB | **112 MB** |

⇒ **It still hits the cap, but for a completely different and far more benign reason.**  It is no longer
collapsing: the energy is settled (ΔE amplitude 2.7e-8), the magnetic order holds without any imposition,
and Δρ has FLOORED at 1.69e-5 against a 1e-5 target — the run's own detector calls it a *"FIT-FLOOR STALL
(Δρ floored, ΔE tiny -- functional/grid)"*, not an oscillation.  That is **A4 territory** (the Δρ/N
convergence gate, `doc/SCFStrategyPlan.md`), which is the third independent thing this session has pointed
at A4.

⚠ **THE HONEST PARITY STANDING, then**: 2736 s against 373 s CPU is **7.3×**, or **3.5× PER ITERATION**
(29.4 s against 8.5 s) since we took 93 iterations to CP2K's 44 — and **0.52× on RAM**.  The per-iteration
gap is bigger than the default row's 2.22× because parity also strips the stream fold, so each collocation
covers ~5× more pairs (2.2 s/call against 0.5 s).  ⇒ **The fold is now the largest single thing the
comparison removes**, which is exactly the qchem-vs-qchem delta rule 3 asks to be reported separately.

⁶ **THE XC-PARITY ROW, and the first MnO row on which qchem beats CP2K on BOTH axes** — 246 s against
373 s CPU and 105 MB against 217 MB.  It is the default recipe with ONE deviation removed: the Becke mesh,
so XC runs the way CP2K runs it (`XC_PairQuadrature` — ρ via the GPW collocation, \f$v_{xc}\f$ pointwise,
\f$H_{xc}\f$ via the exact transpose).  ⚠ **It is NOT the parity row**: the imposition, the low-rank ρ and
the stream fold are all still on.  What it establishes is that the parity ROUTE is no longer what stands in
the way — before 2026-08-28 the same configuration cost 1805 s and 4.5 GB, because a polarized run could not
reach that route at all.

### ★★★ WHAT DELETING THE CACHE DID TO THIS TABLE — read the RAM column and the MnO row together

| system | CPU before → after | peak RSS before → after |
|---|---|---|
| Si Γ | 7.4 → **2.4 s** (3.1× faster) | 267 → **28 MB** (9.5×) |
| Si 2×2×2 Γ-centred | 9.6 → 8.9 s | 269 → **30 MB** (9.0×) |
| Si 2×2×2 shifted MP | 15.6 → 17.5 s | 269 → **31 MB** (8.7×) |
| NaF SR2 Γ | 94.5 → **39.7 s** (2.4×) | 577 → **54 MB** (10.7×) |
| NaF SR2 2×2×2 | 112 → **84.5 s** | 590 → **68 MB** (8.7×) |
| NaF full-SR Γ | 219 → **43.4 s** (5.0×) | 3090 → **58 MB** (**53×**) |
| **MnO AFM-II Γ (imposed, VA)** | 663 → 976 → 837 → 678 → 620 → **584 s** (**0.88× — FASTER**) | 1323 → **491 MB** (2.7×) |

**Two rows now BEAT CP2K on CPU outright** — Si Γ at 0.48× and NaF full-SR at 0.43× — and **every** qchem
row is now well under CP2K's RAM (28–68 MB against 148–186 MB on the small cells), which is the first time
that has been true.  The NaF full-SR row is the headline: 3090 MB and 219 s CPU became 58 MB and 43 s, and
it was the row whose 3 GB used to force `scripts/memsafe`.

⛔ **AND THE MnO ROW GOT SLOWER — 1.26×, and it is not noise.**  That row is a long IMPOSED run (14+17
iterations) where the two box-walk buckets were **477 s of 976 s CPU**, so the cache's measured 2.91× on
those buckets translated almost exactly into the 1.47× the row lost when it went.  ⇒ **The cache was still buying real
CPU on the one row that matters most**, and the case for deleting it rests on the RAM axis and on the two
latent defects it was hiding, not on a free lunch.  ⚠ The 663 s "before" is the 2026-08-19 banked row on an
older binary, so treat the MnO delta as indicative; the directly-measured, same-binary A/B is the
2.91×-on-the-buckets figure in `doc/CollocationRewritePlan.md` step 7.
✅ **AND THE ROW HAS SINCE GONE PAST WHERE THE CACHE LEFT IT — 976 → 837 → 678 → 620 → 584 s, against
the 663 s the 4 GB cache used to buy, on 2.7× less RAM.**  Four changes did it, none needing a re-bank —
the box-walk buckets went 477 → 344 (the `template<int LP>` dispatch) → 192 (the **collocation memo depth
fix**) → 123 (the **gather memo**) → **91 s** (the **exp-table recurrence**).  The first three are
bit-identical by construction; the fourth is not, and moved nothing anybody pins (below).  ⇒ Against CP2K
this row now stands at **1.57× CPU** (was 2.24×), **2.22× CPU per ITERATION** (was 3.2×) — and **0.88× on
WALL, i.e. faster in wall clock.**
### ★★★ THE COLLOCATION MEMO HAD DEPTH 1, AND A POLARIZED RUN ALTERNATES TWO DENSITIES (2026-08-28)

Found by a CALL CENSUS — bucketing the four closure sites so the ledger reports a per-site call count.  On
the MnO row: **64 closure calls, 4 memo hits — 6%.**  `sameD` compared against the LAST collocation only,
so \f$D_\uparrow\f$ and \f$D_\downarrow\f$ evicted each other every single call.  The unpolarized runs the
memo was written against never showed it (Si Γ collocates 1.45×/iteration; MnO was collocating ~20×).

Keeping four older (D, ρ) pairs — same EXACT match rule, so a replay is bit-identical — on the 3-iteration
probe:

| | collocations | memo hits | bucket |
|---|---|---|---|
| depth 1 (as shipped) | 60 | 4 (6%) | 36.7 s |
| **depth 5** | **16** | **48 (75%)** | **9.8 s (3.74×)** |

16 misses is the true number of DISTINCT densities, so depth 5 catches every repeat.  On the full row:
collocate calls **368 → 120**, its bucket **221 → 71 s**, `Etot` unmoved at −61.40297551, +10 MB of RSS.
⚠ This is not an acceleration CP2K lacks — it is removing OUR OWN redundancy; CP2K collocates ρ once per
step.  Knob `GPW_COLLOC_MEMO` (0 restores depth 1).

### ★★★ AND THE SAME CENSUS ON THE GATHER SIDE: 53% OF THE INTEGRATES WERE EXACT DUPLICATES (2026-08-28)

`GPW_INTEGRATE_CENSUS=1` classifies each gather against a short history of (field, screen) hashes.  On the
MnO probe: **30 gathers, 14 distinct, 16 with V AND screen bit-identical** — and not one "same V, widened
screen", so the repeats were pure redundancy, not a screening artefact.

The molecular `IntegrateMemo` could not catch them, for two *correct* reasons: screened calls bypass it (its
key is \f$V\f$ alone while the active set moves with \f$D\f$), and FOLDED calls bypass it (its `nb` records
\f$b\f$ without the orbit multiplicity and the replay never runs `fillImages`).  An imposed polarized run is
both, so it never memoized at all.  ⇒ The fix caches the FINISHED \f$h\f$ on \f$(V_L,\ \text{screen})\f$
in `GPW_Evaluator` — per k-block by construction, which is what makes it exact whatever route produced
\f$h\f$.  Gathers **184 → 74**, bucket **121 → 50 s**, `Etot` unmoved.

⚠ **AND THE PROFILE HAS RE-ORDERED — the box walk is no longer the biggest block on this row.**  ⚠ "This
row" means the DEFAULT row (`CP2K_COMPAT=0`), which runs Becke; the parity row is a different measurement
and is dealt with below.

| block | s (of 367 s wall) | share |
|---|---|---|
| **Becke XC mesh** (ρ sampling 74 + Φ tables 42 + H_xc 23 + mesh build 17 + 3) | **158** | **43%** |
| collocate + integrate-back | 123 | 33% |
| unaccounted (LA, mixing, FFT) | ~75 | 20% |
| local-PP, 1E sums, closures | ~11 | 3% |

⛔ **AND MEASURING IT INVERTED WHAT THAT SHARE SEEMS TO IMPLY — BECKE IS A 3× NET WIN, NOT AN ADVERTISEMENT.**
"43% of the row" does NOT mean "43% to be had by turning it off".  The uniform grid this system needs is
**571787 points against Becke's 48320** — 11.8× more — which is exactly why the `Auto` selector picks Becke
here.  Measured 2026-08-28 with `QCHEM_BECKE_XC=0` and everything else at default (the imposition KEPT, so
the run still converges — full `CP2K_COMPAT` does not, see below):

| MnO AFM-II Γ, VA, imposed | Becke (default) | **uniform XC** |
|---|---|---|
| Etot | −61.40297551 | −61.40295935 (1.6e-5 Ha — the quadrature difference) |
| iterations | 17 | 13 |
| **CPU** | **584 s** | **1805 s (3.09× SLOWER)** |
| wall | 5m28.1s | 30m10.9s |
| **peak RSS** | **491 MB** | **4494 MB (9.2×)** |
| the XC buckets | 158 s | **1592 s** (ρ sampling 874 + Φ tables 488 + H_xc 230) |

⇒ **The problem was not that Becke is an unfair advantage; it was that our UNIFORM route — the parity
route — was ~10× dearer than Becke on this cell.**  ✅ **FIXED 2026-08-28, and it was a CONFLATION, not an
algorithm.**  `XC_PairQuadrature` — ρ via the GPW collocation, \f$v_{xc}\f$ pointwise, \f$H_{xc}\f$ via the
exact transpose, i.e. **CP2K's algorithm, no Φ table anywhere** — already existed and Si Γ already used it.
But `VxcFit::Auto` read `becke || polarized`, and that `polarized` half forced EVERY polarized run onto the
Φ table whatever its grid, because `XC_PairQuadrature::RhoPol` threw.  Nothing spin-specific was ever in the
way: `applyRaw` takes a \f$D\f$, so a channel is one call with that channel's density, and the adjoint
needed no change at all.

| MnO AFM-II Γ, VA, imposed, `QCHEM_BECKE_XC=0` | Φ-table route | **collocation (pair) route** |
|---|---|---|
| **CPU** | 1805 s | **246.5 s — 7.3× faster** |
| wall | 30m10.9s | **4m05.4s** |
| **peak RSS** | 4494 MB | **105 MB — 43× smaller** |
| the XC buckets | 1592 s | **5.2 s** |
| Etot / iterations | −61.40295935 / 13 | −61.40358773 / 25 |

★ **And on that route qchem BEATS CP2K on both axes — 246 s against 373 s CPU, 105 MB against 217 MB.**
⚠ It is still not a parity ROW: the imposition, the low-rank ρ and the stream fold are all still on.  What
it says is that the parity ROUTE is no longer the thing standing in the way.
⚠ The two uniform-XC energies differ by 6.3e-4 Ha because they are different discretisations of the same
quadrature; the pair route is the variational one (\f$H_{xc}=\partial E_{xc}/\partial D\f$ to machine
precision, gate `GPW.RawXCConsistencyFD`).

⇒ The atom-centred XC quadrature is the largest single cost of the DEFAULT row, and it is a **qchem-only
algorithm** — CP2K runs XC on the uniform grid.  ✅ **DECLARED 2026-08-28** as `QCHEM_BECKE_XC` (`RunPolicy::BeckeXC`), so
it is on the deviation table above and `CP2K_COMPAT=1` routes XC to a **basis-sized** uniform grid — the
sizing stays in `qcMesh::ResolveXCMesh`, because handing back a bare `cellKind` would leave `nUniform`'s
basis-blind default of 20 in charge, which is the under-resolution that selector exists to prevent.
⚠ The DEFAULT path is untouched by the change (with `BeckeXC()` true the new overload delegates to the
unchanged two-argument resolver), and that is measured, not asserted: Si Γ reproduces −7.115067844 to all
10 s.f. and `ctest -j8` is 793/793.  What DOES move is every `CP2K_COMPAT=1` row — of which this table has
exactly one, already marked ⚠ STALE.
⇒ `raster` and `cutoffFactor` remain the last two typed options outside the policy (N5).

### ✅ THE EXP-TABLE RECURRENCE — the one place the kernel did not follow CP2K (2026-08-28)

The table build was \f$3n^2\f$ scalar `std::exp` calls where CP2K uses a recurrence.  The exponent is
LINEAR in the inner index of each table (only the cross term \f$2e_ae_bh_{ab}\f$ moves), so a row is one
seed and \f$n\f$ multiplies — **2 exps per row instead of \f$n\f$**.  Seeded at the LARGEST entry and
walked downward, which is the 2026-08-26 underflow rule: a recurrence seeded in the tail can start below
the underflow floor and stay zero through entries that matter.

| | kernel at 32³ | box walk, Si Γ | NaF SR2 Γ | MnO row |
|---|---|---|---|---|
| direct `exp` | 111.8 µs | 0.290 s | 1.71 s | 123 s |
| **recurrence** | **72.9 µs (1.53×)** | **0.204 s** | **1.23 s** | **91 s** |

⚠ It is **NOT bit-identical** — a product of \f$n\f$ rounded factors is not the rounded product — so it was
built default-OFF and defaulted ON only on evidence: against a naive exact reference the contraction goes
1e-15 → **7e-15 relative, flat in box size** (it is \f$n\varepsilon\f$), still **4× better than the WALK's
3e-14**; `ctest -j8` is **793/793 on both settings**; and all three anchors above are unchanged **to all 10
printed s.f.**  ⇒ Anchor-moving in principle, moved nothing in practice.  `GPW_EXP_RECURRENCE=0` is the A/B.

⚠ **AND IT ONLY HALVED THE TABLE TERM, not eliminated it** (89 → 45 µs): with the exps gone the build is
bound by writing \f$n^2\f$ entries, which no algorithm removes — the table has to exist.  The kernel's
\f$O(N^2)\f$ share is 80% → **62%**.

⇒ **And the next KERNEL lever is still NOT the batching.**  Before the recurrence, **84% of the contraction
kernel was the three 2-D Mathieu `exp` tables** (\f$3n^2\f$ scalar `std::exp` calls), worth up to ~35% of this
row.  CP2K builds those tables by RECURRENCE where we call `exp`, and doing the same is the one route that
keeps this table apples-to-apples — a vectorised `exp` would be a qchem-only acceleration under rule 3
above, declared on the deviation line and switched OFF for every head-to-head row, i.e. speed we could not
quote here (user, 2026-08-28).

**Energies moved by ≤ 1e-6 Ha everywhere** — Si Γ −6.6e-7, Si 2×2×2 +1.9e-8, shifted MP 0 (10 s.f.),
NaF SR2 Γ +6.5e-9, NaF SR2 2×2×2 −3.1e-7, NaF full-SR +3.2e-6, MnO −7.0e-7 — which is the anchor re-bank
A1+A7 predicted (doc/OpenWork.md, the anchor-moving sprint) and it does not change any verdict in the Δ
column.

### ★ THE `CP2K_COMPAT=1` ROW — WHAT THE DEFAULT ROW ABOVE IT WAS HIDING (2026-08-26)

**First, the reassuring half: the deviations are pure ACCELERATIONS, not physics.**  Turning all four off
moves the MnO total by **3e-8 Ha** (−61.40297621 → −61.40297618, agreeing to 10 s.f.) on a different
trajectory (10+14 iterations against 14+17).  That is the property `CP2K_COMPAT` most needed to demonstrate
about itself, and it is now measured rather than asserted.

**Then the uncomfortable half.**  Stripped of our own accelerations the standing against CP2K is
**2.5× CPU and 23× RAM**, not the 1.8×/6.1× the default row reports — and 1.79× per ITERATION (38.4 s
against 21.4 s), since the compat run happened to need fewer of them.  The default row is not wrong, but it
is qchem-with-accelerations against CP2K-plain, and it should never be quoted as an algorithm-to-algorithm
comparison.

**⚠ AND `CP2K_COMPAT=1` IS STILL NOT PARITY** — the switch covers four deviations and at least two more
matter (doc/OpenWork.md N5 carries the full list):
- ~~**The pair-stream CACHE**~~ ✅ **DELETED 2026-08-27** (`doc/CollocationRewritePlan.md` step 7), so this
  deviation no longer exists: like CP2K, qchem now re-evaluates the orbital pairs every iteration and keeps
  only a ~0.2–0.4 MB (shell pair, offset) TASK LIST.  The history below is kept because it is what turned
  the campaign toward making the on-the-fly evaluation fast instead of caching harder — and because its
  last row is the one that finally made the cache indefensible.  **What ACTUALLY closed it**: after the
  contraction kernel, the cache bought 2.91× on the two box-walk buckets and ~1.1–1.5× on a whole run,
  against 25× the RAM on the unfolded probe.  ⚠ Read the MnO row above before quoting this as a pure win.
  Original measurement, zeroing the budgets (`GPW_STREAM_BUDGET_PTS=0`, all 8778 pairs on-the-fly) on the
  MnO 3-iteration probe:

  | MnO, per SCF iteration | CPU/iter | peak RSS |
  |---|---|---|
  | CP2K (44 steps, 373 s CPU; its own log prints 8.4–8.5 s/step) | **8.5 s** | **217 MB** |
  | qchem, cache ON  | ~31 s (3.6×) | 4218 MB (19×) |
  | qchem, cache OFF — **as measured 2026-08-25, ESTIMATED** | ~853 s (100×) | 353 MB (1.6×) |
  | qchem, cache OFF — **MEASURED from the ledger, before the box-walk work** | **573 s (67×)** | **166 MB (0.8×)** |
  | qchem, cache OFF — **MEASURED, after it (2026-08-26)** | **260 s (31×)** | **166 MB (0.8×)** |
  | qchem, cache OFF + contracted kernel — **MEASURED 2026-08-27, and this is now the DEFAULT** | **35.5 s (4.2×)** | **155 MB (0.7×)** |

  ⇒ **The cache is not an advantage we hold over CP2K; it is a workaround for an on-the-fly box walk that
  was 67× off theirs, bought with 4 GB.**  Without it our RAM is BETTER than CP2K's (166 vs 217 MB), so
  the target was never "cache more cleverly" or "trade RAM for CPU" — it was **make the on-the-fly pair
  evaluation fast**.  That is now half done: **2.21× measured on this very probe** (doc/OpenWork.md, the
  box-walk section), which is the "even 2 or 3× would pay for itself" the user asked for.
  ⚠ **THE 853 FIGURE WAS AN INSTRUMENT ARTEFACT, NOT A MEASUREMENT.**  It was (2740 s CPU for 3
  iterations) minus an ESTIMATED ~182 s setup, because `RunMnO` drives `SolidCalculation`, which opens no
  report run — so the benchmark's most expensive row was the ONE campaign run with no timing ledger.
  `e8339cf2` gives the arm the same `GpwReport` bracket every other driver holds; the ledger's exclusive
  buckets then sum to the wall clock (1201.3 s of 1202.1 s) and nothing is subtracted.  **The estimate was
  1.5× pessimistic** — hence "67×", not "100×".  There is still no policy hook for the cache, only the raw
  budget knob.
  ⚠ Provenance for both measured rows: `MNO_SKIP_FM=1 GPW_MNO_NMAX=2 GPW_REPORT=1
  GPW_STREAM_BUDGET_PTS=0 GPW_STREAM_BUDGET_PTS_F32=0`, AFM arm, **symmetry FREE and NO fold active**
  (`[fold] collocation streams (T3 pairs): NONE`), `GPW_OMP_THREADS=1`, BLAS pinned to 1, measured 103%
  CPU — i.e. serial, unfolded, and the two binaries differed ONLY in the box-walk diff.  Trajectory
  identical both sides (iters, lastΔρ, m_stag, Eee, site moment); `Efinal` moves 2e-8 Ha.
- ~~**`imposeSymmetry` ITSELF.**~~  ✅ **WIRED 2026-08-26** — `CP2K_COMPAT=1` now implies
  `imposeSymmetry=0` (knob `QCHEM_IMPOSE_SYMMETRY`).

  ⛔ **AND THE RE-TAKEN ROW DOES NOT EXIST: AT TRUE PARITY OUR MnO RECIPE DOES NOT CONVERGE.**  Measured
  2026-08-26, the banked recipe under `CP2K_COMPAT=1` with the imposition vetoed: stage 1 caps at 80
  iterations at −60.431, stage 2 caps at 80 more at **−57.620** — 3.8 Ha short of the −61.40297618 the
  imposed compat run reaches in 24 iterations (2555 s CPU, 42m57s, 4.76 GB).  **So there is currently no
  honest MnO row to put in this table**, and that absence is the finding: the imposed star-average was
  buying CONVERGENCE, not accuracy and not the magnetic basin.
  ★ **The AFM ORDER SURVIVED the free run** (m_stag 0.66/0.59, integrated site moment 4.781 → 4.222 e), so
  the imposition was NOT what held the basin — which is the opposite of what was expected, and it moves
  the question from symmetry to the MIXER.  Note where CP2K stands on the same cell: 44 steps at 8.5 s
  each with Broyden α=0.2 / NBUFFER 8 / MAX_SCF 200, against our α=0.45 / PulayDepth 0 / 80.  ⇒ A fair
  MnO row needs a CP2K-like mixing recipe first; taking one before that would be comparing their converged
  answer against our iteration cap.
  Verified locally 2026-08-26: CP2K does NO symmetry work in these decks.
  The 1129-line `bench_MnO_AFM2_VA_cp2k.log` contains **zero** occurrences of "irrep", "symmetry" or
  "point group" — QuickStep keeps K and P as DBCSR sparse ATOM-BLOCK matrices over the full AO basis and
  diagonalizes the lot (blocking by atom-pair SPARSITY, not by irrep; no SALC blocking, no k-block
  splitting).  The one symmetry knob that exists is BZ-side and defaults OFF: our own Si 2×2×2 log reads
  `BRILLOUIN| K-Point point group symmetrization  OFF` and lists all 8 k-points.
  ⇒ Our imposed MnO row folds the BZ, star-averages ρ every iteration, uses the site-adapted invariant XC
  mesh (~2×) and folds the collocation streams (5.2× on pairs).  **The CP2K row does none of it.**  These
  rows compare a SYMMETRY-exploiting code against a SPARSITY-exploiting one on a small, high-symmetry
  cell — the regime that most favours us, and a 100-water box would invert it.  (CP2K's design centre is
  large disordered systems, where every one of {G}, {k}, {r} has a group of order 1, so the folding payoff
  is exactly 1× for what they build for.  That reading of WHY is inference; the WHAT above is from the logs.)

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
