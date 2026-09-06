# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**This is an INSTRUMENT, not a report.**  A row here is a claim that two codes did the same work on the
same hardware; the sections below are what makes that claim checkable.
📖 **The reasoning, the wrong turns and every superseded number are in `doc/BenchmarkHistory.md`** (split
out 2026-09-04).  Come back here for what is TRUE NOW; go there for why.

**Read in this order:** §1 the process → §2 what parity actually means → §3 the rules → §4 the commands →
**§5a the BIN-1 TABLE** → §5 the whole-run rows → **§5f why the parity row is 2.05×**.
★ **§5a is the one table to point at for per-iteration CPU** — it holds nothing else, by request
(2026-09-05); the deltas and open questions that used to crowd it are in §5e.

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
| **1** | **per-iteration CPU** | ★ **PER SCF ITERATION WE ARE AHEAD ON 7 OF 9 ROWS** (§5a, 2026-09-05, both codes pinned serial): MnO ALL DEFAULTS **0.82×**, MnO FM **0.77×**, `QCHEM_BECKE_XC=0` **0.45×**, Si **0.14× / 0.66× / 0.56×**, NaF full-SR **0.18×**.  The two losses are NaF SR2 **1.45×** (a Becke cost) and — the one that counts — `CP2K_COMPAT=1` **1.72×** (was 2.05× before §5f's lever A).  ⇒ **NOT closed** — but on the SAME ALGORITHM (our fixed-point stage against CP2K's diagonalise-and-mix decks) it is **1.13×**, three gathers against their two, and the third gather is §5f lever B, behind N4.  The 1.72× two-stage figure adds a GDM line search CP2K's runs have no counterpart for |
| **2** | **init / pre-iteration** | ⛔ **PROMOTED, and now measured serially: MnO's setup is 184 s = 47% of the default run against CP2K's 8.1 s (23×)**, of which **136.6 s is two Becke mesh builds** (§5a).  The earlier "~57 s of 328 s" was a THREADED ledger bucket — the build's partition loop is `omp parallel for`.  ✅ With `QCHEM_BECKE_XC=0` our setup is **1.76 s and beats CP2K's 8.1 s** — so bin 2 is a Becke-mesh question, exclusively |
| **3** | **peak RAM** | ✅ **solved, and we win**: ~476 MB defaults, **113–132 MB on the parity routes against CP2K's 217 MB** |
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
| **`[t=12.34 s]` on every section heading, fold and log line** | the run clock, `report::RunElapsed()` (added 2026-09-06, `doc/OpenWork.md` Step 0c) |

⚠ **A STAMP IS A MOMENT, NOT A DURATION — AND A GAP BELONGS TO WHAT RAN *BEFORE* IT, NOT TO THE LINE THAT
CARRIES IT.**  The stamped line is the END of the silent stretch above it.  So in
`scf ▸ siteMoments  [t=95.93 s]` after a stage boundary at 60.17 s, the 35.8 s is the stage rebuild that
FINISHED at 95.93 — it says nothing about the cost of the site moments, which are **1.7 ms per call, 0.053 s
over the whole run** (`scf: order probe`, bucketed 2026-09-06 for exactly this reason).  Attributing a gap
to the item that closes it is the one wrong inference this instrument makes easy; the LEDGER is what prices
an item, the stamps only localise WHEN.

**READ THE STAMPS BEFORE THE LEDGER.**  The `timing` table says WHAT the run spent; the stamps say WHEN,
and the GAPS BETWEEN CONSECUTIVE STAMPS are exactly the time no bucket is charging.  That is how the last
25 s of unbucketed MnO time was located (a 35.8 s silent stretch between two anneal stages) without adding
a bucket first — cheaper than guessing where to put the next probe, and it works on a run you have already
taken.

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

## 5a. ★★★ THE BIN-1 TABLE — per-ITERATION CPU, qchem vs CP2K, one thread each

**ONE TABLE, ON PURPOSE** (user, 2026-09-05: *"the important comparison table was the only table in
section 5a, but now the section is flooded with other tables so I cannot easily point to it"*).  Every A/B
delta, ledger split and open question that used to sit here is in **§5e**.  This is the bin-1 instrument
and the only table to quote for bin 1.

**THREAD STATE — SERIAL ON BOTH SIDES, AND MEASURED SO** (rule 3b): qchem `OMP_NUM_THREADS=1
GPW_OMP_THREADS=1`, measured **99% CPU** on every row; CP2K `OMP_NUM_THREADS=1`, measured 97–99%.
**Taken 2026-09-05**, one box (14 GB, 16 cores), qchem at `9f4f4ae2` built `-O3 -march=native`, CP2K 2025.2,
commands copied from §4.  ★ **Every row RE-TAKEN after §5f's lever A** (the Hartree energy stopped building
a matrix): the parity row fell 2.05× → 1.72×, the Si 8-k rows 0.84×/0.68× → 0.66×/0.56×, MnO defaults
0.90× → 0.82×, `BECKE_XC=0` 0.53× → 0.45×.  ⚠ **The qchem `setup` column is much larger than the 09-04 entry's ~57 s for MnO
because the Becke mesh build THREADS** (`src/Structure/Imp/UnitCell.C:273` — the partition loop is
`#pragma omp parallel for` over quadrature points, and `GPW_OMP_THREADS=1` pins it): 68.3 s serial per
build against the 16.7 s that ledger read with threads free.  The serial figure is the one that belongs
beside a serial CP2K row (`doc/BenchmarkHistory.md` §9).

**HOW TO READ IT.**  `CPU` is whole-run user+sys.  `setup` is that run's own pre-SCF work — qchem: the sum
of the ledger's `setup:` buckets; CP2K: total CPU minus the sum of its printed per-step times.  So

- **`s/it (SCF)` = (CPU − setup) / iterations is the BIN-1 number**, and the `×` column is its ratio;
- **`setup` is the BIN-2 number**, in the same row, on the same run.

⇒ **Compare the SCF columns.**  The total column is kept only because it is what the whole-run table (§5)
divides — and the two disagree by 1.8× on MnO precisely because 44% of that run is setup.

| row | span / k | q iters | q CPU | q setup | **q s/it (SCF)** | c steps | c CPU | c setup | **c s/it (SCF)** | **BIN 1 ×** | q s/it (total) | total × |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Si Γ ᵇ✓ | SIPP_SR, 1 k | 17 | 1.10 s | 0.16 s | **0.055** | 12 | 5.03 s | 0.4 s | 0.383 | **0.14×** ✅ | 0.065 | 0.15× |
| Si 2×2×2 Γ-centred ᵇ✓ | SIPP_SR, 8 k | 16 | 4.46 s | 0.35 s | **0.253** | 13 | 5.60 s | 0.6 s | 0.385 | **0.66×** ✅ | 0.275 | 0.64× |
| Si 2×2×2 shifted MP ᵇ✓ | SIPP_SR, 8 k | 14 | 4.26 s | 1.15 s | **0.218** | 14 | 5.91 s | 0.5 s | 0.386 | **0.56×** ✅ | 0.300 | 0.71× |
| NaF SR2 Γ ᵇ✓ | LOWQ_SR2, 1 k | 23 | 23.0 s | **9.82 s** | **0.573** | 16 | 7.18 s | 0.9 s | 0.394 | **1.45×** ⛔ | 1.000 | 2.23× |
| NaF full-SR Γ ᵇ✓ | LOWQ_SR, 1 k | 30 | 30.0 s | **11.17 s** | **0.628** | 27 | 101.85 s | 6.9 s | 3.519 | **0.18×** ✅ | 1.000 | 0.27× |
| **MnO AFM-II — ALL DEFAULTS** (5 of §2's 7 deviations active ᵃ) ᵇ⚠ | VA, 1 k | 14+17 = **31** | 395.6 s | **184.3 s** | **6.817** | 44 | 372.9 s | 8.1 s | 8.291 | **0.82×** ✅ | 12.76 | 1.51× |
| MnO AFM-II, `QCHEM_BECKE_XC=0` (4 of 7 — ONE rung down, NOT parity ᵃ) ᵇ⚠ | VA, 1 k | 14+25 = **39** | 147.0 s | 1.81 s | **3.723** | 44 | 372.9 s | 8.1 s | 8.291 | **0.45×** ✅ | 3.769 | 0.44× |
| ⚠ MnO AFM-II, `CP2K_COMPAT=1` probe (0 of 7 — `AT PARITY` as far as is KNOWN ᵃ) ᵇ⚠ | VA, 1 k | 10+10 = **20**, CAPPED | 287.5 s | 1.89 s | **14.28** | 44 | 372.9 s | 8.1 s | 8.291 | **1.72×** ⛔ | 14.38 | 1.70× |
| ★ **…the same row's FIXED-POINT stage alone — THE LIKE-FOR-LIKE NUMBER** ᵇ✓ | VA, 1 k | marginal | 37.4 s / 4 it | cancels | **9.35** | 44 | 372.9 s | 8.1 s | 8.291 | **1.13×** | 9.35 | 1.13× |
| MnO **FM** — ALL DEFAULTS ᵇ⚠ | VA, 1 k | 18+15 = **33** | 398.0 s | **184.8 s** | **6.460** | 22 | 192.4 s | 8.7 s | 8.350 | **0.77×** ✅ | 12.06 | 1.38× |

ᵇ **✓ = BOTH CODES RUN THE SAME KIND OF SCF STEP ON THIS ROW; ⚠ = THEY DO NOT.**  A per-iteration ratio is
only a comparison when the iteration is the same thing on both sides, so this was CHECKED per row, in each
run's own trace, rather than assumed:

| row | qchem's step | CP2K's step | |
|---|---|---|---|
| Si Γ, Si 2×2×2 (both) | DIIS + linear mixing on a diagonalise | `DIIS/Diag.` + `P_Mix` | ✓ |
| NaF SR2 Γ, NaF full-SR Γ | DIIS + Kerker on a diagonalise — the ladder's GDM rung is **NOT ENGAGEABLE** (Fermi-smeared occupations sit outside the integer-occupation manifold GDM rotates, and it says so) | `Broy./Diag.` | ✓ |
| the four MnO rows | **two stages**: `Ladder` (fixed point) then `accel: GDM` — DIRECT MINIMISATION with a geodesic line search | `Broy./Diag.` | ⚠ |

⇒ **Only the MnO rows are affected, and only through their second stage.**  CP2K's counterpart for a direct
minimiser is OT — and (checked in the decks and in every log's update-method column) **none of the
benchmarked decks run `&OT` either**: they all diagonalise and mix.  The one deck with an `&OT` section is
the unbenchmarked `naf_gpw.inp`, and its comment says why (diagonalisation diverged on that diffuse basis).
So on these rows the mismatch is not GDM-vs-OT, it is **minimiser vs mixer**, and stage 2's line-search
trial densities land in the probe's per-iteration average with nothing on the other side to match them.

The ★ row isolates the comparable half by DIFFERENCING `GPW_MNO_NMAX` (a single-stage `MNO_ANNEAL=5e-3` run
at N=6 against N=2, so setup and seed cancel and what is left is the marginal cost of one iteration):
**3 gathers + 2 collocations per iteration against CP2K's 2 + 2**, and at our per-call rate that is 9.35 s
against 8.29.  ⇒ **On the same algorithm we are 1.13×, and the whole residual is the third gather** —
\f$V_H\f$ gathered separately from \f$v_{xc}^\sigma\f$ (§5f lever B).

ᵃ **counted off each run's own banner**, not from memory (rule 3b).  Defaults:
`QCHEM_DM_LOWRANK=on* GPW_STREAM_FOLD=on* QCHEM_MIX_RHO_M=off GPW_XC_DM_SOURCE=off
QCHEM_IMPOSE_SYMMETRY=on* QCHEM_BECKE_XC=on* GPW_DAWARE_SCREEN=on*` (`*` = differs from CP2K, five of
them); the middle row is the same with `QCHEM_BECKE_XC=off(stated)`; the bottom row prints
`CP2K_COMPAT=1 -> AT PARITY` with every flag off.

⇒ **PER SCF ITERATION WE ARE AHEAD OF CP2K ON SEVEN OF THE NINE ROWS**, including MnO on its default
route (0.82×) and its FM arm (0.77×).  Two rows are behind, and they are behind for two different reasons:

⛔ **NaF SR2 Γ — 1.45×, and it is a Becke cost, not a GPW one.**  43% of that row's CPU is setup.  §5e.
✅ Its DRIVER is comparable (footnote ᵇ): both codes diagonalise and mix there — the ladder's GDM rung is
not engageable under Fermi smearing — so unlike the MnO rows, this 1.45× is a like-for-like step cost and
nothing about it is waiting on OT.

⛔ **THE PARITY ROW, AND WHAT IT ACTUALLY MEASURES.**  `CP2K_COMPAT=1` removes the stream fold (5.2× on
MnO's pair count) and the low-rank ρ as well, so this is the honest algorithm-to-algorithm row and the
other MnO rows are not.  It reads **1.72×** — but the two-stage probe runs a MINIMISER (GDM) that CP2K's
decks do not, so read the two rows separately (footnote ᵇ):
- **on the same algorithm — fixed point, diagonalise and mix, both codes — we are 1.13×**, and the entire
  residual is the third gather (§5f lever B, blocked behind N4);
- the **1.72×** adds our stage-2 GDM line search, whose trial densities have no counterpart in a CP2K run
  that is mixing rather than minimising.  ⇒ It is not a like-for-like step cost, and chasing it as one
  would be optimising against a comparison that does not exist.
⚠ Two further caveats, both real: `GPW_MNO_NMAX=10` CAPS the probe at 20 iterations with a different stage
mix, so rule 3d applies; and **`CP2K_COMPAT=1` is our best KNOWN parity, not proven parity** — §2's
deviation list has grown every time anyone has looked (4 → 7 items).
⇒ **BIN 1 IS NOT CLOSED, but it is now one gather wide on the comparable algorithm.**

★★★ **AND THE TABLE SIZES BIN 2, MEASURED.**  MnO's serial setup is **184.3 s — 47% of the default run —
against CP2K's 8.1 s (23×)**, and **136.6 s of it is TWO Becke mesh builds** (68.3 s each, one per anneal
stage; the rest is 43.1 s of XC-mesh Φ tables).  Both vanish with `QCHEM_BECKE_XC=0`, where our setup is
**1.76 s and BEATS CP2K's 8.1 s**.  ⇒ Bin 2 is a Becke-mesh question, exclusively.

⛔ **BUT IT IS A GRID-SIZE QUESTION BEFORE IT IS A CODE QUESTION** (user, 2026-09-06): *"I am reluctant to
work on that until identify the proper becke grid size (Nradial=40, Nangular=29 is big) that yields the
same accuracy … as the uniform 20³ mesh.  Becke is always more expensive to setup, no way around that.  It
is parallel which helps."*  The build scales with the point count, and `nR=40, degree=29` sets the entire
Becke side of every cost figure on this page — so the recipe, not the loop, is what to attack first.
⇒ **The calibration is the gate**, and the harness for it exists (`BeckeLadder` + `UniformXCProbe` in
`IntegrationTests/GPW_SCF_UT.C`, scoring \f$\Delta E_{xc}\f$ and \f$\max|\Delta V_{xc}|\f$ on one frozen
density); what is missing is putting the uniform probe on the ladder so the crossover can be read off.
`doc/OpenWork.md`'s next action carries it, including the two standing rules the V2.6a work earned (a
frozen density understates the self-consistent shift on a METAL; Al's angular convergence is
non-monotonic).
⚠ **Two things survive the gate whatever it says**: the build is `omp parallel for`
(`src/Structure/Imp/UnitCell.C:273`, ~4× — the user's *"it is parallel which helps"*), and **there are two
builds for the same cell**, one per anneal stage, which is a redundancy independent of grid size.

⚠ **WHAT MOVED, AND WHY — rule 3a's check, honestly reported.**  **Iteration counts**: removing the
gather's density screen (`a7561e92`) is not bit-identical, so the small rows took a different SCF path —
Si Γ 11 → 17, Si 2×2×2 Γ-centred 7 → 16, shifted MP 16 → 14, NaF SR2 29 → 23.  **Every MnO row still
reproduces its banked iteration count exactly** (31 / 39 / 20 / 18+15).
**Energies**: they now differ in their last digits from the banked values, by two roundoff-scale effects
stacked — the unscreened gather, and lever A's G-space \f$E_H\f$ (§5f), which sums the same integral in a
different order.  Si Γ reproduces to all 10 printed digits; the other Si rows to 2–7 nHa; NaF to 2.7–6.3
µHa; **MnO to 2.2–2.6e-7 Ha** (−61.40297529 against the banked −61.40297551), i.e. 4e-9 relative on a 61 Ha
total and four orders below the 99.65 mHa the row is actually measuring.  ⚠ Anchors pinned tighter than
1e-7 on a periodic total energy will need re-banking; none in the suite are (806/806 green).
The MnO **FM** row is a different case: it had not been re-taken since 08-19 (footnote ⁵), so this is its
first measurement on the current code — 398 s CPU and 481 MB against the stale row's 2321 s and 4947 MB.

---

## 5. THE ROWS — single thread

⚠ **THE WHOLE-RUN TABLE BELOW IS NOT THE BIN-1 INSTRUMENT — §5a ABOVE IS** (user, 2026-09-04: *"the
previous table in section 5 is not very informative"*).  Its `wall` and `CPU` columns are whole-run totals
over runs whose iteration counts differ by 3×, and they charge each run's setup to its per-step cost.
**Read it for two things: the ENERGY column (Δ vs CP2K, the accuracy claim) and peak RSS (bin 3).**
For runtime go to §5a — same runs, decomposed.

Energies in Ha.  **Both columns measured on this box (14 GB, 16 cores) through `scripts/bench`** — the CP2K
side is no longer banked prose: `apt`'s CP2K 2025.2 reproduces every banked 2026.1 deck value to the printed
digits (`doc/CP2KBuild.md`), so both codes are measured under one wrapper.  Provenance per row is in
*How each row was produced*.

> ✅ **BOTH SIDES ARE NOW SERIAL AND MEASURED SO** (2026-09-05): qchem `OMP_NUM_THREADS=1
> GPW_OMP_THREADS=1` at 99% CPU on every re-taken row, CP2K `OMP_NUM_THREADS=1` at 97–99%.  Wall and CPU
> therefore agree to ~1% and the `CPU ×` column is finally a like-for-like ratio.
> ⚠ **The history this replaces, because it will be quoted**: rows before 09-05 had `GPW_OMP_THREADS` unset
> while blaze threaded the BLAS regardless — measured 115–239% CPU — so their wall column FLATTERED qchem
> by the number of cores it happened to take.  **A thread-count knob is not a measurement**; an earlier cut
> of this table printed "1 thr" on the qchem rows on the strength of the unset knob.

**★ EVERY ROW BUT TWO RE-TAKEN 2026-09-05, SERIAL, ON ONE BINARY** (`a3045163`, `-march=native`) — the same
runs §5a decomposes, so the two tables cannot disagree.  The `taken` column says which rows are not from
that cut.

| row | last taken | why it may have moved since |
|---|---|---|
| Si Γ, Si 2×2×2 (both), NaF SR2 Γ, NaF full-SR Γ, MnO AFM-II defaults, MnO `BECKE_XC=0`, **MnO FM** | **09-05, serial, current** | — |
| NaF SR2 2×2×2 | 08-28 | no CP2K counterpart (footnote ²), so it is not a head-to-head row; ⚠ cheap to re-take |
| **MnO `CP2K_COMPAT=1`** (the uncapped 93-step row) | **08-28** | ~45 min to re-take; the 09-05 cut used the documented `GPW_MNO_NMAX=10` probe instead (§5a, §5b).  The old row's *"does not converge"* verdict is RETRACTED — see ⁷ |

⇒ **The `CP2K_COMPAT=1` row is now real rather than aspirational**, because the XC pair route made the
parity ROUTE affordable — see footnote ⁷ (§5d).  CP2K column untouched throughout.

| system | k-mesh | span | qchem Etot | CP2K Etot | Δ (qchem−CP2K) | wall q / c | **CPU q / c** | **CPU ×** | peak RSS q / c |
|---|---|---|---|---|---|---|---|---|---|
| Si (FCC) | Γ | SIPP_SR | −7.115067844 | −7.115057882 | **−10.0 µHa** | **1.1 s** / 5.2 s | **1.1** / 5.0 s | **0.22×** | **30** / 148 MB |
| Si (FCC) | 2×2×2 Γ-centred | SIPP_SR | −7.778472833 | −7.778457865 | **−15.0 µHa** | 4.5 s / 5.8 s | **4.4** / 5.6 s | **0.79×** | **35** / 153 MB |
| Si (FCC) | 2×2×2 shifted MP | SIPP_SR | −7.868473429 ¹ | −7.867436530 | **−1.04 mHa** | 4.3 s / 6.1 s | **4.2** / 6.0 s | **0.71×** | **34** / 153 MB |
| NaF (rocksalt) | Γ | LOWQ_SR2 (both) | −24.4303428161 | −24.431213375 | **+0.871 mHa** | 23.2 s / 7.4 s | **23.0** / 7.2 s | **3.2×** | **55** / 173 MB |
| NaF (rocksalt) | 2×2×2 Γ-centred | LOWQ_SR2 | −24.5468834873 | — ² | | 1m07.8s / — | 84.5 s / — | — | **68** / — MB |
| NaF (rocksalt) | Γ | LOWQ_SR (full) | −24.4309445921 | −24.432293467 | **+1.349 mHa** | **30.2 s** / 1m42s | **30.0** / 102 s | **0.29×** | **58** / 186 MB |
| **MnO AFM-II — ALL DEFAULTS** ⁴ | Γ | **VA (N=118)** | −61.40297529 ⁴ | −61.303325178 | **−99.65 mHa** | **6m38s** / 6m14s | **396** / 373 s | **1.06×** | **476** / 217 MB |
| **MnO AFM-II, `QCHEM_BECKE_XC=0`** ⁶ | Γ | **VA (N=118)** | −61.40358753 | −61.303325178 | −100.26 mHa | **2m28s** / 6m14s | **147** / 373 s | **0.39×** | **113** / 217 MB |
| **MnO FM — ALL DEFAULTS** ⁵ | Γ | **VA (N=118)** | −61.44158219 ⁵ | −61.304782531 | **−136.80 mHa** | **6m40s** / 3m13s | **398** / 192 s | **2.07×** | **481** / 217 MB |
| **MnO AFM-II, `CP2K_COMPAT=1`** ⁷ | Γ | **VA (N=118)** | −61.39789688 ⁷ | −61.303325178 | −94.57 mHa | 45m40s / 6m14s | **2736** / 373 s | **7.3×** | **112** / 217 MB |
| MnO AFM-II | 2×2×2 (`MNO_KMESH=2`) | VA | ❓ | ❓ | | ❓ | ❓ | | ❓ |


### 5d. Footnotes to the table

Compact here; the full stories are in `doc/BenchmarkHistory.md` at the section named after each.

- **¹** Si 2×2×2 shifted MP (−1.04 mHa) — the residual after a **D-aware integrate-back SCREEN defect** was
  fixed 2026-08-19.  It is the suite's ONLY fractional-k SCF coverage (every other k is TRIM, where the
  defect is structurally invisible), which is why it had rotted to −3.7351 while DISABLED.  *(history: the
  screen was reading \f$\mathrm{Re}[D\overline{e^{ikR}}]\f$ where it needed \f$|D|\f$.)*
- **²** No CP2K counterpart: CP2K's NaF decks carry no `&KPOINTS` section, i.e. they are Γ, while the qchem
  test defaults to 2×2×2.
- **⁴** MnO ALL DEFAULTS — **re-measured four times** (08-19/20, 09-04, 09-05), same command, same box;
  every cut agrees on the energy to 4e-8 and on the iteration count (14+17 = 31) exactly.
- **⁵** ✅ **RE-TAKEN 2026-09-05** (it had been the only pre-Step-2 measurement left): serial, VA recipe,
  18+15 = 33 iterations, \f$E_{tot}\f$ −61.44158232 against the 08-19 row's −61.441583060 (7e-7, the
  gather-screen removal).  Its cost columns fell 5.6× in CPU and 10× in RAM; per SCF iteration it is
  **0.83× CP2K** (§5a).
- **⁶** ⚠ **ONE of §2's seven deviations off, NOT parity** (user, 2026-09-05) — it makes XC run CP2K's way
  and leaves everything else ours.  It is the first MnO row on which qchem beat CP2K on both axes and is
  still the standout: **0.53× per SCF iteration**, 107 MB against 217 MB, and a setup of 1.76 s against
  CP2K's 8.1 s.
- **⁷** The parity row: an earlier *"at true parity our recipe does not converge"* verdict was **RETRACTED**
  (08-28).  It still hits the iteration cap, but for a far more benign reason — 5.1 mHa short, not 3.8 Ha,
  with the AFM order surviving both stages.
- **⁸** The 93-iteration parity figure is pre-fix and its stage mix differs — **do not** read 29.4 → 19.3 as
  a 1.5×; per-iteration is not comparable across iteration caps (rule 3d).

### 5b. The three MnO rows are ONE system and ONE recipe

★ **THE THREE `MnO AFM-II` ROWS ARE ONE SYSTEM AND ONE RECIPE, with progressively more of OUR deviations
switched off.**  Read them as a ladder, not as three experiments — the only thing changing is which of §2's
**seven** declared deviations are on.  ⚠ `QCHEM_BECKE_XC=0` is ONE of the seven: it is a rung on the way to
parity, **not** parity (user, 2026-09-05: *"there is much more to parity than just BECKE_XC=0"*).

| the row | what is off | CPU (serial, 09-05) | RSS |
|---|---|---|---|
| **ALL DEFAULTS** | nothing — every deviation at its qchem default (Becke XC, imposition, low-rank ρ, stream fold) | 396 s (31 it) | 476 MB |
| **`QCHEM_BECKE_XC=0`** | ONE of the seven: the Becke XC mesh — so XC runs CP2K's way, everything else still ours | **147 s (39 it)** | **113 MB** |
| **`CP2K_COMPAT=1`** | ALL SEVEN KNOWN — the honest algorithm-to-algorithm row | 288 s (20 it, CAPPED probe) | 132 MB |

⚠ The middle row is the fastest because Becke's atom-centred quadrature is expensive machinery — and on the
default route it is 44% of the run before SCF even starts (§5a).  The bottom row costs the most per step
because parity also removes the STREAM FOLD (5.2× on MnO's pair count) and the low-rank ρ; its 341 s is a
`GPW_MNO_NMAX=10` probe, not a converged run (the uncapped row was 2736 s over 93 capped steps).
⇒ **Always say WHICH of the three** — "the MnO row" has meant all three of these in the space of one day.
⚠ And **`CP2K_COMPAT=1` is our best KNOWN parity, not proven parity**: §2's list started at four items and
is at seven, having grown every time anyone looked (user, 2026-09-05: *"CP2K_COMPAT=1 is supposed to be
parity, but we keep finding new differences"*).  A new find moves the bottom row, in either direction.


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

## 5e. Bin-1 working notes — the deltas and the open questions BEHIND §5a's table

⚠ **§5a holds ONE table on purpose.**  Everything that used to sit beside it and made it hard to point at
lives here.  Nothing in this section is the bin-1 instrument; §5a is.

### What removing the gather density screen cost and bought (2026-09-04)

★★★ **THE UNSCREENED GATHER (2026-09-04) — A REAL SPLIT BY k-MESH, NOT A UNIFORM WIN.**  Removing the
density screen from the integrate-back (see §6, and doc/OpenWork.md for why it is also a CORRECTNESS fix)
unlocks the k-independent memo, so one real-space sweep serves every k-block.  Where there are no k-blocks
to share it across, it is pure added width:

| GATHER seconds PER ITERATION (ledger bucket ÷ iterations) | screened | **unscreened** | |
|---|---|---|---|
| Si Γ (1 k) | 0.00693 | 0.00784 | **+13%** ⛔ |
| Si 2×2×2 Γ-centred (**8 k**) | 0.196 | **0.112** | **−43%** ✅ |
| Si 2×2×2 shifted MP (**8 k**) | 0.187 | **0.111** | **−41%** ✅ |

and on the MnO rows — all Γ, and all with IDENTICAL iteration counts before and after, so total CPU is a
fair comparison there:

| MnO, VA recipe, serial, NATIVE | iters | screened | **unscreened** | |
|---|---|---|---|---|
| ALL DEFAULTS | 31 | 386.1 s | 408.6 s | **+5.8%** ⛔ |
| `QCHEM_BECKE_XC=0` | 39 | 155.7 s | 168.2 s | **+8.1%** ⛔ |
| `CP2K_COMPAT=1` probe | 20 | 338.9 s | 346.1 s | **+2.1%** ⛔ |

⇒ **IT SAVES ~40% OF GATHER TIME ON AN 8-k RUN AND COSTS 2–8% AT Γ.**  Every MnO row is Γ, so the flagship
numbers get slightly worse while the multi-k rows get materially better.  **KEPT**: it is first a
CORRECTNESS fix — the diagonal-seed defect and the stream fold's orbit-invariance, both in
doc/OpenWork.md — and the Γ cost is the price of not having a self-fulfilling truncation in the Fock.
⚠ It is NOT a `CP2K_COMPAT` matter: the screen is gone on every route, not just the parity one.

### ⚠ THE BECKE MESH BUILD THREADS — which collides with §7a's reading (open, 2026-09-05)

§5a's setup column is SERIAL, and 68.3 s of it per build is the Becke mesh.  That loop is
`#pragma omp parallel for` over quadrature points (`src/Structure/Imp/UnitCell.C:273`), so it is NOT
structurally serial — yet §7a's NaF threading run inferred a ~9 s serial residual and matched it to the
setup buckets, of which 7.0 s IS the mesh build.  Both cannot be the whole story.
⇒ **Resolve it while profiling bin 2** (doc/OpenWork.md's named next action): measure the mesh build's own
speed-up curve on MnO and on NaF separately.  Candidate explanation — NaF is a 2-atom cell whose point
count is too small to amortise the region, so it threads on MnO and does not on NaF.

★ **A note on units, since the two columns are subtracted**: the ledger's buckets are WALL seconds and
`CPU` is user+sys.  Every §5a row measured 97–99% CPU, so the subtraction is sound there; it would NOT be
on a threaded row, which is one more reason §5a is a serial-only table.

### The k-SCALING gap — the one structural bin-1 defect the table shows

Same system, same span, only the k-mesh changing (SCF-only seconds per iteration, §5a's numbers):

| | qchem s/iter | CP2K s/iter |
|---|---|---|
| Γ | **0.051** | 0.383 |
| 8 k-points (Γ-centred) | **0.323** | 0.385 |
| **rise** | **6.3×** | **1.01×** |

We start **7.5× ahead** of CP2K at Γ and keep only **1.2×** of it at 8 k — still a win, but the whole
shape of the small-cell story is in that rise.  Cause and the two refuted fixes:
doc/OpenWork.md's k-scaling entry (the gather memo is bypassed whenever a density screen is passed, and 32
of 51 gathers on the 8-k run are the SAME FIELD).  ⚠ Part of their flatness is not an optimisation to copy
but a fold we do not do: CP2K folds a shifted MP mesh 8→4 by time reversal (§6b item 1).

### NaF's loss is a BECKE cost, not a GPW one

`NaF SR2 Γ` is the one small row where we are behind per iteration, and its ledger says why —
Becke mesh build **7.00 s (30% of the run)** + XC-mesh ρ sampling **6.55 s (28%)** + hamiltonian ctor
**2.02 s (9%)**, against a box walk (gather + collocate) of **0.80 s (3%)**.
⇒ On the system where Becke should be at its best, **Becke IS the runtime**, and 42% of that row's CPU is
pre-SCF setup — i.e. it is mostly a bin-2 row wearing a bin-1 label.  That makes §6's parked
\f$\lVert V_{xc}-V_{xc}^{fit}\rVert\f$ study the deciding measurement for this row too.

### The box walk is ~95% of an MnO parity run

109 s of 114.6 s, of which the GATHER is 74% (81.5 s over 65 calls, against 27.9 s over 22).
⇒ Amdahl leaves nothing outside the walk worth touching on that route, and any per-call win is worth ~3×
more on the gather side than on the collocate side.

### STANDING PROBE for bin 1

(user, 2026-08-28: *"we just cut off at ~10 or so iterations, just to get a decent average"*):
`CP2K_COMPAT=1 GPW_MNO_NMAX=10` — ~6 minutes, and quote per-call beside per-iteration.

### ✅ THE PER-STEP COMPARISON — ANSWERED 2026-09-05, see §5f (kept for the question it asked)

★ **RESOLVED BY COUNTING**: a fixed-point iteration issues **4 gathers + 2 collocations**, and both of the
readings below turned out to be true of a different part of it — the term assembly asks twice for Hartree
(§5f lever A) *and* the GDM line search does build extra densities once it engages (lever C, and it is
bin-4 currency).  CP2K's own counters say it does 2 and 2.  The original entry, which framed the question:

The 09-04 ledger shows **~9 distinct KS-field integrations per SCF iteration** on a 2-channel system where
the physics needs ~3 (one \f$V_H\f$ gather + one \f$V_{xc}\f$ per channel); `GPW_INTEGRATE_CENSUS=1` says
all of them are genuinely NEW fields, so it is NOT redundancy a memo could remove.  Two readings, and they
call for opposite responses:
- **the accelerator is buying bin 4 with bin 1 work** — if the Ladder line-search builds \f$H\f$ several
  times per "iteration", then our iteration is not CP2K's step and per-iteration is the WRONG metric; the
  honest one is total \f$H\f$-builds (or wall) to convergence;
- **or the term assembly asks more often than it needs to**, in which case it is a ~3× on the dominant
  bucket and would put the run BELOW CP2K.

⇒ **Distinguish them before optimising either way**: count \f$H\f$-builds per SCF step directly, and read
CP2K's own per-step \f$H\f$ count out of its log.  Cheap, and it decides whether bin 1 is finished.

---

## 5f. ★★★ WHERE THE 2.05× IS: **CALL COUNT, NOT KERNEL** (2026-09-05)

**Our kernel is already FASTER than CP2K's, per call.**  Read the two codes' own counters side by side —
CP2K's `T I M I N G` block (44 SCF steps) against our ledger (the parity probe, 20 iterations):

| | calls per SCF step | s per call | who wins the call |
|---|---|---|---|
| CP2K `integrate_v_rspace` | **2.00** (88/44) | 2.133 | |
| qchem gather | **4.0** (marginal, fixed-point) | **1.882** | **qchem, 0.88×** |
| CP2K `calculate_rho_elec` | **2.05** (90/44) | 2.029 | |
| qchem collocate | **2.0** (marginal, fixed-point) | **1.889** | **qchem, 0.93×** |

Both codes spend ~99% of the run in these two routines (CP2K 370 s of 374; qchem 334 s of 341).
⇒ **The whole 2.05× is that we call the gather TWICE as often as CP2K.**  Nothing in the kernel is behind.

**HOW THE COUNTS WERE TAKEN — differencing, not instrumentation** (`GPW_MNO_NMAX` = n against n+4, so the
setup/seed offset cancels and what is left is the marginal cost of one iteration):

| stage | gathers / iteration | collocations / iteration |
|---|---|---|
| Ladder (fixed point) | **4.0** | **2.0** |
| GDM (direct min) | 5.0 | 2.0, plus the LINE SEARCH's trials once it engages |

★ **AND THE LEDGER NAMES THE FOUR.**  A fixed-point iteration closes 2 `h ball` fields and 2 `h raw` fields.
On this route the raw adjoint is XC's (one field per spin) and the ball adjoint is Hartree's, so:

> **THE HARTREE MATRIX IS BUILT TWICE PER ITERATION** — once for the FOCK at \f$\rho_{mix}\f$
> (`Vee_Hartree::MakeMatrixT`) and once for the ENERGY at \f$\rho_{new}\f$
> (`Vee_Hartree::GetEnergy` → `0.5*cd->DM_Contract(this,cd)`).  Two different densities, so the gather memo
> cannot catch it.  ✅ Corroborated: with a pass-through mixer (\f$\alpha=1\f$, no Kerker) the two densities
> coincide on some iterations and the count falls 4.0 → 3.5.

**⇒ THREE LEVERS WERE ON THE TABLE.  ONE LANDED, ONE IS REFUTED, ONE REMAINS** (2026-09-05).  The
per-iteration budget was 4×1.882 + 2×1.889 = **11.3 s** against CP2K's 2×2.133 + 2×2.029 = **8.3 s**:

| | lever | verdict |
|---|---|---|
| **A** | **\f$E_H\f$ WITHOUT A MATRIX** — \f$E_H=\tfrac12\mathrm{Tr}(DV_H)=\tfrac12\Omega\sum_{\Delta G}\overline{\tilde\rho}V_H\f$, a G-space pairing on the fit ball.  Not an approximation: the gather is the exact adjoint of the collocation, so Parseval makes the two expressions equal by construction.  ★ **AND ONE MAP IS ENOUGH** — \f$V_H=k\tilde\rho\f$ with \f$k\f$ real, so \f$\overline{\tilde\rho}V_H=|V_H|^2/k\f$ and \f$\tilde\rho\f$ is never fetched (a second fetch costs a second IBZ star-average of the whole map: 16 ms/call on Si Γ) | ✅ **LANDED `9f4f4ae2`.  Parity probe 341.4 → 287.5 s CPU (−15.8%), gathers 95 → 66, 2.05× → 1.72×.**  Si 2×2×2 rows −19.5% / −11.4% (the same fix memoizes \f$V_H\f$ across irrep blocks).  806/806, every iteration count and \f$E_{tot}\f$ reproduced |
| **B** | **ONE GATHER PER SPIN: \f$\langle i|V_H+v_{xc}^\sigma|j\rangle\f$** — CP2K's `sum_up_and_integrate`, and the `CompositeExFunctional` argument one level up | ⛔ **REFUTED ON MEASUREMENT — the two fields are not on the same discretization.**  See below |
| **C** | the GDM LINE SEARCH's trial densities — **42 of the probe's 82 collocations**, i.e. ~28% of the run after A | ⏸ **PARKED AS A TODO UNTIL WE HAVE OT** (user, 2026-09-06).  See below |

⛔ **WHY B IS REFUTED (2026-09-05).**  Summing the two potentials needs them on ONE representation, and they
are deliberately on two: **\f$V_H\f$ is a BALL field** (band-limited, Poisson in G, gathered by the ball
adjoint's per-level \f$\{G\}\f$ restriction) while **\f$v_{xc}\f$ is a RAW RASTER field** (pointwise
nonlinear, gathered by the raw adjoint's per-level spectral BOX truncation).  Measured directly — build the
Hartree block both ways on Si Γ — they differ by **6e-5 relative**, so routing \f$V_H\f$ through the raw
adjoint is a different operator, not a re-association: it breaks the adjointness with the collocation that
made \f$\tilde\rho\f$, which is exactly the property lever A's Parseval identity rests on.
⇒ **B is available only if BOTH fields ride the ball route** — i.e. only if XC gives up the raw
\f$\rho_{DM}\ge0\f$ feed, which is the collapse-basin fix (`doc/GPWPlan` 0.5(f2)).  ⇒ **It becomes
available if N4 lands** (`doc/OpenWork.md`: make everything robust with \f$V_{xc}[\rho\ge0]\f$), and not
before.  Filed there, not here.

⏸ **WHY C IS PARKED, NOT PURSUED (user, 2026-09-06).**  *"CP2K uses the OT method instead of GDM so we
can't do proper parity timings against CP2K anyway."*  The line search is DIRECT MINIMISATION; CP2K's
counterpart for that is OT, which we do not have — and (checked in the decks and logs, §5a footnote ᵇ) the
benchmarked CP2K runs are not even using OT: they diagonalise and mix.  ⇒ There is no comparable step to
optimise C against.  **Trial-density cost becomes a real question when OT exists and can be timed against
GDM, minimiser to minimiser; until then it is a TODO, not a lever.**

⇒ **WHERE THAT LEAVES BIN 1.**  Measure the stage CP2K's decks actually match — our FIXED-POINT stage —
and it is **3 gathers + 2 collocations per iteration against CP2K's 2 + 2, i.e. 9.35 s against 8.29 =
1.13×** (§5a's fixed-point row, differenced 2026-09-06).  The term assembly is at its floor: 1 Hartree +
2 XC, and CP2K reaches 2 only by summing \f$V_H\f$ into \f$v_{xc}^\sigma\f$ before integrating.
⇒ **The whole remaining like-for-like gap is lever B, and lever B is behind N4.**

---

## 6. qchem ACCELERATIONS NOT IN CP2K — the qchem-vs-qchem deltas

These are what §2's knobs BUY.  They are excluded from a head-to-head row by rule 3c and reported here
instead, which is the more useful statement.

| acceleration | measured worth | provenance |
|---|---|---|
| Becke XC mesh (`QCHEM_BECKE_XC`) | ⚠ **NEGATIVE on MnO, and worse than it looked**: 1.8× on the SCF iteration (6.82 vs 3.72 s/iter) and **102× on setup** (184.3 vs 1.81 s — two mesh builds at 68.3 s each).  It buys atom-centred accuracy, not speed | §5a |
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
| 3 | **`sum_up_and_integrate` — ONE integrate per spin.**  CP2K adds \f$V_H\f$ and \f$v_{xc}^\sigma\f$ into one real-space potential and integrates it ONCE per spin: 88 `integrate_v_rspace` calls over 44 steps.  We call the gather 4× per iteration (2 Hartree + 2 XC).  ⚠ This is the whole of the parity row's 2.05×, since our per-call gather is 0.88× theirs | the parity row (§5f) | 09-05 |
| 2 | **k-INDEPENDENT per-step cost** — their per-iteration time is flat from Γ to 8 k (0.383 → 0.385 s SCF-only) where ours rises **6.3×** (0.051 → 0.323).  Not a "feature" so much as a consequence of collocating the k-summed density once per step | most of the Γ lead (§5a, §5e) | 09-04 |

⚠ Item 1 is a REAL algorithmic advantage they hold, not a deviation to switch off; item 2 is a gap of ours
with a diagnosed cause.  ⇒ **Neither belongs in `CP2K_COMPAT`** — that switch turns OUR accelerations off,
and turning theirs off is not available to us.

⚠ CP2K also runs `Wavefunction type COMPLEX` on BOTH Si k-rows, including the all-TRIM Γ-centred one where
we use REAL blocks.  So on that row we hold an advantage they decline to take — and are still only at
parity per iteration.

---

## 7. THE ROWS — threaded (12 cores), AND THE OMP GAP ANALYSIS

✅ **FILLED 2026-09-06, BOTH SIDES.**  The cross-code half was the thing missing since 09-04; it is here
now, and it inverts the question that motivated it.  Same binaries and decks as §5a, 12 cores each,
\f$E_{tot}\f$ agreeing to all printed digits with the serial runs on every row.

### 7a. Our side: 12 threads against our own serial

| row | serial wall | **12-thread wall** | **speedup** | 12t CPU / serial | RSS ser / 12t |
|---|---|---|---|---|---|
| NaF SR2 Γ | 23.0 s | **10.96 s** | **2.10×** | 1.41× | 55 / 73 MB |
| MnO `QCHEM_BECKE_XC=0` | 147.0 s | **49.5 s** | **2.97×** | 1.50× | 113 / 197 MB |
| MnO `CP2K_COMPAT=1` probe | 287.5 s | **61.6 s** | **4.67×** | 1.58× | 132 / 210 MB |
| MnO ALL DEFAULTS | 395.6 s | **128.4 s** | **3.08×** | 1.51× | 476 / 563 MB |

★ **The parity route threads BEST (4.67×)** — it is nothing but the box walk, which is our OpenMP region.
The Becke routes carry more serial residue, which 7c breaks down.

### 7b. ★★★ CP2K's side — AND THE HEADLINE IS THAT ITS PARALLELISM DOES NOT HELP HERE

| CP2K row | serial wall | 12 OMP threads | 12 MPI ranks | best speedup |
|---|---|---|---|---|
| Si Γ | 5.18 s | 5.49 s (**0.94×**) | — | **none** |
| NaF SR2 Γ | 7.37 s | 8.95 s (**0.82×**) | — | **none** |
| MnO AFM-II VA | 374.2 s | 342.7 s (**1.09×**) | **259.6 s (1.44×)** | **1.44×** |

⚠ **These are CP2K's own thread counts, not a knob we set and hoped**: its banner reports
`GLOBAL| Number of threads for this process 12`, and the run still measured **196% CPU** — ~2 cores of
work spread over 12 threads.  On the two small cells it is *slower* than serial while billing 2× the CPU.
The MPI arm is the honest one for CP2K (its design centre), and 12 ranks buy **1.44× for 8.3× the CPU**
(3107 s against 374 s serial).

★ **AND ITS OWN TIMING BLOCK SAYS WHERE**: on MnO, the two GPW hot routines — which are **98% of the run**
— barely move between 1 and 12 threads:

| CP2K routine (SELF time) | serial | 12 threads | speedup |
|---|---|---|---|
| `grid_collocate_task_list` | 182.6 s | 167.9 s | **1.09×** |
| `grid_integrate_task_list` | 187.7 s | 167.2 s | **1.12×** |

⇒ **It is precisely the GPW collocate/integrate route whose OMP axis does not scale here** — the same two
routines whose per-call cost we beat by 0.88×/0.93× serially (§5f).  Two readings, and this measurement
does not separate them: (a) the route's OpenMP was never a priority — CP2K's research centre is hundreds
of small molecules in a box, where the parallel axis that matters is MPI over molecules/atoms, not OMP
inside a solid's task list (user, 2026-09-06); (b) a 4-atom cell simply has too few tasks per level to
spread.  ▶ **The discriminator is one deck we do not have**: the same run on a supercell (say 2×2×2 MnO,
32 atoms).  If its grid routines start scaling, it was size; if they do not, it was the code.  ⚠ Until
then, quote the measurement, not either explanation.

⇒ **CROSS-CODE, AT 12 CORES, ON MnO** — per SCF step, since the iteration counts differ (rule 3d):
**qchem 4.14 s/iteration against CP2K's best 5.90 s (0.70×)**, and we get there on **596 s of CPU against
their 3107 s (0.19×)**.  ⚠ Rule 3c applies: that qchem row runs our accelerations.  ⚠ And the standing §2
caveat is the whole explanation — this is a 4-atom, high-symmetry cell, the regime that most favours a
symmetry-exploiting code and least favours CP2K's distribute-over-atoms design.  **On a 100-water box the
verdict would invert**, and none of this says anything about many-node scaling.

⇒ **SO: "ARE WE MISSING PARALLEL OPPORTUNITIES CP2K EXPLOITS?"  NOT ON THIS CLASS OF SYSTEM.**  We
out-scale CP2K on both of its axes here.  The opportunities we are missing are our OWN, and 7c names them.

### 7c. Where OUR threaded time actually goes — bucket by bucket

MnO ALL DEFAULTS, the same run in both columns (serial 395.6 s → 12 threads 128.4 s):

| bucket | serial | 12 threads | speedup | |
|---|---|---|---|---|
| setup: Becke mesh build | 138.2 s | 16.9 s | **8.2×** | ✅ |
| scf: XC-mesh ρ sampling (matrix-free) | 74.9 s | 11.3 s | **6.6×** | ✅ |
| scf: collocate (pair scatter) | 44.8 s | 8.3 s | **5.4×** | ✅ |
| setup: XC-mesh Φ tables | 44.3 s | 6.8 s | **6.5×** | ✅ |
| scf: integrate-back (pair gather) | 16.4 s | 2.7 s | **6.2×** | ✅ |
| **scf: XC-mesh quadrature H_xc** | 11.2 s | **9.2 s** | **1.21×** | ⛔ |
| **scf: E_H V_H field build** | 7.6 s | **8.5 s** | **0.89×** | ⛔ |
| scf/setup: the FFT closures, local-PP short | ~3 s | ~3.4 s | ~0.9× | ⛔ |
| **everything not in a bucket** | ~52 s | **~59 s** | **~0.9×** | ⛔ |

⚠ **THE LAST ROW IS SUPERSEDED — IT IS NOW 0.03 s** (2026-09-06, `doc/ParallelAndOraclePlan.md` 1.1 + 1.1(a)).
Instrumenting it did not find "diagonalisation, orthogonalisation, mixing, the fit solves": the SCF's whole
linear-algebra side is ~5 s of a serial run and the diagonalisation is **0.038 s**.  It found **the
Hamiltonian being built once per ANNEAL STAGE** — `SolidCalculation::BuildStage` re-`Factory`s it and had no
bucket at all.  On the same 12-thread MnO row (121.0 s wall, 39 buckets summing to 120.98 s):

| bucket | GPW_OMP unset | 12 threads | speedup | |
|---|---|---|---|---|
| `setup: hamiltonian ctor` (exclusive of the two below) | 51.4 s `[×2, 25.7/call]` | **46.4 s** `[×2, 23.2/call]` | **1.11×** | ⛔ |
| ⤷ `setup: XC-mesh Φ tables` | 54.1 s `[×2]` | 6.2 s `[×2]` | **8.7×** | ✅ |
| ⤷ `setup: becke mesh build` | 15.4 s `[×2]` | 15.8 s `[×2]` | **1.0×** | ⚠ see below |
| **Hamiltonian construction, all in** | **121.0 s = 39%** | **68.4 s = 56%** | 1.77× | ⛔ |
| the three residue buckets together | 0.36 s | 0.36 s | — | — |
| **everything not in a bucket** | **0.05 s** | **0.03 s** | — | ✅ |

Same build, same recipe, `Etot=-61.40297529` on both arms to all printed digits; 313.9 s wall / 532 s CPU /
**169%** unset, 121.0 s / 569.5 s CPU / 469% at 12 threads.  39 buckets summing to 313.71 of 313.76 s and
120.98 of 121.01 s respectively.

★ **SO THE SHAPE OF THIS ROW IS SETUP, NOT SCF — 69.5 s of setup against 51.5 s of SCF at 12 threads.**
That is bin 2 ("Init (pre-iteration) time") in the gap-close priority order, and it is now the biggest
single lever on the MnO wall — bigger than either pair loop.  **The half that does not thread is the ctor's
own 1.11×**; its children scale fine, which is precisely why it hid.  ▶ Next: bucket the fit-basis half of
that 23.2 s/call, and establish whether an anneal stage must rebuild the Hamiltonian at all (~34 s if not).

★★★ **AND THEN THE CTOR WAS OPENED UP, AND THE ROW MOVED — 2026-09-06, `ParallelAndOraclePlan.md` 1.1(b).**
Six more buckets found that the 23.75 s/call was **two calls to one bad index**: the run folds the same
~97k mesh points twice per Hamiltonian (the orbit-consistency filter in `UnitCell::CreateIntegrationMesh`,
then `FoldMesh` in `GPW_IBS::CreateXCQuadrature`), 9.6 + 9.5 s/call.  `TorusIndex` bucketed EVERY mesh on a
constant 64³ grid — 0.37 points per bucket on average, which is the wrong statistic for a clustered
atom-centred radial mesh.  Grid made as fine as the tolerance allows + centre-bucket-first probe:

| | GPW_OMP unset | 12 threads |
|---|---|---|
| orbit-consistency fold | 9.61 → **0.163** s/call | 9.61 → **0.193** s/call |
| `FoldMesh` | 9.47 → **0.151** s/call | 9.47 → **0.190** s/call |
| Hamiltonian construction, all in | — | 34.4 → **15.5** s/call |
| **WHOLE RUN** | 313.9 → **237.2 s** (**1.32×**) | 120.8 → **83.2 s** (**1.45×**) |

⚠ The "GPW_OMP unset" column is a MIXED arm, not a serial one (see the protocol note below) — both sides
of its 1.32× were taken the same way, so the ratio is sound, but do not read it as a serial figure.
| setup / SCF split | 71.8 / 165.2 s | 32.1 / 51.1 s |
| peak RSS | 496 → 506 MB | 562 → 561 MB |

`Etot = -61.40297529` on every arm before and after, to all printed digits — it is a data structure, not a
numerical method.  814/814 green.  The largest setup bucket is now the Becke mesh build (8.15 s/call), then
the site-adapted angular sets (4.03 s/call, never measured before).

⚠ **THE MnO ROWS IN §5a AND §7a ARE STALE, AND ARE NOT RE-TAKEN PER INCREMENT** (user, 2026-09-06).  A row
is a CROSS-CODE artefact: it costs a CP2K arm as well as ours, and re-taking one after every speedup both
burns the box and lets each increment mask the last.  Same discipline as the anchor-moving sprint
(`doc/OpenWork.md` item **S**): **let the wins accumulate, then re-bank the rows ONCE**, at the end of
Phase 1 — after 1.2 (BLAS-mode arm) and 1.3 (the 2×6), which move them again.  Until then read the ratios
in §5a/§7a as a floor on our side and quote 1.1(b)'s numbers for anything that turns on the MnO wall.

✅ **THE BECKE MESH BUILD "COLLISION" IS RESOLVED, AND THE BANKED ROW WAS RIGHT — MY ARM WAS WRONG**
(2026-09-06).  The flag raised here said the mesh build measured the same in both arms (15.4 vs 15.8 s)
where §7c banks 138.2 → 16.9 s (8.2×).  One run at `GPW_OMP_THREADS=1` settles it:

| arm | CPU% | `setup: becke mesh build` |
|---|---|---|
| `GPW_OMP_THREADS=1` | **99%** | **68.74 s/call** |
| `GPW_OMP_THREADS` unset | 188% | 7.10 s/call |
| `GPW_OMP_THREADS=12` | 469% | 8.15 s/call |

⇒ The mesh build threads at **8.4×**, §7c's 8.2× stands, and §5e's open question is closed.

★★★ **THE REAL FINDING IS A PROTOCOL DEFECT: `GPW_OMP_THREADS` UNSET IS NOT A SERIAL ARM.**  The Becke
build is *"parallel by DEFAULT under QCHEM_OPENMP; `GPW_OMP_THREADS`, when set, is honoured as the thread
CAP"* (`UnitCell.C`) — so UNSET means **all cores** for that loop while the GPW pair loops stay serial.
An "unset" row is therefore a MIXED arm, not a serial one, and its 169–188% CPU says so out loud.  This is
§4's own *"a knob is not a measurement"* rule biting the person who wrote it down.
▶ **A serial row MUST set `GPW_OMP_THREADS=1`** (99% CPU is the check).  Any row in this file whose serial
arm reads much above 100% CPU was taken the mixed way and should be re-taken at the Phase-1 re-bank.

**MnO ALL DEFAULTS, post-1.1(b), `Etot=-61.40297529` on all three arms:**

| arm | wall | vs serial |
|---|---|---|
| `GPW_OMP_THREADS=1` (true serial) | **361.2 s** | 1.00× |
| unset (mixed) | 237.2 s | 1.52× |
| `=12` | **83.2 s** | **4.34×** |

★ **THE BECKE MESH BUILD THREADS AT 8.2× — §5e's open question is CLOSED, and the answer is "it threads
fine on MnO".**  Whatever holds NaF back (7d) is not the loop being serial.

⛔ **THE THREE THINGS THAT DO NOT SCALE, in priority order:**
1. **~59 s of UNBUCKETED work** — the largest single block of a threaded run, and it is invisible because
   nothing times it: diagonalisation, orthogonalisation, mixing, the fit solves, the SCF bookkeeping.
   ▶ **The first action is an instrument, not an optimisation**: bucket the SCF's non-GPW work.
2. **The XC-mesh quadrature \f$H_{xc}=\Phi^\dagger\,\mathrm{diag}(w\,v)\,\Phi\f$ at 1.21×** — and this
   one is a DELIBERATE TRADE, not an oversight (user, 2026-09-06; the reasoning is written into
   `DeltaFit_IBS::AdjointT`).  **Our parallelism lives ABOVE the linear algebra** — per k-block / irrep /
   spin — with BLAS pinned to one thread (`qchem::PinBlasToOneThread`), precisely to avoid OMP nesting.
   The measurement behind it: one dispatched whole-matrix `zgemm` runs **34.1 GFlop/s against 1.87 for
   ANY blocked or viewed form**, so hand-blocking to spread over threads loses 13× to gain 8×.
   ⇒ **The bucket's 1.21× is the PRICE of that policy, and the thing to question is the WIDTH OF THE LEVEL
   ABOVE, not the policy.**  On MnO at Γ that level is 2 spins × 1 k-block = **2-way**, so ten of twelve
   cores have nothing to do in this bucket by construction.  On the 8-k Si rows the same policy has 16-way
   width and nothing is left on the table.
   ▶ **AND THERE IS A CHEAP TEST WORTH RUNNING**: the "one dispatched zgemm" arm needs
   `QCHEM_BLAZE_BLAS=ON`, which was correctly refused on 2026-08-15 because the system BLAS was then
   NETLIB.  ⚠ It is not any more — `libblas`/`liblapack` now resolve to `openblas-pthread` — so the 34
   GFlop/s path is available again, **single-threaded, composing with the pin and with the levels above**.
   That is a rebuild and one run.
3. **The \f$V_H\f$ field build at 0.89×** — that bucket is 6.6% of the threaded wall and it is mostly
   `SymmetrizeGMap`, the IBZ star-average over the point group (48 ops), which is a serial walk over a
   \f$\{G\}\f$ map.  It got LOUDER, not quieter, when §5f's lever A removed the gathers around it.

### 7d. The Amdahl residual, and the NaF anomaly that is still open

Solve \f$S+P/12=t_{12}\f$ against \f$S+P=t_1\f$:

| row | Amdahl \f$S\f$ | as % of serial | what 7c attributes it to |
|---|---|---|---|
| MnO ALL DEFAULTS | **104 s** | 26% | ⚠ **re-attributed 2026-09-06**: the ~59 s "unbucketed" is the **Hamiltonian ctor's 46 s exclusive half** (×2, once per anneal stage) + ~21 s of non-scaling buckets + imbalance |
| MnO `BECKE_XC=0` | 32 s | 22% | same shape, no mesh build |
| NaF SR2 Γ | 9 s | 40% | ⚠ **matches its SERIAL setup buckets (9.8 s) to ~7%** |

⚠ **NaF REMAINS THE ANOMALY.**  Its inferred serial fraction still equals its setup, yet the very same
setup buckets thread at 6–8× on MnO.  The likeliest reading is now SIZE, not structure: NaF is a 2-atom
cell whose mesh build is 7 s and whose Φ tables are 0.35 s — too little work per thread to amortise the
regions.  ▶ Cheap to settle: the same bucket table as 7c, taken on NaF.

⚠ **CPU inflation is 1.41–1.58×** on our side even with the barrier spin gone — thread management, load
imbalance, memory bandwidth.  CP2K's is 1.8× threaded and **8.3× under MPI**, so this is not a gap against
them; it is a cost of the last speedup increment.

---

## 8. WHAT IS STILL MISSING

★ **AS OF 2026-09-04, in the user's priority order** (`doc/OpenWork.md` carries the same list as the
forward queue):
0. ✅ **THREADED ROWS — DONE 2026-09-06, both sides (§7).**  CP2K at 12 threads was the missing piece; it
   turns out its parallelism does not help on these cells (0.82–1.09× on OMP, 1.44× on 12 MPI ranks), so
   the remaining threading work is entirely OUR OWN: ~59 s of unbucketed serial SCF work, a dense XC GEMM
   on blaze's serial kernels, and the IBZ star-average.
1. **CLOSE BIN 1 AND BIN 2 SINGLE-THREADED** — per-iteration CPU and pre-SCF setup, against CP2K.  ✅ §5a
   is FILLED as of 2026-09-05 (nine rows, both codes serial, setup split out): bin 1 is ahead on seven of
   nine and comes down to the `CP2K_COMPAT=1` row at 2.05×; bin 2 is the Becke mesh build.  The whole-run
   table is NOT informative (mixed thread states, 3× different iteration counts) and is kept only for
   \f$E_{tot}\f$ and RSS.
2. ✅ **THE BUSY-WAIT BARRIER** — fixed, and everything re-run at 12 threads (§7).  Cause was identified 2026-09-04:
   nothing in the tree sets `KMP_BLOCKTIME` or `OMP_WAIT_POLICY`, and LLVM's libomp spins **200 ms** after
   every parallel region — with per-shell-pair regions that is mostly spin.  It matches the inflation
   already recorded: **663 s threaded CPU against 500 s serial for the same work**.
3. ✅ **The MnO FM row** (footnote ⁵) — re-taken 2026-09-05: 411 s CPU / 474 MB, 0.83× CP2K per SCF step.
4. ✅ **NaF's CP2K iteration counts** — read off its own logs (SR2 16 steps, full-SR 27); both rows now have
   a per-iteration column in §5a.
5. **An attribution A/B** — the 09-04 deltas are "current vs banked" across everything since 08-28, with
   no parent-commit A/B on the VA recipe.

⏸ **PARKED, deliberately** (user, 2026-09-04: *"defocusing"*): the
\f$\lVert V_{xc}-V_{xc}^{fit}\rVert\f$ fit-quality study.  It stays worth doing one day — every
Becke-vs-uniform cost number here is taken at UNKNOWN-EQUAL accuracy, so §6's "Becke is a negative
acceleration" is a COST statement and not a verdict — but it is not what bins 1 and 2 need.




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