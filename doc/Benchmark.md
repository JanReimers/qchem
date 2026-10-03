# The qchem ↔ CP2K head-to-head — runtime, RAM, energy

**This is an INSTRUMENT, not a report.**  A row here is a claim that two codes did the same work on the
same hardware; the sections below are what makes that claim checkable.  Live content only: the rules, the
run commands, the one table (§5a), the open levers (§10).  **Everything else — the §5–§9 narrative, every
wrong turn and superseded number — is archived verbatim in `doc/Records/BenchmarkHistory.md`, section
"Benchmark.md verbatim as of 2026-10-01"** (headings `## N.` inside it).

**Read in this order:** §1 the process → §2 what parity means → §3 the rules → §4 the commands →
**§5a the BIN-1 TABLE** → §10 open levers.  **Copy the run command from §4 / `scripts/retake5a`; never reconstruct it.**

---
## 1. THE SYSTEMATIC PROCESS — single thread first, then threads; and the four bins

### 1a. Single-thread parity FIRST, then threads — in that order (user)

A threaded comparison against a serial code measures algorithm AND parallel efficiency at once.  (1) Get the
SINGLE-THREAD time in line first: `OMP_NUM_THREADS=1 QCHEM_OPENMP_THREADS=1` (pins our OpenMP regions and the BLAS).
(2) THEN look for OMP-shaped gaps.  ⚠ Our OpenMP threads BUSY-WAIT at the barrier, so a threaded run bills far
more CPU than it uses (294 s serial billed ~590 s at 16 threads) — the CPU column overstates us wherever we thread.

### 1b. The four gap-close bins (user, 2026-08-28) — priority order

> *"1) Per iteration CPU time, 2) Init (pre iteration) time, 3) top RAM usage (mostly solved I think),
> 4) very roughly match the # of iterations.  I think items that fall into bin 4 can just be documented."*

| bin | the axis | where it stands (MnO AFM-II VA, 2026-09-04) |
|---|---|---|
| **1** | **per-iteration CPU** | ★ **PER SCF ITERATION WE ARE AHEAD ON 7 OF 9 ROWS** (§5a, 2026-09-05, both codes pinned serial): MnO ALL DEFAULTS **0.82×**, MnO FM **0.77×**, `QCHEM_BECKE_XC=0` **0.45×**, Si **0.14× / 0.66× / 0.56×**, NaF full-SR **0.18×**.  The two losses are NaF SR2 **1.45×** (a Becke cost) and — the one that counts — `CP2K_COMPAT=1` **1.72×** (was 2.05× before §5f's lever A).  ⇒ **NOT closed** — but on the SAME ALGORITHM (our fixed-point stage against CP2K's diagonalise-and-mix decks) it is **1.13×**, three gathers against their two, and the third gather is §5f lever B, behind N4.  The 1.72× two-stage figure adds a GDM line search CP2K's runs have no counterpart for |
| **2** | **init / pre-iteration** | ⛔ **PROMOTED, and now measured serially: MnO's setup is 184 s = 47% of the default run against CP2K's 8.1 s (23×)**, of which **136.6 s is two Becke mesh builds** (§5a).  The earlier "~57 s of 328 s" was a THREADED ledger bucket — the build's partition loop is `omp parallel for`.  ✅ With `QCHEM_BECKE_XC=0` our setup is **1.76 s and beats CP2K's 8.1 s** — so bin 2 is a Becke-mesh question, exclusively |
| **3** | **peak RAM** | ✅ **solved, and we win**: ~476 MB defaults, **113–132 MB on the parity routes against CP2K's 217 MB** |
| **4** | **iteration count** | ⚠ **RE-JUDGED 2026-09-20 under rule 3f** — the old "31 (defaults) / capped (parity) against CP2K's 44" compared a 1e-5 G-space residual with CP2K's 1e-6 max\|ΔP\|, i.e. two different questions.  On the SAME measure, threshold and loop shape (one density-side history, no Fock DIIS, no MOM): **imposed 22 / 24, free (`CP2K_COMPAT=1`) 33 / 37 against CP2K's 44 / 104** (U=0 / U=4 eV, §5 rows ⁸).  ⇒ CLOSED in our favour; what remains is bin 1's third gather |

The live tracker for these is `doc/OpenWork.md` (accuracy/features) and §10 below (perf levers); this file holds the measurements behind them.

---

## 2. WHAT `CP2K_COMPAT=1` ENCOMPASSES — ⚠ AN EMERGING LIST, NOT A FINISHED ONE

**Every deviation found so far is an ACCELERATION, not physics**: turning them all off moves the MnO total
by **3e-8 Ha** (agreeing to 10 s.f.).  That is the property the switch most needed to demonstrate about
itself, and it is measured rather than asserted (history §2).

⚠ **THIS LIST HAS GROWN EVERY TIME SOMEONE LOOKED.**  It started at four items; it is eight (the eighth, +U's
form, is the first that is physics rather than an acceleration — inert unless a Hubbard manifold is on the run,
so every banked row below is untouched by it).  Assume it is still incomplete — a parity row is only as honest as the last thing we noticed we were doing differently.

| # | knob | what qchem does that CP2K does not | found |
|---|---|---|---|
| 1 | `QCHEM_DM_LOWRANK` | factored/low-rank ρ (\f$D=LL^\dagger\f$) — a SINGLES route; CP2K collocates PAIRS | 08-25 |
| 2 | `GPW_STREAM_FOLD` | orbit fold on the collocation pair streams (5.2× on MnO's pair count) | 08-25 |
| 3 | `QCHEM_MIX_RHO_M` | (ρ,m) mixing channels instead of (up,dn) | 08-25 |
| 4 | `QCHEM_XC_DM_SOURCE` | \f$V_{xc}\f$ fed ρ[D] wholesale instead of ρ_mix | 08-25 |
| 5 | `QCHEM_IMPOSE_SYMMETRY` | space-group imposition: BZ fold + ρ star-average + site-adapted XC mesh.  ⚠ **OVERRULES the caller** | 08-26 |
| 6 | `QCHEM_BECKE_XC` | atom-centred (Becke) XC quadrature instead of the uniform grid.  ⚠ **OVERRULES the caller**; it was **43% of the MnO row** | 08-28 |
| 7 | `GPW_DAWARE_SCREEN` | D-aware collocation box tolerance \f$\varepsilon/|c_{ij}|\f$ instead of flat \f$\varepsilon\f$ | 09-04 |
| 8 | `QCHEM_U_EIGEN` | **the first PHYSICS deviation, not an acceleration** (user ruled it onto this list 2026-09-20): DFT+U evaluated on the Löwdin block's EIGENVALUES (Dudarev, rotationally invariant) where CP2K keeps only its DIAGONAL POPULATIONS (`dft_plus_u.F`, `IF (isgf == jsgf)`).  Inert on a run without a Hubbard manifold; on MnO VA at U=4 eV it is 0.08 vs 0.61 Ha of E_U | 09-20 |

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
[<label> run] symmetry: IMPOSED (Shubnikov from the decoration);  threads: OMP_NUM_THREADS=1 QCHEM_OPENMP_THREADS=1 (BLAS pinned to 1)
[<label> run] CP2K_COMPAT=0 -> DEVIATING;  QCHEM_DM_LOWRANK=on*  GPW_STREAM_FOLD=on*  QCHEM_MIX_RHO_M=off  QCHEM_XC_DM_SOURCE=off   [* = differs from CP2K]
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

### 3f. An ITERATION COUNT is comparable only on the SAME CONVERGENCE MEASURE at the SAME THRESHOLD (user, 2026-09-20)

*"Anytime I look at top and CP2K finishes many minutes before ITMain, my antenna goes up."*  It should: the two
codes were converging DIFFERENT QUANTITIES.  CP2K's `EPS_SCF` is **max|P_out − P_in| over the AO density-matrix
elements**, per spin, un-normalised (`qs_scf_loop_utils.F`, `self_consistency_check`), typically 1e-6 (MnO deck)
or 1e-7.  Our `MinΔρ` on a Kerker/Pulay recipe is the MIXER's residual **max|ρ̃_out(G) − ρ̃_in(G)|** in G-space
(the largest coefficient is N/Ω ≈ 0.08 on MnO, so 1e-5 there is ~1e-4 relative — **~100× looser** than the deck),
and on a linear-D recipe it is a third thing (Σ_blocks‖ΔD‖_F/N_e).  Measured on Si Γ (2026-09-20): the Kerker gate
"converges" in **5** iterations at 1e-3 on the residual (E = −7.114894); CP2K's measure at 1e-6 needs **27**
with Pulay(8) and lands on −7.115068 — 0.17 mHa lower, the true anchor.  On MnO +U the honest count was not
"55 vs 104".  ⇒ `SCFParams::Measure::MaxΔD` puts CP2K's measure on our loop (successive D_out, max over blocks
and spins); **a row that quotes an iteration count, a wall, or a "converges in" sets `Δρmeasure=MaxΔD`,
`MinΔρ=EPS_SCF` and `NMaxIter=MAX_SCF` off the deck** — and states it.  ⚠ Every iteration-count column in §5
written before this rule compares a looser qchem criterion with a tighter CP2K one; re-take before quoting.

---

## 4. HOW TO PRODUCE A ROW — the same wrapper for both codes

⚠ **THE `GPW_SCF.*` FILTERS IN THIS SECTION NAME RETIRED TESTS** (the 2026-09-15 test-suite renaming, rule 3a's
own failure mode).  The live rows and their current filters are in **`scripts/retake5a`** — `GPW_Si.Γ_CP2K`,
`GPW_Si.k222_CP2K`, `GPW_Si.k222s_Imp_CP2K`, `GPW_NaF.Γ_Imp_Anchor`, and `gpwprobe mno` for MnO — which also
sets the rule-3f knobs (`GPW_MEASURE=maxdd GPW_EPS GPW_NMAX GPW_ACC=null GPW_MOM=0 GPW_KERKER_G0=1 GPW_PULAY=8`,
all read by the harness's `EnvOverrides`).  Run `scripts/retake5a` to reproduce the rule-3f table in §5a.

```bash
scripts/bench "Si Gamma qchem" -- build/Release/IntegrationTests/ITMain --gtest_filter=GPW_SCF.SiliconGammaConverges
scripts/bench "Si Gamma cp2k"  -- cp2k -i IntegrationTests/CP2K/si_fcc_gpw.inp
```
**Peak RAM is measured from OUTSIDE, identically for both codes** (`scripts/bench` → `/usr/bin/time -v` "Maximum
resident set size" = `VmHWM`); on a qchem row it also prints qchem's own `VmHWM` as a cross-check.  `PEAK RSS` is
a PROCESS watermark: one config per process (`MNO_SKIP_FM` / `MNO_SKIP_AFM` for MnO).  `QCHEM_OPENMP_THREADS` governs
the GPW pair loops only; BLAS routing is `QCHEM_BLAZE_BLAS`.  qchem has no `GLOBAL| Number of threads` banner, so
`scripts/bench` reports measured CPU% and CPU seconds instead of trusting a knob.  Every GPW run prints Etot at
10 s.f., wall + per-bucket ledger (`GPW_REPORT=1`), PEAK RSS, `[fold]` lines, and `[t=…s]` run-clock stamps.

⚠ **A STAMP IS A MOMENT, NOT A DURATION — a gap belongs to what ran BEFORE the stamped line.**  Read the stamps
before the ledger (the gaps between stamps are the time no bucket charges — that located 25 s of unbucketed MnO
time); the LEDGER prices an item, the stamps only localise WHEN.  A number's precision must not depend on which
other diagnostics are on (`std::defaultfloat` restores the format flag, not the precision).

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
| `[<label>] site moments (Becke-partitioned …) [e]: 0:+… 1:-… net=…` (polarized, atom-centred XC mesh) | the end-of-run result line, always (was `QCHEM_SITE_MOMENTS=1`, retired 2026-10-03); the json `scf.siteMoments` holds the last iteration's |
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
the watermark of the whole process.  `QCHEM_OPENMP_THREADS` governs the GPW **pair loops only** and says nothing
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

Verify each side reproduces its own history before reading a Δ: the CP2K decks against `doc/Records/CP2Kresults.md`
(all five re-validated 2026-08-19, `doc/Records/CP2KBuild.md`), and the qchem runs against the tests' own anchors.


---

## 5a. ★★★ THE BIN-1 TABLE — per-ITERATION CPU, qchem vs CP2K, one thread each

**ONE TABLE, ON PURPOSE** (user, 2026-09-05: *"the important comparison table was the only table in
section 5a, but now the section is flooded with other tables so I cannot easily point to it"*).  Every A/B
delta, ledger split and open question that used to sit here is in **§5e**.  This is the bin-1 instrument
and the only table to quote for bin 1.

**THREAD STATE — SERIAL ON BOTH SIDES, AND MEASURED SO** (rule 3b): qchem `OMP_NUM_THREADS=1
QCHEM_OPENMP_THREADS=1`, measured **99% CPU** on every row; CP2K `OMP_NUM_THREADS=1`, measured 97–99%.
**Taken 2026-09-05**, one box (14 GB, 16 cores), qchem at `9f4f4ae2` built `-O3 -march=native`, CP2K 2025.2,
commands copied from §4.  ★ **Every row RE-TAKEN after §5f's lever A** (the Hartree energy stopped building
a matrix): the parity row fell 2.05× → 1.72×, the Si 8-k rows 0.84×/0.68× → 0.66×/0.56×, MnO defaults
0.90× → 0.82×, `BECKE_XC=0` 0.53× → 0.45×.  ⚠ **The qchem `setup` column is much larger than the 09-04 entry's ~57 s for MnO
because the Becke mesh build THREADS** (`src/Structure/Imp/UnitCell.C:273` — the partition loop is
`#pragma omp parallel for` over quadrature points, and `QCHEM_OPENMP_THREADS=1` pins it): 68.3 s serial per
build against the 16.7 s that ledger read with threads free.  The serial figure is the one that belongs
beside a serial CP2K row (`doc/Records/BenchmarkHistory.md` §9).

**HOW TO READ IT.**  `CPU` is whole-run user+sys.  `setup` is that run's own pre-SCF work — qchem: the sum
of the ledger's `setup:` buckets; CP2K: total CPU minus the sum of its printed per-step times.  So

- **`s/it (SCF)` = (CPU − setup) / iterations is the BIN-1 number**, and the `×` column is its ratio;
- **`setup` is the BIN-2 number**, in the same row, on the same run.

⇒ **Compare the SCF columns.**  The total column is kept only because it is what the whole-run table (§5)
divides — and the two disagree by 1.8× on MnO precisely because 44% of that run is setup.

⚠ **THE `q iters` / `c steps` COLUMNS BELOW ARE RULE-3f-STALE (2026-09-20)**: every qchem count was taken at
1e-3–1e-5 on the MIXER RESIDUAL (a G-space quantity), every CP2K count at `EPS_SCF` on max|ΔP| — not the same
measure, and on MnO ~100× looser on our side.  The **s/it** columns survive (a per-call cost does not care when
the loop stops); the counts and the `total ×` column do not — read them as "how many steps that recipe took to
ITS OWN criterion", never as a convergence-rate comparison.  The re-taken, like-for-like counts are the two
rule-3f MnO rows in §5 (footnote ⁸: 33 / 37 vs 44 / 104) and the Si example in rule 3f itself (5 vs 27 to reach
the same anchor).  A full re-take of this table on `Measure::MaxΔD` is queued in §10 below.

| row | span / k | q iters ⚠3f | q CPU | q setup | **q s/it (SCF)** | c steps ⚠3f | c CPU | c setup | **c s/it (SCF)** | **BIN 1 ×** | q s/it (total) | total × ⚠3f |
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

**★ THE RULE-3f RE-TAKE (2026-09-20, `scripts/retake5a` + the two `gpwprobe mno` rows of footnote ⁸) — the SAME
convergence measure (max|ΔD| / max|ΔP|) at the DECK'S threshold, the deck's loop shape (one density-side history:
Kerker G0=1 + Pulay 8; nothing on the Fock side; no MOM), `CP2K_COMPAT=1` (no Becke mesh, no imposition), both
codes serial.**  Whole-run wall and RSS here are like-for-like BECAUSE the loops now stop at the same place;
CP2K's numbers are the banked ones (its counts were already at `EPS_SCF`).  The one stated difference: the Si
decks mix P directly (α 0.4) where our stable no-Fock-side route is Kerker + Pulay (§10 row "Linear D-mixing
DIVERGES").

| row | EPS / MAX_SCF | **q iters** | q wall | q RSS | **c steps** | c wall | c RSS | iters q/c | qchem E | CP2K E |
|---|---|---|---|---|---|---|---|---|---|---|
| Si Γ (free) | 1e-7 / 60 | **16** | 2.7 s | 39 MB | 12 | 5.0 s | 148 MB | 1.33× | −7.115067447 | −7.115057882 |
| Si 2×2×2 Γ-centred (free) | 1e-7 / 60 | **15** | 8.6 s | 54 MB | 13 | 5.6 s | 153 MB | 1.15× | −7.778472674 | −7.778457865 |
| Si 2×2×2 shifted MP | 1e-7 / 60 | **14** | 10.6 s | 55 MB | 14 | 5.9 s | 153 MB | 1.00× | −7.867452508 | −7.867436530 |
| NaF SR2 Γ | 1e-6 / 200 | **16** | 7.8 s | 84 MB | 16 | 7.2 s | 173 MB | 1.00× | −24.43039482 | −24.431213375 |
| MnO AFM-II (free) | 1e-6 / 200 | **33** | 5m42s | 262 MB | 44 | 6m14s | 217 MB | 0.75× | −61.41154311 | −61.303325178 |
| MnO AFM-II +U 4 eV (free) | 1e-6 / 200 | **37** | 6m23s | 263 MB | 104 | 15m46s | 217 MB | 0.36× | −60.81349836 | −60.68597088 |

⇒ **Iteration counts are AT PARITY on Si/NaF (1.0–1.3×) and in our favour on MnO (0.75× / 0.36×)**; the Si Γ
1.33× is 16 vs 12 on a 3-second run.  ⚠ **The Si 2×2×2 shifted-MP energy is −7.867452508 here, 16 µHa from
CP2K's −7.867436530, where the table above has −7.868473429 (−1.04 mHa).  Checked the same day: it is the
IMPOSITION, not the criterion or the grid** — with `QCHEM_IMPOSE_SYMMETRY=1` every criterion and both XC grids
give −7.868473; free, −7.867452.  The symmetry fold of the shifted (k=±¼, non-TRIM) mesh lowers the energy by
1.02 mHa: an OPEN DEFECT, `doc/OpenWork.md` §3 (footnote ¹ is now history).  The other energies moved by
< 1 µHa; on MnO the 108 mHa absolute offset (§4a) is unchanged.

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
`QCHEM_DM_LOWRANK=on* GPW_STREAM_FOLD=on* QCHEM_MIX_RHO_M=off QCHEM_XC_DM_SOURCE=off
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


## 5. THE ROWS — single thread  (ARCHIVED — whole-run table, footnotes, 5b/5c/5d)
Not the bin-1 instrument.  Verbatim in `doc/Records/BenchmarkHistory.md` § "Benchmark.md verbatim as of 2026-10-01", `## 5.`.  The MnO ⁸ rows (rule 3f) are summarised in §10 "Where it stands".

## 5e. Bin-1 working notes  (ARCHIVED)
Unscreened-gather split by k-mesh, Becke-build threading, the k-scaling gap, NaF's Becke cost, box walk ~95% of MnO.  The open parts are §10 rows.  Archive `## 5e.`.

## 5f. WHERE THE 2.05× IS: CALL COUNT, NOT KERNEL  (ARCHIVED, 2026-09-05)
Conclusion kept: the parity gap is gather CALL COUNT, not kernel speed; lever A landed (2.05× → 1.72×), lever B (one gather per spin) REFUTED as cheap and sits behind N4 (§10), like-for-like 1.13×; lever C parked until OT.  Evidence in the archive, `## 5f.`.

## 6 / 6b. qchem accelerations not in CP2K / CP2K accelerations we lack  (ARCHIVED)
The qchem-vs-qchem deltas §2's knobs buy; archive `## 6.` / `## 6b.`.

## 7. THE ROWS — threaded (12 cores)  (ARCHIVED; 7a–7d)
Headline kept: CP2K's parallelism does not help on these cells (0.82–1.09× on OMP, 1.44× on 12 MPI ranks); remaining threading work is ours (residual → §10 last row).
### 7c. Where OUR threaded time goes — bucket by bucket  (ARCHIVED)
MnO defaults serial 395.6 s → 12 threads 128.4 s; Becke build 8.2×, XC ρ sampling 6.6×; the unthreaded buckets are the Amdahl residual.  Table in the archive, `### 7c.`.

## 8. WHAT IS STILL MISSING / 9. STANDING OBSERVATIONS  (ARCHIVED)
Superseded by §10 and `doc/OpenWork.md`.  ⚠ Lesson kept: any MnO runtime row from a FREE run folds nothing (`[fold]` prints NONE) and must say so; `CP2K_COMPAT=1` is the best KNOWN parity, not proven.

---

## 10. STANDING MEASUREMENTS / OPEN PERF LEVERS (moved from `OpenWork.md` §4b, 2026-10-01)

Verbatim archive of the old tracker: `doc/Records/OpenWork_History5.md`.  Each row: what is open · next concrete action · record.

**Where it stands (numbers from §5a and the rule-3f table):** per SCF iteration serial qchem is AHEAD of CP2K on seven of nine §5a rows; the like-for-like parity row (`CP2K_COMPAT=1`, fixed-point stage) is **1.13×**, 1.72× for the capped two-stage probe, and its whole residual is ONE gather (lever B, behind N4).  Rule-3f whole-run, same measure/threshold/loop: MnO AFM-II **5m42s vs 6m14s wall** (33 vs 44 iterations); Si 2×2×2 8.6 s vs 5.6 s.  Peak RAM on the parity routes 113–132 MB vs CP2K's 217 MB (whole-run MnO 262 vs 217 MB).  Threaded (12 cores) 4.14 s/step vs CP2K 5.90 (they run 0.82–1.09× their own serial).  History: the MnO per-iteration cost started at **67× CP2K's** (573 s vs 8.5 s, `BenchmarkHistory.md` §cache; the earlier "100×" was an instrument artefact).  No QE timing exists.

**Parked levers:** B (one gather per spin) behind N4; C (GDM trial densities) behind OT (OpenWork §2).

| row | what is open | next concrete action · record |
|---|---|---|
| **Linear D-mixing DIVERGES on free Si with no Fock accelerator** (found 2026-09-20 by `scripts/retake5a`) | `GPW_Si.Γ_CP2K` under `CP2K_COMPAT=1 GPW_ACC=null` (plain linear D-mixing, `ProductionGates` α=0.30): descends to −7.11505 by iteration 19, then diverges (E −7.0886 at 60, [F,D] 2e-3 → 0.2).  The adaptive relax raised α 0.30 → 0.45 at iteration 2 and never re-damped.  CP2K's `DIRECT_P_MIXING` at α=0.4 converges the same cell to 1e-7 in 12 steps; our Kerker G0=1 alone takes 54, Kerker + Pulay(8) 13 (E −7.115067447, the anchor).  DIIS on the Fock side had been masking it | two questions: (a) why does the adaptive controller (V1.18) not re-damp on a rising energy — a rising-E, rising-[F,D] trajectory is exactly its trigger; (b) is direct P mixing at fixed α=0.4 stable on our loop (a `GPW_ALPHA` knob + adaptive off)?  Cheap: seconds per run.  Not a physics question; until answered the no-Fock-side route on our side is Kerker + Pulay · `doc/Benchmark.md` rule 3f, `scripts/retake5a` |
| **Re-take Benchmark §5a on CP2K's convergence measure** (rule 3f, 2026-09-20) | every `q iters` / `c steps` / `total ×` cell in §5a compares a 1e-3–1e-5 mixer residual with CP2K's `EPS_SCF` on max\|ΔP\| — two different questions; the s/it columns survive, the counts and totals do not | re-run each §5a row with `Δρmeasure=MaxΔD`, `MinΔρ=EPS_SCF`, `NMaxIter=MAX_SCF` off its deck, the deck-shaped loop (Null accelerator, Pulay 8, no MOM) where the deck mixes Broyden/Pulay, and `CP2K_COMPAT=1`; one sitting, ~1 h serial (Si rows seconds, NaF minutes, MnO ×2 the bulk).  The two MnO rows are DONE (§5 ⁸); the Si Γ example is in rule 3f · `doc/Benchmark.md` rule 3f |
| **Size the Becke grid** (was item 1 / bin 2) | MnO's setup is 184 s = 47% of the default run (CP2K 8.1 s), 136.6 s of it TWO Becke mesh builds — but the RECIPE (nR=40, degree 29) is **over-generous on every system**: 3.5× on Si/Al, ~25× on NaF/Mn (`BeckeLadder`, scored by \f$\max|\Delta V_{xc}(i,j)|\f$ — score the MATRIX, the energy cancels error the operator does not) | ▶ **a POLICY CALL for the user**: at an absolute \f$\max|\Delta V_{xc}|\le\f$ 1e-4 the answer is nR=40, GL-17…21 (2–3× cheaper than production).  ⚠ no default flips on ladder evidence alone — Al is non-monotonic and a frozen ladder understates the self-consistent shift on a metal ⇒ a converged A/B on Al first.  Two items survive whatever the calibration says: (a) why TWO builds for one cell (each anneal stage builds one); (b) the build threads at 8.2× on MnO — NaF's serial fraction is still its setup (likeliest SIZE).  Then the per-element radial scaling, the coarse-end routing calibration, Becke+IBZ (the real-space star-average is untested on this route) · History4 "1. BIN 2" |
| **BM(3) — the stock Lebedev rule is already site-invariant** | Lebedev-29 (302 dirs) is EXACTLY invariant under Si's \f$T_d\f$ site group (0 unmatched of 7248) because Lebedev rules are octahedral orbits; W2b nonetheless builds a site-adapted rule at 886 dirs/atom — a **2.9× sitting unclaimed** | test the STOCK rule for site invariance first and reuse it when it passes (the test must stay EXACT); keep W2b as the fallback; measure on MnO before spending anything — it vanishes as site symmetry drops · History4 row BM |
| **The ρ̃-mixed sampling bucket** (the actual per-iteration XC lever) | `FourierMixCD.C:65` samples \f$\rho(r)=\sum_G\tilde\rho(G)e^{iGr}\f$ by DIRECT SUMMATION at every mesh point, every Kerker/Pulay iteration: 35.0 s / 6 iterations serial against the DM GEMM's 1.70 s.  ⚠ It is NOT the low-rank D bucket (that GEMM is nearly bypassed on ρ̃-mixed recipes — the 7–8× rank win is real but bites only on DM-backed routes) and NOT Φ-sparsity (Φ is 48% dense on MnO; per-atom batching = 1.05×) | attribute the cost INSIDE that bucket before optimising: an FFT to a coarse uniform grid + interpolation to the mesh (O(N log N + npts)) is legitimate only for the cusp-free CORRECTION, not the full ρ; an adaptive G-ball keyed to \f$|\tilde\delta|\f$ is cheap late and exact at convergence.  N4 decides what is sampled · History4 "Vxc MUST BE FED THE DM ρ(r)" |
| **The k-scaling gap** (cross-k gather memo) | Si per-iteration cost rises 6.9× from Γ to 8 k-points; CP2K's rises 1.03× — we start 6.7× ahead and spend it on k.  Cause: the per-offset reductions \f$B_{ij}(n)\f$ are k-INDEPENDENT but `IntegrateMemo` is bypassed whenever a density screen is passed (32 of 51 gathers on Si 2×2×2 are the SAME FIELD).  ⛔ two attempts refuted (History4): withholding the screen perturbs the SEED trajectory; flooring vanishing weights breaks the fold's orbit invariance | the trajectory-exact route: memoize \f$B\f$ over the UNION of active sets and let each block contract only its own; raise `kMaxIntegrateMemos` (4) with it.  Pin 19 already removed the gather's D-screen, which is the other half.  CP2K also folds 8 → 4 by TIME REVERSAL on the shifted mesh — a separate 2× · History4 "THE k-SCALING GAP" |
| **Step 2's remainder — the per-iteration G-space folds** | the {G}-star fold is wired at two STATIC sites; the per-iteration consumers (ρ̃, the Poisson multiply, the V_xc gathers, G_ERI3 columns, seed structure factors) are UNFOLDED: 12–24× on MnO's magnetic group, 48× cubic, unclaimed; **T3.4b** multi-k per-block arming of the pair-stream fold (Γ-only is armed by default) | extend the ball fold to the per-iteration sites (the FFT itself does not fold trivially); T3.4b = union-of-reps stream caches or the star-summed joint scatter.  Space-group collocation reduction route (b) IMPOSES the symmetry — say so · History4 "Step 2", `doc/Records/SymmetryUpgradePlan.md` |
| **\f$B_{ij}(R)\f$ k-independent 1E memo** | "keep k out of the key"; payoff only on multi-k | time with row KP (OpenWork §2) · `GPWPlan1.md` item 5 |
| **Becke partition, the loop itself** | `-march=native` alone measured 1.13× with the loop still SCALAR; `BeckeImage` is array-of-structs with a data-dependent `P>0` exit that saves only 2.3% of work | vectorise first (chunk the exit, SoA); the `norm()` table (~2 MB gather vs a 15-cycle hardware sqrt) is DOUBTFUL and points the opposite way — re-measure after vectorising, do not assume.  Both are small beside V1.22 · History4 "Becke partition, what is LEFT" |
| **`Eval`/`EvalGradient` still run their own image loop** | the per-point callers (KB quadrature) were not migrated to the `LatticeSum1E` point-set seam because the seam re-derives its offset list per call; so the code's *"ONE remaining explicit image list"* is narrowed, not retired | cache the offsets on the seam side, then migrate; it retires the last place GPW enumerates images itself · History4 "THE Φ BLOCH POINT-SUM SEAM" |
| **Φ-sparsity for LARGE cells** | ⛔ refuted on MnO (48% dense) — the item's own caveat "the win grows with cell size" was the whole story | re-measure with `GPW_PHI_SPARSITY=1` on a battery-scale supercell BEFORE building anything |
| **The singles-vs-pairs census** | the factored density turns pairs into singles (8778 → 118 on MnO) but the comparable count is (i,R) singles vs (i,j,R) pairs at the same ε, weighted by box volume — pairs SCREEN far harder (Gaussian product theorem), which is why CP2K collocates pairs | build the singles-side census over the orbital reach before any singles-route work; the crossover is system-dependent and may favour BOTH routes · History4 "THE FACTORED FORM CHANGES THE OBJECT" |
| **Parity deviation #8** | `CP2K_COMPAT=1` is our best KNOWN parity, not proven; the deviation list grew 4 → 7 every time anyone looked | finding the next one is part of any parity claim; the like-for-like row is a `GPW_MNO_NMAX`-capped probe — step one is an uncapped or equal-cap comparison · `doc/Benchmark.md` §2, §5a |
| **Residual 1.33–1.45× CPU inflation at 12 threads** | load imbalance / memory bandwidth; ~4.1 s left in Phase 1, deliberately not taken | revisit on a 2×2×2 LiMn₂O₄ supercell where the buckets reshuffle — optimising the tail of a development cell is fitting to the wrong system |
