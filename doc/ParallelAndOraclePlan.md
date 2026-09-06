# Parallel work, and the second oracle — a sequenced plan (cut 2026-09-06)

**Why this file exists.**  Three threads converged in one session and they have a natural order that is not
obvious from any of them alone: (1) the threaded tables are filled and they say the remaining OMP gap is
OURS, not a CP2K feature we lack; (2) `DFT+U` is wanted soon and needs an oracle; (3) the abstract
interfaces have accumulated implementation detail and want a cleanup pass before a new capability lands on
them.  This is the order, with what each phase unblocks and what it costs.

📖 Measurements behind every claim: `doc/Benchmark.md` §5f (the call-count analysis), §7 (both codes
threaded).  The live tracker stays `doc/OpenWork.md`; **this file does not duplicate it** — the last
section says exactly which of its open items are folded in here and which are not.

---

## THE STANDING CONSTRAINTS (read before running anything)

- **One heavy run at a time**, and heavy sweeps through `scripts/memsafe` (14 GB box).  `ctest -j8`, never
  -j16.  No concurrent builds — a background ninja beside the user's own corrupts `.ninja_deps`.
- **COPY the run command from `doc/Benchmark.md` §4**; a large discrepancy against a probe's documented
  cost is a configuration difference, not a speedup.
- **Source builds are the DEFAULT for every comparison code** (user, 2026-09-06): *"we always want to look
  into the source any way … and with our own build we can possibly get better control and profiling."*
  ⇒ apt availability no longer decides which codes are candidates — ABINIT, SIESTA and Elk are back on the
  list.  Toolchain in place: `gfortran` and `mpif90` (OpenMPI) installed; **flang-21 at
  `/opt/LLVM-21.1.6-Linux-X64/bin`** if a code wants LLVM Fortran.  Build each with symbols
  (`-g -fno-omit-frame-pointer`) so `perf` can read it — that is half the point of building it ourselves.

---

## PHASE 1 — THE OMP GAP THAT IS OURS  ·  no new dependencies

§7 settled that CP2K has no threading trick we are missing on this class of system (its own GPW hot
routines run 1.09×/1.12× on 12 threads).  What §7c found instead is three non-scaling blocks of our own.

### 1.1 ✅ DONE 2026-09-06 (`22b7f0b9`) — AND IT IS NOT THE LINEAR ALGEBRA

Ten buckets added (the SCF step's three phases, the fill, the direct-min step, the line search's per-trial
`MoveOrbitals`, the total energy, the mix, the density install; plus the facade's four construction
phases).  `report::Timed` nests exclusively, so the existing GPW buckets read as children.

**MnO ALL DEFAULTS, serial 393 s / 12 threads 122 s:**

| what the plan suspected | measured |
|---|---|
| diagonalisation | **0.038 s TOTAL** over the whole run |
| fill orbitals / accelerator projections | 0.11 s / 0.002 s |
| total energy (all terms) | 0.84 s |
| Fock assembly overhead (its term buckets excluded) | 1.44 s |
| density mix | 2.47 s |
| **all of it together** | **~5 s of 393 s** |

⇒ **The SCF's linear algebra is not where the time is**, and the instrument was the cheapest possible way
to find that out — a session spent threading the eigensolve would have bought 0.04 s.

★ **WHAT IT DID FIND: `setup: hamiltonian ctor` = 23.9 s**, EXCLUSIVE of its children (the Becke mesh build
16.1 s threaded, the Φ tables 6.2 s).  It runs once, it does not thread, and it is now the largest named
block outside the box walk.
⚠ **Unbucketed: 49 s → 25 s.**  What is left is inside `Iterate` but outside the per-iteration buckets.

⚠ Note for anyone reading another row's ledger: the FACADE path's `seed + ortho` bucket is 0.001 s because
it defers the first Fock into the SCF loop; the test harness's same-named bucket on `RunGpw` CONTAINS that
first Fock.  Same label, different content — check which harness produced the row.

*(The brief this answered: ~59 s of a 128 s threaded run in no bucket at all — larger than every known
non-scaling bucket combined — with the SCF's non-GPW work as the suspect list.  Instrument, re-run, and let
the result choose what comes next.  It did, and it chose neither 1.2 nor 1.3.)*

### 1.1(a) ✅ DONE 2026-09-06 — THE RESIDUAL 25 s WAS THE **SECOND** HAMILTONIAN, AND IT IS NOT IN `Iterate` AT ALL

The brief was "bracket `Iterate`'s loop body and the facade's `Converge` to name the residual 25 s".  Both
brackets are in, and **both came back ~zero** — so the premise they were built on ("what is left is inside
`Iterate`") was wrong, and the instrument said so in one run:

| the residue bucket | charged |
|---|---|
| `scf: iterate (residue — loop body outside the named buckets)` | **0.33 s** (×2 stages) |
| `scf: converge (residue — final density + m(r) extraction)` | **0.0009 s** |
| `setup: facade ctor (residue — decisions between the named buckets)` | **0.011 s** |

✅ **FIXED 2026-09-06 by `doc/CleanupCandidates.md` R2.22 (`0210cfb9`): it no longer does.**  The rebuild
was never a physics requirement — `~tSCFIterator` deleted the Hamiltonian its own caller had built, so the
facade had no way to keep one.  With the iterator holding it non-owning, ONE Hamiltonian serves the whole
schedule: **MnO 12-thread wall 83.2 → 67.6 s (1.23×), `Etot=-61.40297529` bit-identical, 814/814**, and the
ledger's `setup: hamiltonian ctor` / `setup: becke mesh build` show one call each.  The accelerator is still
per stage (stale Pulay/DIIS, and the type changes Ladder→GDM) — cheap, and never the cost.

★ **WHERE IT ACTUALLY WAS: `SolidCalculation::BuildStage` REBUILDS THE WHOLE HAMILTONIAN FOR EVERY ANNEAL
STAGE, and that call had no bucket.**  1.1 read the ledger's single `setup: hamiltonian ctor` entry and
concluded "it runs once".  It runs **once per stage**, and the MnO recipe (`MNO_ANNEAL="5e-3,0"`) has two —
so stage 1's rebuild was invisible, and it is the 25 s.  With `BuildStage` charging the SAME bucket, the
ledger now reports it as a per-call price instead of two unrelated rows.

**MnO ALL DEFAULTS, BOTH ARMS, same build, `Etot=-61.40297529` on both to all printed digits** —
`GPW_OMP_THREADS=12`: 121.0 s wall / 569.5 s CPU / 469% / 562 MB.  `GPW_OMP_THREADS` unset: 313.9 s wall /
532 s CPU / **169%** / 496 MB (⚠ "unset" is not "serial" — blaze still took 1.7 cores; §4's own pin):

| bucket | GPW_OMP unset | 12 threads | ratio | what 1.1 believed |
|---|---|---|---|---|
| `setup: hamiltonian ctor` (**exclusive**) | 51.43 s `[×2, 25.7/call]` | **46.36 s** `[×2, 23.2/call]` | **1.11×** ⛔ | 23.9 s, "runs once" |
| ⤷ `setup: XC-mesh Φ tables` | 54.10 s `[×2]` | 6.25 s `[×2]` | **8.7×** ✅ | 6.2 s |
| ⤷ `setup: becke mesh build` | 15.43 s `[×2]` | 15.77 s `[×2]` | **1.0×** ⚠ | 16.1 s |
| **Hamiltonian construction, all in** | **121.0 s = 39%** | **68.4 s = 56%** | 1.77× | ~46 s = 38% |
| all three residue buckets together | 0.36 s | 0.36 s | — | (did not exist) |
| **everything not in a bucket** | **0.05 s** | **0.03 s** | — | 25 s |
| *(39 buckets, summing to)* | *313.71 of 313.76 s* | *120.98 of 121.01 s* | | |

★ **THE CTOR'S EXCLUSIVE HALF IS THE ONE THAT DOES NOT THREAD (1.11×)** — its two children do (Φ tables
8.7×), which is exactly why it went unnoticed: the bucket beside it scales beautifully.  At 12 threads it
is 38% of the wall on its own.

✅ **ONE CROSS-ARM READING DID NOT MATCH THE BANKED TABLE — FLAGGED, THEN RESOLVED AGAINST ME.**  The Becke
mesh build measured 15.4 s unset against 15.8 s at 12 threads (no change), where `doc/Benchmark.md` §7c
banks 138.2 → 16.9 s (8.2×).  One `GPW_OMP_THREADS=1` run settled it: **68.74 s/call serial against
8.15 s/call at 12 threads = 8.4×**.  §7c was right and my arm was wrong — **`GPW_OMP_THREADS` UNSET IS NOT
A SERIAL ARM**: the Becke build is parallel by default and reads the variable only as a thread CAP, so
"unset" runs that loop on all cores while the GPW pair loops stay serial (169–188% CPU, which is the tell).
A serial row must set `GPW_OMP_THREADS=1` and check for 99% CPU.  Recorded in `doc/Benchmark.md` §7c as a
protocol rule, because it invalidates the serial arm of any row taken the mixed way.

⇒ **THE LEDGER IS NOW A PARTITION OF THE RUN**, in both arms, which is the property that makes it an
argument rather than a list: "everything not in a bucket" can no longer be the largest block in the table,
so a future session cannot be sent chasing one.  And the split it reveals is the headline for bin 2:
**setup 69.5 s against SCF 51.5 s at 12 threads** (122.8 vs 190.9 unset) — threaded, this run spends MORE
time building the Hamiltonian than converging it.

⛔ **ONE SUSPECT REFUTED IN PASSING.**  The direct-min line search calls `itsHamiltonian->GetTotalEnergy`
DIRECTLY, bypassing the iterator's bucketed `TotalEnergy()` helper — up to 12 full density builds + energy
evaluations per GDM iteration that nothing had ever measured.  It looked like the residual's obvious home.
It is **1.21 s over 44 calls** (0.027 s/call).  Bucketed now, and cheap.

▶ **NEXT — (b)**: open up the **46.4 s exclusive** Hamiltonian ctor.  ✅ **DONE — see 1.1(b) below.**

★ **FOLDED IN AND DONE — `doc/OpenWork.md` item 5 (Step 0c), "the instruments report WHAT, not WHEN"**:
a timestamp per report item would localise the residual without adding a single bucket, because the GAPS
BETWEEN SECTIONS are exactly the unbucketed time.  ✅ **Built with (a), and it earned its keep on the first
run** — every console heading now carries the run clock (`grids ▸ becke  [t=11.50 s]`), and the same MnO
log reads:

```
[MnO AFM-II Gamma stage 1/2] … iters=14 …          ← last stage-1 line, t = 60.17 s
scf ▸ siteMoments  [t=95.93 s]                     ← first stage-2 line
```

**35.8 seconds in which the run printed nothing** — the stage-2 Hamiltonian rebuild, named by the STAMPS
alone before any bucket was read.  (The same reading at the head of the run: `grids ▸ becke` at 11.50 s,
first SCF item at 31.64 s = the 20 s of stage-0 ctor.)  ⇒ The claim in item 5 was right, and this is the
run that demonstrates it.

### 1.1(b) ✅ DONE 2026-09-06 — THE CTOR WAS ONE BAD INDEX, AND THE RUN IS **1.45× FASTER**

Same method as 1.1: bucket before optimising.  Six buckets took the ctor's 23.75 s/call apart in one run,
and the answer was not distributed at all — it was two calls to the same routine.

| phase of the ctor | s/call | was it suspected? |
|---|---|---|
| **XC mesh orbit-consistency fold** (`UnitCell::CreateIntegrationMesh`) | **9.61 s** | no — it had no bucket |
| **XC mesh orbit fold** (`GPW_IBS::CreateXCQuadrature` → `FoldMesh`) | **9.47 s** | no |
| Becke mesh build | 8.15 s | already bucketed |
| site-adapted angular sets (W2b) | 4.03 s | no |
| XC-mesh Φ tables | 3.14 s | already bucketed |
| Hartree CD fit basis · IonIon Ewald · term ctors · PP models (GTH) · `DeltaFit_IBS` ctor · quadrature copy | **< 0.001 s each** | all of them |

★ **THE RUN FOLDS THE SAME ~97k MESH POINTS TWICE PER HAMILTONIAN** — once in `UnitCell` to build the
orbit-consistency keep-mask, once in `GPW_IBS` to build the fold it keeps (the code says so: *"the caller
rebuilds its orbit fold FROM the finished mesh"*).  **19.1 s per call, 38.2 s of a 121 s run.**

★★★ **BUT DE-DUPLICATING THEM IS THE WRONG FIX, BECAUSE ONE FOLD WAS ALREADY 50× TOO SLOW.**
`TorusIndex` (`src/Symmetry/Lattice_3D/Imp/Fold.C`) matched images through a bucket grid whose resolution
was a **constant 64 per axis**: its ctor started at 64 and only ever SHRANK, and the shrink condition
(`1/nb <= 2*tol`) is never true at the mesh tolerance 1e-8.  Average occupancy looked perfect — 97256
points in 64³ = 262144 buckets is 0.37 per bucket — **and the average is the wrong statistic**: an
atom-centred RADIAL mesh is clustered, the inner shells put thousands of points within one bucket-edge of a
nucleus, and every query near a nucleus scanned all of them.  The buckets are a sparse `unordered_map`, so
a far finer grid costs no memory it does not use.

**FIXED**: the grid is now as fine as the tolerance allows (capped at 2^20/axis so the key stays in int64),
and `Find` probes the query's OWN bucket first — these meshes are op-invariant by construction, so
\f$Wp+\tau\f$ is another mesh point to ~1e-15 and the centre bucket hits essentially every time, turning
27 hash lookups into 1.  The full 3×3×3 neighbourhood remains the correctness path for a point near a
bucket face.

| | before | after | |
|---|---|---|---|
| orbit-consistency fold | 9.61 s/call | **0.193 s/call** | **50×** |
| `FoldMesh` | 9.47 s/call | **0.190 s/call** | **50×** |
| Hamiltonian construction, all in | 34.4 s/call | **15.5 s/call** | 2.2× |
| **MnO ALL DEFAULTS, 12 threads, WHOLE RUN** | 120.8 s | **83.2 s** | **1.45×** |

✅ `Etot = -61.40297529` — **bit-identical to all printed digits**, as it must be: this is a data structure,
not a numerical method.  813/813 green.  Pinned by
`SymmetrizeMesh.TorusFoldIsIndependentOfTheBucketGridOnAClusteredMesh`, which folds a five-decade clustered
set at three tolerances (hence three different grids) and asserts the orbits are the same — the property
that makes any future index change safe.

⚠ **THE SECOND FOLD IS STILL THERE, AND IT IS NO LONGER WORTH REMOVING** — 0.19 s/call.  Recorded so nobody
re-derives it: the duplication is real, it is now 0.3% of the ctor, and a correctness-preserving merge
(the two folds differ in input — pre- vs post-filter point list — and in tolerance) would buy 0.4 s a run.
**Do not spend a session on it.**

▶ **WHAT (b) LEAVES OPEN, now that the ctor is 15.5 s/call:**
1. ✅ **The Becke mesh build is the largest setup bucket (8.15 s/call at 12 threads) and it THREADS at 8.4×**
   — settled above, §7c stands.  Serially it is 68.74 s/call, i.e. **38% of a true-serial 361 s run**, so it
   is the top target for 1.2/1.3's serial arm rather than for threaded work.
2. **Site-adapted angular sets, 4.03 s/call, and they do NOT thread** (4.05 serial / 4.04 unset / 4.03 at
   12 threads — flat across every arm).  Never measured before, never suspected; now ~5% of the threaded
   wall and rising as everything around it gets faster.
3. **Why it is built TWICE at all.**  ⚠ **ANSWERED, AND IT IS OWNERSHIP, NOT PHYSICS**: `tSCFIterator`'s
   destructor does `delete itsHamiltonian` — it deletes an object it did not create — so `BuildStage` MUST
   hand each stage a fresh one or the previous iterator's destructor takes it down.  The test harness says
   so in a comment: *"Fresh Hamiltonian + accelerator per stage (the iterator OWNS + deletes them; a kT
   change must not carry stale DIIS history across the re-seed)"* — and the stated physics reason (stale
   DIIS history) applies to the **accelerator**, which genuinely must be fresh.  The Hamiltonian is a pure
   function of (structure, basis, species, functional, xcMesh, vxcFit), none of which change between
   stages, and it rides along only because of the `delete`.  Sharing it would now save **~15.5 s of an
   83 s run (19%)**.  ★ This is `doc/CleanupCandidates.md` material as much as a perf item — CLAUDE.md
   says `delete` should be rare or non-existent — so it is filed there; see **R2.22**.

### 1.2 THE BLAS-MODE SERIAL ARM  ·  `-DQCHEM_BLAZE_BLAS=ON`, **pin kept**

The prize is SERIAL, not threaded: the two Φ-table GEMMs are **86 s of a 396 s serial run** (ρ sampling
74.9 + H_xc 11.2), and `DeltaFit_IBS::AdjointT` already measured the two paths — one dispatched
whole-matrix `zgemm` at **34.1 GFlop/s against 1.87 for any blocked or viewed form**.  That dispatch was
correctly refused on 2026-08-15 because the system BLAS was then netlib; `libblas`/`liblapack` on this box
now resolve to **openblas-pthread**, so the fast path is available again.
⚠ **Keep `qchem::PinBlasToOneThread`** — this arm is about a faster serial GEMM, not about threaded BLAS.
⚠ **CMake gotcha**: `set(... CACHE ...)` does not override an existing cache entry — use a fresh tree or
pass `-DQCHEM_BLAZE_BLAS=ON` to `cmake` on the existing one explicitly.
**Accept on**: serial CPU down on the MnO default + `BECKE_XC=0` rows, 806/806 green, and the anchor delta
RECORDED (it will not be bit-identical — a different summation order, same class as §5f lever A).
★ **THIS MOVES ANCHORS ⇒ it belongs to `doc/OpenWork.md` item S, the anchor-moving sprint** (A2–A6 still
open, to be re-banked in ONE window so they do not mask each other).  Do not re-bank it alone.

### 1.3 NESTED THREADING — the 2×6 (user, 2026-09-06)

Only if 1.2 wins.  Today our parallelism lives ABOVE the linear algebra (per k-block / irrep / spin) with
BLAS pinned to one thread, and at Γ with 2 spins **that level is 2 wide** — ten of twelve cores idle in the
XC quadrature bucket by construction.  ⇒ A THREAD BUDGET in `qchem.Parallel`:
`outer_width × blas_threads ≈ cores`, with `outer_width` = the width the caller actually has
(n_spin × n_kblocks), so Γ/2-spin gets 2×6 and an 8-k run keeps 12×1.
★ **The nesting hazard is NOT the one the tree hit before**: OpenBLAS-pthread runs its own pool, so
outer-OMP × inner-pthread is not an OMP nested region (`BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0` stays as
is).  The real risk is plain oversubscription, which is what the budget exists to prevent.
**Measure** 2×6 against 12×1 and 1×12 on MnO Γ.  **Exit criterion for Phase 1**: MnO Γ at 12 cores,
**3.08× → ≥5×**.

---

## PHASE 2 — THE SIZE QUESTION  ·  one deck each

### 2.1 OUR OWN SCALING CURVE — Si supercells, 2 → 4 → 8 → 16 atoms

Cheap, and it is the honest answer to "do we have OMP gaps at PRODUCTION size": our 3.08× is measured on a
4-atom cell, and a small cell starves threads on both sides of this comparison.  If our speedup climbs with
size, Phase 1's exit criterion should be read at the large end.

### 2.2 THEN CP2K's 32-ATOM MnO SUPERCELL (user)

Settles the one caveat §7b could not: is CP2K's 1.09× a route that was never parallelised (their research
centre is hundreds of small molecules, where MPI over molecules is the axis that pays) or simply too few
tasks to spread?  Cheap on their side — their regime — and we need not run ours.  It closes a doc caveat
rather than changing our plans, hence second.

---

## PHASE 2.5 — THE SOLID/OOD CLEANUP CAMPAIGN  ·  before any new capability lands (user, 2026-09-06)

> *"Implementation details keep creeping into our abstract interfaces."*

A new capability (+U) will be built ON these faces, so it goes after the pass, not before.  The campaign is
not open-ended — it is the already-recorded debt, plus the two smells this session surfaced:

| what | where it is written down |
|---|---|
| `doc/CleanupCandidates.md` R1 (correctness-adjacent) and R2 (mechanical hygiene) | that file, and `doc/OpenWork.md` item 7 |
| **V1 — the interface-design questions**, which is exactly the "implementation detail in an abstract face" class | `CleanupCandidates.md` V1 |
| **`FIT_SF_Ortho` — separate the METRIC axis into faces** (specced 2026-08-23, not built): `OverlapDiagonal` sits on the metric-NEUTRAL face, so `Fit_IBS` answers in the wrong normalisation.  ⚠ Acceptance: must NOT become a `dynamic_cast` type switch | `doc/OpenWork.md` item 6 |
| **The identity-question smell**: `IsPolarized()` / `IsRelativistic()` on the Hamiltonian faces — the `IsSlaterBasisSet()` shape the user has already ruled against once (`IsGeometryOnly()`, 2026-09-04).  Test: does the answer SELECT WHICH CODE RUNS? | `CleanupCandidates.md` (the Hamiltonian-faces note) |
| ★ **NEW — the GRID × FIT-BASIS audit** (user, 2026-09-06): the integration grid and the fit basis are ORTHOGONAL axes, and which pairings are worth running is high-level POLICY.  Audit the faces for anywhere a grid choice IMPLIES a basis choice (or vice versa) — "Becke" must not mean "delta" in any signature | this file; `IntegrationTests/GPW_SCF_UT.C` V2.8 |
| **The `dynamic_cast` survey** — CLAUDE.md's standing TODO: casts must be abstract→abstract, never abstract→concrete, and a failing one should throw with context | `doc/FittingCleanupPlan.md` item C |
| **V1.32** — de-template the finite `IrrepCD` leaf | `doc/OpenWork.md` item 7 |

⚠ **Do the cleanup with the anchor sprint in mind**: anything that changes a summation order joins item S's
single re-bank window (see 1.2).

---

## PHASE 3 — DFT+U, ORACLE FIRST

### 3.1 ✅ THE ORACLE ALREADY EXISTS AND IS ALREADY VALIDATED — CP2K supports DFT+U

`&DFT_PLUS_U` per `&KIND` with `U_MINUS_J`, and `PLUS_U_METHOD MULLIKEN | LOWDIN` at the DFT level
(verified in the installed CP2K 2025.2 input reference).  ⇒ **No new package is needed to start +U.**
Produce the oracle row BEFORE writing our +U: MnO AFM-II, the deck we already trust, one `U_MINUS_J` on the
Mn kind.  This repo's habit is that a capability lands with an anchor; here the anchor is nearly free.

### 3.2 IMPLEMENT +U WITH THE PROJECTOR FLAVOUR MATCHED DELIBERATELY

A cross-code +U comparison means nothing unless the PROJECTOR definition matches — Mulliken and Löwdin
give different occupation matrices for the same density.  CP2K offers both, which is convenient: pick one,
implement that one, and **declare the choice on the `RunPolicy` deviation line** beside the other seven, so
a run states which +U it is running.
★ **`doc/OpenWork.md` N3 (charge and spin need separate preconditioning) belongs WITH this**, not after it:
+U on an antiferromagnet is a spin-channel-sensitive capability, and today's mixer takes charge medicine in
the magnetisation channel.  Landing +U on top of that mixer is how a +U bug and a mixing bug become
indistinguishable.

### 3.3 A SECOND +U ORACLE ONLY IF 3.1/3.2 DISAGREE INEXPLICABLY

Different basis, different projector, independent implementation.  See Phase 4 for the candidate.

---

## PHASE 4 — THE SECOND CODE  ·  when there is a question for it

★ **The question it must answer is already open**: `doc/OpenWork.md` item 4 / Step 5 — **MnO's −99.65 mHa
against CP2K, operator not yet named**.  One oracle cannot tell us whether that gap is ours or theirs; a
second, with a different basis and a different symmetry treatment, can.  That — not curiosity — is the
trigger for this phase.

| candidate | GPW-like | symmetry folding | DFT+U | build | verdict |
|---|---|---|---|---|---|
| **Quantum ESPRESSO** | ✗ (PW/PAW) | ✓✓ full space-group IBZ reduction — the only candidate that might beat us at criterion 2 | ✓✓ mature (and DFT+U+V) | source; gfortran + OpenMPI + ScaLAPACK all present | **the accuracy + symmetry oracle**; NOT a timing peer |
| **GPAW** | ~ (PAW; LCAO + real-space grid mode is the closest in spirit) | ~ | ✓ | source (Python + C) | the closest thing to a GPW timing peer among easy builds |
| **Elk** | ✗ (all-electron LAPW) | ✓✓ | ✓ | source | all-electron TRUTH for a small cell; slow, and that is fine for an oracle |
| **ABINIT** | ✗ (PW/PAW) | ✓✓ | ✓ | source (no apt candidate here) | overlaps QE; pick one |
| **SIESTA** | ~ (NAO + real-space grid — arguably the closest architecture to ours) | ~ weak historically | ✓ | source | interesting for architecture comparison, weakest on criterion 2 |

**Recommendation**: **QE first** (it answers Step 5 and criterion 2 at once), **Elk** if Step 5 needs an
all-electron arbiter, GPAW/SIESTA only if we want a second *timing* peer for the GPW-like route.

---

## WHAT IS DELIBERATELY NOT IN THIS PLAN

| not now | why, and what would change it |
|---|---|
| the Becke MESH BUILD (bin 2's 136.6 s) | gated on the grid-size calibration — and the V2.8 ladder says the recipe is 3.5–25× over-generous, so the SIZE question comes first.  ⇒ the threshold is a policy call (`doc/OpenWork.md` item 1) |
| §5f **lever B** (one gather per spin) | REFUTED as a free change: ball vs raw discretisation, 6e-5 relative.  Unblocked only by **N4** |
| §5f **lever C** (GDM's trial densities) | parked behind **OT** — no comparable step in a CP2K run that is mixing rather than minimising |
| chasing CP2K's OMP | settled by §7; their grid routines do not scale here either |

---

## WHAT WAS FOLDED IN FROM `doc/OpenWork.md`, AND WHAT WAS NOT

**Folded in** (they now have a place in a sequence rather than a standing list):

| item | where it lands here | why |
|---|---|---|
| **S** — the anchor-moving sprint (A2–A6) | **1.2** | the BLAS arm moves anchors; it must share the ONE re-bank window |
| **5** — Step 0c, instruments report WHAT not WHEN | **1.1** ✅ **DONE 2026-09-06** | same mechanism, same file, one job — and it located 1.1(a)'s answer before any bucket was read |
| **6** — `FIT_SF_Ortho` metric axis | **2.5** | it is precisely "implementation detail in an abstract face" |
| **7** — continuous cleanup (`CleanupCandidates` R1/R2/V1, V1.32) | **2.5** | the campaign IS this item, given a window |
| **N3** — charge/spin mixing channels | **3.2** | +U on an AFM is spin-channel-sensitive; land the mixer split with it |
| **4 / Step 5** — MnO accuracy, name the operator | **4** | it is the QUESTION that justifies a second oracle |

**Not folded in, and why**: **N4** (the \f$V_{xc}[\rho\ge0]\f$ track — it gates lever B and is its own
plan), **OT** (gates lever C; a separate build), **N1's coverage gap** (`RunGpw`/`RunGpwAnnealed` bypass the
facade, so the detectors never run on the benchmark/ladder/scaling runs this plan is full of — ⚠ worth
doing before Phase 2's runs are trusted, but it is a tracker item, not a phase), **2** (benchmark protocol —
substantially satisfied by §5a/§7's per-row thread state), **8** (the 136-function span — unrelated).
