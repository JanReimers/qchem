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

▶ **NEXT (start here)**: (a) bracket `SCFIterator::Iterate`'s loop body and the facade's `Converge` to name
the residual 25 s; (b) then open up the 23.9 s Hamiltonian ctor — it is one call, it is serial, and at 12
threads it is 20% of the wall.
⚠ Note for anyone reading another row's ledger: the FACADE path's `seed + ortho` bucket is 0.001 s because
it defers the first Fock into the SCF loop; the test harness's same-named bucket on `RunGpw` CONTAINS that
first Fock.  Same label, different content — check which harness produced the row.

### 1.1 (original brief) — INSTRUMENT THE UNBUCKETED WORK — do this first, and do not optimise before it

**~59 s of a 128 s threaded MnO run is in no bucket at all**: diagonalisation, orthogonalisation, mixing,
the fit solves, SCF bookkeeping.  That is larger than every known non-scaling bucket combined (~21 s), and
nothing times it.  ⇒ Add report buckets around the SCF's non-GPW work, re-run the 12-thread MnO row, and
let the result choose between 1.2, 1.3 and something not yet on this list.
★ **Fold in `doc/OpenWork.md` item 5 (Step 0c) here** — "the instruments report WHAT, not WHEN" is the same
mechanism and the same file; adding a timestamp per report item while adding buckets is one job, not two.

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
| **5** — Step 0c, instruments report WHAT not WHEN | **1.1** | same mechanism, same file, one job |
| **6** — `FIT_SF_Ortho` metric axis | **2.5** | it is precisely "implementation detail in an abstract face" |
| **7** — continuous cleanup (`CleanupCandidates` R1/R2/V1, V1.32) | **2.5** | the campaign IS this item, given a window |
| **N3** — charge/spin mixing channels | **3.2** | +U on an AFM is spin-channel-sensitive; land the mixer split with it |
| **4 / Step 5** — MnO accuracy, name the operator | **4** | it is the QUESTION that justifies a second oracle |

**Not folded in, and why**: **N4** (the \f$V_{xc}[\rho\ge0]\f$ track — it gates lever B and is its own
plan), **OT** (gates lever C; a separate build), **N1's coverage gap** (`RunGpw`/`RunGpwAnnealed` bypass the
facade, so the detectors never run on the benchmark/ladder/scaling runs this plan is full of — ⚠ worth
doing before Phase 2's runs are trusted, but it is a tracker item, not a phase), **2** (benchmark protocol —
substantially satisfied by §5a/§7's per-row thread state), **8** (the 136-function span — unrelated).
