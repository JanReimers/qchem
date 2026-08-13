# Cleanup candidates (running list for the SOLID/OOD/CleanCode session)

Things noticed in passing while adding features — flagged here instead of fixed inline, so the
refactoring session can batch them.  (User keeps the master list; merge freely.)

Reorganized 2026-08-05 (Claude review pass): everything below the SOLID section is grouped by
status — **READY / VERIFY / DECIDED-ELSEWHERE / DONE** — with per-item verdicts from a 4-sweep
code verification.  User prose retained (typos fixed only); user replies of 2026-08-05 folded in.

Split 2026-08-09: the closed record moved to **`doc/CleanupHistory.md`** (nothing trimmed — each closed
item leaves a one-line stub here and its full text there, so `see R2.8`-style references still resolve).
Read the history file's preamble before deciding to prune it; the rejected attempts have repeatedly been
worth more than the landed ones.

---

# START HERE (handoff, 2026-08-10)

## Just landed: **V1.31** (`627a4ff9`), **V1.27**'s live half (`a1e1f9bb`), **R2.18**+**R2.19**
(`86c5b24d`), **R2.9** (`268473b9`), **R1.7** (`26af31b6`).  Full records → doc/CleanupHistory.md.
**Merged both ways with the MnO clone and pushed** (`2e587ee6` on origin/main); qchem6 is synced.

## Next task: pick from the open list below — nothing is half-done

- **V1.27's remaining half 🔶 DESIGN RULED (user 2026-08-10), not started.**  *"I suspect we are
  ultimately going to need all 4 combos of {Molecular,Solid}x{PP,Non-PP}SCFIterator.  This should be done
  with mixins."*  That retires the "rename" framing: there is no NAME that fixes a MISSING DIMENSION.
  See the item under READY/R2 for the column analysis that shapes the mixin split.
- **`SCFParams::MinVirial = 1e-13`** — documented *"effectively off"* while the test is
  `error < MinVirial`, so the default is maximally ON.  Comment or default is wrong.  Affects ALL-ELECTRON
  runs (where the gate is real), so it is a decision about what the default gate SHOULD be, not a typo fix.
- Then the R2/V backlog as before.

Four things worth carrying forward, because they generalize past their own items:

- **R1.7 — an ISP split can pay off on DIP grounds.**  Once `MakeDirect`/`MakeExchange` were off the
  client-facing face, a grep showed NOTHING outside qcBasisSet had ever named an `ERI4` — so the substrate
  went into an `Internal.` module and two libraries stopped importing the 4-index type.  When an item says
  "this interface promises more than its clients consume", check whether the surplus also crosses a LIBRARY
  boundary; if it does, the ISP fix and the DIP fix are one edit.
- **R2.9(ii) — when an item offers two fixes, check whether it named the right DEFECT.**  The item said
  "returns a reference to shared scratch — document or return by value".  By-value turned out to be
  impossible (it OVERRIDES a reference-returning virtual) and documenting understated the problem: the real
  defect was that the two implementations of ONE interface made DIFFERENT reference-lifetime promises, and
  a caller holding the base pointer could not tell which it had.  Fixing the asymmetry was cheaper than
  either offered option.  **Re-derive the defect before picking from an item's menu** — the menu was written
  earlier, with less in front of it.

## A convention now worth stating once, because three items in a row turned on it

**`GetXxx()` is the CACHED accessor and returns a REFERENCE; `MakeXxx()` is the uncached compute and
returns BY VALUE** (user, 2026-08-10).  A caller needing an owned copy asks the `Make` half — it never asks
`Get` to change its return type.  Consequences that already paid out:
- R2.9(ii) considered "return by value" on `GetMatrix`: wrong half of the pair.
- R2.18: qcHamiltonian was spelling `Make` as `CalculateMatrix`/`CalcMatrix` — two names, one role.  Renamed.
- R2.19: a class with a `Make` that only forwarded turned out to need no `Make` at all.  **The test is
  useful in general: if your `MakeXxx` does not COMPUTE anything, you are not an implementer, you are a
  forwarder — override `Get` and hand back what you are forwarding to.**
**Visibility is NOT part of the convention** and is deliberately unsettled — see R2.18's user ruling.
Protected was the original state; the public ones are drift, each from a case where "purely internal"
proved wrong.  Do not standardise either way yet.

### STANDING PRACTICE (user, 2026-08-10: *"note why — yes exactly!!"*)

**When a `MakeXxx` has to become public, say WHY at the declaration**: which client needed the by-value
form, and why the cached `Get` would not do.  One line, at the point of the change.

The reason this is worth a rule rather than a habit: the encapsulation policy is deliberately being left to
EMERGE from refactoring, and a policy can only emerge from evidence.  A visibility change with no recorded
reason destroys exactly the evidence the eventual decision needs — after five of them nobody can tell
whether the pattern is "fitters need raw integrals", "tests need to bypass the cache", or five unrelated
accidents.  The rule generalises past `Make`: **any deliberate loosening of encapsulation should carry its
reason, because the decision to tighten it again later can only be made from those reasons.**

- **V1.27 — name a capability query for what the CLIENT consumes, not for the CAUSE.**  `IsVirialValid()`,
  not `IsPseudopotential()`: PPs are not the only thing that breaks the virial theorem, so the cause-named
  version would have been correct today and wrong at the first non-Coulombic term that is not a PP.  Third
  instance of the same lesson (R2.13 "Becke", R2.17 "SiteAdaptedBecke", R1.7's ERI4 face).
- **V1.31 — when a cache needs an awkward key, suspect the LOOP, not the key.**  `SymFockCache` needed an
  elementwise density compare because of WHERE it sat; the position was forced by a pair loop that should
  not have been running for that basis at all.  Removing the loop deleted the cache, its staleness test and
  its incomplete key together.  **A cache that is hard to invalidate is often a cache that should not
  exist.**

## ✅ V1.31 DONE `627a4ff9` — full record → doc/CleanupHistory.md.  (analysis kept below)

## The analysis that produced it — and the diagram that changed the answer (2026-08-10)

The user asked for a flow diagram of one Fock build, on the hunch that *"this whole thing is just designed
wrong.  We are somehow caching the wrong thing in the wrong place."*  **The hunch was right, and the chain
shows it in one read:**
- `ContractAll` ALREADY runs the whole sweep exactly once per density serial.
- Below it, the composite's pair loop reaches the SALC decorator ~N² times — a loop built for ERI4 PAIR
  BLOCKS, which the SALC path does not have at all (R1.7).
- `SymFockCache` exists only to stop that loop from doing ~N² whole-molecule AO builds.  It collapses
  N² → N.
- **But \f$J\f$ is LINEAR in \f$D\f$:** \f$\sum_C J_{AO}(O_C D_C O_C^\mathsf{T}) =
  J_{AO}(\sum_C O_C D_C O_C^\mathsf{T})\f$.  Those N builds ARE one build.
⇒ It caches a **partial** AO Fock (J from ONE irrep block's density), at the **basis** level, inside a loop
that should not be iterating for this path.  Wrong thing, wrong place — and the elementwise D compare is
FORCED by that position: down there the only thing in hand is a matrix, which is why no version is
reachable.  The staleness question was a symptom, not the disease.
**So the fix is not a better staleness test and not a rehoming: give the SALC path its OWN
`AccumulateDirectAll` — sum the AO densities, build ONCE, slice N times.** The memo then has nothing to
memoize and deletes itself, taking the missing-`Ocd` key defect with it.  Needs a design call on where the
SALC path branches (a capability question on the face R1.7 created is the obvious candidate), which is why
it is still not started.

*(the earlier analysis, still accurate, follows)*

## Earlier: V1.31 blocked on a ruling — do not start it as written

Attempted 2026-08-10; the settled design does not survive the call path.  **Two facts kill it** (full
evidence in the item): `SymFockCache` is an INTRA-SWEEP memo, not a cross-iteration cache — the term's own
version guard means the sweep runs once per density serial, so the cache can never hit across iterations;
and a version counter would BREAK the polarized SALC path, because `tPolarized_CD::AccumulateDirectAll`
runs Up then Down inside ONE sweep against the same key while `tPolarized_CD::Version()` forwards to Up.
The elementwise D compare is the only thing separating the two channels.  `M_Sym.water_HF_polarized` is
live and green on exactly that path.  "Move it to the term" is separately not executable: the memo is an
INTERMEDIATE AO build inside `dm->AccumulateDirectAll(X)`, on the far side of the density's virtual API,
at a granularity no term can see.
**The question for the user is now: SCOPE the memo to one sweep (staleness stops being a question), or go
further and exploit linearity — \f$\sum_C J_{AO}(O_C D_C O_C^\mathsf{T}) = J_{AO}(\sum_C O_C D_C
O_C^\mathsf{T})\f$ — for ONE AO build per sweep instead of N, which deletes the memo outright.  Both need
a sweep boundary that does not exist yet, so they are one question.**

**Also queued, unruled:** **R1.9** — molecular `BasisSetID()` prints hex.  Small, but it invalidates a
diagnostic and a cache-key claim that `SymFockCache`'s own comment makes.

## Coordination — the MnO campaign runs in the qchem6 clone, in parallel

**DO NOT TOUCH** (their working set): `src/SCFAccelerator`, `src/SCFIterator`, `src/WaveFunction`,
`IntegrationTests/GPW_SCF_UT.C`, `src/BasisSet/Lattice_3D/`.  If a task needs one of them, STOP and ask
rather than reaching in — a mid-flight collision there costs them a multi-hour run.
Theirs by assignment: **V1.30** (urgent — makes imposition opt-in), **V1.24(i)+(iii)**, **V1.28/V1.29**
(Shubnikov + the SSB workflow).

## Box discipline — 14 GB, SHARED

`ctest -j8`, never `-j16` (CLAUDE.md).  When an MnO run holds ~6 GB, **only `-j2` is safe** — GPW tests
peak 1–2 GB each, so `-j8` wants 8–16 GB and `-j4` can still overrun.  Compiling at `-j2` is fine
alongside a run; ask before taking the box for a full sweep.

## Build corollary learned 2026-08-09 (not yet in CLAUDE.md)

CLAUDE.md says build `allTests`, not `ITMain`.  The corollary: **a NEW test exe must be added to the root
`CMakeLists.txt` `allTests` DEPENDS list, or it is silently never built or run** — `ctest` then emits a
`<name>_NOT_BUILT (Not Run)` placeholder, which does NOT appear in the "N tests passed" line.  This bit
`UTSCFAccelerator`: 163 lines of accelerator unit tests that had never executed.

## State at handoff

- Branch `solid-cleanup` @ `26af31b6`; **689/689 ctest green** (`-j2`, alongside a ~6 GB MnO run; 38 disabled),
  warning-free build.
- Ahead of `origin/main` (at `77130eff`).  Includes the `allTests` fix the MnO side also needs — worth
  pushing early.
- Landed 2026-08-10: **R1.7**.
- Landed 2026-08-09: **D6**, **V1.26**, **V2.4** (grid selector armed on measured evidence), **V2.6**
  (measured; defaults deliberately unchanged), **Step 4** (the `cHamiltonian` + complex-accelerator public
  doors, and `SolidCalculation` — the named home for above-SCFIterator decisions).

---

## High-level goals (user)

- I am mostly concerned about abstract interfaces.  These are what the client code for any module
  is supposed to see.  This project is an attempt (probably never been done before) to capture one
  set of high-level interfaces that support SCF electronic structure calculations for any
  structure: atoms, molecules and 3D lattices (also 1D polymers, 2D graphene), using any basis
  set.  In other words 95% of the high-level abstract interfaces (charge density, Hamiltonian,
  orbitals, wavefunction, SCF iterator, accelerators) should be structure neutral and basis-set
  neutral.  Even the basis set interfaces (qcBasisSet) are (were) structure neutral.
- I have so far identified 4 libraries that have a structure-neutral face with specific structure
  specializations:
  1. qcStructure: Structure has derived classes for Atom, Molecule, UnitCell.  Lattice_3D and
     ReciprocalLattice have no inheritance relation with Structure.
  2. qcSymmetry: separate folders for Atom/Molecule/Lattice_3D symmetry types.  Spin degrees of
     freedom are orthogonal to spatial symmetry ... but that only works for non-relativistic irreps.
  3. ElectronConfigurations: we need separate classes for non-aufbau (specific # of electrons per
     irrep) filling.  We also have separate classes for Atom, Molecule and Crystal aufbau filling
     ... can these be combined into one general aufbau filler?
  4. qcBasisSet -> qcAtom_BS, qcMolecule_BS, qcLattice_3D_BS
- With the introduction of lattice calculations some structure-specific classes have been creeping
  into the qcChargeDensity and qcHamiltonian libraries.  Also possibly into qcWaveFunction and
  qcOrbitals.  If we can find ways of getting these refactored back to more generic status, that
  would be a big win.

## SOLID principles

1. SRP: Single Responsibility Principle
   - This principle is more about concrete classes than it is about interfaces.
   - We have a group of classes called Evaluators that do the low-level work of evaluating all
     required integrals, op(r)/grad(r) and a couple of other minor things, for a particular basis
     function type.  These have multiple responsibilities but they are mostly built up through a
     network of mixins in order to achieve the end result.  The objective with the Evaluators is to
     have a well-defined list of functions you need to implement in order to have the framework
     just work for an SCF calculation.  So we are not trying to satisfy SOLID at this level of code.
   - Mixins are a very powerful technique for building up a final concrete class.  We can apply the
     SRP to the mixin components, but obviously it does not make sense for the final class.  It is
     by definition multi-responsibility.
2. OCP: Open Closed Principle
   - We are still in the R&D stages for this project so abstract interfaces will be modified.  We
     just need a good reason to do so.
   - For example if we extend a structure-neutral interface in order to get some Lattice_3D feature
     working, we need to ask the question: what does this change mean for atoms and molecules?
3. LSP: Liskov Substitution Principle
   - Virtual functions that default to some sort of "not implemented" behaviour violate the LSP.
4. ISP: Interface Segregation Principle
   - There is always an urge to add getters and setters.  For setters, always ask: is this
     something I can set at construction time?  For getters always ask: what are we going to do
     with the getter data?  Why not ask the owning class to do that task instead (maintain
     encapsulation).
   - I do not claim that all the abstract interfaces in qchem are truly segregated.  A good example
     is the irrep basis set interfaces.  They do three things:
     a. Deliver integral tables (interface mixin from multiple).
     b. Evaluate op(r) and grad(r) (interface mixin from VectorFunction<T>).
     c. Expose the irrep symmetry, as an abstract interface pointer.
   - a & b are built up from many mixin interfaces.  ISP?  I think so.
   - Any time our abstract interfaces get augmented, those new functions should be part of an
     existing responsibility, not introducing a new one.
5. DIP: Dependency Inversion Principle
   - This is possibly the most powerful concept, and often enables adherence to the previous 4
     principles.
   - This applies equally well for C++ classes, C++ modules (DAG enforced by compiler) and
     libraries (DAG enforced by linker).
   - Example: lib A depends on lib B.  B needs a way to send info to A.  Create an abstract
     interface AI (ha ha) in B; classes in A derive from B::AI; when instances of the A classes are
     passed to B (as an AI*), B can then call back into A.
   - src/ChargeDensity/ChargeDensity.C tStatic_CC, tDynamic_CC are a working example that the
     Hamiltonian library classes derive from and pass back into the ChargeDensity library.
   - As well as being a dependency inversion, this probably has a "Gang of Four" pattern name.
     (Claude: the closest GoF name is **Observer** — a callback interface owned by the lower layer;
     in modern terms simply a "callback interface".)

### A heuristic that earned its keep (2026-08-07 session)

**When a candidate "special case" turns out to be SIMPLER than the general case, it is not a special case
-- it is the general case with terms set to zero.  When it needs machinery the general case does not have,
it is real.**

Four times in one session the apparent asymmetry was an artifact of naming or placement, not physics:
- `InsertStandardTerms` -- the axis was never double/dcmplx; `Ham_PP` is `<double>` and was already on the
  other side.  The generalization extended FURTHER than the code assumed (R2.8).
- `Symmetry::SymOp` -- already generic; only its ADDRESS was crystal-specific (R2.17.3).
- `PW_Hartree` -- not "a Hartree term that also does electron-ion", just two terms wearing one name (R2.14).
- a site-adapted mesh for an ATOM -- looked like a degenerate case wanting a stub; is the EASIEST genuine
  case (one orbit, τ=0, A=I, no torus), so the implementation is three lines and correct (R2.17.3).

The control case, so this does not decay into "everything generalizes": **`PW_XC`'s seed exception is
REAL.**  A matrix-free density has no D, so ρ_DM=φᵀDφ does not exist -- no amount of zeroing produces it,
and it needed actual machinery (the route latch, R2.16).  That is what a true special case looks like.

Practical use: when an item asks "does this extend to all structure types?", write out what the candidate
exception would COST.  Cheaper than the general case ⇒ fold it in.  More expensive ⇒ it is genuinely
different, and the capability belongs on the types that have it (the `tSpinResolved_CD` idiom), NOT on the
base with a stub.

There are 6 other principles related to package cohesion and coupling.  My experience is that the
first 5 are the most commonly misunderstood and misapplied.

---

# Worklist (status-organized, 2026-08-05)

Legend:
- **READY** — claim verified against the tree; fix direction clear; executable without a design
  session.  (Correctness-adjacent items front-loaded.)
- **VERIFY** — needs a design discussion, a measurement, or a repro before acting.
- **DECIDED-ELSEWHERE** — the governing decision/record lives in another doc or pin; engage it
  there rather than relitigating here.
- **DONE** — kept for the record.

Per-item status markers (added 2026-08-06), so status is visible where you READ the item rather than
only in the LANDED section below:
- **✅ DONE `<sha>`** — landed and green on branch `solid-cleanup`.
- **⏳** — written and building, but not yet fully verified.
- **🔶** — the DESIGN is settled (usually by a user ruling) but no code has been written yet.
- unmarked — untouched.

Execution venue (user, 2026-08-05): a branch in the second working tree **~/Code/qchem1**, so the
MnO campaign proceeds undisturbed in qchem6.

## READY

### R1 — correctness-adjacent (do these first)

- **R1.0b SHARED-RADIAL (SP/"L"-shell) support in the Gaussian94 reader — USER IDEA 2026-08-06.**
  A Cartesian d shell carries an s-type contaminant r²e^{−αr²}, which is what made the Mn GPW basis
  rank-deficient (λmin 1.15e-07, cond 8.2e7 — see `GPW_SCF.MnAtomInBoxDChannel`).  The cure that landed
  is empirical (drop the s window to 2 functions and let the contaminants span the rest).  The user's
  structural alternative: generate bases at the **valgen** stage with the SAME RADIAL SET SHARED ACROSS
  l — the standard trick Pople calls an **"SP shell"** (a.k.a. **"L shell"**, e.g. 6-31G), generalized
  here to s+d.  CP2K's format encodes it natively (one set with `lmin..lmax`, e.g. `2 0 1 4 4 4`) and the
  MOLOPT families share one exponent set across s, p AND d.
  **BLOCKED ON:** the flagged `PG_Cart::IrrepBasisSet` bug — the Gaussian94 reader MERGES same-exponent
  shells across l, which is why valgen carries the standing "KEEP EXPONENTS DISJOINT ACROSS l" rule.
  Multi-l shells must be representable before shared radials can be emitted.
  **MEASURE, DON'T ASSUME:** sharing exponents does not by itself remove the contaminant redundancy —
  the contaminant is r²e^{−αr²} (n=2) while an s at the same α is e^{−αr²} (n=0), i.e. independent
  functions; the measured cause was the NUMBER of s functions whose span mimics r²e^{−αr²}.  Shared
  radials buy compactness + CP2K-comparable structure; whether they also buy conditioning is an
  experiment (`GPW_SCF.MnAtomInBoxDChannel` + the vet's λmin/cond readout is the instrument).
  **ALSO WORTH FIXING WHILE THERE:** our GPW path is Cartesian-only — SPHERICAL bases throw
  ("not a molecular Gaussian basis (no `Molecule::LatticeSum1E`)", the parked S3b work).  Spherical d
  has no contaminant at all and is what CP2K actually solves in (its log for our own shell list: 55
  Cartesian vs 47 spherical functions), so wiring the spherical lineage into the lattice sums would
  retire this whole class of problem AND make the CP2K comparisons apples-to-apples.

- **R1.0 The "FAKE RADIAL" `op(r)` on atomic bases — USER DIRECTION 2026-08-06, own session.**
  An ATOM irrep block's 3-D face `operator()(rvec3_t)` returns the purely RADIAL chi_i(|r|) with the
  irrep's Y_lm silently omitted.  Any consumer that quadratures it as a real 3-D function gets nonsense
  when the integrand carries angular structure — this produced the occupied-d KB defect (l=0 projector
  leaking into EVERY l block; every l>=1 projector integrating to ~1e-33), invisible for the whole life
  of the PP code because MnO is the first system with OCCUPIED d projectors (fix c2d86ec9;
  doc/SymmetryUpgradePlan.md §7 step 7).
  *Interim fix that landed:* `BasisSet::ImplicitAngular_IBS` (`ImplicitL()`/`RadialValues(r)`) makes the
  fakeness explicit in the type system and `PP_NonLocal` cross-casts to it for a per-l radial assembly —
  that CONTAINS the trap, it does not remove it.
  *The real fix (user), two steps:* **(1)** consumers stop touching `op(r)` directly — code that only
  wants an INTEGRAL asks the basis for it, so the basis owns how its own functions are represented (the
  CLAUDE.md "prefer classes to do/answer high-level operations" bias; kills this whole bug class).
  NAMING (user, 2026-08-06): NOT a new `Integrate` verb — these are **`MakeOverlap` / `MakeOverlap3C`
  OVERLOADS**, joining the existing `Integrals_Overlap` family, with the arity following the call site:
    * \f$\langle\chi_i|f\rangle\f$ → a VECTOR → a `MakeOverlap(f, mesh)` overload
      (`PP_NonLocal`'s projection vector b_i).
    * \f$\langle\chi_i|f|\chi_j\rangle\f$ → a MATRIX → a `MakeOverlap3C(f, mesh)` overload
      (`PP_Local`'s V_loc block — f in the operator/C slot, exactly the existing 3C shape).
  Call sites on the raw-op(r) mesh route today: `PP_Local::CalculateMatrix`
  (`WeightedOverlap(mesh,*bs,VlocField)`), `PP_NonLocal::CalculateMatrix` explicit-angular branch
  (`Overlap(mesh,*bs,BetaYlmField)`), `Fit_IBS` (`Overlap(itsMesh,*this,f)`), plus the XC/fitting paths
  sharing the WeightedOverlap shape.
  **(2)** then make `op(r)` HONEST (real Y_lm linear combinations), nothing depending on the fake any
  more.  Design question to settle FIRST in that session: an atomic block carries l plus an m-LIST and
  the spherical solver keeps m-degeneracy in the OCCUPATIONS, so decide what `op(r)` returns per (i,m)
  before touching the density/XC paths, which currently rely on the radial-only form.
  *On completion:* RETIRE `ImplicitAngular_IBS` and collapse `PP_NonLocal::CalculateMatrixRadial` back
  into the single 3-D route.  Keep `A_PP.PerLKleinmanBylanderOracle` (s/p/d/f + cross-l zeros) and the
  CP2K atom oracles as the refactor invariant — they pin the physics independently of which route
  computes it.

- **R1.1 ✅ DONE `06e23f5d`. `FittedVxcPol::GetEnergy` clobbers `te.Exc`** — `te.Exc = 0.0;` before delegating.  **→ doc/CleanupHistory.md**
- **R1.2 ✅ DONE `06e23f5d` (with one CORRECTION, below). `=` vs `+=` on `EnergyBreakdown`.**  Assigners:.  **→ doc/CleanupHistory.md**
- **R1.3 ✅ DONE `38a1ebd6` — fixed via `Clone()`, not stopgapped. `UnitCell` SLICING copy** — SCFIterator.C:163.  **→ doc/CleanupHistory.md**
- **R1.4 ✅ DONE `72fecf8d` (THROW, not assert — deviation explained). Silent zero `Gradient()` overrides** — FourierMixCD.C:75 and IrrepCD<dcmplx>.  **→ doc/CleanupHistory.md**
- **R1.5 ✅ DONE `72fecf8d`. `tChargeDensity::EvalBatch` duplicates `ScalarFunction::operator()(rvec3vec_t)`.**.  **→ doc/CleanupHistory.md**
- **R1.6 ✅ DONE `06e23f5d`. `Write()` streams raw POINTERS (hex addresses)** — Imp/FittedVxc.C:111 (`os << itsLDAVxc`).  **→ doc/CleanupHistory.md**
- **R1.7 ✅ DONE `26af31b6`. `SymmetryAdapted_IBS::MakeDirect/MakeExchange` return empty `ERI4{}` silently** —
  split into `Orbital_HF_IBS` (contraction) + `Internal.Orbital_ERI4_IBS` (substrate); the substrate is now
  invisible outside qcBasisSet.  **→ doc/CleanupHistory.md**
- **R1.8 ✅ DONE `06e23f5d`. `FittedVee` casts `bs` and dereferences with NO assert** (Imp/FittedVee.C:41-42) — the.  **→ doc/CleanupHistory.md**
- **R1.9 ✅ DONE `3882938e`. Molecular `BasisSetID()` streamed its SEPARATORS as hex addresses.**
  One `#include <sstream>` in `PGData.C`, plus the SURVEY the item asked for: eight more module TUs
  streamed literals without `<ostream>`/`<sstream>` and now include it explicitly — **including the two
  that were getting it transitively via `<iomanip>`, because relying on a transitive include is the same
  implicit-visibility bet that caused the bug.**  Verified by the route that found it (the ID now reads
  `" PG { Primative 1@(0,0,0.117):S ..."`).  692/692 green; the key STRING changes, which is safe (the
  cache is per-process, the dim guards would catch an under-specific key, and the canonical-pair ordering
  only needs to be consistent within a run).
  **The transferable bit:** `<string>` is NOT enough to make `operator<<(ostream&, const char*)` visible
  with this toolchain — a module TU that streams must include `<ostream>`/`<sstream>` ITSELF.  *(original
  text follows)*
  **Molecular `BasisSetID()` streams its SEPARATORS as hex addresses — FOUND 2026-08-10 while
  writing R1.7's diagnostic.**  `PGData::BasisSetID()` (Molecule/Evaluators/PG_Cart_MnD/Imp/PGData.C:31-40)
  reads `os << " PG { " << *radial << "@" << centre << ":" << pol << " " ... << "}"`.  The OBJECTS print
  fine; every STRING LITERAL comes out as `0x7784d915b66e`-style hex.  Observed verbatim in an exception
  message, so this is measured, not inferred.
  **Mechanism (and it is a C++20 MODULES trap, not the R1.6 one):** `operator<<(basic_ostream&, const
  CharT*)` is a FREE FUNCTION TEMPLATE in `<ostream>`, while `ostream::operator<<(const void*)` is a
  MEMBER.  That TU's global module fragment includes only `<vector>` and `<string>`, so the free template
  is not visible and `const char*` falls through to the `const void*` member.  Members survive; free
  operator templates do not.  Expect the fix to be one `#include <ostream>` (or `<sstream>`, which the
  file also uses without including).
  **Why it matters beyond cosmetics:**
  - **The doc comment on `SymFockCache` claims the opposite** — "a stable string identity, not a raw
    pointer -- deterministic, no void*" (SymmetryAdapted_IBS.C:30-31).  That claim is FALSE today for
    every molecular basis, and it is the comment a future reader would trust.
  - The ID is the cache KEY (`DB_Cache.C:58`) and the ERI4 CANONICAL-ORDER discriminator
    (`a->BasisSetID() <= b->BasisSetID()`).  Correctness holds today only because the cache is RAM-only
    and per-process; string-literal addresses are stable within one run.  **Any disk/persisted cache, or
    any cross-run comparison of keys, breaks the moment it is added** — and ASLR means the key differs
    every run.
  - It defeats every diagnostic that quotes the ID, which is exactly how it was found.
  **Do NOT bundle this with a rename or a key-format change**: fixing the include CHANGES the key string.
  Verify with the existing `IntegralsCache` dim-mismatch guards (DB_Cache_RAM.C:211-232), which exist
  precisely to catch a key that is not specific enough.
  **Worth a SURVEY, not just the one file:** any TU whose global module fragment lacks `<ostream>`/
  `<sstream>` but streams literals has the same silent defect.  R1.6 was the same SYMPTOM from a
  different cause (`qchem::op<<` binding `const Streamable&`), so a grep for hex in test output will
  catch both.

### R2 — mechanical hygiene

- **R2.1 ✅ DONE `06e23f5d`. `tDM_CD::DM_ContractBlocks` → pure virtual.**  The asserting default is DEAD — all three.  **→ doc/CleanupHistory.md**
- **R2.2 ✅ DONE `48e25b74`. Collapse `Kinetic` + `PW_Kinetic` → `Kinetic<T>`.**  Both are 0.5×(kinetic matrix);.  **→ doc/CleanupHistory.md**
- **R2.3 ⛔ WITHDRAWN — NOT a free dedup; re-filed as part of V1.1 (verified 2026-08-07).**.  **→ doc/CleanupHistory.md**
- **R2.4 ✅ DONE `38a1ebd6`. Stale-comment/import batch**: Band_DFT_IBS.C header claims PlaneWave_IBS implements it.  **→ doc/CleanupHistory.md**
- **R2.5 ⚗️ HAMILTONIAN HALF DONE (3 of 5 sites); the 2 ChargeDensity sites remain. `exit(-1)` in library code → throw** (5 sites): ~~Imp/LDAVxc.C:33-43 (dies with R2.6),
  Imp/FittedVxcPol.C:49-53, Imp/VxcPol.C:41~~, and `tPolarized_CD::MixIn`/`GetChangeFrom`
  (Imp/ChargeDensity.C:149-167 — also an LSP narrowing: accepts any tDM_CD, requires Polarized).
  Contrast the correct pattern at Imp/HF_HT.C:28-30.  These kill the pybind GUI / test runner.
  Fold into the D7 cast-survey custom-exception work.
  **DONE 2026-08-07:** LDAVxc's two died with the class (R2.6); FittedVxcPol/VxcPol now
  `throw std::runtime_error` naming the Spin::None-on-a-polarized-term mistake.  `qcHamiltonian` is now
  `exit()`-free.  **STILL OPEN: the two `tPolarized_CD` sites in qcChargeDensity** — left deliberately,
  because their fix is not just a throw: the LSP narrowing (the signature accepts any `tDM_CD` but the
  body requires a Polarized) is the actual defect, and that is V1.6/V1.8's seam.
- **R2.6 ✅ DONE 2026-08-07. The `LDAVxc` bundle** — a "Hamiltonian term" whose `CalcMatrix`/`GetEnergy` call.  **→ doc/CleanupHistory.md**
- **R2.7 ✅ DONE 2026-08-07. `FittedCD::Clone()` — delete.**  Pure virtual (FittedCD.C:28) whose SOLE implementation.  **→ doc/CleanupHistory.md**
- **R2.8 ✅ DONE 2026-08-07. `InsertStandardTerms<dcmplx>` = assert(false)** (Imp/HamiltonianImp.C:49-53) — a.  **→ doc/CleanupHistory.md**
- **R2.9 ✅ DONE `268473b9` (all three sub-items). Small Hamiltonian hardening** — (iii) the whole-basis
  latch now THROWS on change (not the asked-for assert: `-DNDEBUG`); (ii) `tDynamic_HT_Imp_NoCache` keys its
  scratch by `Irrep`, so both siblings of one interface finally promise the SAME reference lifetime — the
  actual defect, and neither of the two fixes the item proposed; (i) `XC_GridEngine` is const + `mutable`
  with a `shared_ptr<const>` holder, and its two non-cross-invalidating rho caches are pinned by asserts.
  **→ doc/CleanupHistory.md**
- **R2.10 ✅ DONE 2026-08-07. `Fit_IBS::SetMesh` → ctor parameter.**  Two-phase construction; the construction-time.  **→ doc/CleanupHistory.md**
- **R2.11 ✅ DONE 2026-08-07. `DB_Cache_RAM.C`** — a screenful of `-Winconsistent-missing-override` warnings on every.  **→ doc/CleanupHistory.md**
- **R2.12 ✅ DONE 2026-08-07. `UnmatchedCounts`/fold `tol` defaults** — 1e-8 fractional as a literal in three places.  **→ doc/CleanupHistory.md**
- **R2.13 ✅ DONE 2026-08-07. Becke strings/labels rename in `Delta_*`/`XC_GridEngine`.**  Verified: the classes are.  **→ doc/CleanupHistory.md**
- **R2.14 Hamiltonian term-naming sweep** (user conventions):
  - Enn/Vnn nuclear-nuclear repulsion; Een/Ven electron-nuclear attraction; Eee/Vee
    electron-electron repulsion; Eex/Vex exchange; Ecorr/Vcorr correlation; Exc/Vxc
    exchange-correlation.
  - For solids these don't work perfectly: at G=0 only certain combinations of Enn+Een+Eee are
    finite.  Also PPs are a combination of Ven+Vee + KB projectors; use Ven for PPs, where n is
    understood to be a shielded nucleus.
  - Renames wanted: `PP_Local`→`Ven_Local`, `PP_NonLocal`→`Ven_NonLocal`, `PW_Pseudo`→Ven-flavored
    name (user suggested Ven_PP), `PW_Hartree`→Vee-flavored name (user suggested Vee_PP) — final
    names = user's call; PW_Kinetic needs no rename (dies in R2.2).
  - **USER 2026-08-07: `Ven_PP` / `Vee_PP` both approved in principle.  Also asked whether PW_Pseudo's
    short/long components should be renamed too, assuming PW_Hartree is strictly a PP term and strictly
    ee repulsion.  BOTH ASSUMPTIONS ARE WRONG — verified against PWTerms.C before renaming anything:**
    - **`PW_Pseudo` has NO long-range component.**  Its matrix is `MakeLocalPotentialShort` + (optional)
      `MakeSeparablePotential` (Imp/PWTerms.C:50-53) and its energy is `te.Een` + the SHORT G=0
      alignment.  The LONG-range local part was moved OUT to `PW_Hartree` — the CP2K local-PP split
      (doc/GPWPlan.md 0e-PP): the deep-well erf potential folds into the ONE G-space Poisson solve
      instead of a per-orbital-pair sharp-field sweep.  So there is no short/long pair to rename here;
      the term is entirely short+nonlocal.  `Ven_PP` alone would read as "the whole PP", so prefer
      `Ven_PP_Short` (or `Ven_PP` + a doc line saying the long part lives in the electrostatics term).
    - **`PW_Hartree` is NOT strictly a PP term.**  `itsLocal` may be NULL — a pure all-electron / no-PP
      run gets plain Hartree with no core-charge fold (PWTerms.C:75-89 and the ctor doc).  A `_PP`
      suffix would be false for that configuration.
    - **`PW_Hartree` is NOT strictly ee repulsion.**  It contributes to BOTH energies (Imp/PWTerms.C:
      135-141): `te.Eee += 0.5*(total-eLong)` (electron-electron Hartree, WITH the ½ double-counting)
      and `te.Een += eLong` (electron-ion long-range, NO ½), plus the LONG G=0 alignment into
      `te.E_alphaZ`.  One Poisson solve, two physical energies — exactly the user's own "for solids
      these don't work perfectly" caveat made concrete.
    - **Naming consequence:** `Vee_PP` is misleading twice over.  Candidates that survive the facts:
      `Vee_VenLong` (says both energies, keeps the Vee/Ven vocabulary), or a what-it-IS name like
      `PW_Electrostatics` / `Poisson_PW` (one Poisson solve; the energy SPLIT is then a documented
      property rather than a promise the name makes).  Recommend `Vee_VenLong`, since the term
      genuinely owns two contributions and the name should not hide the second.  **Final call: user's.**
  - **✅ RESOLVED + DONE 2026-08-07 — USER CHOSE THE SPLIT ("clean regardless of blast radius"), which
    makes the naming question evaporate: with V_long its own term, each term carries ONE energy.**
    - `PW_Pseudo`  → **`Ven_PP_Short`**  (static: V_loc,short + KB nonlocal; + the SHORT G=0 alignment)
    - *(new)*      → **`Ven_PP_Long`**   (STATIC: V_loc,long, the Gaussian core charge; + the LONG G=0
      alignment).  It was always density-INDEPENDENT — `MakeLocalPotentialLong(structure, model)` takes
      no density — so it was a static term living inside a dynamic one.
    - `PW_Hartree` → **`Vee_Hartree`**   (dynamic: V_H[ρ] ONLY — the only piece that depends on ρ)
    - The `if (itsLocal)` runtime test is GONE: `Ven_PP_Long`'s ctor THROWS on a null model, and a run
      with no local PP omits the term (the `Ham_PP` `if (sep) Add(PP_NonLocal)` idiom).  The
      `0.5*(total-eLong)` add-then-subtract-back is gone with it — each term contracts its own matrix.
    - **USER NAMING CONVENTION (2026-08-07, general): a term carrying ONE SIDE of a short/long split
      must SAY SO in its name.**  Applies to any future range-split term, not just this one.
    - Anchor check: the energy re-association (`0.5*(Tr(D(V_H+V_long)) − Tr(D V_long))` → `0.5*Tr(D V_H)`)
      was the flagged risk.  667/667 green, all GPW anchors unmoved.
  - `Delta_XC` → e.g. `DeltaFittedVxc` (it IS a FittedVxc with a δ-function fit basis); rename
    along with the V2.1 decision.
  - **`PW_XC` — the LAST `PW_`-prefixed term after the 2026-08-07 renames.  USER ASKED (2026-08-07): is
    "PW" correct, and does it mean the PW ORBITAL basis or the PW FIT basis?  Answer: the FIT basis — and
    the name is wrong the same way "Becke" was in R2.13, i.e. attached to the wrong noun.**  Verified:
    - **fit basis: genuinely must be orthonormal G-space.**  `Fitting::Factory(cFIT_SF_ABS)` ASSERTS
      `bs->isOrtho()` (Imp/FunctionFitter.C:56) and returns a `GriddedScalarFitter`, which OWNS the FFT
      quadrature grid — `Grid()` needs a `G_FieldEvaluator`.  So the fit is a projection (no metric solve)
      AND the quadrature is the FFT on that basis's own raster; both follow from the fit basis.
    - **orbital basis: NOT plane-wave.**  `CalcMatrix` needs only a `Band_FT_IBS` (G-space 3-centre
      tensors), which BOTH `PlaneWave_IBS` and `GPW_IBS` implement.  Decisive: `GPW_IBS::
      CreateVxcFitBasisSet` returns a `PlaneWaveFit_IBS`, so GPW — GAUSSIAN orbitals — feeds a plane-wave
      fit basis to this term.  (The Answered-questions section already recorded "yes, GPW uses it".)
    - **density:** a `FourierDensity`.  Also not a PW-specific requirement.
    - **USER FOLLOW-UP: "if it only requires fit_bs->isOrtho() then `OrthoFittedVxc` is better -- I want to
      be VERY precise about what the term actually needs."**  Right instinct, but checking inverted the
      answer: the term needs TWO INDEPENDENT capabilities and `isOrtho()` is the WEAKER one.
      - `isOrtho()` -- the METRIC axis (the projection IS the fit).  Genuinely general: an orthonormal
        WAVELET basis (BigDFT/Daubechies) satisfies it.  So ortho ≠ PW in general — the user's "only ortho
        fit basis in the universe?" question answers NO.
      - `G_FieldEvaluator` -- the QUADRATURE axis (the FFT raster the fit is sampled on and E_xc integrated
        on).  **RECIPROCAL-SPACE BY INTERFACE, not merely by implementation**: its vocabulary is `ΔG_Map`,
        \f$e^{i(B\Delta m)\cdot r}\f$, `ForwardFFT`, and `GridCoeff(Vt, ivec3_t dm)` keyed by an INTEGER
        reciprocal-index difference.  A wavelet basis could not implement it.  **This is the BINDING
        requirement.**
      - ⇒ `OrthoFittedVxc` would name the EASY half and drop the HARD one; `PWFittedVxc` names the hard one
        and implies the easy one for everything that exists.  Keep PW.
    - **✅ The better fix, done 2026-08-07: the contract was enforced in TWO PLACES.**  `Factory(cFIT_SF_ABS)`
      asserted `isOrtho()` at CONSTRUCTION while `OrthoScalarFitter::FitGrid()` asserted the
      `G_FieldEvaluator` at FIRST GRID USE -- so an ortho-but-not-G-space basis constructed happily and
      tripped later, somewhere else.  Two-phase contract, the same smell as R2.10's `SetMesh`.  Both checks
      now sit in `Factory`, where the object is built, with the two axes named.  The class NAME then only
      has to distinguish siblings; it does not have to carry the contract.
    - **Proposed name, symmetric with the `Delta_XC` line above: `PWFittedVxc`.**  Then the family reads as
      "FittedVxc + WHICH FIT BASIS": `FittedVxc` (Gaussian aux) / `DeltaFittedVxc` (δ-functions) /
      `PWFittedVxc` (plane waves) — and "PW" modifies the noun it is actually true of.
    - **Sequence:** do it WITH the `Delta_XC` rename (V2.1), not before — renaming one of a matched pair
      leaves the family less consistent than it is now.
- **R2.18 ✅ NAMES DONE `86c5b24d`; the encapsulation half deliberately left open (user ruling below).
  The `Make`/`Get` cached-accessor PAIR is the project convention — qcHamiltonian was the one
  library that spelled it differently, and inconsistently with itself (USER, 2026-08-10:
  *"GetMatrix goes through caching ... if there is no cache it calls MakeMatrix() which does return by
  value.  So if you need a return by value override it should be the MakeXXX() call."*).**
  - **The convention, stated:** `GetXxx()` is the CACHED accessor and returns a REFERENCE; `MakeXxx()` is
    the uncached compute and returns BY VALUE.  A caller that wants an owned copy asks the `Make` half —
    it does NOT ask the `Get` half to change its return type.  qcBasisSet follows it everywhere:
    `MakeOverlap`/`Overlap`, `MakeKinetic`/`Kinetic`, `MakeNuclear`/`Nuclear`, `MakeRepulsion3C`/
    `Repulsion3C`, `MakeDirect`/`Direct`, `MakeExchange`/`Exchange`, `MakeCharge`, `MakeNorm`,
    `MakeInvOverlap`, `MakeRestMass` — 11 `Make` verbs, no exceptions.
  - **qcHamiltonian spells the same role two other ways:** `tStatic_HT_Imp::CalculateMatrix` and
    `tDynamic_HT_Imp::CalcMatrix` / `tDynamic_HT_Imp_NoCache::CalcMatrix` (HamiltonianTerm.C:47,80,126).
    Two names for ONE role in one file, and neither is the project verb.  Rename both to `MakeMatrix`.
    Mechanical (compiler finds every override), but it touches every concrete term, so it wants its own
    commit rather than riding on a behaviour change.
  - **✅ USER RULING 2026-08-10 on the two halves — they have OPPOSITE priorities.**
    - **NAMES: high priority, do it.**  *"Yes we should clean up the names in qcHamiltonian ... consistency
      is high priority."*  ✅ DONE `86c5b24d`.
    - **ENCAPSULATION (public vs protected `Make`): low priority, DELIBERATELY LEFT OPEN.**  *"All the
      MakeXXX() functions were originally protected.  For DFT the 3C versions (MakeOverlap3C,
      MakeRepulsion3C) still are.  It seemed like these were purely internal functions ... but that turned
      out to be incorrect in some cases.  Anyway I have no strong policy on this (encapsulation level)
      right now.  Maybe the right policy will emerge as we refactor.  My intuition says that it is a low
      priority decision."*
    - **So the history is the opposite of what the item assumed:** protected was the ORIGINAL state
      everywhere and the public ones are the DRIFT — each one a case where "purely internal" turned out to
      be wrong.  `MakeOverlap3C`/`MakeRepulsion3C` are the surviving originals, and qcHamiltonian's
      `MakeMatrix` is not an outlier at all; it is simply un-drifted.  **Do NOT "fix" the visibility to
      match qcBasisSet** — that would be standardising on the drift.
    - **Do not open this as its own item.**  The policy is expected to EMERGE from refactoring: when a
      `Make` has to go public, note WHY (which client needed the by-value form and why the cached one would
      not do).  Those reasons are the evidence a policy would be made from; the decision is cheap to defer
      and expensive to guess.
  - **ONE STALE COMMENT LEFT BEHIND ON PURPOSE:** `src/SCFIterator/Imp/SCFIterator.C:186` still says
    "Vxc::CalcMatrix".  It is a comment only, and SCFIterator is on the MnO campaign's DO-NOT-TOUCH list,
    so it was NOT edited — reaching into their working set for a comment is not worth a collision.  Sweep it
    when that list is released.
  - Found while writing up R2.9(ii), where "return by value" was considered as a fix to `GetMatrix` — the
    convention says that was the wrong half of the pair to reach for.

- **R2.19 ✅ DONE `86c5b24d`. `FittedVxcPol` copied a matrix its child already owned — its own HF twin
  `VxcPol` had the fix, three files away.**  Found 2026-08-10 following the user's R2.18 remark.
  - `FittedVxcPol::CalcMatrix` (Imp/FittedVxcPol.C:45-73) is a PURE FORWARDER: both branches end in
    `(s==Spin::Up ? itsUpVxc : itsDownVxc)->GetMatrix(...)`, i.e. they return BY VALUE a matrix the
    child's own per-Irrep cache already holds stably.  The `tDynamic_HT_Imp_NoCache` base then stores that
    copy in scratch purely to have something to return a reference to.  One full matrix copy per call, per
    spin, per irrep, per SCF iteration, to satisfy a signature.
  - **`VxcPol` — the polarized HF term, the direct analogue — already does it right** (Imp/VxcPol.C:39-51):
    it overrides `GetMatrix` and RETURNS THE CHILD'S REFERENCE.  Same `Spin::None` throw, same
    `Polarized_CD` cross-cast, no copy, no scratch, no `NoCache` base.  The two polarized wrappers sit in
    one library and disagree; the copying one is the outlier.
  - **Fix:** give `FittedVxcPol` the `VxcPol` shape — override `GetMatrix`, forward the child's reference,
    drop the `tDynamic_HT_Imp_NoCache` base.  The seed fallback (spin-agnostic density → the unpolarized
    child block) forwards a reference just as well.
  - **`FittedVcorrPol` must KEEP `NoCache`** — verified: its `CalcMatrix` genuinely computes (it refits
    `itsVcFitter` per spin and returns `Overlap(dftbs)`), because ONE fitter is shared by both channels, so
    the result cannot be cached across a spin pair.  So this is not "delete NoCache"; it is "stop using it
    where a forwarder was meant".  R2.9(ii)'s Irrep-keyed scratch still earns its place for
    `FittedVcorrPol`, and would become that class's private business alone.
  - Cheap to verify: the polarized molecular DFT tests (`M_DFT*`/`M_Calculation` polarized cases).
  - **What landed `86c5b24d`:** `FittedVxcPol` now overrides `GetMatrix` and returns the child's reference;
    it no longer derives from `tDynamic_HT_Imp_NoCache` and has no `MakeMatrix` at all — it computes
    nothing, so it has nothing to Make.  That is the cleanest confirmation the Get/Make split was the right
    lens: a pure forwarder has a `Get` and no `Make`, and the old code had to invent a `Make` (and a scratch
    slot to hold its result) purely to satisfy the base it should not have had.
  - **Left as-is deliberately:** `FittedVcorrPol` KEEPS `tDynamic_HT_Imp_NoCache` — verified it genuinely
    recomputes (one `itsVcFitter` shared across both channels, refit per spin), so R2.9(ii)'s Irrep-keyed
    scratch still earns its place and is now that one class's private business.

- **R2.20 The oracle helpers are production data living in a TEST module (USER 2026-08-12, deferred).**
  `RelativeError` / `RelativeHFError` / `RelativeDFTError` / `RelativeDHFError` sit in
  `IntegrationTests/TestUtils.C` (module `qchem.Unittests.TestUtils`), but they are thin wrappers over
  `thePeriodicTable()`'s NIST/Dirac reference energies — production data, not test scaffolding.
  - **Consequence today:** `CLIapps/scfrun.C` — a shipped CLI, not a test — imports a *Unittests* module
    for them, and because TestUtils is a per-target module FILE SET (not a library), every consumer
    RECOMPILES it: ITMain, scfrun, and now UTSCFIterator.
  - **Fix:** move the four into a real library beside `qchem.PeriodicTable`.  Then scfrun imports no test
    module at all, and TestUtils shrinks to what is genuinely test-only.
  - **NOT a scfrun facade problem** (the premise this item started from): scfrun already drives
    `AtomCalculation` and `Calculation` directly.  The oracle helpers are its ONLY reach into TestUtils.
  - User: *"could live in PeriodicTable ... but like you said, later."*  Deferred, not rejected.

- **R2.16 Construction-time facts re-asked at RUN time (USER PRINCIPLE, 2026-08-07).**
  **User ruling:** *"I much prefer that the whole Hamiltonian is decided and fixed at construction time.
  The only dynamic aspect is the ChargeDensity that we feed it."*  Survey done while splitting the PW
  electrostatics terms; the sites fall into three groups.
  - **(A) Branch on the DENSITY — legitimate** (the density IS the dynamic input): `XC_GridEngine::
    Rho`/`RhoPol`'s three-way dispatch (a pure capability query, fine); `FittedVxcPol`/`FittedVcorrPol::
    CalcMatrix`'s `dynamic_cast<Polarized_CD*>` seed fallback (a POLICY, not a query — the ρ↑=ρ↓=ρ/2
    collapse; defensible but worth revisiting when the seeds are spin-resolved).
    **One is NOT benign: `PW_XC::RefreshRhoGrid`** sets `itsRhoIsRaw` per iteration from what the density
    answers, and `CalcMatrix` then picks the RAW vs BALL adjoint — its own comment calls BALL
    "NON-variational".  So the FUNCTIONAL BEING MINIMISED can change mid-run, hidden behind the
    density-is-dynamic exemption.  ✅ **ADDRESSED 2026-08-07 (latch + throw); see the analysis below.**
    - **What the two routes ARE** (the user asked; worth writing down once).  Both produce ρ(r) on the XC
      grid, by different routes.  **BALL**: take ρ̃(G) on the finite {G} sphere and inverse-transform.  The
      truncation rings (Gibbs), so ρ goes NEGATIVE on sharp products — tripping the XC ρ>0 guard — and the
      round trip is a projection, so H_xc ≠ ∂E_xc/∂D: **non-variational**.  **RAW**: collocate
      ρ_DM(r)=φᵀDφ directly in real space (D-weighted level densities combined spectrally, zero-padded,
      no ball restriction).  D is PSD ⇒ ρ_DM ≥ 0 pointwise, and `applyRawAdjoint` is the EXACT transpose
      of `applyRaw`, so H_xc = ∂E_xc[ρ_DM]/∂D to machine precision: **variational**.  They are not two ways
      to compute one number — they minimise DIFFERENT functionals.
    - **The capability is construction-time**, as the user hoped: `GPW_Evaluator::Overlap3CTensor` sets
      `applyRaw`/`applyRawAdjoint` UNCONDITIONALLY (Imp/Evaluator.C:393-394) and a plane-wave basis never
      does (GMap.C:182: "plane waves have no raw representation").  Checked: `RasterFields::HartreeOnly`
      only tunes grid sharpness, it does NOT suppress the raw pair.  So route == orbital-basis lineage.
    - **But a PURE construction-time choice is impossible, for a real physical reason: the SEED.**  A
      matrix-free density (iteration 0 — SeedCD/PolarizedSeedCD) has no D, so ρ_DM=φᵀDφ does not exist and
      it can ONLY answer BALL.  That is inherent to seeding, not a design defect, and iteration 0's energy
      is discarded anyway.  Any construction-time route must therefore still exempt the seed.
    - **What landed instead — the property that actually matters.**  `PW_XC` now LATCHES the route on the
      first DENSITY-MATRIX-backed density and THROWS if it ever changes after that; the matrix-free seed is
      explicitly exempt.  So the seed→SCF transition is allowed (and is the only allowed one), while any
      mid-SCF flip — a mixer whose raw shadow "late-activates", a composite where one block lacks raw —
      becomes a loud error instead of a silent functional swap.  That is the guarantee the user's principle
      was really asking for; the remaining freedom is exactly the freedom physics requires.
    - **Still open (smaller now):** a `route` ctor argument would let a caller FORCE ball (an A/B instrument
      — today only the `GPW_XCROUTE` env var reports the route, it cannot select it).  That needs a
      capability question on the neutral `Band_FT_IBS` face so the factory can ask without a concrete cast.
      Worth doing when someone actually wants the A/B.
  - **(B) Construction-time facts re-asked per call — the actual violation.  ✅ ALL CLEAR 2026-08-07.**
    - `if (itsLocal)` in the old `PW_Hartree::LongBlock`/`GetEnergy` — died with the term split.
    - `if (itsSep)` in `Ven_PP_Short::CalculateMatrix` — the KB projectors became their own term,
      **`Ven_PP_NonLocal`** (ctor throws on a null model; a local-only PP omits the term).  This also makes
      the PW lineage mirror the molecular `PP_Local`/`PP_NonLocal` pair it had drifted from.
    - `!isFinite()` in both PP terms' `GetEnergy` — the G=0 alignment coefficient is now computed ONCE in
      the ctor (`itsAlphaZ`, exactly 0 for a finite structure, which IS the correct value — no faking), and
      `GetEnergy` just scales it by the current electron count.  Bonus: `SumFormFactors` no longer re-runs
      every iteration.
    - `IonIon` — **and this one was hiding a real cost.**  `te.Enn = NuclearRepulsion(...)` recomputed the
      full EWALD LATTICE SUM on every `GetEnergy`, i.e. every SCF iteration, to return the same number: E_nn
      depends only on geometry and ion charges, both fixed at construction.  Now evaluated once in the ctor.
      **Also fixes a gap in R1.2**: IonIon was on that item's ASSIGNER list but is absent from the landed
      `06e23f5d` note — its `te.Enn =` survived the `=`→`+=` sweep.  Now `+=`.
    - Left alone deliberately: `IonIon::Write`'s `isFinite()` branch picks a DESCRIPTION string
      ("pair sum" vs "Ewald lattice sum").  Cosmetic, no behaviour rides on it.
  - **(C) The GOOD pattern, already in the tree**: `Ham_PP` (Imp/Hamiltonians.C:111) —
    `if (sep) Add(new PP_NonLocal(...))`.  Absence of a capability means the term is NOT IN THE LIST.
    Decided at construction; zero runtime tests.  This is the target shape for every (B) site.
- **R2.17 `UnitCell::CreateSiteAdaptedBeckeMesh` — the name carried three things it should not
  (USER CRITIQUE 2026-08-07).**  ✅ TWO OF THREE DONE; the third is a design call.
  1. ✅ **"Becke" merely repeated `mp.cellKind`** — and the body ASSERTED that value, i.e. the NAME was
     carrying a precondition the PARAMETER already states.
  2. ✅ **"SiteAdapted" named one of two STRATEGIES**, which is an implementation detail the caller should
     not be choosing.  Both points fixed by making it an overload: `CreateIntegrationMesh(mp, ops)` —
     the presence of \a ops is what distinguishes it, which is what a signature should say.
     **Bigger win than a rename:** the caller (`GPW_IBS`'s XC-quadrature factory) was branching on
     `mp.cellKind` to pick site-adapted-vs-group-average — a BASIS deciding how a STRUCTURE builds its own
     mesh.  That branch moved INTO the overload, so the basis now just asks for "a mesh invariant under
     these ops" and the cell owns how.  Same altitude error as V1.10b's mixer and R2.16's runtime probes.
  3. ✅ **DONE — USER CHOSE (a)+(c) 2026-08-07.  `std::vector<Symmetry::Lattice_3D::SymOp>` is
     structure-specific, but a site-adapted MOLECULAR mesh is equally plausible** (user).
     **(a) the TYPE is promoted:** `SpinAction` + `SymOp` moved out of `Symmetry/Lattice_3D/Fold.C` into a
     new root module `qchem.Symmetry.SymOp` (`src/Symmetry/SymOp.C`, namespace `qchem::Symmetry`), beside
     `Irrep.C`/`Spin.C` — the root holds what all three structure families share.  `Lattice_3D` ALIASES
     both, so every existing `Symmetry::Lattice_3D::SymOp` spelling still compiles unchanged (23 files
     untouched).  `UnitCell::CreateIntegrationMesh(mp, ops)` now takes the neutral spelling, which was the
     whole point: the signature is no longer crystal-specific.
     **(c) the METHOD stays on `UnitCell`:** there is no molecular implementation yet and inventing an
     unused one would be speculative.  When one is wanted, this signature is already the neutral one to
     hoist onto `Structure`.
     **WHEN THAT HOIST HAPPENS — reasoned through 2026-08-07, so it need not be re-derived:**
     - **There is NO free generic default on `Structure`.**  `MakeInvariant` (the group-average route the
       uniform branch uses) calls `Wrap01` on every image (SymmetrizeMesh.C:216) — it folds on the
       FRACTIONAL TORUS.  For a molecule/atom that is simply wrong: a mesh point at 3.7 Bohr would wrap to
       0.7 of a cell that does not exist.  So the base cannot offer "build the plain mesh, then average
       it"; each structure family must implement its own.
     - **Do NOT give `Atom` an ignore-ops-and-warn body** (the shape first proposed).  It is the LSP hole
       this document forbids ("virtual functions that default to some sort of 'not implemented' behaviour")
       and that R1.4 / R1.7 / V1.6 / V1.7 are all instances of.  A warning is also the wrong instrument:
       warnings are for caller ERRORS, and passing ops to an atom is not an error — a symmetry-broken or
       maximally-stretched atom has a real finite point group.  The warning would say "I ignored what you
       correctly asked for".
     - **And it is unnecessary, because `Atom` is the EASIEST genuine case, not a degenerate one.**  One
       centre ⇒ one orbit trivially; τ=0 (a point group fixes the origin); no torus metric; so the ONLY
       thing ops can affect is the angular set.  The whole implementation is `MakeInvariantAngularMesh(ops,
       L)` in place of the default angular quadrature — ~3 lines, and CORRECT.  (Today
       `Atom::CreateIntegrationMesh` is a one-liner onto `MakeMolecularMesh`, which at natom==1 is "just the
       shifted product grid".)  It is strictly less work than the crystal case, which needs orbits, τ, the
       torus stabilizer test and the bond-direction screen.
     - **So two honest shapes, both stub-free:** (1) hoist to `Structure` and let `Atom` implement it for
       real; or (2) do NOT declare it on `Structure` — declare it on the structures that have it, the
       `tSpinResolved_CD` cross-cast-capability idiom.  Weigh (2) seriously: the atomic solver exploits
       sphericity through IRREPS (l,m), not through mesh symmetrisation, so a caller handing ops to an
       `Atom` may never materialise — the same (c) reasoning that deferred the molecular implementation.
     **User notes worth keeping:** τ=0 for molecular point groups is fine ("we are not fighting
     performance or RAM problems with this code") — a point group fixes a point, so it HAS no translation
     part, and a consumer that Cartesianises via A·W·A⁻¹ needs no special case because a molecule's A is I.
     And on atoms: a finite op list genuinely cannot represent a closed-shell atom's continuous O(3)
     symmetry — "we are into Lie groups" — but nothing in the code asks it to.  Discrete ops are exactly
     right for a symmetry-broken/stretched configuration or for a site group inside a crystal; the
     continuous case is served by the `Symmetry::Atom` spherical machinery, which works in (l,m) instead of
     enumerating operations.  That caveat is now recorded on the struct so nobody later tries to enumerate
     O(3) for an atom.
     *(The findings that shaped the choice:)*
     - The builder needs BOTH the linear part and a translation (`op.W`, `op.tau`) — a screw axis or glide
       plane has a nonzero τ that decides which atoms share an orbit.  So a plain
       `std::vector<Matrix3D<double>>` of Cartesian rotations LOSES information the crystal needs; the
       neutral type has to be the (W, τ) pair.
     - Which is exactly `Lattice_3D::SymOp` minus its namespace — and it already works for a molecule:
       τ=0, and the Cartesianisation `A·W·A⁻¹` is the identity when A is (`Molecule`'s A is I).  So the
       STRUCT is already neutral; only its ADDRESS is not.
     - **The decision is therefore a qcSymmetry organisation question, and it is the user's:** the doc's
       own high-level goal says qcSymmetry has "separate folders for Atom/Molecule/Lattice_3D symmetry
       types", so promoting `SymOp` to a neutral home cuts across that taxonomy.  Options: (a) a neutral
       `Symmetry::SymOp` above the three folders, with Lattice_3D aliasing it; (b) leave the type where it
       is and give `Structure` a virtual taking it (qcStructure already depends on qcSymmetry); (c) leave
       as-is until a molecular site-adapted mesh is actually wanted.
     - Note the site-stabilizer TEST also differs (torus metric mod 1 for a crystal, plain distance for a
       molecule) — but that is implementation, and belongs in each override, not in the argument type.
- **R2.15 ✅ INTERFACE DONE 2026-08-07 (`nAngular` → `angularDegree`, behaviour-preserving); ONLY THE
  DEFAULT FLIP REMAINS (decision 6 below, deliberately separate -- it changes every unpinned run).
  What landed:** the field is a DEGREE for every scheme; Lebedev resolves it through `ResolveLebedev`
  (round UP to the cheapest tabulated rule of at least that degree, ANNOUNCING any substitution and
  naming the four deliberately-excluded orders); `LebedevMenu()` is the one place the ladder is written
  down and every audit test now reads it, so a table and a test can no longer disagree; `LebedevAngular`
  is exported so rule-level audits address a RULE while the interface speaks DEGREES.
  **`ClassifyOrbits` answers "how do we capture the same-degree tuples in code?" (user):** it MEASURES
  which high-symmetry directions a rule occupies (⟨100⟩/⟨110⟩/⟨111⟩) from the directions themselves
  rather than annotating per rule -- the same discipline the degree measurement forced, for the same
  reason.  Two facts only visible once measured: rules 12 and 24 occupy NO high-symmetry direction at all
  (pure general orbits, so `angRot` has nothing to steer for them), and the ⟨110⟩ column is sparse and
  non-monotonic up the ladder (50, 146, 170, 194, 434).
  **A latent bug the change fixed:** `GPW_SCF_UT.C` set the knob to `L` under a scheme chosen by
  `GPW_BECKE_ANG` -- GaussLegendre by default, Lebedev on override.  With the calibrated L=29 the override
  asked for a 29-DIRECTION Lebedev rule, which does not exist, so it would have hit `default: assert(false)`.
  The A/B instrument was broken for its main setting; one meaning for both schemes fixes it.
  **Migration hazard worth remembering:** a blanket rename turned one site's "24 directions" into
  "degree 24", which resolves to a 302-point grid -- a 12x change.  Every Lebedev site needed the
  count→degree CONVERSION, not a rename.  That is exactly why renaming the FIELD (rather than
  redefining `nAngular` in place) was the right vehicle: the compiler forced all ~20 sites into view.
  *(original item text and the decision list follow)*
  **Superseded groundwork note: the degree MEASUREMENT test; the interface change + default flip
  still need the decisions below.  `nAngular` → degree-typed angular interface.**
  **USER 2026-08-07: agrees degree should be the canonical knob.**  Four things still to settle, then the
  work is mechanical:
  1. **Resolution policy: round UP.**  The available Lebedev degrees are sparse and irregular
     (0,1,3,3,5,7,8,11,15,17,19,21,23,29,31,35), so a requested degree usually is not in the table.
     Round-up is the only defensible rule -- it guarantees AT LEAST the requested exactness; nearest or
     round-down can silently under-integrate.
  2. **The degree collisions.**  `nDir` 1 and 2 both sit at the bottom (degrees 0 and 1) and 6 and 8 BOTH
     measure degree 3 -- so degree→count is not invertible.  Round-up + cheapest-wins resolves it
     (degree 3 → the 6-point rule), which makes the 8-point rule unreachable.  Nothing uses it; confirm.
  3. **`EulerMaclaren` has NO degree at all -- so `angularDegree` would be a LIE for it.**  Its θ rule is a
     transformed trapezoid (the `m∈{1,2,3}` clustering kills endpoint derivatives -- the Euler-Maclaurin
     trick), which buys asymptotic CONVERGENCE for smooth integrands, never exactness.  MEASURED degree is
     **−1**: it cannot integrate even the constant to 1e-10, because its weights only sum to 4π
     approximately.  The tree already said so in three places (MolecularMeshTests.C:216, GPW_SCF_UT.C:2089,
     MeshPrimitives.C:130).  So the item's premise "a DEGREE for GL/EM" is wrong -- there are THREE
     semantics on one field: Lebedev=COUNT, GL=DEGREE, EM=RESOLUTION.  Decide: exclude EM from the
     degree-typed face, or RETIRE it (nothing selects it in production -- it appears only in its own
     tests, and GL already gives arbitrary L WITH exactness at the same ~L²/2 product-grid cost).
  4. **✅ EM RETIRED 2026-08-07 (user ruling): "we should just retire EulerMaclaren".**  Gone entirely --
     the builder TU, the enum value, the `em_m` parameter, the `,em` field of `MeshParams::ID()` (the RAM
     cache is per-process, so no cross-run key consequence), the factory case, the report name and the
     three tests.  Both survivors now have a real polynomial degree, which is exactly what lets the knob
     become degree-typed: the three-semantics problem collapses to two, and both are degrees.
  5. **"Lebedev will usually be the most efficient and therefore the most important option -- remaining
     decisions should prioritize getting Lebedev as sensible and honest as possible" (user, 2026-08-07).**
     The MEASURED menu, which is what those decisions should be made against:
     | nDir | 1 | 2 | 6 | 8 | 12 | 24 | 30 | 50 | 86 | 110 | 146 | 170 | 194 | 302 | 350 | 434 |
     |------|---|---|---|---|----|----|----|----|----|-----|-----|-----|-----|-----|-----|-----|
     | degree | 0 | 1 | 3 | 3 | 5 | 7 | 8 | 11 | 15 | 17 | 19 | 21 | 23 | 29 | 31 | 35 |
     - **The ladder has GAPS, and they are principled: 9, 13, 25, 27 are missing.**  Those are the orders
       excluded by the generator audit -- 74 (deg 13), 230 (deg 25), 266 (deg 27) carry NEGATIVE WEIGHTS,
       and the deg-9 32-point rule was removed for a weight-sum bug.  So a request for degree 25 must jump
       to 302 (degree 29): a **55% cost jump** over the 194-point rule it skips.
     - **⇒ HONEST means the resolver ANNOUNCES the rounding** when requested ≠ delivered, and says WHY the
       gap exists.  Silently substituting a 55%-more-expensive grid is exactly the kind of invisible
       decision this whole session has been removing.
     - **6 and 8 share degree 3 and are NOT redundant.**  They differ in ORBIT DIRECTION -- 6 puts points
       on the \f$\langle100\rangle\f$ axes, 8 on the \f$\langle111\rangle\f$ body diagonals.  That
       distinction is load-bearing for the site-adapted work (§6a: a special orbit lying along a structure
       axis is the thing `angRot` exists to steer away from), so cheapest-wins resolution must not be sold
       as "8 is dominated".  Document the pair; keep both.
     - Keep the table's `degree` field = the CONSTRUCTED, guaranteed degree (29 for 302, 35 for 434) even
       though the monomial scan over-delivers -- guaranteeing less than you deliver is honest; the reverse
       is not.
  6. **ICOSAHEDRAL RULES ARE CRYSTALLOGRAPHICALLY FORBIDDEN -- and the code already handles it, by
     construction rather than by luck (user observation 2026-08-07).**  The 12-direction rule is the
     icosahedron, whose 5-fold axes NO Bravais lattice has, so it can never be invariant under a crystal
     point group.  Using it to seed an IMPOSED-symmetry mesh would silently break the T2 invariance
     precondition.
     - **Verified it cannot happen:** the imposed path never touches the Lebedev tables.
       `MakeInvariantAngularMesh` builds an invariant set from scratch, from a deterministic FIBONACCI
       SPHERE seed pool chosen precisely because it "never lands on symmetry axes" -- then symmetrises
       under the site group.  So the crystallographic constraint is satisfied structurally.
     - **On the FREE-run path there is no invariance requirement**, and Leb-12's genericity w.r.t. a cubic
       lattice is a FEATURE: it cannot accidentally put a quadrature point on a bond axis.
     - Worth keeping visible because a reader could reasonably assume any Lebedev rule can seed an imposed
       mesh.  It cannot, and the icosahedral rule is the sharpest illustration of why.
  7. **The default flip is SEPARATE -- and it was BLOCKED BEHIND D6, which was not obvious until located.
     ✅ D6 LANDED 2026-08-07, so this is now UNBLOCKED: the flip is a one-line edit to
     `qcMesh::BeckeXCParams` in src/Mesh/Mesh.C, awaiting only the measurement below.**
     ~~The "free-run Becke default" is not a library default at all: it lives in `BeckeXCParams()` in
     **IntegrationTests/GPW_SCF_UT.C:189** (`GPW_BECKE_L`, default 29, GaussLegendre).~~  That WAS D6's
     complaint -- "the de-facto PRODUCTION Becke recipe living in the integration-test harness".  So:
     - flipping it then meant editing a TEST FILE to change production behaviour, which is backwards;
     - and it still moves every pinned GPW anchor, since the recipe feeds `ResolveXCMesh` for all of them.
     **Order was: D6 first (recipe moved to the library beside `MeshParams`), then flip.**
     - **The measurement, per the D8 standing pin** ("fit quality is measured by grid-convergence of
       ρ/property vs a fine reference -- NEVER ΔE_total"): compare Leb-302 against GL-29, both against a
       FINE reference (Leb-434 or GL-35), on a property rather than a total energy.  The instrument now
       works: `GPW_BECKE_ANG=lebedev` was broken until R2.15 gave both schemes one meaning for the knob.
     - The case for the flip is the measured 302 vs 450 directions at equal degree 29 (**33% fewer
       points**) -- but that is the COST side; the measurement is what shows the accuracy side is equal. -- it changes every unpinned run's numbers.  Land the type change
     first (behaviour-preserving), then flip with the measurement.
  **The rename is what makes the migration safe:** `nAngular` → `angularDegree` breaks every designated
  initializer `.nAngular=`, so the compiler forces a visit to each call site; each Lebedev site then
  converts to the degree reproducing today's grid exactly (table lookup ⇒ bit-identical), and GL/EM sites
  are unchanged.
  **✅ GROUNDWORK LANDED: `Mesh_AngularDegree` (4 tests, src/Mesh/tests/MeshPrimitives.C).**  Degree-typing
  makes each rule's degree LOAD-BEARING -- it stops being a comment and becomes the contract -- so the
  degrees are now MEASURED, monomial-by-monomial against the closed-form sphere integral
  (\f$\int x^ay^bz^c d\Omega\f$), with no spherical-harmonic dependency.  It immediately found **two
  understated annotations**: `nDir=6` claimed L=1 but is the classical degree-**3** octahedral rule, and
  `nDir=2` claimed L=0 but is degree **1**.  Under degree-typing both would have pushed callers to a more
  expensive rule than needed.  Corrected in place.
  - Contract is `EXPECT_GE`, not `EQ`: round-up needs "at least D", so over-delivery is not a defect -- and
    a monomial scan CAN exceed the constructed degree, because monomials odd under a rule's octahedral
    symmetry vanish identically on both sides and never discriminate (302 measures 31 vs its stated 29;
    434 measures ≥40 vs 35).
  - Also confirms the item's headline INDEPENDENTLY: Leb-302 and GL-29 both reach degree 29, at 302 vs 450
    directions = **67%**, exactly as claimed.
  - Worth having regardless of R2.15: this file's own header records a 32-direction rule shipped VERBATIM
    from the old library with sum W = 0.971·4π, found and removed by hand.  A degree measurement catches
    that class of defect on the first run.
- **R2.15 (original text) `nAngular` → degree-typed angular interface.**  `nAngular` is a COUNT for Lebedev but a
  DEGREE for GL/EM (and the imposed site-adapted builder consumes it as the degree) — the dual
  semantics BLOCKS flipping the free-run Becke default to the measured-equal Leb-302 (67% of
  GL-29's directions).  Fix: `angularDegree` + per-scheme count resolution; the default flip rides
  along.  (The warn/auto-resolve ergonomics half of this item stays in D5.)

## VERIFY

### V1 — interface-design questions

- **V1.1 `Orbital_DFT_IBS` ⇄ `Band_FT_IBS` merge.**  User (2026-08-05): Orbital_DFT_IBS simply
  specifies what integrals an IBS must supply to support DFT; "Band" and "FT" have no place in that
  specification.  History: Band_FT_IBS once had ~12 extra Fourier/Grid members, twiddled down over
  many sessions to essentially ONE — the main conceptual gap was reluctance to acknowledge that the
  PW representation of ρ is a fit, just a trivial one whose expansion coefficients come in one
  step.  Verified current state: (i) the scalar-type mismatch (Orbital_DFT_IBS<T> templated but
  fit-face args hard-wired REAL, rFIT_*_ABS in Fit_IBS.C:45-46,87-88, vs Band_FT_IBS hard-wired
  dcmplx) — user: trivial to fix, and the FIT_*_ABS<T> templating is already pinned with the fitter
  work; (ii) `ERI3<T>` vs `G_ERI3` return types (G_ERI3 was built as the harmonizing data-structure
  spec); (iii) extras: the using-decl is noise; `CreateXCQuadrature` (new, 2026-08 W1) is
  HOISTABLE to the neutral DFT face (molecules implement it returning their Becke quadrature — a
  unification, not a divergence); leaving exactly ONE genuinely Fourier member,
  `MakeOverlap(f(G))` (Band_FT_IBS.C:77) — itself re-expressible through the fit abstraction
  (fitted-potential coefficients → contraction), removing "FT" from the face.  Engage the pinned
  FACTOR-not-FUSE analysis (fitting-boundary pin + doc/FittingCleanupPlan.md) — the argument-type
  question is settled there; what remains is execution sequencing with the fitter templating.
  **Absorbed from the withdrawn R2.3 (2026-08-07):** the 4-line `Overlap3C`/`Repulsion3C` cache-lookup
  bodies in Imp/Band_FT_IBS.C:15-25 and Imp/Orbital_DFT_IBS.C:10-20 are the SAME code modulo exactly the
  two blockers listed above (return type `G_ERI3` vs `ERI3<T>`, argument `cFIT_*` vs `rFIT_*`) plus
  `theCache<dcmplx>` vs `theCache<T>`.  They are the smallest concrete instance of the merge, so they make
  a good FIRST target once those are decided — and a good litmus test that the decision actually works.
  **USER (2026-08-05): merge would be fantastic; THE one big remaining issue = for molecules we do
  a Dunlap fit for ρ (Coulomb metric, Repulsion integrals, charge-constrained) while for solids we
  don't (orthonormal projection ⇒ metric degenerate) — requires discussion.**  Groundwork for that
  discussion already exists in the fitting-boundary pin: the metric axis is REAL for AO (two
  different solves) and DEGENERATE for FT (projection IS the fit for both metrics), which is why
  the PW fitter implements both metric faces at once — the merge discussion is "how does the
  merged face express the metric choice without naming it", not "which metric wins".
- **V1.1b 🔶 ANALYSIS DONE, awaiting the user's re-read of the paper. The `Eee = 2·EeeFit − EeeFitFit` expression is DUNLAP-SPECIFIC — it is part of V1.1's
  metric discussion, not a free-standing formula (user, 2026-08-05; user wants to re-read the
  paper).**  Verified conventions: `GetSelfRepulsion()`=½⟨ρ̃|ρ̃⟩ (Imp/FittedCDImp.C:55) and
  `eeeFit`=½⟨ρ|ρ̃⟩, so the expression is exactly Dunlap's ROBUST form
  \f$E_J\approx\langle\rho\tilde\rho\rangle-\tfrac12\langle\tilde\rho\tilde\rho\rangle\f$, whose
  error is SECOND order in (ρ−ρ̃) **only because the fit is done in the COULOMB metric** — the
  stationarity condition ⟨ρ−ρ̃|c⟩_Coulomb=0 is what annihilates the first-order term.  The ctor
  says so out loud: "Charge-CONSTRAINED Coulomb-metric density fit (Dunlap-Connolly-Sabin 1979)"
  (Imp/FittedCDImp.C:21-23).
  - **Under an OVERLAP-metric fit the robustness justification evaporates** (stationarity is then
    ⟨ρ−ρ̃|c⟩_overlap=0, which does NOT kill the Coulomb first-order error) — so a "straight rho fit"
    needs a different energy expression.  User's instinct confirmed.
  - **Under an EXACT/orthonormal fit it degenerates CORRECTLY**: ρ̃=ρ ⇒ ⟨ρρ̃⟩=⟨ρ̃ρ̃⟩=⟨ρρ⟩ ⇒
    2A−B = ½⟨ρρ⟩ = the exact Hartree.  So the PW/GPW path is safe as-is; it is an overlap-metric
    GAUSSIAN fit that would break.
    **CORRECTION (user, 2026-08-07): "EXACT/orthonormal" conflates two INDEPENDENT axes, and the slash
    hides an assumption that is FALSE for GPW.**  Orthonormality of the fit basis means S=I, so the normal
    equations collapse to c=⟨f|ρ⟩ — that is about COST and CONDITIONING.  Whether ρ̃=ρ is a separate
    question, about whether ρ lies in span{G}, and orthonormality says nothing about it.  Split by lineage:
    - **PW orbitals: genuinely exact.**  ρ=ψ*ψ is band-limited to the difference set {G_i−G_j}, and
      `PlaneWave_IBS::CreateCDFitBasisSet` builds the CD fit basis at the 4× cutoff that covers it
      (`Ecut()*std::max(4.0, mp.relCutoff)`, Imp/PlaneWave_IBS.C:112-121, whose comment states the rule).
      So ρ̃=ρ holds by construction, and the degeneracy argument above is sound.
    - **GPW (Gaussian orbitals): NOT exact.**  A Gaussian product has infinite bandwidth, so a finite {G}
      ball truncates it — which is exactly what `ReportGridCharge`/`GridChargeLost` MEASURES (CP2K's
      "Electronic density on regular grids" line).  ρ̃≠ρ there, so "the PW/GPW path is safe as-is" is
      overstated for the GPW half: its Hartree is a genuine approximation, exact only as the density
      cutoff grows.  (It does NOT use the 2A−B form, so no robustness claim is being violated — but the
      reason it is safe is "it never uses the Dunlap expression", not "its fit is exact".)
    - **On "rank 2 represented by rank 1" (user's intuition):** the represented object is ρ(r), the
      DIAGONAL ρ(r,r) — genuinely a function of one point, so no rank is being squeezed.  D→ρ̃ IS
      many-to-one (it sums D_ab over each difference G_b−G_a), so D is not recoverable from ρ̃; but
      Hartree and LDA need only the diagonal.  The full ρ(r,r′) would be needed for exact exchange — and
      consistently, `IrrepCD<dcmplx>` NA-asserts on the HF accumulators.  The rank collapse is real; it
      just costs nothing for the functionals this path supports.  **This is exactly V1.1's "molecules Dunlap-fit, solids don't"
    issue seen from the ENERGY side** — whatever the merged face does about the metric, this
    expression moves with it.
  - **An overlap-metric fit path ALREADY EXISTS**: `NumericCD` (the SAD seed) overrides
    `GetUnconstrainedFit` with c₀=S⁻¹⟨f|ρ⟩ (NumericCD.C:48-51).  Today that is harmless and the
    invariant holds — but only BY ACCIDENT OF TYPING: `FittedVee::GetEnergy` and
    `tSCFIterator::TotalEnergy` take `tDM_CD*`, and the matrix-free seeds are `tChargeDensity` but
    NOT `tDM_CD`, so a seed can never reach the energy expression (it only feeds `CalcMatrix`,
    where the metric does not matter — J(ρ̃) is just the fitted density's potential).  Nothing
    states or enforces this; widen the energy path to `tChargeDensity`, or make a matrix-free
    density a `tDM_CD`, and the formula silently loses its second-order property.
  - Action: name the invariant where it lives (the fitter/energy pair), and fold the "which metric
    ⇒ which energy expression" rule into V1.16's explicit metric-strategy face.  Also ties the
    third leg of the CD taxonomy (D2's seeds) to a real correctness boundary.
- **V1.2 🔶 FEASIBILITY PROBE DONE + user APPROVED attempting it (appendix below). `Orbital_PP_IBS` — a structure-neutral PP-integral face (dependency INVERSION).**  User
  framing (2026-08-05): PPs require certain NEW TYPES of integrals from the IBS; the question is
  whether there is a structure-neutral way to ask for them without spilling PP details — if yes, we
  can break the qcBasisSet(qcLattice_BS)→qcPseudopotential dependence (today PlaneWave_IBS/GPW_IBS
  implement `Pseudopotential::Integrals_Pseudo<dcmplx>` whose args are PP types).  Candidate
  neutral primitives: (1) ⟨i|V|j⟩ for a species-attached local radial field (lattice-summed; the
  long/short split as a range parameter, not a PP concept); (2) projector brackets ⟨i|β_lm⟩ for
  species-attached radial×Y_lm functions.  qcPseudopotential then calls the face from ABOVE —
  the dependency edge inverts.  Passes the pseudo-wall pin (these ARE new integral types).
  Molecular PP (term-side quadrature today, by design) could optionally adopt the face later.
  Related dead surface: `Integrals_Pseudo::MakeLocalPotential` (unsplit matrix) is documented
  unit-test-only.  **USER (2026-08-05): approved to attempt — try it and see if we hit any
  roadblocks.**  (Feasibility probe of the actual `LocalPotential`/`SeparablePotential` payloads
  and the DAG: see the probe notes appended below when available.)
- **V1.3 ✅ MECHANISM DONE `72fecf8d` (both ε-adapters deleted; via `GetEMatrix`, NOT `DM_ContractBlocks` — see "WHAT
  LANDED INSTEAD" below).  The quadrature-term face (second list) is still open.  `FittedEpsXc`/`FittedVxc` simplification.**  Physics answer (user asked 2026-08-05):
  yes, genuinely different matrices — v_xc = δE_xc/δρ = ε_xc + ρ·∂ε_xc/∂ρ, so
  D·⟨i|v_xc|j⟩ = ∫v_xc·ρ ≠ ∫ε_xc·ρ = E_xc (Slater exchange: factor 4/3 — the retired ¾-virial
  shortcut, FittingCleanupPlan §I.1).  A class merge is impossible: both need `GetMatrix` with the
  SAME signature returning DIFFERENT matrices, and `DM_Contract` dispatches through that one
  signature.  The over-complication is real: FittedEpsXc is a full Dynamic_CC with its own fitter,
  version stamp, and cached matrix, to deliver ONE scalar per cycle — and it has already been
  cloned (`FittedEpsCPol`, Imp/FittedVcorrPol.C:78-98).
  **USER DESIGN DIRECTION (2026-08-05): a framework RE-ALIGNMENT for everything ex/corr/xc, in the
  Hamiltonian library, structure-neutral (not PW/GPW-specific).**  All other terms assume
  E_zz = D·V_zz — energy computed from the (possibly cached) matrix V_zz.  The ex/corr/xc terms
  break that assumption (E_xc = D·⟨ε̃⟩ ≠ D·⟨v_xc⟩).  So: a new interface (extension of
  `tDynamic_HT`?) that forces **`GetVMatrix()` and `GetEMatrix()`** and sorts out which gets used
  where (V → Fock assembly; E-matrix → the energy contraction).  This dissolves the DM_Contract
  signature collision (no second Dynamic_CC object needed), subsumes FittedEpsXc AND its
  FittedEpsCPol clone into one seam, and gives GGA/corr a home from day one.  The quadrature
  routes (Delta_XC/PW_XC compute E by grid integral, NO E-matrix) are handled by the user's
  companion design (2026-08-05): **a second variant/alternate of the `tDynamic_HT` face tailored
  for quadrature-evaluated terms; the Hamiltonian keeps SEPARATE LISTS of each term type and
  loops through them — perfect LSP in action** (no term is ever forced to fake a matrix it
  doesn't have; new term types are easy to add as new lists).

  **IMPLEMENTATION PLAN (grounded 2026-08-06, compute-free pass over the term hierarchy).**  Two
  facts found while reading it change the shape of this item:
  1. **The "separate lists per term type" design is ALREADY IN PRODUCTION.**  `tHamiltonianImp` keeps
     THREE lists — `itsSHTs` / `itsDHTs` / `itsHF_HTs` — behind three faces (`tStatic_HT`,
     `tDynamic_HT`, `tDynamic_HF_HT`) with three different `GetMatrix` arities, and both the matrix
     and energy loops simply walk each list in turn (Imp/HamiltonianImp.C:63-82).  So adding a fourth
     kind is idiomatic here, not novel — it is the user's own pattern, already load-bearing.
  2. **A term whose energy is NOT D·V already exists, and it is the precedent to copy.**
     `tDynamic_HF_HT`'s own doc: it "is deliberately NOT a `tDynamic_CC`: its energy comes from its
     OWN cached blocks (`DM_ContractBlocks`), not a per-irrep `GetMatrix` round-trip."  That is
     exactly the E≠D·V escape the xc family needs, and it requires NO new machinery.
  - ~~**Mechanism (small, and it deletes code):** `FittedVxc` builds its ⟨i|ε̃_xc|j⟩ blocks into a map
    and takes E_xc from `cd->DM_ContractBlocks(...)`, exactly as `Vee`/`Vxc` already do.~~
    **THIS ROUTE IS BLOCKED — found on contact, 2026-08-07.**  `Vee`/`Vxc` can build a whole-system block
    map only because they are `tDynamic_HF_HT`, whose `GetMatrix` RECEIVES `wholeBasis` (the composite
    basis — the cross-irrep view).  `FittedVxc` is a `tDynamic_HT`: its 3-arg `GetMatrix` gets ONE irrep's
    basis and never the composite, so it has no way to ENUMERATE the irrep blocks a `DM_ContractBlocks`
    map must contain.  Faking it (latch a `map<Irrep,const odftbs_t*>` from successive `CalcMatrix` calls,
    à la `Dynamic_HF_HT_Imp::itsWholeBasis`) would add exactly the raw-pointer latching R2.9(iii) already
    flags, AND is unsound in timing: `GetEnergy` runs on the density AFTER the step (SCFIterator.C:279
    contracts ρ_out while the Fock was built from ρ_in), so cached ε blocks would be fit against the wrong
    density.  Plus the string→Irrep key change below.  Enumerating the irrep blocks is precisely the job
    `DM_Contract` already does — that callback IS the density telling the term its own leaf bases.
  - **WHAT LANDED INSTEAD — the user's own V/E face, used as the mechanism.**  The collision was never
    `DM_Contract`; it was one SPELLING.  `tDynamic_HT` already IS-A `tDynamic_CC` (Hamiltonian.C:56), and
    both faces named the method `GetMatrix`, so one term could expose only ONE matrix — hence a second
    object to carry the other.  So `tDynamic_CC`'s method is now **`GetEMatrix`** (ChargeDensity.C), with
    `tDynamic_HT::GetEMatrix` defaulting to `GetMatrix` (E = D·V, true for every term but the xc family)
    and `IrrepCD::DM_Contract` calling `GetEMatrix`.  `FittedVxc` and `FittedVcorrPol` now override it
    with their ε fit — so **`FittedEpsXc` AND its clone `FittedEpsCPol` are both DELETED**, and the two
    ε fitters became plain members beside the v fitters (same fit basis, so the 3-centre setup is still
    shared).  Zero new state, no key change, no basis latching, and no numbers move: `DM_Contract` calls
    exactly the same fit+`Overlap` on exactly the same `(bs, spin, cd)` as the adapters did.
    - The V half deliberately KEPT the name `GetMatrix` rather than becoming `GetVMatrix`: that rename is
      a large mechanical sweep (`tStatic_HT`/`tDynamic_HT`/`tDynamic_HF_HT`/`tHamiltonian` + every term
      and call site) with no behavioural content.  Do it as a standalone naming commit if wanted — the
      V/E DISTINCTION is already in the type system and documented on both faces.
  - **Still outstanding (the declarative half, unchanged):** the quadrature terms (`Delta_XC`, `PW_XC`)
    own a mesh, integrate E directly, and have NO E-matrix — they must never be forced to fabricate one.
    That is the SECOND term list / separate face, and it is untouched by this commit.
  - ~~**Snag to fix en route:** `DM_ContractBlocks` is still keyed by `std::string` BasisSetID~~ — moot for
    V1.3 now (nothing here touches `DM_ContractBlocks`), but the observation stands on its own and is worth
    keeping: V1.4 moved `DM_RhoAtPoints` to `Irrep`, and `tHT_Common`'s term cache was ALREADY
    `std::map<Irrep,hmat_t<T>>` (HamiltonianTerm.C:23).  So Irrep is the established key on the term
    side and BasisSetID is the odd one out; extend V1.4 to `DM_ContractBlocks` in its own pass.
    (This also retro-justifies V1.4: the term caches had been Irrep-keyed all along.)
- **V1.4 ✅ DONE `80fc2ae8`. `DM_RhoAtPoints` Phi key → Irrep (USER RULING 2026-08-05).**.  **→ doc/CleanupHistory.md**
- **V1.5 `G_FieldEvaluator` 3-face ISP split.**  Verified: 11 public methods with near-1:1
  method↔consumer mapping — `EvalField*`/`ForwardFFT`/`GridCoeff`/`FieldCoeffs` →
  OrthoFunctionFitter only; `RhoOnGrid`/`Integral`/`EmitGridReport` → PWTerms only;
  `ApplySpectralFilter` → DensityMixer only (Kerker); `MakeFourierDensity` → SeedCD only.  Three
  disjoint responsibilities (evaluate-a-ΔG_Map, FFT quadrature, one-off analytic services).
  Caveat: deliberately promoted to the "full PW grid-engine DIP seam" in FittingCleanupPlan §K —
  the split must engage that decision.  Correction to the old bullet: PWVxcField is not a 4th
  consumer (helper inside PWTerms; the Hamiltonian reaches the engine via `FunctionFitter::Grid()`,
  FunctionFitter.C:139).  Related: XC_GridEngine's `Lattice_3D::Fold` + dcmplx dependency bars any
  molecular reuse of the quadrature engine — same design session.
- **V1.6 `tDM_CD::Accumulate*` — face split, NOT pure-virtual.**  Verified override matrix: every
  concrete family relies on the default for exactly 2 of the 4 (IrrepCD lacks `*All`;
  Composite/Polarized lack `*Both`) — pure-virtual just forces 6 new asserting stubs.  The `*Both`
  pair is an internal Composite↔leaf collaboration protocol (only called from
  Imp/CompositeCD.C:48,59) — no business on the public face.  Shape: a whole-system face and a
  pair-partner leaf face (the `tSpinResolved_CD` cross-cast idiom).  NDEBUG hazard on record:
  these void assert-only bodies are silent NO-OPS in Release — a bare IrrepCD through
  `Vee::AccumulateAll` yields a zeroed J and a silently wrong Fock.
- **V1.7 The periodic trio (`GetFourierDensity`/`GetRhoOnGrid`/`GetRepulsion3C`) — 9 asserting
  stubs, the largest LSP block in qcChargeDensity.**  Re-declared + NA-asserted on
  Polarized/Composite/IrrepCD for BOTH T (Imp/ChargeDensity.C:97,120,138; Imp/CompositeCD.C:197,
  226,249; Imp/IrrepCD.C:267,286,303).  The correct mechanism ALREADY EXISTS in the same file —
  `FourierDensityBase<T>` gives dcmplx the capability and double an empty base — but the derived
  classes re-declare outside it and assert.  Fix: declarations live only on the dcmplx side
  (if-constexpr-guarded definitions or a `tPeriodic_CD` mixin).
- **V1.8 `IrrepCD`↔`IrrepCD` concrete same-class casts in the hot path** (Imp/IrrepCD.C:84,98,
  218,227: `Accumulate*Both`/`MixIn`/`GetChangeFrom` take abstract `tDM_CD&` and narrow to the
  concrete leaf to touch `itsDensityMatrix`; the in-file comment names "the IrrepCD↔IrrepCD
  idiom").  Abstract→concrete, the pattern the project rule forbids; also makes MixIn
  unimplementable for any future leaf.  Wants a double-dispatch primitive or an abstract
  density-block face.  (Design with V1.6 — same seam.)
- **V1.9 ✅ DONE `38a1ebd6`. `Structure`→concrete-`UnitCell` down-casts in 4 libraries**.  **→ doc/CleanupHistory.md**
- **V1.10 Two abstract→CONCRETE basis casts in src/** — Imp/SymmetryAdapted_IBS.C:109,118
  (Orbital_HF_IBS* → concrete SymmetryAdapted_IBS, solely to reach `itsO`) and
  Internal/Imp/Orbital_DHF_IBS.C:89,109 (Orbital_ERI4_IBS& → Orbital_RKB_HF_IBS_Imp&).  Both are
  "give me your private state" reaches — promote the needed answer to an abstract question on the
  face.  (Unit-test exemption does not apply; these are src/.)
  *(NOT in this list, and deliberately so: `Orbital_ERI4_IBS::Substrate` added by R1.7 is an
  abstract→ABSTRACT cross-cast — the sanctioned direction — and it THROWS naming both bases.)*
- **V1.10b ✅ DONE (see LANDED). Mixer.  **→ doc/CleanupHistory.md**
- **V1.11 Occupation seam: two-phase SCF-WF construction + the `TakeElectrons*` family.**
  (i) `SCFWaveFunction::Init/SetMOM/SetSmearing/AdoptMOMReference/ReleaseMOMReference`
  (SCFWaveFunction.C:39-76): run-config known at construction time delivered by post-ctor setters;
  every concrete WF defends against the un-configured state.  (ii) Five `TakeElectrons*` virtuals
  on `TOrbitals<T>` (Orbitals.C:128-152) with a dual-meaning tuple slot (the doc comment itself
  warns "here the double is −TS, NOT leftover electrons"); every new occupation policy so far added
  a virtual to the abstract face (OCP).  One `Fill(const OccupationPolicy&)` returning a named
  struct; the policy object also feeds the WF ctor.  Cross-ref: this IS SCFStrategyPlan's
  "occupation seam" — design it there.
- **V1.12 `EnergyBreakdown` — 13 public data members (OCP+SRP).**  Every new term family edits the
  struct + totals + `op+=` + Display; `GridChargeLost` is a GPW health DIAGNOSTIC ("not an energy",
  its own comment) and `MinusTS` is WF-side entropy — both riding the energy value object; also
  lattice-only `E_alphaZ` in the neutral struct.  Candidate: keyed contributions + a small fixed
  set of roles for the totals; move GridChargeLost to the run report/IterationTrace (which already
  carries it).
- **V1.13 ✅ DONE 2026-08-07 — executed as the compiler-verified DELETION R2.6 made possible.  **→ doc/CleanupHistory.md**
- **V1.14 Report-emission creep on neutral faces + global bool toggles** — `EmitBasisUsage`
  (WaveFunction.C:48, defaulted no-op), `EmitRadialReport` (IrrepBasisSet.C:68), `EmitGridReport`
  (G_FieldEvaluator.C:60, PURE — forces every implementor), plus function-local-static
  `bool& ReportBandGap()`/`ReportGridCharge()` process-globals that leak state between tests (the
  SCFIterator comment admits it).  Fix: a reporter/visitor that PULLS; toggles on SCFParams.
- **V1.15 `tBasisSet<T>::Create*FitBasisSet` defaults** — the generic body hard-codes
  `Orbital_DFT_IBS<double>` regardless of T (only the explicit dcmplx specializations save it) and
  derefs an unguarded iterator (a 1E/HF-only basis ⇒ null ⇒ UB in Release).  Also
  `CreateXCQuadrature`'s default body is byte-identical in two places (Imp/BasisSet.C:42 ==
  Band_FT_IBS.C:53).  Decide: pure, or a documented "this basis does not fit" contract; hoist the
  shared default.  (Interacts with the V1.1 CreateXCQuadrature hoist.)
- **V1.16 `ProjectedDensity_AO::GetRepulsion3C` asserting default** (Fitting/Imp/FunctionFitter.C:
  28-35) — the overlap-vs-Coulomb metric choice is encoded as "override the OTHER method and this
  one becomes poison", an implicit unenforceable pairing (NumericCD must just remember).  Make the
  metric an explicit refinement face.
- **V1.17 `tWaveFunction::GetSpinDensity()` returns null as the unpolarized answer**
  (WaveFunction.C:39) — a capability half the hierarchy lacks, on the base, every client
  null-checking a raw pointer.  Correct idiom one library over: `tSpinResolved_CD` as a cross-cast
  face.  Move to a `SpinResolvedWF` face — also aligns with the spin-native-is-primary bias (the
  polarized WF is the primary type, not a special case bolted on via nullable getter).
- **V1.18 `FourierMixCD` tell-don't-ask + `MakeDensityMixer` ISP.**  `RhoTilde()` hands out the
  raw ΔG_Map and PulayMixer runs the whole DIIS algebra outside the density
  (DensityMixer.C:183-233); `SetRawRho` + external `RasterKerker` is a get/compute/set straddle —
  the mixing algebra wants to be density-side ops or a dedicated mixable-density face.  Also
  `MakeDensityMixer` takes `const tDM_CD*` but uses only GetTotalCharge + a FourierDensity cast
  (DensityMixer.C:312-320) — excludes the matrix-free seeds from seeding the mixer BY TYPE, not
  intent.
- **V1.19 Seed assembly: give `Structure` the question.**  Seed code reads concrete `Atom` public
  fields (Imp/Seed.C:102-104 `a->itsZ`,`a->itsR`; Imp/SeedCD.C:91 `itsSpinFlip`) and clones the
  UnitCell into (unflipped, flipped) groups because `MakeFourierDensity(st, formFactor(Z,g2))` is
  species-keyed — the flip-group sub-cell duplication is the SYMPTOM; the neutral face yielding
  raw atoms is the CAUSE.  A `ForEachSite(fn(Z,R,flip))` / per-atom form-factor visitor (or the
  per-atom-index `G_FieldEvaluator` overload — weigh against the pseudo-wall pin) removes the
  field access, the UnitCell casts, AND the sub-cell duplication.  Also `MakeSeedDensity`'s closed
  enum switch has assert(false)+nullptr arms (Imp/Seed.C:152,158) — in Release an unsupported
  (strategy × T) cell returns null and SCFIterator silently runs a core guess;
  compile-time-over-runtime says make unsupported combos unconstructible or throw.
- **V1.20 `SymmetryAdapted_IBS` → Internal: decide the library-family rule first.**  Verified: all
  consumers inside the src/BasisSet/ tree BUT in a different CMake target (qcMolecule_BS).  Moving
  to `.Internal.` creates a cross-LIBRARY internal import under the CLAUDE.md rule.  Precedent the
  rule is already bent: `qchem.BasisSet.Internal.GMap` imported by src/Fitting + src/ChargeDensity.
  Decide: is the qcBasisSet* family "one library" for Internal purposes?
- **V1.21 `BandStructure.C`: promote or demote.**  Confirmed test-only (sole import =
  tests/BandStructureUT.C:13); worse, tests/PlaneWaveUT.C:59 defines its OWN local `SolveBands`
  instead of importing.  Either promote (band plots are on the viz roadmap) or demote into the
  test tree; either way kill the duplicate.
- **V1.22 `MakePeriodicBeckeMesh` ε-tail drops vs orbit consistency (W2c find).**  The builder's
  borderline drop decisions (`<eps` screens + `w>0` keep) are per-point and bit-sensitive, so the
  site-adapted caller post-filters orbit-incomplete points (`CreateSiteAdaptedBeckeMesh`).
  Cleaner: make the drop decision ONCE per representative (angular dir × radial shell) and apply
  it to the whole atom orbit inside the builder — removes the second fold pass + the filter.
- **V1.23 `Symmetry::Lattice_3D::DirectOf`** — currently unused after the CreateXCQuadrature move
  (GPW uses its native direct ops).  Keep (documents the U=Wᵀ convention; T3 will want it) or fold
  its doc into `ReciprocalOp` and drop — decide at the refactor session.
- **V1.24 `GDMParams::FDMax` naming + the fallback commit (2026-08-03 imposed×GDM investigation).**
  **(ii) ✅ DONE 2026-08-09 by the MnO dev — CLOSED.**  The fallback commit is fixed.  Diagnosed here
  2026-08-03 and independently re-derived from the MnO trace six days later, using this item's own two
  listed reproducers — the worklist described the bug before the campaign hit it.
  - **REFINEMENT, because the item's own prescription was half wrong:** it said "hold position on fallback".
    Holding position gives Δρ=0, which MnO's convergence gate reads as CONVERGENCE — a false success, the
    worst failure mode available.  The landed fix DEGRADES TO THE MIXED STEP instead, which is what needed
    the new `RejectStep` contract to be safe (a rejected step is the CALLER's judgement, since only the
    caller can evaluate the energy).  **Generalisable: "do nothing" is not a safe fallback when the
    convergence test reads CHANGE — no-change and converged are indistinguishable to it.**
  **(i) STILL OPEN, now owned by the MnO dev** (`GDMParams` became public in `SCFAccelerator.C` with the
  2026-08-09 merge, which is their file): `FDMax` reads as a step size but is the ENGAGEMENT gate on ‖[F′,D′]‖ (the geodesic step size
  is the quadratic-model `itsStdef` capped by `Trust` radians); rename to something like
  `EngageBelowFD` and consider making the norm intensive (per-element or per-electron) so one value
  transfers across basis sizes.  (ii) `DirectMinStep`'s 12-backtrack line search COMMITS a tiny
  non-descent step on fallback (SCFIterator.C ~L424) — the measured uphill leak on imposed NaF-SR
  (+23–56 mHa over 100 iterations along projector-curved diffuse directions).  Fix: hold position
  on fallback (or accept `best` only within a noise floor); pair with soft-direction
  preconditioning — the 1/(ε_a−ε_i) diagonal Hessian blows up the step exactly along the
  near-degenerate diffuse modes.  **(iii) SOFT-DIRECTION PRECONDITIONING is the remaining half and is still
  open** (MnO dev, 2026-08-09): the new precondition check makes GDM DECLINE when the OCCUPATION is wrong,
  and does nothing for a legitimately soft direction — the other failure.  Reproducers: DISABLED_ImposedGDMProbe_SiDiamondIBZ (healthy),
  DISABLED_NaFImposedGDMSmearProbe (pathological, NAFGDM_* knobs); GPW_GDMTRACE=1 shows
  DESCENT/FALLBACK per step.
- **V1.25 Minor CD-interface trims** — non-const `Polarized_CD::GetChargeDensity(Spin)` overload
  has no external consumer (removable); `tSpinDensity` holds two raw non-owning `tDM_CD*` with
  unmanaged lifetime.

- **V1.26 Uniform-vs-Becke: a SMOOTHNESS question that reduces to a COST CROSSOVER — so `Auto` becomes a
  SELECTOR, an explicit choice is honoured, and only the strictly-dominated choice is warned about.
  USER RULING 2026-08-07 (in two parts), arriving right after D6 landed.**
  > "Deciding between Uniform and Becke grids should be based on overall smoothness (PPs for sure, maybe
  > PPs and orbital basis functions).  PWs with ultra soft PPs can 'get away with' Uniform.  GPW with any
  > high exponents (like F with alpha_max=40Ha) have to (in practice) use Becke.  But there will be a
  > mostly grey area in between.  So this must be a user decision at the highest CalculateSolid level.
  > What the software needs to do is warn the user if a Uniform grid is too coarse to handle the sharpest
  > basis function (squared) or handle the sharpest PP."

  Verified against the tree before writing anything down:

  1. **The physics is quantitatively RIGHT, and it is a COST crossover, not an accuracy trade-off.**  Using
     the code's OWN Nyquist mapping (`UnitCell::CreateIntegrationMesh`, \f$n=\lceil 2a\sqrt{2E}/\pi\rceil\f$)
     against the code's OWN density-scale floor (\f$E=\texttt{cutoffFactor}\cdot\alpha_{\max}\f$, cutoffFactor=2):
     | case | a (bohr) | E=2α_max | n/axis | n³ points |
     |------|---------|----------|--------|-----------|
     | Si diamond, sipp α_max=2   | 10.26 |  4 | 19 |   6,859 |
     | Si diamond, α_max=8        | 10.26 | 16 | 37 |  50,653 |
     | NaF rocksalt, F α_max=40   |  8.70 | 80 | 71 | 357,911 |
     | MnO rhombohedral, α_max=40 |  9.30 | 80 | 75 | 421,875 |
     Becke at the production recipe is **49,384 points** for the 2-atom Si cell (measured this session).  So
     at α_max=40 the uniform grid needs ~7x MORE points than Becke — the user's "have to, in practice, use
     Becke" is not a preference, it is the cheaper grid by a large factor.  Conversely at α_max=2 uniform
     wants 6,859 points, ~7x FEWER than Becke.  **The crossover is computable from α_max BEFORE any grid is
     built**, so the warning can carry the cost comparison, not just a complaint.
  2. **`MeshParams::nUniform`'s default of 20 is BASIS-BLIND, and that is the live hazard.**  Inverting the
     same formula, n=20 resolves only \f$\alpha_{\max}\approx2.3\f$ (a=10.26) to \f$3.3\f$ (a=8.7).  Si/sipp
     (α_max=2) *just* squeaks under it — which is precisely why the uniform XC route looks fine on Si and on
     nothing else.  Any consumer that reaches the uniform mesh without setting `eCut` gets that silently.
  3. **The exact idiom the user is asking for ALREADY EXISTS one grid over** — `GPW_Evaluator`'s ctor
     (Evaluator.C:277-281) warns on `cerr` when an EXPLICIT `densityEcut` is below
     `cutoffFactor*alpha_max`, names α_max and the floor, and **honours the explicit value anyway**
     ("we don't hide it, but we don't silently override the explicit choice either").  That is the user's
     ruling already implemented for the density/collocation grid.  The new work is to apply the same
     idiom to the two floors it does NOT cover (below), not to invent a mechanism.
  4. **Both sharpness sources are already available, in ONE common currency — a Gaussian exponent.**
     - basis: `Lattice_3D::MaxExponent()` (already the density floor's input).
     - PP: `LocalPotential_Gaussian::ShortRangeGaussian(Z)` returns \f$c\,r^{2n}e^{-\alpha r^2}\f$ terms with
       \f$\alpha=1/2r_{loc}^2\f$ — an ABSTRACT capability face (optional, reached by the sanctioned
       abstract→abstract cross-cast; `MultiSpecies_LocalPotential` forwards it).  So no new getter and no
       `r_loc` leak into the abstract interface is needed.  A model with no closed-Gaussian short part simply
       cannot be checked — say so rather than pretend (the LSP discipline this session has been applying).
     - Sharpest local PPs in the SHIPPED GTH database, α=1/(2r_loc²): **Na q9 15.4**, Ne 13.9, Mg 13.5,
       He 12.5, **F 10.5**, O 8.2, Mn q15 3.75, Si 2.6, Na q1 0.64.  Two orders of magnitude of spread, so
       "the sharpest PP" is a real discriminator and not a rounding effect.
  5. **The PP floor is genuinely MISSING, which is the strongest argument for this item.**  The integrand
     is \f$\langle\chi_i|V_{short}|\chi_j\rangle\f$, exponent \f$2\alpha_{\max}+\alpha_{pp}\f$ — but
     `GPW_Evaluator::PPMeshParams()` sets `mp.eCut = densityEcut` (Evaluator.C:817), i.e. the DENSITY floor
     \f$2\alpha_{\max}\f$ with **no \f$\alpha_{pp}\f$ term at all**.  So the user's two criteria are not
     redundant: for a soft basis with a semicore PP (Na q9, α_pp=15.4) the PP binds and the current floor
     misses it entirely; for F (2α_max=80 vs α_pp=10.5) the basis binds.  *Scope note:* that uniform mesh
     has ONE consumer today — the KB-projector grid fallback (Evaluator.C:1119); the analytic
     `MakeOverlap` route above it returns early, so the exposure is the fallback path.
  6. **The facade LAYER exists; only the name did not** (user clarification 2026-08-07: "I just mean the
     equivalent of `AtomCalculation` for solids... whatever code exists above SCFIterator in the DAG,
     `RunGPW`?").  It is **`RunGpw(lat, mol, GpwOptions)`** — the solid twin of `Calculation`/`AtomCalculation`,
     with `GpwOptions` as its `CalcOptions` — and it lives in the anonymous namespace of
     IntegrationTests/GPW_SCF_UT.C.  (`RunGPW` is a positional convenience wrapper over it.)  So this is
     D6's disease one level up: D6 moved the Becke RECIPE into the library, but the OPTIONS STRUCT and the
     DRIVER a user would set `cellKind` on are still in the test harness.  Their home is `src/Calculation/`,
     beside the two existing facades.

  **✅ POLICY SETTLED — USER RULING 2026-08-07, in the cost-crossover framing:** *"in principle we can
  auto-select uniform vs Becke based on the n³ vs n_Becke if the user leaves the default at auto.  Otherwise,
  only warn the user if they are forcing the least efficient option."*  So `UnitCellKind::Auto` SURVIVES and
  gains a real job — it becomes a COST SELECTOR rather than the unconditional `Auto`→Becke that D6 landed —
  and an explicit `Uniform`/`Becke` is honoured, with a diagnostic only when it is the losing choice.
  R2.15 item 7 (GL-29 vs Leb-302) is unaffected: that is an angular-rule choice *within* Becke.

  **The selector is cheap and needs no grid build** — both costs are closed-form:
  - uniform: \f$n^3\f$ with \f$n=\lceil 2a\sqrt{2E}/\pi\rceil\f$, i.e. \f$\propto a^3\alpha_{\max}^{3/2}\f$ —
    grows with cell volume AND basis sharpness.
  - Becke: \f$n_{atoms}\times n_{radial}\times n_{dirs}(\text{scheme, degree})\f$ — **independent of
    \f$\alpha_{\max}\f$**, \f$\propto n_{atoms}\f$.  \f$n_{dirs}\f$ is exact up front: GL degree L gives
    \f$(L+1)^2/2\f$ (=450 at L=29, from `GaussLegendreAngular`'s \f$n_\theta=(L{+}1)/2\f$, \f$n_\phi=L{+}1\f$),
    Lebedev degree 29 gives 302.  Verified against a live run: free Si = 2x40x450 = 36,000; the IMPOSED
    site-adapted mesh measured 886 dirs/atom -> 70,880 before the ~30% tail drop -> 49,384 built (so the
    estimate is a conservative over-count, and the selector must use the run's actual imposed/free mode).
  - The two different scalings are what make the decision ROBUST far out and genuinely grey near the crossover
    — which is the user's own "mostly grey area in between", now quantified rather than felt.

  **THREE REFINEMENTS (Claude, accepted into the design pending user objection) — all one asymmetry:**
  - **(a) Cost at Nyquist parity systematically FAVOURS uniform, so a bare comparison over-chooses it in exactly
    the grey zone.**  The Nyquist \f$n\f$ sizes the grid to resolve the DENSITY (band-limited at
    \f$2\alpha_{\max}\f$).  The XC integrand is \f$v_{xc}(\rho)\propto\rho^{1/3}\f$ — pointwise-nonlinear, NOT
    band-limited (this is why `relCutoff` exists at all, and why the plan calls Becke "the near-ideal grid for
    the one pointwise-nonlinear sharp-at-core term").  So the uniform side of the comparison is a LOWER BOUND
    on what it actually needs, while Becke's radial clustering handles \f$\rho^{1/3}\f$ at the core natively.
    The two sides are not equal-accuracy.
  - **(b) ⇒ a MARGIN, not a tie-break.**  Choose uniform only when it is cheaper by a clear factor (~2x), else
    Becke.  That makes the default SAFE in the grey zone and encodes the user's own phrasing: "PWs with ultra
    soft PPs can *get away with* uniform" is a margin statement, not a coin flip.
  - **(c) ⇒ the diagnostic is ASYMMETRIC.**  Forced Uniform when Becke is cheaper is STRICTLY DOMINATED (more
    expensive AND less accurate) → **warn**.  Forced Becke when uniform is cheaper is merely wasteful, never
    wrong → **info line, not a warning**.  The options are not symmetric, so the diagnostics must not be.
    *This also disposes of the warning-noise worry:* both deliberate A/B tests are Si/sipp, where uniform is
    the cheap side (6,859 vs 36,000) — `DeltaFitUniformGridMatchesPWFit_SiGamma` forces the CHEAPER grid
    (silent) and `BeckeXCMatchesUniformXC_SiGamma` forces the SAFE one (info only).  Neither becomes noise.

  **Where \f$\alpha_{pp}\f$ lands — the cost framing answers the question that was open.**  It is an INPUT TO
  THE SELECTOR, not a separate warning: it raises the uniform grid's requirement (integrand exponent
  \f$2\alpha_{\max}+\alpha_{pp}\f$) and so pushes the crossover toward Becke.  It therefore changes behaviour
  ONLY where the user left `Auto`, which is by definition "you choose for me" — neither a silent behaviour
  change nor a mere diagnostic.  **Keep separate and still open:** `PPMeshParams()` sizing its OWN mesh at
  \f$2\alpha_{\max}\f$ with no \f$\alpha_{pp}\f$ term is an under-resolution that survives whichever grid the
  selector picks (point 5 above).

  **Execution order:**
  1. ✅ **DONE.** Name the Nyquist mapping ONCE — \f$n(a,E)\f$ is a literal inside `CreateIntegrationMesh` and every
     consumer of this item would otherwise duplicate it (the R2.12 "name it once" case) — plus its inverse
     \f$E(a,n)\f$, which is what turns a bare `nUniform` into a comparable cutoff.
     `qcMesh::UniformDivisions` / `UniformCutoff`, in `qchem.Mesh` beside the two `MeshParams` fields they
     relate; `CreateIntegrationMesh` calls the first.  Pinned against the original literal by test.
  2. ✅ **DONE.** A cost-estimate pair (uniform \f$n^3\f$, Becke \f$n_{atoms}n_r n_{dirs}\f$) beside them, so the selector
     and the diagnostic share ONE cost model rather than two that can drift apart.
  3. ✅ **DONE (selector DISARMED — see below).** The selector + the asymmetric diagnostic.
  4. Promote `GpwOptions` + `RunGpw` to `src/Calculation/` so `cellKind` is a documented facade knob (point 6).
     Steps 1-3 did not depend on step 4.

  **✅ LANDED 2026-08-07 (steps 1-3).  681/681 ctest green, warning-free build, ZERO anchors moved.**
  New module **`qchem.Mesh.XCPolicy`** (src/Mesh/XCPolicy.C) + 11 tests (src/Mesh/tests/XCPolicyTests.C).
  - **Policy split OUT of `qchem.Mesh`, which D6 had merged in.**  `qchem.Mesh` is now the geometry-free VALUE
    layer (Mesh, MeshParams, the Nyquist arithmetic) — no I/O, no environment, no opinions; `qchem.Mesh.XCPolicy`
    holds the recipe, the cost model, the selector and the diagnostics.  That is the split D6 should arguably
    have made on day one (the `<cstdlib>`/`getenv` in the core value module was the tell), and it costs one
    import at the two consumers.
  - **The two sharpness sources are read off ABSTRACT capability faces, as predicted, with no new virtuals:**
    `BasisSet::Molecule::LatticeSum1E::MaxExponent()` and
    `Pseudopotential::LocalPotential_Gaussian::ShortRangeGaussian(Z)` (→ \f$\alpha=1/2r_{loc}^2\f$), both by
    sanctioned abstract→abstract cross-cast, gathered by a facade-level `GatherSharpness` helper that states
    FACTS and decides nothing.  It sits beside `RunGpw` and moves with it in step 4.
  - **⛔ THE SELECTOR IS DELIBERATELY DISARMED (`Auto` still → Becke; arm with `GPW_XCGRID_AUTOSELECT=1`).**
    Armed, it would have flipped **Si to UNIFORM** — undoing the 2026-08-01 Becke-default flip, which was
    itself the product of a measurement, and moving every pinned GPW anchor.  Refinement (a) says exactly why
    the model would be wrong to trust here (it sizes for the band-limited density, not for the nonlinear
    \f$v_{xc}\f$), and `kUniformMargin=2.0` is a GUESS.  **Per the D8 pin, a grid choice is measured, not
    reasoned about** — so the verdict is computed and ANNOUNCED on every run and acted on by none.
  - **The announce line IS the calibration instrument, and it already earned its keep** — first harvest,
    from unmodified pinned tests:
    | test | α_max | α_pp | uniform | Becke | verdict |
    |------|-------|------|---------|-------|---------|
    | `SiliconGammaConverges`      |  2 | 2.58 |     4,913 | 72,000 (imposed) | UNIFORM |
    | `NaPseudoAtomInBoxDoublet`   |  2 | 0.64 |    32,768 | 36,000 (imposed) | BECKE   |
    | `O2TripletInBoxMatchesFinite`|  2 | 8.15 |   132,651 | 72,000 (imposed) | BECKE   |
    | `MnAtomInBoxDChannel`        | 36 | —    | 1,860,867 | 18,000 (free)    | BECKE   |
    - **O2 is the vindication of carrying \f$\alpha_{pp}\f$ separately** (point 5): SAME basis as Si
      (α_max=2), opposite verdict, and the ONLY difference is oxygen's PP (α_pp=8.15 vs 2.58) raising the
      required cutoff from 6.6 to 12.2 Ha.  Without the PP term this run would have been scored as "soft"
      and picked uniform.  The two criteria are independent in practice, not just in principle.
    - **Mn is the user's α_max case, quantified: 103x.**  Sharpness moves the uniform side by two orders of
      magnitude while the Becke side does not move at all — the crossover is real and the far field is not
      close.
    - **The Si verdict is the one to distrust, and it is exactly where the margin bites** — 4,913 vs 72,000
      is a 15x claim in favour of a grid the project MEASURED to be the worse one.  Either the margin is far
      larger than 2, or (likelier) point-count is the wrong currency near the crossover.  That is the
      calibration question, now backed by data instead of intuition.
  - **The asymmetric diagnostic behaves as designed on the real suite** — verified by running it over all 20
    GPW tests, not assumed: `DeltaFitUniformGridMatchesPWFit_SiGamma` (explicit Uniform, cheaper AND
    adequately resolved at nUniform=20 ⇒ 9.4 Ha vs 6.6 required) is **silent**; `BeckeXC_IBZ_SiDiamond`
    (explicit Becke where uniform is cheaper) gets the **info note**, not a warning.  Total noise added to the
    suite: **two lines, on one test** — no deliberate A/B drowned, which was the R2.11 risk flagged up front.
  - **⚠️ AND THAT ONE TEST IS A COUNTEREXAMPLE TO REFINEMENT (c)'s RATIONALE — found by the diagnostic, on its
    first run over the suite.**  `SiPseudoAtomInBoxMatchesFinite` (a=16 box, nUniform=20) draws BOTH warnings,
    and the resolution one is a TRUE POSITIVE (dr=0.8 bohr resolves 1.9 Ha where the system needs 6.6).  But
    the test pins Uniform DELIBERATELY, and the file already says why: **Becke is not unconditionally safe.**
    A freely-rotating DEGENERATE density (its half-filled 3p atom) gets orientation-DEPENDENT error on Becke's
    FIXED-AXIS angular grid — v_xc is nonlinear, so an anisotropic ρ's error rotates with it — turning an
    energy-neutral rotation into a **~Ha-scale oscillation**.  Smearing fixes it (fractional occupation
    restores the symmetric density), so the exception is narrow — but it is real, and it is the one case where
    the "Becke is never wrong, only wasteful" asymmetry INVERTS.
    - **Consequence, applied:** the warning now states the COST FACT and names the exception instead of
      prescribing Becke.  Its earlier wording ("prefer UnitCellKind::Becke") was advice this very test
      documents as harmful.
    - **Worth keeping visible:** this is the second time in two sessions that a diagnostic's FIRST run over
      the suite corrected the reasoning that motivated it.  The cheap lesson is to run a new warning across
      everything before believing its premise, not just to check it is quiet.

  **STILL OPEN after this landing:**
  - **Calibrate `kUniformMargin`, then arm.**  The D8-compliant measurement (grid-convergence of ρ vs a fine
    reference — never ΔE_total) on the Si case, which is where the model and the earlier measurement disagree.
    `XCPolicy.AutoStillResolvesToBeckeWhileTheSelectorIsDisarmed` is the guard that must be deliberately
    updated when it is armed.
  - **`PPMeshParams()` still sizes its own mesh at \f$2\alpha_{\max}\f$ with no \f$\alpha_{pp}\f$** (point 5).
    Independent of the selector; unchanged by this work.
  - **The Becke recipe's own defaults are load-bearing for the selector** (user, 2026-08-07: "there is a lot
    riding on the defaults for nRadial and nDirs(degree)... the degree can be determined from point symmetry
    of the atom site, but I have the impression that you can often get away with much lower degrees than the
    point symmetry dictates").  They set the ENTIRE Becke side of the comparison, so an over-generous recipe
    biases the selector toward uniform — the same direction as refinement (a), i.e. the two errors COMPOUND
    rather than cancel.  Filed as V2.6; the warning is on `BeckeXCParams`' doc comment so a reader meets it
    where the numbers live.

- **V1.27 ✅ THE LIVE HALF DONE `a1e1f9bb` (the virial gate/column now turns itself off via
  `IsVirialValid()`); the ITERATOR RENAME itself is still open.  `MolecularSCFIterator` /
  `SolidSCFIterator` are named for the STRUCTURE and discriminate on the PSEUDOPOTENTIAL.  USER 2026-08-09.**
  > "They were originally about what iteration columns to display.  MolecularIterator actually displays
  > columns that make sense for non-PP calculations and SolidIterator displays columns that make sense for
  > PP calculations.  They were originally misnamed simply because most Molecule runs were non-PP and all of
  > our solid runs are PP."
  - **The code already states the PP reason, in `SolidSCFIterator`'s own doc:** the virial is DROPPED because
    "GTH local + KB projectors break the Coulombic-homogeneity assumption behind 2+V/K".  That is a property
    of the POTENTIAL, not of the lattice — a molecular PP run breaks it identically.
  - **✅ USER RULING 2026-08-10 — the remaining half is a DECOMPOSITION, not a rename.**
    *"I suspect we are ultimately going to need all 4 combos of {Molecular,Solid}x{PP,Non-PP}SCFIterator.
    This should be done with mixins."*
    - **Why a rename could never have worked** (this was the confusion in the earlier framing): the two
      classes are not two variants of one thing.  `MolecularSCFIterator` is `tSCFIterator<double>` and
      `SolidSCFIterator` is `tSCFIterator<dcmplx>` — so the Molecular/Solid axis IS the matrix-element-type
      axis, and PP-ness is a SECOND, independent axis that has no representation at all.  `MolecularSCFIterator`
      is moreover an EMPTY subclass whose only job is to name the molecular path.
    - **The two mixin axes are NOT symmetric — measured, not assumed.**  Of the four ways
      `SolidSCFIterator`'s display differs from the base, only ONE is about pseudopotentials:
      | difference | driven by |
      |---|---|
      | no virial column/gate | **PP-ness** |
      | `ΔE/E` gates instead of `Δ[F,D]` | grid/collocation (non-variational SCF) |
      | `ρ_lost/N` (grid-charge leak) | grid/collocation |
      | gap is a PERMANENT column | periodic (near-gapless flapping is a solid pathology) |
      So the axes are **PP-ness** (virial) and **grid/periodic-ness** (the other three) — genuinely
      orthogonal, which is exactly the case mixins are for, and exactly why one inheritance chain could not
      express it.
    - **Note the PP axis is ALREADY runtime-adaptive** via `IsVirialValid()` (the landed half), so the
      mixins are mostly needed for the grid/periodic columns.  A molecular PP run today gets the right
      VIRIAL behaviour with no new class; what it cannot get is a column set that is neither the
      molecular nor the solid one.
    - **The 4th combo is not hypothetical:** an all-electron PERIODIC run (the parked APW/LAPW tests) is
      Solid × Non-PP, and would want the virial back WITHOUT the collocation columns.
  - **⚠️⚠️ BIGGER THAN THE ITEM SAYS — INVESTIGATED 2026-08-10, NEEDS A USER RULING BEFORE EXECUTION.
    The virial is not only a DISPLAY column; it is a CONVERGENCE GATE.**
    - `Imp/SCFIterator.C:336`: `itsConverged = ... && fabs(eb.GetVirial()+idealVirial) < ipar.MinVirial;`
      — a hard `&&`.  So a PP run is gated on a quantity the tree itself documents as meaningless for it.
    - **And the codebase has been working around it BY HAND at every PP call site, nine times:**
      `ValenceBasisGen.C:59,92,113` (`p.MinVirial = 1e30` + the comment *"the virial theorem does NOT hold
      under a pseudopotential"*) and `IntegrationTests/A_PP.C:43,116,134,174,203,308` (`1e10`/`1e30`, *"virial
      off (N/A to PP)"*).  **That is exactly the R2.16 anti-pattern**: a construction-time fact ("this run
      uses pseudopotentials") re-supplied by every caller as a magic threshold, where FORGETTING it is
      silent.  The `Calculation` facade's own `{.pseudopotential=true}` path does NOT set it.
    - **The default is entangled and should be looked at in the same pass:** `SCFParams::MinVirial = 1e-13`
      with the comment *"1e-13 => effectively off; the textbook -V/K=2 virial is not gated for molecules"*.
      But the test is `error < MinVirial`, so a SMALLER threshold is STRICTER — 1e-13 makes the clause
      essentially unsatisfiable, i.e. the gate is maximally ON, not off.  Either the comment or the default
      is wrong.  Every serious caller overrides it (A_HF_dfPin uses 4e-2 for a genuine all-electron pin),
      which is why this has not bitten.
    - **PROPOSED SHAPE (not implemented — the OCP rule says an abstract-interface addition wants a reason
      agreed first).**  Follow the pattern already in the tree rather than inventing one: `Add()` already
      computes `itsIsPolarized`/`itsIsRelativistic` by OR-ing the TERMS' own `IsPolarized()`/
      `IsRelativistic()`, and the iterator already derives `idealVirial` from
      `itsHamiltonian->IsRelativistic()` (Imp/SCFIterator.C:238).  So: add `IsPseudopotential()` to the term
      face (default false; true on `PP_Local`, `PP_NonLocal`, `Ven_PP_Short`, `Ven_PP_NonLocal`,
      `Ven_PP_Long`), OR it up in `tHamiltonianImp::Add`, and let the iterator ask ONCE — dropping both the
      virial COLUMN and the virial CLAUSE.  The nine hand-set thresholds then become redundant.
    - **✅ RULED AND DONE `a1e1f9bb` — and the ruling improved the NAME, which improved the item.**
      *"I suspect PPs are not the only thing in the whole electronic structure universe that breaks the
      virial theorem.  So if we make a new function it should be `IsVirialValid()` ... instead of
      `IsPseudopotential()`.  So then it seems clear `IsVirialValid()==false` should result in no virial
      gate in the SCF iterator loop."*  (The `idealVirial`-returns-NaN alternative was raised and rejected
      by the user as a smell — rightly: it would push a sentinel onto every consumer of `idealVirial`.)
    - **Name the PROPERTY THE CLIENT CONSUMES, not the CAUSE.**  `IsPseudopotential()` would have been
      correct today and wrong at the first non-Coulombic term that is not a PP (an external or model
      potential, a finite field, a cutoff Coulomb).  This is R1.7's lesson a third time (after R2.13's
      "Becke" and R2.17's "SiteAdaptedBecke"), and it is now cheap to state: **when adding a capability
      query, write down what the caller will DO with the answer; if the name does not match that sentence,
      it is naming the implementation.**
    - **The fold is AND, not OR** — the one implementation subtlety.  `IsPolarized`/`IsRelativistic` are
      OR-ed in `Add()` (one term is enough to make the Hamiltonian polarized); validity is conjunctive
      (the virial holds only if EVERY term is Coulombic).  Same shape, opposite operator.
    - **Relativistic is deliberately NOT this flag.**  The Dirac virial is still VALID, it just has ideal
      ratio 1 instead of 2 — which `IsRelativistic()` already selects.  Folding it in would have repeated
      the wrong-noun error the new name exists to avoid.
    - The display methods ask the Hamiltonian directly rather than threading a flag through the virtual
      signature, so no sentinel stands in for "no virial column".
    - **The nine hand-set `MinVirial` thresholds are now redundant** (ValenceBasisGen ×3, A_PP.C ×6).
      Left in place deliberately — harmless, and their comments document the physics — but a later sweep
      can drop them, and any NEW PP caller needs nothing.
    - **STILL OPEN, deliberately separate:** `SCFParams::MinVirial = 1e-13` with the comment *"effectively
      off"*.  The test is `error < MinVirial`, so smaller is STRICTER — the default is maximally ON, not
      off.  Either the comment or the default is wrong.  Untouched here because it affects ALL-ELECTRON
      runs (where the gate is real) and V1.27 was about PP runs; fixing it is a separate decision about
      what the default gate SHOULD be.
  - **⚠️ AND IT IS A LIVE DEFECT, not just a name.**  `Calculation` supports `{.pseudopotential=true}` (the
    `sipp` molecular PP runs) and uses `MolecularSCFIterator`, which inherits the base layout — so **a
    molecular pseudopotential run displays a virial column the tree itself documents as invalid for it.**
    Verified: `src/Calculation/Imp/Calculation.C:29` and `Imp/AtomCalculation.C:26` both take
    `MolecularSCFIterator`.
  - **`MolecularSCFIterator` is an EMPTY subclass** — `using tSCFIterator<double>::tSCFIterator;` and nothing
    else.  It has no molecular behaviour at all; it is `tSCFIterator<double>` under another name, which is
    why the misnaming cost nothing until a molecular PP run existed.
  - **`SolidSCFIterator` conflates TWO axes, and only one of them is structural:**
    - `CreateMixer` (V1.10b) — genuinely PERIODIC: Kerker/Pulay are G-space and need a lattice.
    - the column set — mostly PSEUDOPOTENTIAL (the virial), partly METHOD (ΔE gates instead of Δρ because a
      collocation SCF is non-variational).  The frontier-gap column is the one plausibly solid-specific
      (metals), though a near-degenerate molecular open shell wants it too.
  - **This is R2.8's finding again, one library over** (see the LANDED note): there, `double` vs `dcmplx`
    looked like the discriminator for `InsertStandardTerms` and the real axis was BARE vs PSEUDISED nuclei,
    with `Ham_PP` sitting on the `<double>` side to prove it.  Same accidental correlation here —
    dcmplx ↔ solid ↔ periodic ↔ PP all coincide in today's test matrix, so any of them "works" as a
    selector until a molecular PP run asks for the PP columns.
  - **Fix direction (not yet designed):** the column set is a TRACE POLICY, chosen from what the run
    actually is (pseudised? variational? gapless?), not from the iterator's type.  Likely a small
    `TraceColumns` value the facade supplies — which is exactly the kind of above-SCFIterator decision
    `SolidCalculation` (Step 4) exists to own, so sequence it after Step 4.

- **V1.28 ⚠️ IMPOSING SYMMETRY ON AN AFM STRUCTURE WOULD DESTROY THE AFM ORDER — the density star-average is
  spin-blind, structurally.  Flagged 2026-08-09, unprompted, because it is directly in the MnO path:** the
  stated plan is "get AFM working with no imposeSymmetry, then start imposing symmetry (Shubnikov groups) and
  cut the RAM substantially".  Step two walks into this.
  - **What is verified (read-only, no runs):**
    1. `Symmetry::SymOp` ALREADY carries `SpinAction sigma` — documented as "Shubnikov spin action (σ), §4
       tier 4a", with a good non-collinear caveat.  So the vocabulary exists.
    2. But `SpinAction::Flip` is constructed in **exactly one place in the whole tree, a unit test**
       (`src/Symmetry/tests/L_Fold.C:213`).  Nothing in `src/` ever produces a spin-flipping op.
    3. And `ReciprocalOp` — which is what the density fold actually consumes (`tComposite_CD::itsPointOps` is
       a `std::vector<ReciprocalOp>`) — is `{U, tau}` with **no σ field at all**
       (`src/Symmetry/Lattice_3D/SpaceGroup.C:48`).  So even if detection produced Flip ops, the reciprocal
       path could not carry them.  `GPW_Evaluator::RecipSymOps()` drops σ on the floor by construction:
       `rops.push_back({Transpose(op.W), op.tau})`.
    4. The ops come from `lat.GetSpaceGroup()` — detected from ATOM POSITIONS, which have no spin.
    5. Each spin channel is its own `tComposite_CD<dcmplx>` and star-averages under the SAME op set; the
       polarized total is a plain ↑+↓ sum (`Imp/ChargeDensity.C:83`).
  - **⇒ The consequence (physics, inferred from the above):** the CHEMICAL space group of rocksalt MnO
    contains operations mapping the Mn↑ sublattice onto the Mn↓ sublattice.  Star-averaging ρ↑ under those
    forces it to be invariant under an operation that exchanges the sublattices — i.e. it **averages the two
    magnetic sublattices together and collapses the AFM order toward the nonmagnetic solution.**  The run
    would not crash; it would quietly converge to the wrong state, which is the expensive failure mode.
  - **This is the SAME BUG CLASS main just fixed one component over.**  `041ddff3` ("Solve the MnO AFM
    collapse: the rho-tilde density mixers are SPIN-BLIND") fixed spin-blindness in the MIXER.  The SYMMETRY
    FOLD is spin-blind too — and worse, blind by TYPE rather than by oversight, since `ReciprocalOp` has
    nowhere to put σ.  The mixer fix is the precedent for what this needs.
  - **Fix direction:** give `ReciprocalOp` a σ (it is the one type on the path that lacks what `SymOp`
    already has), have the fold apply it by swapping channels when σ=Flip, and derive the op set from the
    MAGNETIC (Shubnikov) group rather than the chemical one — an op that exchanges sublattices belongs in the
    group only when paired with a spin flip.  Until then, **`imposeSymmetry` and a polarized AFM density are
    mutually exclusive and should say so**: the cheap interim is a hard throw when
    `imposeSymmetry && dynamic_cast<const tSpinResolved_CD<T>*>(seed)` with a non-trivial op set, in the
    V1.10b "fail loudly in the R&D phase" spirit — far better than a silently demagnetised run.
  - **TRIGGER: the moment MnO AFM converges free and `imposeSymmetry` is turned on.**  Do the interim throw
    before that switch is flipped, not after.
  - **🔔 TRIGGER FIRED 2026-08-09** (MnO dev): AFM-II converges free at E=−60.92 with a properly staggered
    moment (+0.506/−0.529, m_net −0.02).  MnO is safe *today* only because `RunMnO` sets
    `imposeSymmetry=false` by hand — i.e. safe by call-site discipline, which is exactly the fragility V1.30
    flags about the `true` default.
  - **⚠️ THE PROPOSED GUARD IS TOO STRICT — measured 2026-08-09, before implementing it.**  "Throw when
    `imposeSymmetry && spin-resolved`" would break FOUR PASSING TESTS, all legitimate:
    `Polarized{,Seed}SingletMatchesUnpolarizedSiGamma` (multiplicity 1), `O2TripletInBoxMatchesFinite` (3),
    `NaPseudoAtomInBoxDoublet` (2) — every one of them polarized AND `imposeSymmetry=true` via the default.
    The ζ=0 singlets have ρ↑=ρ↓ so every chemical-group op is a symmetry of BOTH channels; O₂ and Na have
    NON-STAGGERED moments, so the ops map ↑ sites to ↑ sites.  Imposing is correct in all four.
  - **⇒ THE CORRECT DISCRIMINATOR IS STAGGERED-vs-NOT, and it is measurable: the ops are a symmetry of the
    CHARGE but not of the SPIN.**  Equivalently \f$m=\rho_\uparrow-\rho_\downarrow\f$ must be invariant
    under the op set — ζ=0 gives \f$m\equiv0\f$ (trivially invariant), FM gives an \f$m\f$ that follows the
    atoms (invariant), AFM gives a staggered \f$m\f$ (NOT invariant).  Only the last is caught.  That
    condition is also the definition of "this needs a Shubnikov rather than a chemical group", so the interim
    guard and the eventual fix share ONE criterion — which is the argument for measuring rather than
    proxying.
  - **✅ FORK RESOLVED 2026-08-09 (MnO dev's ruling) — none of (a)/(b)/(c), and the reasoning changes the
    item.**  Two of the three points land outright and one collides with a measurement:
    - **(2) ACCEPTED — the raw map needs no accessor.**  It is the same density built with an EMPTY op set,
      and `tComposite_CD` takes its ops at CONSTRUCTION.  That is exactly why `ReportSymmetryFound` can
      already measure a per-op defect on a free run.  The A/B is constructible from the public face; I had
      been treating the construction as fixed and only the query as variable.
    - **(3) ACCEPTED, and decisive — it is this project's own prior decision.**  §3
      (SymmetryUpgradePlan.md:444) pins "impose-on-assert as the default, the release-audit bundled with any
      imposition, and the diagnostic on FREE runs — imposition is never silent by construction."  Widening
      `FourierDensity` to measure a PRE-symmetrization defect INSIDE an imposed run is the shape §3 rejected.
      **(b) is withdrawn.**  If raw access is ever wanted it should be a QUESTION — `SymmetryDefect(ops)` —
      not a raw-map getter, per CLAUDE.md's `IrrepCD`-has-no-`GetDensityMatrix()` exemplar.
    - **(1) CONTRADICTED BY MEASUREMENT.**  `imposeSymmetry ∧ spin-resolved ∧ |ops|>1` is satisfied TODAY by
      the four correct tests listed above.  And it cannot be tightened with more configuration terms, because
      **"staggered" is not expressible in the configuration**: seeds are keyed PER-ELEMENT
      (`SeedCD::itsScaleByZ`), so MnO's two Mn sites are indistinguishable to the library, and no
      `SpinPattern`/per-site-moment concept exists anywhere in `src/` — MnO's staggered seed is built in
      qchem6's `RunMnO`.  A type+config predicate can therefore only OVER-fire.
  - **⇒ THE DURABLE ANSWER FALLS OUT OF (3) RATHER THAN FIGHTING IT.**  §3 already says the defect diagnostic
    belongs on the FREE run — and `ReportSymmetryFound` already MEASURES it per-op there.  So the check is not
    "throw when imposing on a polarized density"; it is **"do not impose an op the free run's own defect
    measurement reports as broken."**  The measurement exists; it is simply never fed forward.  That is the
    same criterion as the eventual Shubnikov work (an op belongs in the group only if the density respects
    it), and it is precisely step 5 of the user's SSB workflow (V1.29).
  - **⇒ AND THE IMMEDIATE MnO PROTECTION IS V1.30, NOT A NEW GUARD.**  Flipping
    `GpwOptions::imposeSymmetry` to `false` makes imposition OPT-IN, so nobody can flip the switch by
    accident and `RunMnO`'s hand-opt-out stops being load-bearing.  That removes the hazard this interim
    throw was invented to cover, with no false positives and no new interface.  **V1.30 supersedes the
    interim throw as the urgent item.**  (Owned by the MnO dev — `GpwOptions` is in their file set, and the
    flip needs the gates that relied on `true` made explicit plus a full sweep, so it is not a one-liner.)
  - *(superseded fork, kept for the record:)*  **IMPLEMENTATION FORK.**  The measurement must see the RAW per-channel map: each
    `tComposite_CD` star-averages internally, so by the time `tPolarized_CD` sees its channels the change has
    already happened and the defect measures ≈0.  Options:
    (a) add a raw accessor + `PointOps()` to `tComposite_CD` and do the check in `tPolarized_CD`, needing an
        abstract→concrete cast to reach it (the pattern this codebase flags);
    (b) widen the `FourierDensity` capability face with a raw accessor — clean casts, but widens an abstract
        interface for a guard;
    (c) measure at the `SymmetrizeGMap` call site inside `tComposite_CD`, which has raw + ops in hand for
        free — but a single channel cannot see the cross-channel signature (charge-invariant, spin-not), and
        a large defect there is BENIGN for an unpolarized run, where projecting a broken seed is the whole
        point of imposition.
    Recommend **(b)**: the raw density is a legitimate question about a `FourierDensity` (V1.26's selector and
    the §3 defect diagnostic both want it too), and it is the only option needing no concrete cast.

- **V1.29 The spontaneous-symmetry-breaking DISCOVERY workflow — it is an established method, and step 4
  needs the electronic Hessian rather than noise.  USER 2026-08-09 sketched: (1) imposed run, (2) converge,
  (3) save ρ, (4) reseed a FREE run to see where the orbitals want to move, (5) infer a subgroup, (6)
  re-impose on it.  Asked whether this is known in the literature.**
  - **It is, under two names in two communities.**  Quantum chemistry: **SCF stability analysis** (Seeger &
    Pople 1977), whose classic case is the RHF→UHF *triplet instability*; when the broken solution is the
    deliverable it is **broken-symmetry DFT** (Noodleman).  Solid state: the electronic analogue of
    **soft-mode following**, with step 5 being **isotropy-subgroup analysis** from Landau theory — tooled by
    ISOTROPY/ISODISTORT (Stokes & Hatch) and the Bilbao Crystallographic Server (AMPLIMODES; MAXMAGN /
    k-SUBGROUPSMAG for the magnetic case).  The magnetic version is **magnetic representation analysis**,
    which returns the Shubnikov subgroup directly.  *(Citations from recall — verify before quoting.)*
  - **⚠️ STEP 4 AS SKETCHED CANNOT WORK, for exactly the reason the user suspected.**  A symmetric converged
    solution is a stationary point of the energy in the FULL space, not merely within the symmetric manifold.
    Reseed a free run with it and the gradient is zero by symmetry — the SCF sits still.  Noise does move it,
    but as a random walk: slow, irreproducible, and silent about WHICH mode is unstable.
  - **⇒ The replacement is the ELECTRONIC HESSIAN (orbital-rotation / stability matrix), and it subsumes
    steps 4 AND 5.**  Its negative eigenvalues ARE the symmetry-breaking instabilities; each eigenvector
    transforms as an irrep of the parent group, which is precisely the order-parameter label an
    isotropy-subgroup lookup consumes.  So: how many instabilities, which ones, and the subgroup — mechanically,
    with no noise and no saddle-point ambiguity.
  - **Affordable in practice:** stability codes never form the matrix — Davidson for the lowest few
    eigenvalues.  And the SPIN-FLIP block alone (the triplet instability) is much cheaper and is exactly the
    AFM question, so the magnetic case is the cheap case.
  - **USER QUESTION 2026-08-09: "you have to turn off imposeSymmetry to get the Hessian, so it will be RAM
    expensive?"  Half right, and the better answer inverts the cost.**
    - Correct that the symmetric manifold cannot show you the instability: that is *why* the symmetric
      solution looks stationary — `imposeSymmetry` keeps exactly the TOTALLY-SYMMETRIC block.
    - But the Hessian **BLOCK-DIAGONALISES BY IRREP of the parent group** (orbital rotations decompose into
      irreps; the Hessian has no elements between different irreps).  So you never need the full unrestricted
      Hessian — sweep the NON-symmetric blocks ONE AT A TIME, each a fraction of the rotation space.  The
      symmetry machinery becomes the tool rather than the obstacle.
    - **And the block you find the negative eigenvalue in IS the order-parameter irrep** — step 5's label
      falls out of the bookkeeping instead of being inferred from a drifting density.
    - **RAM is not the binding cost.**  Davidson needs Hessian-VECTOR products (≈ one Fock build each) and a
      handful of trial vectors: footprint ≈ (a few) × orbital space, NOT \f$N_{rot}^2\f$.  The cost is CPU.
      ⇒ this step does NOT have to wait on the Shubnikov RAM savings.
    - For AFM the relevant block is the SPIN-FLIP one, which factorises separately from the spatial irreps in
      the non-relativistic collinear tier — the cheap corner, not the expensive one.
  - **Reusable:** the stability matrix IS the RPA/TDDFT (A,B) matrix.  If excited states ever land, this is
    the same machinery — worth knowing before designing it as a one-off.
  - **Do not over-build it for MnO:** AFM-II is known experimentally, so the known magnetic structure can
    simply be imposed.  The discovery workflow earns its keep on materials whose order is UNKNOWN.
  - Depends on V1.28 (σ on `ReciprocalOp` + a Shubnikov op set) for step 6 to be expressible at all.

- **V1.30 ✅ APPEARS ALREADY FIXED — item is STALE (verified against the tree 2026-08-10, not re-derived
  from the doc).**  `GpwOptions::imposeSymmetry` is now `= false` with the comment "DEFAULT off (full mesh,
  free run)" (Lattice_3D/BasisSet.C:72), and `SolidCalcOptions::imposeSymmetry` is `= false`
  (Calculation/SolidCalculation.C:103) — so the two facades AGREE and both default off, which is exactly
  what this item asked for.  Landed on the MnO side (it was theirs by assignment) and the worklist entry
  was never closed.  **Left marked rather than deleted so the next MnO session does not spend time
  re-doing it; confirm and close at the merge.**
  *(original text follows)*
  **`GpwOptions::imposeSymmetry` defaults to `true` while its own comment says "OPT-IN".  Found
  2026-08-09 by the Step 4 facade gate, which resolved a different Becke cost than the driver did.**
  - The field is `bool imposeSymmetry = true;` and the comment three lines down reads *"OPT-IN per the §3 pin
    (an imposed default would also ~2x the suite's XC grids)"*.  One of them is wrong, and the comment is the
    one carrying a reason.
  - **How it surfaced:** `SolidCalculationMatchesTheSiAnchor` printed `Becke 36000 pts [..., free]` where the
    driver's Si run prints `Becke 72000 pts [..., imposed]` — the 2x the comment predicts.  The energies still
    agree (the selector routes Si to uniform either way), so the gate passed; the discrepancy showed up only
    in the announce line.  Worth noting that the gate caught a defaults divergence it was not written to look
    for, which is the argument for having the facade announce its decisions at all.
  - **⚠️ Now actively hazardous, not cosmetic.**  Per V1.28 an imposed run star-averages each spin channel
    under the CHEMICAL space group, whose sublattice-exchanging ops would average an AFM structure's magnetic
    sublattices together.  Default-ON means the MnO work inherits that silently unless every call site
    remembers to turn it off.  **Default-OFF is the safer way to be wrong**: an imposition you did not ask for
    is invisible in the result, a missing one only costs time.
  - `SolidCalcOptions` deliberately defaults it to `false` (the comment's stated intent) and says so at the
    field, so the two facades diverge ON PURPOSE rather than by drift.  Decide which is right and align them.

- **V1.31 ✅ DONE `627a4ff9`. `SymFockCache` deleted; the SALC path builds ONE whole-AO Fock and slices it.**
  The memo was caching a partial AO Fock at the basis level inside a loop that should not have been
  iterating; the fix removed the loop, not the staleness test.  **→ doc/CleanupHistory.md**
  *(the full analysis, including the retracted first draft and the refuted ruling, follows -- it is the
  part worth reading)*
  **This item was FILED WRONG on 2026-08-10 and corrected the same day by the user — the correction is
  the more useful half, so it is kept.**
  - **What the first draft claimed:** "caching has escaped `DB_Cache_RAM` — three caches, three
    invalidation disciplines", with `SymFockCache` and `itsJKs` described as caching "AO J/K", and
    `DB_Cache_RAM`'s entries described as "immortal by design".
  - **Both halves of that were wrong, and the user caught it from the description alone:**
    - **`SymFockCache` and `itsJKs` do NOT hold J/K TABLES.**  They hold the CONTRACTED matrix
      \f$J[D]\f$ — `BuildAOFock` runs `raw->AccumulateDirect(M,Dao,raw)`, so `M` is an
      nAO×nAO Fock CONTRIBUTION, not an `ERI4`.  `SymFockCache::Entry` also stores a copy of `D`
      itself, but only as a staleness token.  Both are therefore charge-density-DEPENDENT, and by the
      user's own rule — *"anything charge density dependent is obviously not [a cache candidate]"* —
      they correctly do NOT belong in `DB_Cache_RAM`.  **There is no missing delegation here.**
    - **`DB_Cache_RAM` entries are NOT immortal.**  `Jac`/`Kab` (the ERI4s) are timestamped on every
      access and LRU-evicted by `RunGarbageCollector` once `itsTotalRAM > itsMaxRAM`, protecting the
      entry just inserted (Imp/DB_Cache_RAM.C:354-365).  So the framework already is what the user
      described: the home for the huge J/K tables, WITH garbage collection.  (What never happens is
      clearing between RUNS — a different property, and the one `EmitReport`'s comment is about.
      Non-ERI4 entries are small and are not GC'd.)
  - **So the actual, much smaller item:**
    1. `SymFockCache` (SymmetryAdapted_IBS.C:32-42, Imp:45-66) memoizes \f$J_{AO}[D_{cd}]\f$ per
       cd-irrep, invalidating by an ELEMENTWISE compare against a stored copy of D.
    2. `Dynamic_HF_HT_Imp::itsJKs` (Hamiltonian/Internal/Terms.C:163-167) memoizes the per-irrep
       \f$J[D]\f$ / \f$K[D]\f$ blocks, invalidating by a VERSION COUNTER (`itsCD_Version`).
    These are the same class of quantity — "the Fock contribution for the current density" — memoized
    twice, one library apart, and they do not agree on how staleness is decided.  **Look first at
    whether 1 can simply be deleted and the SALC path routed through 2.**  If it cannot, they should
    at least share the version-counter discipline: the elementwise compare costs an O(n²) scan plus a
    full stored copy of D, per irrep, per lookup.
  - **⛔ THE RULING DOES NOT APPLY HERE — REFUTED 2026-08-10 BY READING THE CALL PATH.  Read this BEFORE
    the ruling below, which is preserved because the reasoning is right in general and wrong for this
    site.**  Two facts, both checked against the tree:
    - **(1) `SymFockCache` is an INTRA-SWEEP memo, not a cross-iteration cache.**  Its only drivers are
      `Vee::AccumulateAll` and `Vxc::AccumulateAll`, each called ONLY from
      `Dynamic_HF_HT_Imp::ContractAll`, which is version-guarded and therefore runs the whole scatter
      sweep EXACTLY ONCE per density serial.  So across iterations the cache can never hit — a new
      density is a new sweep.  What it actually buys is WITHIN one sweep: the composite visits ~N²/2
      canonical irrep pairs and each SALC `AccumulateDirect` would rebuild the AO Fock, so the memo turns
      ~N² AO builds into N.  **Its cross-iteration invalidation logic is answering a question that cannot
      arise.**
    - **(2) A version counter would BREAK the polarized SALC path — and there is a live green test that
      would catch it.**  `tPolarized_CD::AccumulateDirectAll` runs the Up composite then the Down
      composite INSIDE ONE SWEEP.  Both reference the same SALC basis, so both hit the same cache under
      the same key (`cd->BasisSetID()` is content-based: raw id + "[label]").  And
      `tPolarized_CD::Version()` FORWARDS TO ITS UP CHILD (ChargeDensity.C:278) — so Up and Down are
      indistinguishable by version.  A version guard would serve the Down pass the Up channel's
      \f$J_{AO}\f$: a silently wrong polarized Fock.  **The elementwise D compare is the ONLY thing
      separating them**, which makes it load-bearing rather than lazy.  `M_Sym.water_HF_polarized` (and
      `water_DFT_polarized`) exercise exactly this and are green today.
    - **⇒ And "move the memo to the term" is not executable either**, for a reason independent of the
      above: the memo sits at a granularity the term cannot see.  The term memoizes FINAL per-irrep blocks
      (`itsJKs`); these are INTERMEDIATE AO builds produced deep inside `dm->AccumulateDirectAll(X)`, on
      the far side of the density's virtual API.  There is no hook at which a term could supply or scope
      them.
    - **What the finding actually points at — SCOPE it, do not version it.**  If the memo's life is one
      sweep, staleness is not a question to answer but one to ELIMINATE: nothing needs comparing or
      versioning if the thing is discarded at the sweep boundary.  Both defects below dissolve with it
      (the missing `Ocd` in the key cannot collide inside a single sweep either).
    - **And a better fix may sit behind that: ONE AO build per sweep instead of N.**  \f$J\f$ is LINEAR
      in \f$D\f$, so \f$\sum_C J_{AO}(O_C D_C O_C^\mathsf{T}) = J_{AO}(\sum_C O_C D_C O_C^\mathsf{T})\f$
      — sum the AO densities first and build once.  That removes the memo entirely rather than rehoming
      it.  It needs a place to accumulate across the sweep, i.e. the same missing sweep boundary, so the
      two questions are one question.  **NEEDS A USER RULING; do not guess it.**

  - **~~✅ USER RULING 2026-08-10: the VERSION COUNTER is the right discipline — *"I think holding a
    CD_Version is better (cleaner) than comparing D element by element."*  This settles the item, because
    of where the version LIVES.**
    - `Version()` is on `rChargeDensity`, drawn from `ChargeDensity::NextDensityVersion()` — a documented
      program-wide monotonic clock whose comment already records a bug from getting this wrong (a serial
      colliding across density KINDS made a dynamic term reuse the iter-0 seed Fock for the iter-1
      density, silently breaking the SCF).  So this is not merely tidier; it is the discipline the
      codebase has already paid for once.
    - **But `SymFockCache` cannot adopt it where it lives.**  `qcChargeDensity` links `qcBasisSet`
      (ChargeDensity/CMakeLists.txt:38), so qcBasisSet is BELOW it in the DAG and cannot name a charge
      density; and the contraction face hands the decorator only the raw matrix `Dcd`, never the density
      object.  The elementwise compare is therefore not laziness — it is the only staleness test
      available at that altitude.
    - **⇒ The ruling forces the fix to be option (c): DELETE `SymFockCache` and let the term own the
      memo,** in `qcHamiltonian`, where `cd->Version()` is visible and `itsJKs` already does exactly
      this.  The alternative — threading a version token through `Accumulate*` — would push a caching
      concern back onto the contraction face R1.7 just finished cleaning, and would invert the
      dependency besides.  Do NOT take it.
    - Cross-check before deleting: `SymFockCache` is optional today (the ctor's `shared_ptr` defaults to
      null and the decorator builds directly when it is absent — "fine for tests"), so the delete has a
      working no-cache path to fall back on while the term-level route is wired up.~~
    *(struck through: the reasoning above is sound as a general preference and the DAG argument is correct
    — it is the premise "this is a cross-iteration cache whose staleness test should be a version" that
    the call path refutes.  Kept in full per this doc's own rule that refuted prescriptions are worth more
    than landed ones.)*
  - **One real defect, independent of the above: `SymFockCache`'s key is INCOMPLETE.**  The entry
    depends on `raw` and `Ocd`; the key is `cd->BasisSetID()` alone.
    `SymmetryAdapted_IBS::BasisSetID()` is `raw id + "[label]"`, so `raw` is encoded by luck; `Ocd` is
    not encoded at all.  Two SALC transforms of the same raw basis (different tolerance, different
    group) sharing one cache — which the ctor's `shared_ptr` parameter permits — collide silently.
    That is exactly the failure `DB_Cache_RAM`'s dim guards (Imp/DB_Cache_RAM.C:211-232) were added to
    catch, and a privately-rolled cache gets no such guard.  **This is the strongest argument for the
    user's delegate-all-caching principle in this file** — not RAM, not GC, but the fact that a
    hand-rolled cache silently opts out of the key-completeness checks.
  - **Worth recording as the GOOD case: the SALC decorator is otherwise a model citizen.**  Everything
    density-INDEPENDENT it exposes already goes through `theCache<T>()` — the 1-electron blocks via the
    inherited cached accessors, and the transformed 3-centre tensors via `Overlap3C`/`Repulsion3C`
    (Imp/Orbital_DFT_IBS.C:10-20), whose cache-miss hooks are its own `MakeXxx3C`.  `SymFockCache` is
    its ONE hand-rolled cache, and it sits precisely at the density-dependent boundary.
  - **And on "SALC-transformed J/K are good cache candidates" (user, 2026-08-10): there are none to
    cache.**  R1.7 established that the SALC path has no per-irrep-pair `ERI4` at all — it builds the
    whole-AO Fock and slices, so the only ERI4s in existence are the RAW ones, which are already in the
    framework (with GC).  If per-irrep transformed 4-index blocks are ever wanted, that is a decision to
    REVISIT R1.7, not a cache-delegation task; note it would trade the AO build for an
    \f$O(N^4)\f$ four-index transform per irrep pair, which is why the AO-slice route was chosen.


### V2 — measurements / sweeps

- **V2.1 Spin-native collapse: `Delta_XC` ×2 vs `Delta_XC_Pol`+`Delta_VcorrPol`.**  Per the
  spin-native-is-primary bias the unpol pair should become the ζ=0 collapse: the polarized pair at
  ρ↑=ρ↓=ρ/2 IS the unpolarized answer (singlet-collapse gate proves it to 0.12 mHa).  Collapsing
  halves the class count at ~2× XC pointwise work for closed shells (the Ham_PP trade).  Decide
  with a perf measure.  Falls out for free if it lands: XC_GridEngine's second (scalar) rho cache
  dies (pair only; total = up+dn).  Supporting duplication catalog (verified near-verbatim):
  Delta_XC vs Delta_XC_Pol CalcMatrix/GetEnergy differ only in Rho vs RhoPol; VxcPol vs
  FittedVxcPol are both "own two channel terms, cast to Polarized_CD, dispatch on spin" — a
  generic pair-dispatcher should span the HF pair too (the user's original "generic pol Vxc" ask).
  - **USER 2026-08-08, restating the target:** *"All Vxc's should be spin native (working with Spin::None
    when required).  Ideally there should be only one PolarizedVxc that simply stores two abstract Vxc
    pointers and does the obvious function forwarding.  (Same for Vex, I forget how Vcorr works for
    polarized CDs.)"*
  - **Answering the Vcorr question, because it decides the mechanism: EXCHANGE is channel-separable,
    CORRELATION IS NOT.**  The tree already states it at `VWN_Correlation.C:47` — *"v_c^sigma COUPLES both
    channels (through r_s and zeta) -- it is NOT a function of rho_sigma alone."*  \f$\epsilon_c(r_s,\zeta)\f$
    interpolates paramagnetic→ferromagnetic through \f$f(\zeta)\f$ and the spin stiffness (VWN Eq. 4.4), so
    \f$v_c^\sigma\f$ depends on BOTH densities.  Exchange escapes this only because of the spin-scaling
    relation, which is why `SlaterExchange` can be channel-native at all.
    ⇒ **"Two Vxc pointers + obvious forwarding" works for Vex and CANNOT work for Vcorr.**
  - **The spin-native half of the ask is ALREADY DONE, one layer below where the item was looking.**
    `SlaterExchange` carries `itsSpin` and halves ρ only for `Spin::None` (the closed-shell collapse);
    `VWN_Correlation`'s PRIMARY face is the two-channel `GetVc(rup,rdn,s)`, with the scalar face documented
    as the ζ=0 collapse and byte-identical to the historical code.  The FUNCTIONALS are spin-native with
    Spin::None as the special case — exactly the stated bias.  The duplication is at the TERM layer.
  - **⇒ THE SHARPER TARGET: 3 term classes → 1, not 3 → 2.**  `DeltaFittedVxc`, `DeltaFittedVxcPol` and
    `DeltaFittedVcorrPol` are structurally identical — same two bases, same `{functional, engine}` members,
    same four methods — differing only in which functional face they hold and how `CalcMatrix` feeds it.
    Give every XC functional ONE spin-native face, \c double GetV(double rup, double rdn, const Spin&):
    exchange implements it ignoring the other channel, correlation uses both, and unpolarized is
    \c GetV(ρ/2,ρ/2,Spin::None).  Then ONE term class holds ONE functional and covers all three cases, and
    correlation is expressible — which the two-pointer sketch cannot manage.
    This is the session's own heuristic (line ~83): the candidate special case (unpolarized) is the general
    case with a term set to zero (ζ=0), so it is not a special case at all.
  - **⚠️ CROSS-CHECK AGAINST THE NON-COLLINEAR GOAL (user 2026-08-09: "ultimately I want non-collinear AFM
    order to be the default assumption in code, not the exception").**  The proposed face
    \c GetV(rup,rdn,Spin) is HALF future-proof and half not:
    - the two DENSITY arguments are fine — \f$(\rho_\uparrow,\rho_\downarrow)\f$ and \f$(\rho,m)\f$ are a
      linear change of variables, and non-collinear LSDA evaluates the collinear functional along the LOCAL
      quantization axis after diagonalising the 2×2 spin density, so it still needs exactly two magnitudes.
    - the \c Spin ARGUMENT is the collinear assumption in disguise.  Non-collinear needs a DIRECTION
      (\f$\hat m(r)\f$), not an Up/Down label; the answer is a 2×2 potential matrix, not a scalar tagged
      with a channel.
    - `Symmetry::SymOp`'s own doc already issues this warning for the symmetry side — *"consumers must not
      bake 'two scalar channels' into interfaces beyond the collinear tier"*.  It applies verbatim here.
    ⇒ Land V2.1's collapse (it is a strict improvement today), but spell the face so the SPIN argument is the
    replaceable part.  Do NOT let \c Spin leak into the XC functional's contract as a permanent parameter.
  - **USER RULING 2026-08-09 on the REPRESENTATION: SU(2)/matrix, not SO(3)/direction.**  *"'direction m̂(r)'
    is the SO(3) way of thinking about it.  I suspect the SU(2) language is more natural in software.  We just
    need to complexify the up/down coefficients."*  Agreed, with one correction to the mechanics:
    - the DIAGONALS stay REAL (they are densities); what you gain is ONE COMPLEX OFF-DIAGONAL.
      \f$\rho=\tfrac12(n\,\mathbb{1}+\mathbf{m}\cdot\boldsymbol\sigma)\f$ with \f$\rho_{\uparrow\uparrow},
      \rho_{\downarrow\downarrow}\in\mathbb{R}\f$ and \f$\rho_{\uparrow\downarrow}\in\mathbb{C}\f$,
      \f$\rho_{\downarrow\uparrow}=\rho_{\uparrow\downarrow}^*\f$ — 2 reals + 1 complex = **4 real DOF**,
      matching \f$(n,\mathbf{m})\f$ exactly.  \f$m_z=\rho_{\uparrow\uparrow}-\rho_{\downarrow\downarrow}\f$,
      \f$m_x=2\mathrm{Re}\,\rho_{\uparrow\downarrow}\f$, \f$m_y=\mp2\mathrm{Im}\,\rho_{\uparrow\downarrow}\f$.
      So the face generalises \c (rup,rdn) → \c (rup,rdn,rud) with \c rud complex, NOT by complexifying the
      diagonals.
    - **Why the matrix form is the better CONTRACT** (four independent reasons, worth keeping): v_xc becomes a
      2×2 Hermitian matrix entering the Fock matrix directly in spin space, with no rotation to track;
      symmetry acts as \f$\rho\to U\rho U^\dagger\f$, \f$U\in SU(2)\f$, which COMPOSES (an enum does not);
      the collinear tier is exactly the diagonal case \f$\rho_{\uparrow\downarrow}=0\f$, i.e. the same
      collapse pattern the project already prefers; and there is no direction to parameterise, hence no branch
      or gimbal handling.  LSDA still diagonalises the 2×2 locally
      (\f$n_\pm=\tfrac12(n\pm|\mathbf{m}|)\f$) and applies the collinear functional along the local axis —
      but that is INSIDE the functional, not in the contract.
    - **Lands directly on `SpinAction`:** \c {None,Flip} is the two-element subgroup of SU(2) (Flip being one
      specific element).  The general Shubnikov op is \f$\{W|\tau\}\times SU(2)\times\f$ an optional
      ANTIUNITARY factor (time reversal) — that antiunitary piece is what makes it a MAGNETIC group rather
      than a spin-space group, and it is the part an \c SU(2)-matrix-only generalisation would miss.  Relevant
      to V1.28's σ-on-`ReciprocalOp` work: size the field for {rotation + antiunitary flag}, not just a flip.
  - Still to decide with a measurement, unchanged: the ~2x XC pointwise cost for closed shells (the Ham_PP
    trade).  And the free consequence still stands — `XC_GridEngine`'s second (scalar) rho cache dies.
- **V2.2 GPW default seed policy.**  `GpwOptions.seed` defaults to `Uniform`, which the Na-doublet
  campaign showed has a STABLE wrong basin for electron-sparse systems (lone-electron doublet
  converged 72 mHa high with every health metric green).  The molecular facade already defaults
  DFT to SAD.  Candidate: default GPW to `IonicSAD` (SAD-family), Uniform = explicit opt-in —
  needs a suite sweep since every pinned GPW anchor re-seeds.
- **V2.3 Polarized PLANE-WAVE Vxc fit route** — `Ham_PW_DFT` polarized currently THROWS for
  `VxcFit::PlaneWave`: per-channel PW_XC needs per-spin rho-grid caches (PW_XC's `itsRhoGrid` is
  keyed on `cd->Version()` alone, which a polarized density aliases across channels — the trap the
  engine's RhoPol pair-cache fixes).  Design note in the throw message.

- **V2.4 ✅ DONE 2026-08-08 — margin validated, selector ARMED.**  Converged-run A/B on both systems the.  **→ doc/CleanupHistory.md**
- **V2.5 `PPMeshParams()` sizes its uniform mesh with no \f$\alpha_{pp}\f$ term.**  `mp.eCut=densityEcut`
  \f$=C\alpha_{\max}\f$, but its integrand is \f$\langle\chi_i|V_{short}|\chi_j\rangle\f$ with exponent
  \f$2\alpha_{\max}+\alpha_{pp}\f$.  Independent of V1.26's selector (which already accounts for
  \f$\alpha_{pp}\f$ in its CHOICE); this is the mesh sizing itself.  One consumer today — the KB-projector
  grid fallback (GPW Evaluator.C:1119) — so the exposure is bounded, but the floor is simply missing.
  Raising it moves grids, hence anchors: measure first (D8), same instrument as V2.4.

- **V2.6 Are the Becke recipe's `nRadial`/`angularDegree` defaults over-generous?  USER 2026-08-07:**
  *"There is a lot riding on the defaults for nRadial and nDirs(degree) for the Becke grid.  The degree can be
  determined from point symmetry of the atom site, but I have the impression that you can often 'get away
  with' much lower degrees than the point symmetry dictates."*
  - **Why it is now load-bearing beyond accuracy:** since V1.26 those two numbers set the ENTIRE Becke side of
    the Uniform-vs-Becke cost comparison (\f$n_{atoms}n_{radial}n_{dirs}\f$ — nothing else enters).  An
    over-generous recipe therefore biases the selector TOWARD uniform.  That is the SAME direction as
    refinement (a)'s bias, so the two errors **compound rather than cancel** — which is a second, independent
    reason the Si verdict should not be trusted, and a reason to do this measurement BEFORE V2.4's.
  - The degree-vs-site-symmetry point is the §6a/W2b question from the other end: the site-adapted builder
    derives an angular rule FROM the site group, so "what the point group dictates" is already computable —
    the open question is the gap between that and what the integrand actually needs.  The instrument exists
    (`Mesh_AngularDegree` measures a rule's degree monomial-by-monomial; the GPW Becke gate measures dExc/dVxc
    against a fine reference), so this is a sweep, not a design problem.
  - **✅ MEASURED 2026-08-07 on THREE systems — `GPW_SCF.DISABLED_BeckeRecipeLadder_{SiGamma,NaF,MnSextet}`
    (hand-run).**  Method, D8-compliant by construction: converge ONCE, then quadrature the SAME frozen
    density on a ladder of meshes against a fine reference in the same family (nR=100, GL-41, same
    `mhl_alpha`).  No SCF re-runs and no ΔE_total — scored on E_xc and the V_xc MATRIX.  Freezing ρ is what
    makes it a measurement of the QUADRATURE rather than of the SCF.  Free meshes throughout (the ladder
    measures the RULE, not the symmetry fold).

    **`max|dVxc|` — the binding metric (the Becke gate's tolerance is 1e−3):**

    | ANGULAR (nR=40) | Si covalent | NaF ionic | Mn open-shell d |   | RADIAL (GL-29) | Si | NaF | Mn |
    |---|---|---|---|---|---|---|---|---|
    | GL-5  | 1.4e−2 | 2.4e−2 | 1.3e−3 | | nR=10 | 5.0e−2 | 2.5e−1 | 3.0e−1 |
    | GL-7  | 5.6e−3 | 1.1e−2 | 6.0e−4 | | nR=15 | 7.9e−3 | 4.0e−2 | 7.1e−2 |
    | GL-9  | 1.4e−3 | 4.2e−3 | **1.1e−5** | | nR=20 | 3.5e−3 | 3.8e−3 | 1.2e−2 |
    | GL-11 | **9.5e−4** | 3.0e−3 | 6.8e−6 | | nR=25 | **9.0e−4** | 1.2e−3 | 1.5e−3 |
    | GL-15 | 8.0e−5 | **7.6e−4** | 1.4e−6 | | nR=30 | 1.5e−4 | **3.2e−4** | **2.0e−4** |
    | GL-17 | 5.8e−5 | 2.2e−4 | 1.3e−6 | | nR=40 | 1.7e−5 | 7.9e−5 | 1.3e−6 |
    | GL-29 (prod) | 1.7e−5 | 7.9e−5 | 1.3e−6 | | nR=60 | 3.7e−6 | 2.1e−5 | 2.7e−7 |

  - **VERDICT — the two axes are not alike, and only ONE is over-generous.**
    - **ANGULAR: over-generous by ~2.8x.**  Degree needed: Mn 9, Si 11–15, NaF 15–17.  **Recommend degree 17**
      (162 directions vs 450), leaving the worst case (NaF) at 2.2e−4, a 4.5x margin.  **Degree 15 is NOT
      recommended** even though it suffices on Si — NaF sits at 7.6e−4, inside tolerance by only 1.3x.  That
      gap between the one-system and three-system answers is precisely the over-fit this item existed to
      catch: a default tuned on Si alone would have shipped 15.
    - **RADIAL: 40 is right, and all three agree.**  nR=30 is the first rung inside tolerance everywhere,
      nR=20 is 4–12x out, nR=25 is marginal on two of three.  On Si the E_xc approach is also NON-MONOTONIC
      (nR=20 worse than nR=15).
  - **⚠️ THE BONDING-CHARACTER PREDICTION IS REFUTED, AND SO WAS MY EXPLANATION OF THE RADIAL AXIS.**  Both
    corrections come from the same source: I framed the angular requirement as a property of the SITE's own
    density (which harmonics its point group allows), and it is not.
    - **What actually drives it: the Becke PARTITION between DISSIMILAR neighbours.**  Mn is a SINGLE atom in
      a box — no interatomic partition at all — and is trivially easy (degree 9), despite being nominally the
      "most aspherical" system.  (It is also not aspherical: high-spin d⁵ is the one d configuration that is
      spherically symmetric, so that test does not probe asphericity even in principle — a design error in my
      choice of system, worth stating.)  NaF — two ions of very different size and sharpness — is the
      HARDEST.  Si, whose partition is between IDENTICAL atoms, sits between.  So the ordering is by
      PARTITION-SURFACE difficulty, not by site asphericity, and **ionic is the least forgiving, not the
      most**.  The tree already recorded the mechanism (`src/Structure/tests/MolecularMeshTests.C`: "the fuzzy
      Voronoi switching shell is angular-quadrature limited"); the site-symmetry framing simply ignored it.
    - **Still untested: genuine on-site asphericity.**  None of the three has it (Si and NaF are closed-shell,
      Mn's d⁵ is spherical).  A real probe would be a crystal-field-split d (MnO) or the O₂ π* triplet, which
      is already an enabled test.  Until then "does asphericity matter?" is open — the measurement so far says
      only that it is not what separates these three.
  - **⚠️ CORRECTION to the earlier "no cheap a-priori radial diagnostic is possible" claim (commit f514c407).**
    That claim blamed the Becke PARTITION for the node-count heuristic's failure.  **Mn refutes it**: a single
    atom has no interatomic partition and its radial axis is still binding (nR=10 → 3.0e−1).  The real reason
    is simpler and the fix is available: MHL clusters as \f$r\propto x^m\f$, so "9 nodes inside the feature"
    can be 9 nodes bunched at tiny \a r with a gap exactly at the density peak — the heuristic counts nodes
    where it should measure SPACING.  **The quantity that does track the measurement, across all three
    systems, is \f$r_{peak}/\Delta r(r_{peak})\gtrsim3\f$** for the density peak of
    \f$r^2e^{-2\alpha_{\max}r^2}\f$:

    | | nR=20 | nR=25 | nR=30 | nR=40 | nR=60 | first rung in tolerance |
    |---|---|---|---|---|---|---|
    | Si (α_max=2, peak 0.354) | 2.4 | — | **3.4** | 4.4 | 6.5 | nR=30 |
    | Mn (α_max=36, peak 0.083) | 1.3 | — | 1.9 | **3.0** | 4.0 | nR=40 |

    So a basis-driven `nRadial` warning looked shippable — it just needed local spacing at the peak, not a
    node count in a window.  Filed as **V2.7**.
    **⚠️ BUT THE THRESHOLD IS REFUTED BY Al TOO (checked 2026-08-07, against its MEASURED α_max=4 rather than
    a guess).**  The criterion gives 3.2 at nR=30, which the "≳3" rule passes — while measurement puts Al at
    nR=30 at 9.3e−3, 9x out of tolerance, needing nR=40 (3.8).  So the threshold fitted to Si and Mn does not
    cover the metal, on the RADIAL axis either.  V2.7 survives as a SHAPE (local spacing at the peak is the
    right quantity — it explains why the node count fails) but its threshold must be calibrated on a metal,
    which means ≳3.8 and only 3 supporting points.  **Third time in this item that an insulator-fitted grid
    rule broke on Al** — the angular recommendation, the degenerate-shell assumption, and now this.  The
    transferable rule: calibrate a grid criterion on a simple metal, or do not ship it as a global default.
    **🔔 2026-08-09 — the third system arrived, and it is the right one.**  MnO run 24 ends on a FIT-FLOOR
    STALL at −60.134: the fit/grid floor now bounds the answer, not the SCF.  So MnO is a system where the
    radial criterion is not academic, and it is neither an insulator ladder point nor a simple metal — an
    open-shell AFM oxide with a sharp O and a semicore-ish Mn.  **Calibrate V2.7 on {Si, Mn-atom, Al, MnO}
    rather than the original two.**  It also makes **V2.5** (`PPMeshParams`' missing \f$\alpha_{pp}\f$
    floor) testable for the first time on a system that is genuinely grid-bound.
  - **CONSEQUENCE FOR V1.26/V2.4, running OPPOSITE to this item's original prediction.**  The item argued an
    over-generous recipe biases the selector toward uniform.  That was right — so FIXING it makes Becke
    **more** competitive: GL-29→GL-17 is 450→162 directions, dropping the Becke side 2.8x.  Si: 36,000 →
    12,960 against uniform's 4,913, i.e. from 7x dearer to 2.6x.  **V2.6a must land before V2.4**, or the
    margin gets fitted to a Becke cost about to change by 2.8x.

- **V2.6a ⛔ ATTEMPTED AND REJECTED 2026-08-07 — flip `angularDegree` 29 → 17.**  Made the one-line change,.  **→ doc/CleanupHistory.md**
- **V2.7 A basis-driven `nRadial` adequacy warning, on the \f$r_{peak}/\Delta r\f$ criterion.**  The
  V2.6 correction above makes this shippable where the first attempt was not.  It is the radial sibling of
  V1.26's `UniformDivisions` bridge: MHL's \f$r(x)=\alpha(x/(1-x))^m\f$ gives \f$\Delta r\f$ in closed
  form, and \f$\alpha_{\max}\f$ is already gathered for the selector (`XCMeshSharpness`), so the check
  costs nothing new.  Calibrate the threshold on more than two systems first — 3 is fitted to Si and Mn.

### V3 — repro / campaign bugs (Spin-SAD, 2026-08-04)

- **V3.1 Polarized ATOMIC solver fails on an empty minority channel** — `AtomCalculation(11, 0,
  {.pseudopotential=true, .valence=1, .pol=Polarized})` (Na q1 doublet: nUp=1, nDown=0) dies with
  "Invalid setup of symmetric matrix" before the first Fock.  Same failure class as the tier-4b
  finding (`cSCFAcceleratorDIIS::GetNProj` segfault on an empty EC).  Worked around for the Na
  library entry (1 valence electron ⇒ the pair is EXACT without an SCF, hand-constructed in
  atomic_valence_densities.json), but any future fully-polarized one-electron species hits it
  again — the valgen `--spin` path should either fix the empty-channel atom SCF or synthesize the
  exact pair itself when nDown==0.
- **V3.2 Unoccupied-l shells break the atom SCF** — the same Na valgen run ALSO failed when the
  recipe carried the (unoccupied) p window from valence_lowq_sr (`--shell 1:2:0.09:0.3`): "Invalid
  setup of symmetric matrix" even unpolarized-adjacent.  GenerateValenceBasis documents "higher-l
  polarization shells ride along un-validated (as intended)" — but at least the polarized
  GenerateSeedDensity path chokes on them.  Reproduce + fix or document the restriction.

### V4 — watch triggers (act when the trigger appears)

- **V4.1 `CollocMemo` dual duty** — replay memo (exact-D level densities) + adjoint D-screen in
  one struct with different lifetimes (last-D vs union-Dscr).  Fine today; split when a THIRD
  consumer appears.
- **V4.2 `SolveSPD`/NNLS in `SymmetrizeMesh.C`** — module-private dense Cholesky + Lawson-Hanson
  NNLS beside the mesh code; promote to qcMath when a SECOND consumer of small dense LS/NNLS
  appears (Blaze has no NNLS).

## DECIDED-ELSEWHERE

- **D1 `Band_DFT_IBS`** — verified fully dead in code (zero imports/derivations/casts; only the
  CMake entry).  Governing record: FittingCleanupPlan §D deliberately KEPT the module as the
  intended future-GPW `<double>` interface.  NEW FACT flagged for re-decision: GPW has since
  landed (`GPW_IBS.C`) and does NOT implement it.  (Stale comments → R2.4.)
- **D2 "Why is a FourierCD different than a tChargeDensity?"** — governed by the
  pw-fitting-uniform-interface pin (no "Fourier" in abstract faces; PW fitting looks identical to
  molecular fitting) — same campaign, don't re-derive here.  User context preserved: at a high
  level Fourier/PW is an implementation detail; `G_ERI3` (src/BasisSet/Internal/GMap.C) is the
  multi-purpose data structure harmonizing the return-type differences — it ends up being a spec
  for the needed data-structure flexibility; high-level code should not have to worry about G or
  no G, PW, Fourier.  Verified corrections to the old bullet: the class is `FourierDensity` — a
  standalone cross-cast MIXIN, not a tChargeDensity; SCFIterator.C is NOT a consumer (stale
  imports → R2.4); missed consumers = DensityMixer (3 cast sites), tPolarized_CD ↑+↓ forwarding,
  tComposite_CD k-sum forwarding, + 3 direct implementors (SeedCD, PolarizedSeedCD, FourierMixCD).
  Taxonomy note: the seeds (NumericCD/SeedCD/PolarizedSeedCD) are tChargeDensity but neither
  DM_CD nor FittedCD — a third leg beside the doc's two primary variants (DM_CD = ERI-capable
  first stop after any SCF iteration; FittedCD = transient DFT fit: Create/DoFit/GetRepulsion),
  orthogonal aspects: Polarized/Unpolarized, Composite (multiple irreps — is this really
  different than Polarized, since irreps can hold spin?), version tracking, mixing algorithms,
  spin density (plotting).
- **D3 `FIT_SF_ABS::SymmetrizeRaster`** — (SRP) symmetry op on a fit-basis interface.  Removable
  once a τ-acting DIRECT-raster fold variant of `FoldGrid` exists (T3 groundwork): precompute the
  raster fold, ctor-inject into `tComposite_CD` beside `itsPointOps`, delete the virtual (uniform
  route's voxel-shift moves out of `PlaneWaveFit_IBS`).
- **D4 The basis-as-policy-carrier virtual family** — `GetReciprocalPointOps` /
  `GetDetectedReciprocalOps` on `tBasisSet` + `SetSymmetryOps` on the GPW evaluator +
  `GPWParams.imposeSymmetry` (renamed from `reduceBZ` 2026-08-02): all exist because the run's §3
  symmetry policy rides the basis.  Verified sole consumers: CompositeWF.C:169 (handing ops to
  tComposite_CD's stream fold; empty for atoms/molecules — lattice policy riding the neutral
  base) and GPW_SCF_UT.C:310 (src/-dead diagnostic getter).  Governing record: the §3
  `SymmetryPolicy` object (SymmetryUpgradePlan / Lattice_3D-owned space group) — when plumbed,
  these collapse to consumers reading ONE context and the per-factory bool dies.
- **D5 Becke `MeshParams` ergonomics under `imposeSymmetry`** — `mp.nAngular` IS honored but
  `mp.angular` (the GL/Lebedev scheme) is silently REPLACED by the site-adapted rule (announced,
  not warned).  Wanted: `nAngular<=0` → auto-resolve the calibrated default; console warn on
  scheme override; a real error (not a bare assert) when the requested L is unachievable for a
  low-symmetry site (C1/Cs seed-pool exhaustion).  Lands with the `SymmetryPolicy`/facade pass.
  (The degree-typed `angularDegree` interface half is executable now → R2.15.)
- **D6 ✅ DONE 2026-08-07. `BeckeXCParams()` lives in the TEST file + `ResolveXCMesh` (test driver)** — the.  **→ doc/CleanupHistory.md**
- **D7 The `dynamic_cast` survey = FittingCleanupPlan §C** (the one surviving item there; the
  "I want more" vs "what are you" criterion is written there).  Run §C as part of THIS session —
  the cast findings above (V1.8, V1.9, V1.10) are its seed list; give survivors the custom
  exceptions (subsumes R2.5's throw work).
- **D8 Standing pin governing every fit-touching item here**: fit quality is measured by
  grid-convergence of ρ/property vs a fine reference — NEVER ΔE_total (fits are non-variational).
- **D9 The ρ̃ mixer still FUSES the preconditioner with the extrapolator** (opened 2026-08-10 by the
  joint-history fix, SymmetryUpgradePlan §7 step 7).  `PulayMixer` owns BOTH the Kerker filter and
  the history — `ApplyJoint` takes its own Kerker step — and `PolarizedDensityMixer` is still a
  `tDensityMixer` rather than the filter stage of someone else's step.  The end state, argued in
  the plan, is two objects:
      `residual → channel-basis preconditioner (may differ per channel) → joint extrapolator (one B, one c)`
  which is the VASP/QE/CP2K architecture (CP2K's `BETA 1.5` IS a Kerker preconditioner in front of
  ONE Broyden history).  What already landed: the history is joint (`tFieldExtrapolator` +
  `MixJointly`), so the fusion can no longer split a history — this item is the SHAPE, not a
  defect.  Doing it retires the `tDensityMixer` inheritance on `PolarizedDensityMixer` (which is
  what made "pure forwarding, without knowing which leaf it holds" look reasonable), makes Broyden
  a drop-in beside Pulay, and generalises to the 2×2 spin density matrix of the non-collinear /
  Shubnikov work, where "which channel combinations get which filter" survives and "two
  independent leaves" does not.  Re-measure run 11's `(ρ,m)` verdict when it lands: it was recorded
  REFUTED at `PulayDepth=0`, i.e. with no extrapolator for the preconditioner to shape.
- **D10 `DisplayEigen` builds rows and renders them in ONE pass, so nothing is assertable** (opened
  2026-08-10; user: *"printing those tables has been a constant source of bugs"*).  The bug list is
  the evidence: doubly-empty levels dropped below the frontier (fixed 08-08), `setprecision(0)`
  rounding a smeared 0.996 to "1/1" (same day), and an ABSENT channel's ε filled from the other
  channel so a row read as a level empty at an energy where the opposite spin is occupied, with
  ϵ↑−ϵ↓ printing exactly 0.00000000 (fixed 08-10 -- it manufactured MnO's spin-up "hole").  Every
  one of these is a pure-function property of the ROWS, and every one was found by reading a run
  log.  Extract the row build -- `(occ↑, ε↑, label, occ↓, ε↓, Δ)` per level, from the two
  `EnergyLevels` -- as a pure function, and the next one is a unit test instead.
  **Second, deeper defect the same fix should retire**: rows are paired by `(n, sym)`, but `n`
  indexes a DEGENERATE GROUP, and the grouping differs between spin channels once ϵ↑≠ϵ↓ -- which is
  why the MnO table skips indices and runs them out of order, and why levels go missing from one
  channel at all.  Pair by energy/character instead.
- **D11 The occupation RULE is configured AFTER the seed fill** (found 2026-08-10 by the shared-μ
  work; the cause of a live charge-losing bug, patched at the fill).  `tSCFIterator`'s constructor
  runs the seed fill (`itsWaveFunction->Init(...)`), while `SetMOM`/`SetSmearing` are called in
  `Iterate` -- a later call.  So **iteration 0 is filled under a different occupation rule than
  every subsequent iteration**, and nothing says so.  It was benign only as long as every fill path
  was integer-per-block; the moment a μ-SOLVING path existed it produced wrong electron counts,
  because at kT=0 the Fermi count is a staircase in μ and the target falls between steps (Al
  global-μ metal seed: Σw·n = 2.25 vs Ntot=3, on main, for as long as that path has existed; a
  shared-spin Mn sextet: 6↑/0↓ vs Ntot=7, and the ρ̃ mixers were then constructed on those counts).
  The guarding `assert(itsSmearingkT>0.0)` inside the fill is compiled out in Release, so it was
  silent.  **Patched** by gating the μ path on `kT>0` -- correct in its own right, since a seed has
  no self-consistent spectrum for a reservoir to redistribute over -- but the ORDERING is the real
  defect: an object should not perform its first fill before being told how to fill.  Fix by passing
  the occupation rule (SmearingkT / MOM) at CONSTRUCTION, or by moving the seed fill into `Iterate`.
  Then restore a hard failure (not a bare assert) on a μ fill with kT=0.
- **D12 `GPW_Evaluator::Eval` truncates by RADIUS, not by MAGNITUDE** (2026-08-11; the standing pin is
  [[feedback_no_cut_lattice_sums]] -- "THERE IS NO CUT": an ε-converged series with a magnitude
  screen, no radius in ANY interface).  Two sites, both in
  `src/BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C`:
  1. `itsMaxReach = sqrt(-log(1e-10)/MinExponent())` -- ONE global radius from the most DIFFUSE
     exponent, then applied to every function and every image.  ε-derived, so not a hard cut, but it
     is the worst-case exponent charged to all functions.
  2. `BuildImages(cell, max(2*maxReach + 2*maxCellEdge, maxReach + 2*cellRad), ...)` -- a GEOMETRIC
     bound that must be re-derived per cell SHAPE.  **It has already been wrong here**: the `max()`
     of two formulas exists because the historical one "under-enumerated" for the oblique MnO
     rhombohedral cell.  A patch, not a derivation -- the next unusual cell can break it again.
  **This is an INCONSISTENCY, not a missing capability**: the analytic 1E/V_local lattice sums in the
  same run already screen per-pair on magnitude and report the reach as an OUTPUT
  (`[lattice sums] eps=1e-10 (GPW_SCREEN_EPS) pair reach=30.3485 au = 381 cells`).  So one evaluator
  carries two truncation schemes, and only one of them obeys the pin.  Fix: screen each (function,
  image) term on its own contribution magnitude -- per-pair, like the analytic path -- and let any
  radius fall out as a reported consequence.  Then there is no formula left to get wrong.
  NB investigated 2026-08-11 as a candidate for the MnO sublattice defect and REFUTED as the cause
  (the enumeration is generous for this cell: maxReach ≈ 15.2 bohr against a = 8.4); this item is the
  DESIGN debt, which stands on its own.
