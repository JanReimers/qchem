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

**Harvested again 2026-09-08** (user: *"done and todo items interleaved ... hard to read and assess"*).
Fourteen items had been closed IN PLACE and never moved, holding **1165 lines of closed record among the
open ones**.  They are in `doc/CleanupHistory.md` under *"HARVEST 2026-09-08"*, verbatim, with the usual
stubs here.  ★ **Where a closed item had a SURVIVING REMAINDER, the remainder stayed in this file** — it is
open work — and only the closed record moved; the stub says so in its first line.  Verified line-for-line:
zero non-blank lines exist in neither file.

▶ **THE STANDING RULE this makes explicit:** an item's DONE verdict belongs in `CleanupHistory.md` the day
it is written, not "eventually".  A ✅ that stays in the worklist is indistinguishable from work.

---

# START HERE (handoff, 2026-08-17, updated after the V1.11 landing)

## The state that matters: ALL SIX RealComplexPlan prerequisites are DONE

`doc/RealComplexPlan.md` §7 lists what must land before the real/complex type refactor starts, because
those items sit on the very faces it restructures.  Status:

| prereq | state |
|---|---|
| **V1.6** `tDM_CD::Accumulate*` face split | ✅ `2d0f6982` |
| **V1.7** the periodic trio's 9 asserting overrides | ✅ `2d0f6982` |
| **V1.8** `IrrepCD`↔`IrrepCD` concrete casts | ✅ `2d0f6982` (evaporated with V1.6) |
| **V1.10** abstract→concrete basis casts | ✅ `2d0f6982` |
| **V1.1** basis-face merge | ✅ `d49db261`+`4b37221d`+`2bfb83b4` (2026-08-16) |
| **V1.5** `G_FieldEvaluator` ISP split | ✅ `f18a6ee9`+`9ebaebdb` (2026-08-16; §K blocker had dissolved) |
| **V1.11** occupation seam | ✅ `43bbebad`..`2398dd07` (2026-08-17, five increments) — **the LAST one, DONE** |

**What V1.1's landing bought (2026-08-16):** ONE structure-neutral DFT basis face.  The 3-centre tensor
is `Projector3<T>` (one struct; dense / delta-support / matrix-free are REALIZATIONS, a property of the
producing basis, not of the scalar type); `Band_FT_IBS` is deleted and the lattice lineage is
`Orbital_DFT_IBS<dcmplx,dcmplx>`; **`Orbital_DFT_IBS<double,dcmplx>` — a real TRIM k-block sharing the
run's complex G-space fit basis — is now a spelling the MnO real-TRIM work can instantiate**, which is
what that session is waiting on.  The cached accessors key `theCache<TFit>()` (the tensor's own type
axis), so the TRIM case caches correctly from day one.

**What the four 2d0f6982 items bought:** `tDM_CD<T>` has shed four exact-exchange methods and nine
periodic denials.  The complex path now carries NO exact-exchange machinery and the real path none of the
reciprocal-space kind -- so the plan's "widening points" are two small faces (`tHF_System_CD`,
`tHF_Pair_CD`) instead of denials scattered across three density families.  A live `-DNDEBUG` hazard died
with them: `Vee`/`Vxc` used to hit an assert-only `void` body -- a silent NO-OP in the shipped build,
i.e. a zeroed J and a wrong Fock -- and now throw.

## Next: the RealComplexPlan itself is UNBLOCKED — and MnO real-TRIM with it

Every §7 prerequisite has landed.  The type refactor's staging (RealComplexPlan §5) can start at Step 1
(`IsReal()` queries), and the MnO session's real-TRIM work (which ruled 2026-08-16/17 that it WAITS on
V1.1 + V1.11) has everything it needs: `Orbital_DFT_IBS<double,dcmplx>` is a live spelling and the WF/
occupation seam carries no state the variant-child restructuring would entangle.
- **V1.5 landed 2026-08-16** (see the item / doc/CleanupHistory.md).  The session's status audit of
  FittingCleanupPlan found it more finished than it knew: **H** essentially landed via V1.26/V2.4,
  **I.2** landed (`relCutoff` live), **§K's `cvec_t` sub-bullet is overtaken** (`ProjectedScalar_G` no
  longer exists) — so **K's core (the fit-{G} densification) is UNBLOCKED but DEFERRED (user ruling
  2026-08-17, same reason as V1.22: deliberately non-bit-identical on grids every GPW run consumes,
  while the real-TRIM campaign runs against pinned numbers — wrong week; ready at the next campaign
  breakpoint)**.  **I.1's residual**:
  the `GetEpsXc()=0.75*GetVxc()` base default (ExchangeFunctional.C:34) — exact for Dirac exchange, a
  silent-wrong inherited default the day a GGA functional forgets to override; fold into whatever touches
  the functionals first.

## Coordination — CHANGED 2026-08-16 (user); UPDATED 2026-08-17 (real-TRIM running, concurrent cleanup)

**The DO-NOT-TOUCH on the MnO working set is LIFTED** — the MnO session is now WAITING on this branch's
V1.1 work to implement real TRIM irreps at the special k-points ((0,0,0), (½,½,½), …).  User verbatim:
*"DO NOT TOUCH is now please -TOUCH."*  Worth pushing early for that reason.

**2026-08-17: the real-TRIM session is RUNNING (off main, RealComplexPlan §5).**  Its working set, by
staging step — a concurrent cleanup session should stay OUT of these:
- Step 1 (now): `src/Symmetry` (`BlochQN::IsReal`), the `IrrepBasisSet` face, **every Hamiltonian term
  file** (`PreservesReal()` added per term), the run report, the GPW evaluator (`IsTRIM`).
- Step 2 (next): `src/ChargeDensity` composites (`tComposite_CD`/`IrrepCD` child slot), `src/SCFAccelerator`.
- Step 3+: `src/BasisSet/Lattice_3D`, `src/WaveFunction` (the variant child).  Plus `GPW_SCF_UT.C` runs.
Consequently ALSO parked for now (they live in those files): V1.12, V1.17, V2.1, V2.3, R2.14's remaining
renames, the I.1 residual; R1.0 stays "own session" (user), R2.21 stays deferred (user); **§K deferred
(user 2026-08-17, the V1.22 reason: non-bit-identical while the campaign runs — not just "wants the box")**.

**Concurrent-cleanup assignment (2026-08-17): V3.1 + V3.2 (the atomic-solver bugs — atom path only),
R2.20 (oracle helpers out of the test module), then the qcMesh batch R2.15 (+V2.7 riding along;
V1.22 optional, it is bit-SENSITIVE).**  All verification beside the TRIM runs goes through
`scripts/memsafe -H 6G ctest -j2` (box discipline).

## Also open, unruled, and independent of the type work

- **V1.27's remaining half 🔶 DESIGN RULED, not started** -- the {Molecular,Solid}x{PP,Non-PP} iterator
  decomposition.  NOTE it collapsed once: the virial axis is a DERIVED BOOL (`IsVirialValid()`), not a
  class axis, so it is 2 classes + 1 bool, and `AccumulateColumns` already landed the column-list part.
  What remains is whether the grid and periodic axes ever need splitting (see the item).
- Then the R2/V backlog as before.

## PROCESS -- two habits that cost real time this session, both recurring

- **Never pattern-edit a file you have not read.**  Four blanket-replace slips: a Python whole-file
  rewrite flipped CRLF→LF (fixed structurally by `.gitattributes`, `de4a4663`); a `\s*delete cd;` regex
  removed three UNRELATED deletes and mangled an if/else in the parallel campaign's file; `self().` swaps
  hit a class's own private helpers and a constructor's member-init list.  The check that works, and is
  cheap: list every match with its ENCLOSING FUNCTION before replacing anything.
- **Never leave a `ctest` sweep running while editing or rebuilding.**  Twice.  The second time produced
  22 failures that were pure relink artifact -- and were nearly reported as a regression.  A sweep is only
  evidence if nothing changed under it.

Six things worth carrying forward, because they generalize past their own items:

- **V1.6/V1.8 — name a partner face for the OPERATION, not the data.**  `CompleteDirectPair(Ji,Jj,Di,bs_i)`
  says "finish this contraction with your block"; an abstract block ACCESSOR would have compiled, passed
  review, and quietly undone the project's own no-`GetDensityMatrix()` tenet.  Naming it for the operation
  is also what made the concrete cast EVAPORATE rather than move: the caller needs the partner's
  cooperation, never its type.
- **V1.6/V1.7 — a capability face is not segregated while the other half still DECLARES the methods.**
  Empty bodies beat asserts (an assert-only `void` is a silent NO-OP under `-DNDEBUG`; an unreachable
  empty body cannot be entered at all), but both are the interface failing.  The shape that actually
  removes them is a real-path-only CRTP mixin selected by `conditional_t` — the idiom `ProjectedDensityBase`
  and `FourierDensityBase` were already using in the same file.


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

**SUPERSEDED 2026-08-16 (see the START HERE section): the DO-NOT-TOUCH is LIFTED — the MnO session is
waiting on this branch's V1.1 work for real TRIM irreps.**  *(The old rule, for the record:)*
~~**DO NOT TOUCH** (their working set): `src/SCFAccelerator`, `src/SCFIterator`, `src/WaveFunction`,
`IntegrationTests/GPW_SCF_UT.C`, `src/BasisSet/Lattice_3D/`.  If a task needs one of them, STOP and ask
rather than reaching in — a mid-flight collision there costs them a multi-hour run.~~
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

  **SCOPED 2026-08-22 (user): STEP (1) ONLY — step (2) is explicitly left out** of the session that does
  this.  Making atomic `op(r)` honest (real Y_lm combinations) is the riskier half, it drags the
  density/XC paths that rely on the radial-only form, and NOTHING in step (1) needs it: step (1) removes
  the consumers, which is what kills the bug class; the fake then sits behind a face nobody calls.

  **WHAT THE XC SEPARATION WORK ADDED TO THIS ITEM (2026-08-22, commit `b5d5d363` + the discussion after).**
  - **The cure already has a working precedent, not just a design.**  `Fit_IBS` ALREADY answers
    operations — `Overlap(Sf)` = \f$\langle f_a|f\rangle\f$ and `Norm()` — and uses `op(r)` only
    INTERNALLY, on its own mesh, to compute them.  So "ask the basis for the integral" is not a new shape
    for fit bases; it is the shape they were already in, with `op(r)` leaking as a second, redundant face.
  - **`DeltaFit_IBS` is step (1) applied to one basis.**  It has no `op(r)` to fall back on (a delta is a
    distribution), so it exposes only operations.  That is why the tension surfaced there first.
  - ⚠ **A NAMING TRAP in the `MakeOverlap` fold-in.**  For a delta basis
    \f$\langle\delta_g|f\rangle = w_g f(r_g)\f$, so its `Sample` (which returns the bare values
    \f$f(r_g)\f$, because the XC functional is applied POINTWISE) is **not** `MakeOverlap` — it is the
    FIT COEFFICIENT vector in the diagonal metric, \f$c_g=\langle\delta_g|f\rangle/\langle\delta_g|
    \delta_g\rangle = f(r_g)\f$.  The `MakeOverlap`/`MakeOverlap3C` ruling covers the
    \f$\langle\chi_i|f\rangle\f$ / \f$\langle\chi_i|f|\chi_j\rangle\f$ integrals; the delta
    route's per-point verb belongs to the FIT vocabulary instead.  Decide which of `Integrate` /
    `SiteIntegrals` / `Sample` / `Symmetrize` (landed in `b5d5d363`, ahead of this fold-in) keep their
    names and which become `MakeOverlap*` overloads — do not blanket-rename on the family name alone.
  - **THE ENABLING MOVE (user): take `VectorFunction<T>` OFF `IrrepBasisSet<T>` and put it on
    `Orbital_1E_IBS<T>`.**  Surveyed: only FOUR sites evaluate a basis pointwise — `TOrbital.C:47`
    (psi(r)), `IrrepCD.C:543` (rho(r)), `PWTerms.C:702` (the Phi table, batched) — all ORBITAL — and
    `FunctionFitterImp.C:27` (the fitted field \f$\sum_a c_a f_a(r)\f$), the ONLY fit-basis one, which
    becomes an operation `EvalFit(c,r)`.  The plane-wave lineage already proves that shape:
    `G_FieldEvaluator::EvalField(coeffs,r)` exists precisely because a \f${G}\f$ basis cannot answer it
    through `op(r)` either.  NB the removal is from the ABSTRACT FACE, not from concrete classes: `Fit_IBS`
    keeps deriving `VectorFunction<double>` explicitly (it needs it for its own `qcMesh::Overlap`), it just
    stops PROMISING evaluation to consumers.  Bonus: `GetVectorSize()` exists only to satisfy
    `VectorFunction` and shadows `GetNumFunctions()` — one of the pair goes.
  - **IT CLOSES THE XC ITEM'S LAST ESCAPE, as an instance rather than a special case.**  The one surviving
    getter after `b5d5d363` is `Quadrature::Mesh()`, needed because `cDM_CD::DM_RhoAtPoints(points,
    Phi-tables)` is initiated by the DENSITY (it owns D, and lives above qcBasisSet).  But
    \f$\Phi_{gi}=\langle\delta_g|\chi_i\rangle/w_g\f$ is a CROSS-BASIS OVERLAP — exactly an
    "ask the basis for the integral" (`Repulsion(const rFIT_CD_ABS&)` is the existing cross-basis
    precedent) — so under step (1) `DM_RhoAtPoints(points, tables)` becomes `DM_RhoAtPoints(basis)`, Phi
    and BOTH contractions come into the basis, and no coordinate leaves anywhere.
  - **ORDER:** the naming fold-in + the `VectorFunction` relocation first; the rho flip then falls out.
    Step (2) stays out.

  **✅ THE RELOCATION LANDED (2026-08-22, step (1)'s enabling move).**  `VectorFunction<T>` is off
  `IrrepBasisSet<T>`; a basis that really can be evaluated says so by deriving the new
  `BasisSet::Evaluatable_IBS<T>` (which also absorbs the `GetVectorSize()`/`GetNumFunctions()`
  duplication).  Who promises it now: every ORBITAL lineage (atom, the three molecular Gaussian ones,
  the PW/GPW mixin, the two DHF halves) plus `Fit_IBS`, which needs its own values to compute its own
  mesh integrals.  Who does NOT: `FIT_SF_ABS` / `FIT_CD_ABS` — the neutral fit faces — and `DeltaFit_IBS`,
  whose two THROWING `op(r)`/`Gradient` overrides are simply GONE: they existed only because the interface
  demanded an answer a distribution cannot give.  The last fit-basis consumer, the AO fitter's
  \f$\sum_a c_a f_a(r)\f$, became `FieldEvaluator<double>::EvalField(c,r)` on the NON-ORTHO faces (the
  real-space sibling of `G_FieldEvaluator::EvalField`, which already had exactly this shape for
  \f$\{G\}\f$).  Only TWO pointwise call sites remain in the whole tree, both orbital and both honest:
  `TOrbital.C:47` (\f$\psi(r)\f$) and `IrrepCD.C:543` (\f$\rho(r)=\phi^\dagger D\phi\f$).  756/756.

  **⇒ THE TARGET IS LISKOV SUBSTITUTABILITY OF ALL FIT BASES (user, 2026-08-22), AND THE ALGEBRA SAYS IT
  COSTS NOTHING.**  The molecular XC term uses exactly TWO operations (`Imp/FittedVxc.C:74-100`): `DoFit`
  of a COMPOSED field (`VxcDensity` = \f$v_{xc}\circ\rho\f$ — the basis samples it on its OWN mesh, so no
  points escape and no "sample" method is needed), and `Overlap(orbitalBasis)` =
  \f$\sum_a c_a\langle O_i|f_a|O_j\rangle\f$.  Put the δ basis through the same two, using its own
  integrals (\f$\langle\delta_g|f\rangle=w_g f_g\f$, \f$\langle\delta_g|\delta_h\rangle=w_g
  \delta_{gh}\f$, \f$\langle\chi_i|\delta_g|\chi_j\rangle=w_g\chi_i(r_g)\chi_j(r_g)\f$):
  - `DoFit(`\f$v_{xc}\circ\rho\f$`)` → \f$c=S^{-1}\langle\delta|v\rangle = v_{xc}(\rho(r_g))\f$ —
    which is what `Sample` + the functional does today.  **`Sample` IS `DoFit`.**
  - `Overlap(orb)` → \f$\sum_g c_g w_g\chi_i\chi_j = \Phi^\dagger\mathrm{diag}(wv)\Phi\f$ — *identically*
    the singles quadrature.
  - \f$E_{xc}=\sum_a e_a\langle\rho|f_a\rangle\f$, and \f$\langle\rho|\delta_g\rangle=w_g\rho(r_g)\f$,
    so **`Integrate` is that same contraction** — and it stays \f$O(n_{pts})\f$, not \f$O(n_{pts}n^2)\f$,
    because the identity \f$\mathrm{Tr}(D\sum_a e_a\langle i|f_a|j\rangle)=\sum_a e_a\langle\rho|f_a
    \rangle\f$ lets δ take the cheap side (its ρ samples are already in hand as the functional's argument).
  - `SiteIntegrals` is NOT part of the XC contract at all — an atomic-moment observable; it belongs to
    whoever owns the partition, not to the fit face.

  **⇒ AND THE TWO PERFORMANCE OBSTACLES CONFORM THE *OTHER* WAY (user ruling, 2026-08-22).**  Both are
  about ρ sampling, and in both the PERIODIC pattern is the better one, so the MOLECULAR route moves:
  - ρ through the **D-GEMM** against a cached Φ, not pointwise through the density's `ScalarFunction`
    ("small improvement" for molecules too, user).  Needs the Φ tables where the density can ask for them
    — i.e. the same `DM_RhoAtPoints(basis)` flip.
  - **sample ρ ONCE and derive both** \f$v_{xc}\f$ and \f$\epsilon_{xc}\f$ from it, as the periodic
    route does, instead of fitting the two fields with two independent samplings (`itsFitter` +
    `itsEpsFitter` in `FittedVxc`).

  ⚠ **`isOrtho()` IS CORRECTLY NAMED — the older DOC was wrong** (user, 2026-08-22; corrected in place).
  It asks ORTHOGONAL (metric DIAGONAL ⇒ no linear SOLVE), not orthoNORMAL (metric = I).  A plane-wave
  basis happens to be orthonormal; a δ basis is orthogonal with \f$\langle\delta_g|\delta_g\rangle=w_g\f$.
  Both answer `true`.  The distinction bites the moment a diagonal-but-not-unit basis is fitted through the
  general path: \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$, and dropping the denominator
  (right for orthonormal, wrong here) is an error of a factor \f$w_g\f$.

  **✅ INCREMENT 1 LANDED: the density asks the basis (`DM_RhoAtPoints(quadrature)`).**  The Φ tables and
  their build moved OUT of the XC quadrature and INTO the δ basis, behind a new `BasisSet::Collocation`
  face (`NumPoints` / `Values(orbitalBlock)` / `Sample(field)`) — because Φ is an integral over the
  basis's own functions, which is the R1.0 rule.  `DM_RhoAtPoints(points, Φmap, ΦRmap)` became
  `DM_RhoAtPoints(q)`: each density block asks `q` for ITS table with ITS OWN scalar, so the mixed-run
  (3c-3) second map argument is gone, the "block not yet tabled ⇒ fall back to pointwise" branch is gone
  (there is no first pass to heal), and with it the whole `ensureBlock` hint — `Rho`/`RhoPol` lost their
  four ensure overloads down to two methods.  756/756, Si gate unmoved at 11/11 iterations.

  **What is left of the coordinate escape: ONE line** — `itsFit->Mesh().Weights()` in the singles
  quadrature's \f$\Phi^\dagger\mathrm{diag}(wv)\Phi\f$ GEMM.  It closes when `Matrix` itself becomes a
  basis operation, i.e. in increment 2 (δ conforms to `DoFit`/`Overlap`), where that contraction is
  `Overlap(orbitalBasis)` on the fit basis and the weights never leave.  The negative-ρ diagnostic already
  stopped touching weights: it asks for its two masses as integrals.

  **✅ INCREMENT 2 (part): the quadrature contraction is the BASIS's, and the metric diagonal is explicit.**
  - \f$\langle\chi_i|v|\chi_j\rangle=\sum_g w_g\overline{\chi_i}v_g\chi_j\f$ moved out of the XC
    strategy and onto the δ basis as `Quadrature(orbitalBlock, v)` (two overloads, mixed-run as `Values`).
    It owns the points, the weights AND the Φ table, so the whole contraction is its operation and a caller
    hands only the bare field.  **That closed the last coordinate escape**: nothing outside the basis reads
    a weight or a point any more.
  - `FIT_SF_ABS<T>::OverlapDiagonal()` (user directive): PW → ones, δ → \f$w_g\f$, Gaussian →
    \f$1/\mathrm{Norm}_a^2\f$.  Pinned by `GPW.OverlapDiagonalPerRepresentation` — the two periodic cases
    against their definitions, with an explicit `EXPECT_NE(1)` on the δ one so the distinction stays
    load-bearing rather than decorative.
  - `OrthoScalarFitter` → `OrthoNormalScalarFitter` in its own
    `src/Fitting/Internal/OrthoNormalFunctionFitter.C`, **named for the metric it assumes**: it takes
    `OverlapDiagonal()=={1,1,...}` and never calls it.  (For a raster that diagonal is npts ones nobody
    should materialise, and the forward FFT already delivers the coefficients.)

  **✅ AND THAT BLOCKER IS GONE — `FunctionFitter_Scalar` LOST ITS TYPE PARAMETER** (user, 2026-08-22).
  `T` appeared in exactly ONE place on that face, `hmat_t<T> Overlap(const robs_t<T>*)`, because the
  scalar fitter's INPUT is an untyped real-space FIELD.  (Contrast `FunctionFitter_Density<T>`, whose
  `DoFit` takes a typed `ProjectedDensity<T>` — it keeps its parameter, and the asymmetry is the tell.)
  So the class-level `T` was standing for TWO independent things: what my fit basis is made of, and what
  scalar the orbital block I contract against uses.  Those coincided — real Gaussians with real orbitals,
  complex plane waves with complex orbitals — until 3c-3 made mixed runs real, which is the same failure
  mode as the fit/grid welding: one parameter for two axes, and it only shows when the axes disagree.
  Split per ISP into a neutral `FunctionFitter_Scalar` (DoFit / ReScale / Write / the evaluatable field)
  plus `FitContraction<U>`, which a fitter declares once per block scalar it can actually serve: the
  molecular Gaussian fit declares `double` only (its basis has no Bloch 3-centre path), a raster fit
  declares `dcmplx`, and a δ fitter will declare both.  No fitter carries an `Overlap` it would have to
  throw from, and consumers cross-cast to the contraction they need.  757/757, molecular and periodic.


  **✅ INCREMENT 2 COMPLETE: δ GOES THROUGH THE FITTER, AND THREE FACES DISSOLVED.**  `Fitting::Factory`
  now returns a `DeltaScalarFitter` for a δ basis, and `XC_SinglesQuadrature::Matrix` is *fit the sampled
  field, then contract against this block* — the same two calls the molecular XC term makes, on either
  representation.  The strategy no longer performs a quadrature; it composes one.  What that let go:

  | removed | why it could go |
  |---|---|
  | `BasisSet::Quadrature` (whole face + module) | its `Mesh()` getter had three consumers, all of which wanted an OPERATION: `NumPoints` / `Sample` / `Integrate` now sit on `FIT_SF_ABS<T>`, answered by all THREE fit bases (Gaussian over its Becke mesh, PW over its raster, δ over whichever mesh it was built on). No mesh changes hands anywhere. |
  | `BasisSet::Collocation` (whole face) | `NumPoints`/`Sample` went up to `FIT_SF_ABS`; `Values` (the Φ tables) went onto `FIT_SF_Delta`, which is the only representation that tabulates. `DM_RhoAtPoints` takes that face directly. |
  | `FIT_SF_Delta::Integrate` | every fit basis can integrate a field sampled at its own points — it was never δ-specific. |
  | `FunctionFitter_Scalar : ScalarFunction<double>` | a δ fit has no value BETWEEN its points (Σc_gδ_g is a distribution).  Same lesson as `IrrepBasisSet` losing `VectorFunction`: the fitters that CAN evaluate derive it themselves; the one consumer (a test) cross-casts. |
  | `XC_PairQuadrature::Mesh()`, `MatrixT`, the field's point-identity check | all downstream of the above. |

  **`FIT_SF_Delta` is down to FOUR members**, and each is now something only a δ representation can
  answer: `Values` ×2 (the tables the density collocates against), `Quadrature(orb,v)` (the weighted
  contraction the fitter drives), `SiteIntegrals` (the atomic partition — still the odd one out, an
  observable rather than part of the XC contract), and `SymmetrizeSpin` (the magnetic pair projection, δ
  being the only representation a polarized run can use).  757/757, Si two-route gate bit-unmoved.

  ⚠ **TWO CORRECTIONS TO THE ABOVE (user, on reading it).**
  - **`Sample` is not an overlap integral.**  It returns `NumPoints()` values, \f$f(r_g)\f$, one per
    POINT; the projection \f$\langle f_a|f\rangle\f$ is `FIT_SF_NonOrtho::Overlap(f)`, one per FUNCTION.
    They have the same length only for a δ basis (functions = points) and even there differ by the weight,
    \f$\langle\delta_g|f\rangle=w_g f(r_g)\f$.  Documented on the face now, with the contrast spelled out.
  - ⛔ **"a fit basis is a family of weight vectors over shared points, so it always HAS points" is WRONG**
    as a general statement, and it was justifying the three quadrature ops on the neutral face.  That is
    the δ picture: there the functions ARE the points and the mesh is constitutive.  A GAUSSIAN auxiliary
    basis is a family of FUNCTIONS, and the mesh it carries is the DEVICE it computes its own integrals
    with — if those integrals were analytic it would have no points at all.  The honest justification, now
    in the code, is the narrower one: *every scalar fit basis in this project integrates NUMERICALLY on a
    quadrature it owns*, which is a property of the three implementations and not a claim about fit bases.
    An analytic-integral fit basis would be the counterexample that pushes these three onto a narrower face.
  - And the related fact, recorded rather than hidden: **on the molecular lineage `NumPoints`/`Sample`/
    `Integrate` have NO caller yet.**  They get one in the ruled molecular conformance (ρ through the
    D-GEMM, sampled once with both functionals applied to the values), which is the next increment.

  **✅ INCREMENT 3 LANDED 2026-08-23 — ONE FIT-BASIS INTERFACE: THE POINT VOCABULARY IS OFF THE FIT FACE.**
  The user's ruling that opened it: *"a DeltaFitBasis is also a family of functions … as such it should
  have exactly the same interface as the Gaussian fit basis set."*  `NumPoints` / `Sample(f)` /
  `Integrate(values)` described a QUADRATURE, and read as correct only because for δ n_functions ==
  n_points, so the wrong accessor returned the right number.  `FIT_SF_ABS<T>` now carries exactly two
  integrals, both **per FUNCTION**, and every representation answers both:

  | | Gaussian (`Fit_IBS`) | δ (`DeltaFit_IBS`) | plane wave (`PlaneWaveFit_IBS`) |
  |---|---|---|---|
  | count | `GetNumFunctions()` (already there — the removed `NumPoints` was its shadow) | n_pts | n_G in the ball |
  | `Overlap(f)` = ⟨f_a\|f⟩ | its existing mesh quadrature, moved UP from `FIT_SF_NonOrtho` | \f$w_g f(r_g)\f$ | \f$\sqrt\Omega\,\tilde f(G_a)\f$ (sample + forward FFT + ball gather) |
  | `Integrals()` = ⟨f_a\|1⟩ | **is** `Charge()` — one quantity, two faces, so it forwards | \f$w_g\f$ | \f$\sqrt\Omega\,\delta_{G,0}\f$ |

  Each fitter then applies its own metric to that one projection, which is where **`OverlapDiagonal()`
  finally got its production consumer**: `BasisSet::OrthogonalFit(basis,f)` = ⟨f_a|f⟩/⟨f_a|f_a⟩ is now the
  δ fitter's whole `DoFit`, and the same one-liner is the default `tDM_CD::DM_RhoAtPoints` (a matrix-free
  density expressing itself over a δ basis).  E_xc is *coefficients · `Integrals()`* inside
  `XC_SinglesQuadrature`, not a weight sum — algebraically the same and, written as the same loop, the
  same summation order.

  **VERIFIED, and the two risks named in the spec both came back clean.**  758/758 (count UP by the new
  gate, not down), and the pinned Si two-route SCF prints **−7.115067665** (pair) / **−7.115059008**
  (singles) at 11/11 iterations, bit-unmoved.
  - The ⚠ BITWISE pin was real but did not bite: δ's fit is now `fl(w_g f_g)/w_g` where it used to be
    `f_g` directly, and that is NOT exact — but the perturbation is ~1 ulp with random sign over ~10⁵
    coefficients, i.e. ~1e-13 Ha against a gate that prints 1e-10.  Measured, not assumed.  (Two things
    were needed to keep it that small: divide on the REAL parts, since `std::complex` (x,0)/(y,0) goes
    through xy/y² and is worse; and keep `qcMesh::Integrate`'s summation order in the new dot product.)
  - **NEW FINDING — the unification claim is TRUE but the plane-wave fitter must NOT use it, and the
    reason is not bit-preservation.**  `OrthoNormalScalarFitter::Overlap(bs)` looks \f$\tilde v\f$ up at
    ORBITAL index differences \f$m_i-m_j\f$, which run to TWICE the orbital ball — outside the fit ball,
    where a projection onto the fit basis's own \f$\{G\}\f$ honestly reads ZERO while the raster's
    `GridCoeff` wraps and returns the resolved value.  Fitting through `Overlap` would therefore delete
    most of \f$v_{xc}\f$, not move a last bit.  This is `G_FieldEvaluator::ProjectField`'s
    truncation-vs-aliasing warning meeting a real call site: for assembly purposes that fitter's "fit
    basis" is the RASTER's \f$\{G\}\f$, not the declared ball, and it now says so by asking the raster
    face for both halves.  `PlaneWaveFit_IBS::Overlap`/`Integrals` are still implemented (the face is the
    face) and pinned by the new gate.
    ⚠ NB this is NOT a name collision with a self-overlap: `Overlap(f)` is the PROJECTION and the metric
    \f$\langle f_a|f_b\rangle\f$ is `Integrals_Overlap::Overlap()`, and a plane-wave fit basis **has no
    metric member at all** — it rides `EPW_Irrep_IBS` (op(r)/Gradient/GetNumFunctions), not the orbital
    `EPW_Orbital1E_IBS` tier that carries `MakeOverlap` (user, 2026-08-23; verified).  Only a class that
    genuinely carries BOTH — `FIT_SF_NonOrtho` and `Fit_IBS` — needs the `using` that un-hides one past
    the other.
  - Consequently `G_RasterTransform` grew `RasterSize()` / `Sample(field)` / `Integral(values)` — the
    three raster-array questions that used to ride the fit face.  Point vocabulary is CORRECT there: that
    face is keyed by voxel and integer reciprocal index by construction, and a plane-wave basis's voxel
    count is a different number from its function count, which is precisely the conflation being removed.
  - **NEW GATE `GPW.FitProjectionAndIntegralsPerRepresentation`** pins the claim rather than asserting it:
    δ's ⟨δ_g|f⟩ and ⟨δ_g|1⟩ against their definitions, PW's ⟨e^{iG}|1⟩ = √Ω δ_{G,0}, and — the
    representation-independent invariant the XC energy actually rides — ∫f == c·⟨f_a|1⟩ on BOTH
    representations for two fields of known cell integral.  That last one is what would catch the
    \f$1/\sqrt\Omega\f$ normalisation slip the plane-wave side is exposed to.
  - **Secondary finding, recorded not fixed:** `Fit_IBS::OverlapDiagonal()` is the ONLY member of the
    Gaussian fit face in the UN-normalised convention (\f$1/\mathrm{Norm}_a^2\f$), while `Overlap(f)`,
    `Charge()`/`Integrals()` and the `InvOverlap()` metric all fold the norm in.  Latent, not live —
    `isOrtho()==false` sends every Gaussian fit to the S⁻¹ solve, which never reads the diagonal — but a
    future general orthogonal fitter pairing the two would be wrong by \f$\mathrm{Norm}_a^2\f$.  Not
    flipped here because the consistent answer is all-ones, which makes it indistinguishable from the
    plane-wave one, and the gate that keeps the orthogonal-vs-orthoNORMAL distinction load-bearing is on
    the periodic pair: the fix wants its own measurement.

  **⇒ NEXT INCREMENT, SPECCED AND NOT BUILT (user ruling, 2026-08-23): SEPARATE THE METRIC AXIS INTO
  FACES, BOTH SIDES TOGETHER.**  The `Fit_IBS::OverlapDiagonal()` convention clash above is a symptom, not
  the disease: `OverlapDiagonal()` sits on the metric-NEUTRAL face, so a basis that never has a diagonal
  metric is obliged to invent an answer — and the answer it invented is in a different normalisation from
  every other member of its own face.  Move the question to where it is always meaningful and the wrong
  answer stops existing.

  **THE SHAPE (decided).**
  - New `FIT_SF_Ortho<T>` carrying `OverlapDiagonal()`; `FIT_SF_NonOrtho` keeps `Norm()`/`InvOverlap()`.
    `FIT_SF_ABS<T>` keeps only what EVERY fit basis answers: `Overlap(f)`, `Integrals()`, `Symmetrize()`.
  - `Fit_IBS` **loses `OverlapDiagonal()` outright** — the un-normalised member is deleted, not corrected,
    which is the only fix that cannot drift again.  `FIT_SF_Delta` and `PlaneWaveFit_IBS` derive the new
    face; `OrthogonalFit(b,f)` takes `const FIT_SF_Ortho<T>&` and its `assert(b.isOrtho())` disappears
    because the parameter type now carries the guarantee.
  - **The CD side moves in the same increment** (user): `FIT_CD_ABS::isOrtho()` has the identical shape and
    the identical Factory pattern, and its `FIT_CD_NonOrtho` Coulomb-metric split is the mirror image.
    Doing one and not the other leaves the two fit faces asymmetric, which is what the ISP work has been
    steadily removing.
  - **NO orthonormal marker face.**  Orthonormal is orthogonal with a unit diagonal, which
    `PlaneWaveFit_IBS` already answers on the Ortho face; a memberless `FIT_SF_OrthoNormal` would be the
    null-object pattern rejected on 2026-08-01 and again on `FIT_SF_Delta`.  `OrthoNormalScalarFitter`
    keeps stating its assumption in its name and its assert.

  **⚠ THE ACCEPTANCE CRITERION, AND IT IS THE POINT OF THE ITEM (user, 2026-08-23).**  *"We can remove
  `isOrtho` — but only if it does not get immediately replaced by SOLID/LSP-violating if statements"* of
  the form `if (dynamic_cast<FIT_SF_NonOrtho*>(fbs)) {non-ortho stuff} else {ortho stuff}`.  That is a type
  switch wearing a cast, and it is worse than the bool it replaced.  So the increment is **rejected** if
  any such branch appears.  Three facts make that a low risk rather than a hope:
  - **MEASURED: all EIGHT `isOrtho()` call sites are `assert`s.  There is not one live branch on it in the
    tree** (`OrthogonalFit`, `DeltaScalarFitter`'s ctor, and six in the four `Fitting::Factory` overloads).
    So retiring it has nothing to replace — each assert is a run-time re-check of a contract the type
    system would carry for free, exactly the R2.10 two-phase-construction smell one level up.
  - **Narrowing a PARAMETER TYPE is not a cast and not a branch.**  `OrthogonalFit(const FIT_SF_Ortho<T>&)`
    is the sanctioned replacement: the caller must already hold the capability, the compiler checks it, and
    no code asks a question at run time.  Prefer this everywhere over "cross-cast, then test".
  - **The ONE genuine branch in the tree is NOT a metric branch and the split does not touch it.**
    `Factory(cFIT_SF_ABS&)` picks `DeltaScalarFitter` vs `OrthoNormalScalarFitter` by
    `dynamic_pointer_cast<cFIT_SF_Delta>` — a "what IS it" cast.  Both of those bases are orthogonal, so
    the metric faces cannot distinguish them: the difference is the REPRESENTATION (δ divides by its
    diagonal and contracts through a Φ table; PW batch-projects by FFT and contracts through
    `Overlap3C` + the raster's aliased `GridCoeff`).  It must therefore be left exactly as it is by this
    increment, at the creation boundary where a factory is allowed to know, and NOT laundered into a
    metric test that would read as principled and be wrong.

  **⇒ AND THE OPEN QUESTION THAT WOULD DELETE THAT LAST BRANCH — needs a ruling before anyone attempts it.**
  The scalar fitter has exactly two operations, and on BOTH representations both are already delegated
  straight back to the basis (`DoFit` → `OrthogonalFit`/the ortho projection; `Overlap(orb)` →
  `FIT_SF_Delta::Quadrature` / `orb.Overlap3C` + the raster).  If those became basis operations outright —
  `Fit(field)` and `Contract(orb, c)` — there would be ONE fitter class, zero casts and zero branches, with
  the polymorphism living where the metric and the representation both actually are.  That is
  "replace conditional with polymorphism" taken to its end, and it needs no library-DAG inversion (the
  operations sit on the basis; the basis does not construct a fitter, which would make qcBasisSet depend on
  qcFitting and close a cycle).  ⚠ But it dissolves `FunctionFitter_Scalar` as a polymorphic type, which is
  a much bigger claim than a metric split — hence: recorded, not assumed.

  ⛔ **RETRACTED, SAME DAY: "row 4 deliberately did not land" — it was solving the wrong problem.**  I
  wrote that replacing `Values` with \f$\langle\delta_g|\chi_i\chi_j\rangle\f$ needs an `Overlap3C(δ)`
  overload on EVERY orbital lineage.  That follows only if the ORBITAL basis receives, which I assumed
  from the spec's suggested mechanism.  The user's correction: *"the delta fit basis should be able to do
  \f$H_{xc}(i,j)=\langle i|\text{delta fit to }V_{xc}|j\rangle\f$ **without** calling
  `orb->Overlap3C(delta fit basis)`.  In principle it just needs `op(r)` from the orbital basis; in
  practice it uses cached \f$\Phi\f$, but that is an implementation detail."*  Verified in the code: the
  contraction's only call into the orbital basis is `(*bs)(sub)`, the batched point-set `op(r)`.  With the
  FIT basis receiving, no orbital lineage is touched at all.

  **✅ SO IT LANDED, 2026-08-23.  758/758, Si two-route gate bit-unmoved (11/11 iterations).**

  **THE RULE THAT ORGANISED IT (user).**  A fit basis is not a public quadrature-integration engine.  It
  has exactly two jobs: **(1)** provide the integrals needed to do the fit; **(2)** provide integrals
  *over itself, not over third parties* so a Hamiltonian term can form \f$H_{ij}\f$/\f$E\f$.  And the
  house names for integrals are fixed: \f$\langle i|j\rangle\f$, \f$\langle i|f|j\rangle\f$,
  \f$\langle i|f\rangle\f$ are all **Overlap**; \f$\langle i|1\rangle\f$ is **Charge**; two-electron is
  **Repulsion**.  Words like grid / raster / quadrature / collocate / values are implementation, and a
  comment on one of these must say which \f$\langle\cdot|\cdot|\cdot\rangle\f$ it returns.

  | was | is | why |
  |---|---|---|
  | `FIT_SF_ABS::Integrals()` | **`Charge()`** | \f$\langle f_a|1\rangle\f$ is a Charge by the house convention; it is the SAME declaration `FIT_CD_NonOrtho` carries, so `Fit_IBS` satisfies both faces with one override (the CD one went by-value to make that legal) |
  | `FIT_SF_Delta::Quadrature(orb,v)` | — | it named the mechanism, and it took \f$v\f$ |
  | `FIT_SF_Delta::Values(orb)` | — | \f$\Phi_{gi}=\chi_i(r_g)\f$ is **not an integral** (the integral is \f$w_g\chi_i(r_g)\f$), so no name fits it: a value table has no business in an interface |
  | both | **`Overlap3C(orb)` ×2** | \f$\langle\chi_i|\delta_g|\chi_j\rangle\f$ — one integral over MY function, rank-3, so delivered as the house `Projector3` |

  **WHERE THE COEFFICIENTS WENT, and it answers the "why is `v` in the signature" objection.**  The
  object that IS a fit is the FITTER (basis + coefficients), and its face already reads right:
  `FitContraction<U>::Overlap(orb)`.  Below it, `DeltaScalarFitter` asks the basis for
  \f$\langle i|f_a|j\rangle\f$ and contracts its own \f$c\f$ into it — `Overlap3C` then contract, the
  *same two calls the Gaussian fitter already makes*, so Liskov substitutability is now literal rather
  than by analogy.  No coefficient vector appears in a basis signature.  Symmetrically, the density
  contracts its own \f$D\f$ through the forward direction of the same object, so \f$\Phi\f$ became
  private and increment 1's *"exactly one of them travels"* is settled the other way: NEITHER travels,
  because the integral is the thing that crosses.

  ⚠ **`Projector3` grew a THIRD contraction, and this is the one genuine wrinkle.**  A density holds
  \f$D\f$ in two forms — full, or the thin pivoted-Cholesky factor \f$L\f$ — and both are contractions of
  the same integral, so `applyRawFactored(L)` joins `applyRaw(D)`.  WHICH form the caller has, and whether
  the rank pays, stayed density-side with the factor and its memo (`FactoredRho`); only the GEMM moved.
  Free side-benefit: the ρ forward and the H adjoint now sit in one file, next to the BLAS-dispatch
  constraint they already shared and had documented twice.

  ⚠ **COVERAGE OF THE MOVED CODE, instrumented and counted rather than assumed.**  The factored branch
  fires exactly ONCE in the whole integration suite, and it is the `double` (real TRIM) instantiation:
  **`ForwardFactoredT<dcmplx>` is never exercised**.  Not a regression — the same body was previously a
  `FactoredRho<PeriodicIrrepCD<dcmplx>>` template that no ENABLED test instantiates through (the MnO
  recipes that would are `DISABLED`) — but it is a real hole under code that just moved, and a green
  suite should not be read as covering it.

  **✅ AND TWO MORE MEMBERS LEFT THE δ FACE, same day (user).**
  - **`Overlap3C` was never δ-specific** — every fit basis has \f$\langle\chi_i|f_a|\chi_j\rangle\f$, and
    `Orbital_DFT_IBS::MakeOverlap3C` already serves a Gaussian AND a plane-wave fit basis, so it was never a
    molecular-vs-lattice split.  It is now on **`FIT_SF_ABS<T>` with the default the user specified,
    `return orb.Overlap3C(*this)`**: Gaussian and PW inherit their existing machinery untouched (same cached
    tensor, verified bit-identical), δ overrides because it needs none of it.  Both other fitters now call
    `fitBasis->Overlap3C(orb)`, so all three representations are one line, and the genuine
    "is this a DFT-capable orbital basis?" cross-cast happens ONCE, in that default, instead of per fitter.
    ⚠ *Same-scalar only*, and it is a type limit rather than an oversight: the orbital-side tensor is keyed
    by the FIT scalar (`Projector3<TFit>`) while a δ table is keyed by the ORBITAL block's, so the mixed
    real-TRIM-on-complex-fit case (3c-3) cannot share one signature.  It stays as an extra δ overload (real
    arithmetic for a real block) and a `NarrowExact` at the raster route's call site.
    ⚠ *Module-cycle note for whoever touches it:* the default's BODY cannot live in the interface —
    `qchem.BasisSet.Orbital_DFT_IBS` imports `qchem.BasisSet.Fit_IBS`, so importing it back would close a
    cycle.  It is a free template `OrbitalOverlap3C<T>` declared in the interface, defined and EXPLICITLY
    INSTANTIATED (two scalars) in the implementation unit, which may import it.
  - **`SymmetrizeSpin` did not need to differ by representation either** (user's question, and the answer
    is no).  It is now on `FIT_SF_ABS<T>` with the default `{Symmetrize(rho); Symmetrize(m);}` — which is
    **bit-identical to the branch δ already ran** whenever the run carried no σ tags.  A Gaussian basis gets
    two no-ops, a raster gets two star-averages, and δ's override shrinks to the one case that genuinely
    differs: the Shubnikov projection, where a Flip op maps \f$\rho_\uparrow\f$ onto
    \f$\rho_\downarrow\f$ so the pair does NOT separate.  The old justification — "δ is the only
    representation a polarized run can use" — was about which representation gets SELECTED by the
    Hamiltonian's `VxcFit::Auto` policy, and a policy fact has no business shaping a basis face.

  ⇒ **`SiteIntegrals` then read as another public quadrature service** (user, 2026-08-23).  Half-fair: its
  argument IS an expansion over this basis's own functions
  (\f$\mu_A=\sum_{g\in A}c_g\langle\delta_g|1\rangle\f$ — a partitioned `Charge`, not a third-party
  integrand), so it is not the \f$\langle i|f|j\rangle\f$-for-arbitrary-\f$f\f$ smell.  What is foreign
  is the SITE structure: partitioning functions into atomic basins is not a fit-basis concept, it is the
  mesh's, and it was here only because the mesh is private.  Cured below.
  ⚠ Naming debt, unchanged: `applyRaw`/`applyRawAdjoint`/`applyRawFactored` are exactly the
  implementation vocabulary the house convention rejects.  Left alone because they are internal and shared
  with the GPW and plane-wave lineages.

  **✅ AND `SiteIntegrals` IS GONE TOO — BY INJECTION, NOT A GETTER (user ruling, 2026-08-23).**  The survey
  that preceded it: `grep` finds exactly ONE shareable handle to a mesh in the whole tree,
  `FitQuadrature::mesh`, and its sole owner was `DeltaFit_IBS`.  Everyone else who wants a unit-cell
  quadrature (the KB projector assembly, `PP_Local`/`PP_NonLocal`, `Fit_IBS`) builds one by value and
  either drops it per call or copies it into itself.  So the atom PARTITION — a general-purpose object with
  nothing to do with \f$v_{xc}\f$ or \f$\rho\f$ fitting — existed at run time in exactly one place: a fit
  basis.  That, and not any argument about the operation, is why the observable was stuck on that face.
  - **The cure is that the CREATOR hands its creation to both collaborators.**  `CreateVxcFitBasisSet` (the
    one factory call that already makes every representation decision) gained an optional out-parameter
    carrying the `shared_ptr<const qcMesh::Mesh>` it built; `Ham_PW_DFT::BuildTerms` passes it straight to
    `MakeVxcTerms` → `MakeXCQuadrature` → `XC_SinglesQuadrature`.  The moments are taken THERE, where
    \f$\rho_\sigma\f$ is already cached, through `qcMesh::SiteIntegrals`.
  - **Why not a `Mesh()` getter on the fit basis:** same face cleanup, but it re-opens the escape the
    2026-08-22 arc spent itself closing, and once it exists anything may use it for anything.  A creator
    injecting downward and an object publishing its own private state are different acts.
  - **Cost, stated rather than hidden:** an invariant became a convention.  \f$\rho\f$ is indexed by the
    δ basis's FUNCTIONS and the site blocks by mesh POINTS; they agreed *by construction* while one object
    owned both.  They now agree because it is literally the same immutable `Mesh` handed to both — pinned
    by a size assert, which is weaker than the type system was.
  - A raster (PlaneWave) fit basis leaves the out-parameter null, so no Becke mesh is built on a route that
    would not use it, and `SiteMoments` answers empty exactly as before.
  - **VERIFIED END TO END, not just green:** a null partition would make the moments silently empty and the
    suite would still pass, so it was instrumented — `PolarizedRunKeepsItsSpin` reports partition non-null,
    `npts==f.size()`, `NSites()==1`, and prints `[site moments] 0:5.00462`.  The same probe is what turned
    up **R1.0d** above.

  **✅ AND THEN `FIT_SF_Delta` WENT ENTIRELY (user: "remove FIT_SF_Delta completely… isOrtho can move up
  to the imp layer").**  Its last member was *"I can serve a REAL orbital block"* — a CAPABILITY, not an
  identity — so it became **`Integrals_Overlap3C<U>`**, the basis-side mirror of
  `Fitting::FitContraction<U>` and split for the same ISP reason: a fit basis's own scalar and the scalar
  of the block it contracts against are INDEPENDENT axes that only look like one until 3c-3 makes them
  differ.  `FIT_SF_ABS<T>` supplies it for its own scalar; `DeltaFit_IBS` additionally declares
  `Integrals_Overlap3C<double>`, which is the whole of what a δ representation adds.  `isOrtho()` moved
  down to the concrete class — "my metric is diagonal" is a fact about an implementation, not a type other
  code needs to name.

  **TWO CONSEQUENCES I DID NOT ANTICIPATE, both improvements:**
  - **The fitter Factory stopped asking "what IS it".**  It selected on
    `dynamic_pointer_cast<cFIT_SF_Delta>`; with that type gone the check HAD to become what the basis can
    DO — raster-backed ⇒ FFT projection + aliased grid lookup, otherwise the general orthogonal fit.  Both
    candidates are orthogonal, so the METRIC cannot separate them; the TRANSFORM face can, and is the
    property actually being used.  This was the one genuine "what is it" branch flagged in the metric-axis
    spec above as the thing a metric split must NOT launder — it is simply gone instead.
  - **A shared helper replaced what would have been two copies.**  The fitter (adjoint) and a density
    block (forward) both need "the 3-centre face for THIS block's scalar", in different libraries: one
    `Overlap3CFace<U>(basis)` — `if constexpr` same-scalar conversion, reference cross-cast otherwise.
    Compile-time dispatch on the scalar, not a run-time type switch.

  **✅ AND Fit_IBS.C IS NOW INTERFACES ONLY (user, 2026-08-23: "FitQuadrature, enum class VxcFit,
  OrthogonalFit do not belong in the interface definition file").**  Correct — no fit-basis face names any
  of them.  338 lines, declaring `Integrals_Overlap3C`, `FIT_CD_ABS`, `FieldEvaluator`, `FIT_CD_NonOrtho`,
  `FIT_SF_ABS`, `FIT_SF_NonOrtho`, `Fit_IBS` and their aliases, and nothing else.

  | moved out | to | why there |
  |---|---|---|
  | `FitQuadrature`, `enum class VxcFit` | **`qchem.BasisSet.Fit_Types`** | the factory VOCABULARY: `VxcFit` goes IN (which representation?), `FitQuadrature` comes OUT (the finished quadrature).  A LEAF module — it needs the mesh and the fold and nothing else — so it sits BELOW both `Orbital_DFT_IBS` and `tBasisSet`, the two that declare `CreateXCQuadrature`, and neither has to import the fit faces to name its own return type. |
  | `OrthogonalFit`, `Overlap3CFace` | **`qchem.BasisSet.FitOperations`** | free ALGORITHMS written in terms of a face, not part of any face's contract.  Each exists only because two consumers in DIFFERENT libraries need the same few lines (qcFitting's scalar fitter and qcChargeDensity's density). |

  **✅ AND `FieldEvaluator<T>::EvalField` WENT — the whole chain was dead (user, 2026-08-24).**
  The user's reasoning, which is the general form of the objection: *"any entity that has access to the fit
  expansion coefficients c also has access to the basis functions f(r), so it can do the field evaluation
  without going through this peculiar interface.  Why is it peculiar?  Because we are exposing the fit
  coefficients — instant bad smell."*
  ⇒ And the deeper reason the alternative is not available either: pulling \f$f_a(r)\f$ UP to the fitter is
  the PRE-2026-08-22 code, and it is what forced \c IrrepBasisSet to promise \c op(r) that an atom block
  (fake radial) and a \f$\delta\f$ basis (a distribution) could not honestly give.  **The motive was never
  abstract purity: `DeltaFit_IBS::op(r)` would be a very expensive and totally useless function** (user).
  So there were exactly three options — pass \a c down (exposes it), pull \f$f_a(r)\f$ up (resurrects the
  fake promise), or DELETE — and only the third exposes nothing.

  **MEASURED BEFORE REMOVING, because the static trace had already been wrong twice this session.**  Both
  bodies instrumented, full suite run: **0 hits across 758 tests**, 0 on a direct `M_DFT.*:M_HF.*` run, and
  nothing in `pybind/`, `viz-demo/` or `CLIapps/`.  The static picture agreed and explains why:
  `EvalField` ← `FitImpBase::operator()` (the ONLY direct caller) ← `FittedCD_Imp::operator()` (already
  annotated *"// No UT coverage"*) ← nobody.  `FittedVee` is the sole production holder of a `FittedCD` and
  uses `DoFit` / `GetRepulsion` / `GetSelfRepulsion` only.
  ⚠ Keep straight which `EvalField`: `G_FieldEvaluator::EvalField(ΔG_Map, r)` is a DIFFERENT function on a
  different class and is genuinely live — it is how a plane-wave fit evaluates its own inverse transform.

  | removed | |
  |---|---|
  | `FieldEvaluator<T>` + its derivations on `FIT_CD_NonOrtho` / `FIT_SF_NonOrtho` | the face whose whole job was to carry \a c across |
  | `Fit_IBS::EvalField` / `EvalFieldGradient` | |
  | `FitImpBase::operator()` / `Gradient`, and its `ScalarFunction<double>` base | the sole caller |
  | `FunctionFitter_Density<T> : ScalarFunction<double>` | the AO fitter could no longer keep the promise, so the FIELD became a capability: `OrthoFunctionFitter` (plane-wave) derives `ScalarFunction` itself — it genuinely can, by inverse transform over its own \f$\{G\}\f$ — and a consumer cross-casts.  Exactly what `FunctionFitter_Scalar` did earlier for the same reason. |
  | `FittedCD : ScalarFunction<double>` + its two forwarding lines | |

  ONE call site changed (`PlaneWaveDFTUT.C:1252`, now a cross-cast — the same one its scalar-fitter sibling
  200 lines earlier already makes).  758/758, Si gate bit-unmoved.
  ⚠ **The capability is cheap to restore and should NOT come back in this shape.**  `Fit_IBS` still derives
  `Evaluatable_IBS<double>` privately (it needs its own \c op(r) for `Norm()` and `Overlap(f)`), so ten
  lines bring it back — but as an operation returning the RESIDUAL, not one taking coefficients.  See
  **R1.0f** for what a GUI actually needs and why \f$v_{xc}(r)\f$ is not a fit-basis question.

  ⇒ **NEXT, and the user's own reading of the obstacle.**  `FIT_SF_ABS<T>::Overlap3C` should take
  `const Orbital_DFT_IBS<T>&` rather than `Orbital_1E_IBS<T>` — the caller must already hold a DFT-capable
  basis, so the `dynamic_cast` inside the default disappears (compile-time over run-time).  **Attempted and
  it does not compile**, for a structural reason worth recording: a cross-module forward declaration makes
  a DISTINCT entity —
  *"error: declaration of 'Orbital_DFT_IBS' in module qchem.BasisSet.Orbital_DFT_IBS follows declaration in
  module qchem.BasisSet.Fit_IBS"* — and a real import is a cycle, because
  `Orbital_DFT_IBS::Overlap3C(const FIT_SF_ABS<TFit>&)` already names the fit face.  ⇒ **The cure is not to
  merge the modules but to move the DEFAULT further toward the implementation** (user): make
  `FIT_SF_ABS<T>::Overlap3C` pure, and put the delegating `orb.Overlap3C(*this)` body in a mixin declared
  where BOTH faces are visible, which the concrete Gaussian and plane-wave fit bases derive.  Then only
  that mixin names `Orbital_DFT_IBS`, no cycle, no cast — and it also retires the explicit-instantiation
  trick the current default needs, which is a good sign it is the right shape.

  **✅ INCREMENTS 8–10 LANDED 2026-08-24 — THE FITTING BOUNDARY.**  The four items the user set on
  2026-08-24 (doc/OpenWork.md), three of them closed outright.  758/758 after each; the pinned Si
  two-route gate printed **−7.115067665** / **−7.115059008** at 11/11 iterations after each — bit-unmoved.

  **8 (items 0 + 1) — the fitter's contraction face is keyed on BOTH scalars and takes the DFT block.**
  `FitContraction<U>::Overlap(const Orbital_1E_IBS<U>*)` became
  `FitContraction<U,TFit>::Overlap(const Orbital_DFT_IBS<U,TFit>&)`.  One change, two reasons:
  - *Both axes*, because doc/RealComplexPlan.md 3c-3 makes them differ in production — a real TRIM block on
    a periodic run contracts against the run's COMPLEX fit basis, and `<U>` alone means `<U,U>`, a face no
    fitter in the tree declares.  Exactly the correction `BasisSet::Integrals_Overlap3C` took in increment
    7, so the two mirror faces now say the same thing in the same words, with **no default for `TFit` on
    either** (`TFit==U` is a fact about the molecular lineage, not a property of the face).
  - *The DFT block*, because what the contraction needs from the block is its 3-centre tier.  THREE casts
    disappeared outright (`FunctionFitterImp`, `OrthoNormalScalarFitter`, `DeltaScalarFitter::Contract`)
    and the survivors moved UP into `XC_SinglesQuadrature::MatrixT` / `XC_PairQuadrature::MatrixT` — the
    right place, since `cDynamic_HT::MakeMatrix` is deliberately method-neutral (a Kinetic term must not
    see DFT types) while an XC strategy is DFT-specific by construction.
  - **The molecular terms needed no new cast at all**: `FittedVxc` / `FittedVcorrPol` already held
    `odftbs_t` and were merely widening it back to the 1E base to make the call.  Their pointer casts
    became reference casts on the way past, which retired one unchecked deref.

  **9 (item 3) — `Symmetrize` / `SymmetrizeSpin` are off the fit face.**  User: *"these look like functions
  that need access to Fit_IBS's internal (none of your business) integration mesh.  They do not belong in a
  fit basis set interface."*  Right twice over: star-averaging a coefficient vector over a crystal point
  group is not a FITTING question, and the only thing the basis contributed was geometry it happened to
  own.  Same cure as `SiteIntegrals`, and the two halves died differently because the operation genuinely
  differs by representation:

  | | what it is | where it went |
  |---|---|---|
  | δ | one line over the free `SymmetrizeValues` / `SymmetrizeValuesSigned`; the basis supplies only the `Fold` | **INJECTED**, exactly as increment 6 injected the mesh: `CreateVxcFitBasisSet`'s out-parameter widened from `shared_ptr<const Mesh>` to the whole `FitQuadrature` (the fold and the Shubnikov tags are its SIBLING FIELDS), and `XC_SinglesQuadrature` calls the free algorithms itself |
  | plane wave | a voxel permutation plus an FFT shift-theorem glide | **`G_RasterTransform`** (default no-op), where every other raster-only operation already lives.  `PlaneWaveFit_IBS`'s body is unchanged and now overrides that face; its one caller, `Composite_Fourier::GetRhoOnGrid`, cross-casts for the raster — the same "I want more" ask that route is built on, and it is holding a raster array by construction |

  `FIT_SF_ABS<T>` is now four integrals and `isOrtho()`, with **no operation on it that needs a mesh**.
  The grey/magnetic/free three-way stayed ONE method rather than becoming a call-site branch, and its
  no-tags path is still `Symmetrize(rho); Symmetrize(m)` — bit-identical to the deleted base default.
  ⚠ **The MnO Shubnikov gate caught the omission on the first run**, which is the argument for injection
  over a getter working as intended: a probe that built the basis over one bundle and the strategy over
  another failed loudly.  The fix is a probe helper (`SinglesEngineOver`) that hands ONE bundle to both
  collaborators, exactly as the production factory does.

  **10 (item 2) — `Overlap3CFace` deleted, `OrthogonalFit` moved to qcFitting.**
  - `Overlap3CFace` was one reference cross-cast wrapped in an `if constexpr`, saving two lines at three
    call sites — and after increments 7/8 the `if constexpr` arm bought nothing: the same-scalar case is a
    derived-to-virtual-base cast that cannot fail, and every call site runs once per Fock block, not in a
    loop.  Each site now spells its own cast, beside the `Orbital_DFT_IBS` cast it already made.
  - `OrthogonalFit` moved `src/BasisSet/FitOperations.C` → `src/Fitting/FitOperations.C`, module
    `qchem.BasisSet.FitOperations` → `qchem.Fitting.FitOperations`; qcBasisSet loses the module entirely.
    User: *"OrthogonalFit belongs in the qcFitting library.  Client code needs to use that framework, not
    dodge around it."*  Performing a fit is what that library is for, and a fit algorithm sitting in the
    BASIS library was an invitation to use it instead of a fitter.

  ⇒ **THE ONE THING LEFT, AND IT NEEDS A RULING — the second half of that sentence (R1.0g).**  Two callers
  still use the algorithm INSTEAD of a `FunctionFitter_Scalar`: `XC_SinglesQuadrature`'s matrix-free
  branches and `tDM_CD::DM_RhoAtPoints`'s default.  Both want the COEFFICIENT VECTOR, and the fitter face
  deliberately has no accessor for one — adding `Coefficients()` is precisely the smell deleted in
  increment 6.  The spec's own ⚠ concedes the coefficients must reach the term regardless: \f$v_{xc}\f$ is
  pointwise NONLINEAR, so for δ they ARE \f$\rho(r_g)\f$ and the functional is applied to them directly.
  **So the open question is not how to hide them but WHAT TYPE carries them** — and a bare `rvec_t` is the
  weakest possible answer.  The shape of the answer is visible: the object that would carry them is *an
  EXPANSION over a fit basis*, offering `Integrate()` (\f$c\cdot\langle f_a|1\rangle\f$), `Map(functional)`
  and `Contract(orb)` — which is what `XC_Quadrature` and `FunctionFitter_Scalar` jointly are today, split
  across two objects with the coefficient array passed between them as a raw vector.  Introducing it
  reshapes the term/quadrature/fitter triangle, so it is a design ruling, not a refactor.  Being in
  qcFitting at least makes the bypass legible as what it is: framework algorithm, no framework object.

  **✅ RULED AND BUILT SAME DAY (increment 11) — and the answer was cheaper than the question.**  The user:
  *"DeltaScalarFitter needs to hold a shallow copy of \f$\Phi_{gi}=\chi_i(r_g)\f$."*  Two facts made it
  work with no new machinery: **`Projector3` already IS that shallow copy** (three closures over the fit
  basis's address-stable table — ~100 bytes against tens of MB, so the table itself never leaves), and the
  **fitter has RUN lifetime** where the density does not (`XC_SinglesQuadrature`'s ctor builds it once;
  `SCFIterator.C:343` builds a fresh density every iteration, which is why a density-side Φ cache was never
  an option).  New face `Fitting::ScalarProjector` — `Overlap3C(block)` ×2 for a matrix-backed density to
  contract its own \f$D\f$/\f$L\f$ into, `Project(field)` for a matrix-free one, `NumCoefficients()` so a
  composite can size its accumulator — and `tDM_CD::DM_RhoAtPoints` takes it instead of the fit basis.
  `OrthogonalFit` is down to TWO callers, both inside qcFitting: **no client bypass remains.**

  **And the coefficient question dissolved rather than being solved**: a PROJECTION returns coefficients, a
  FIT stores them, and \f$\rho\f$ only ever needed the first — so `Project()` returns `rvec_t` and no
  `Coefficients()` getter is needed anywhere.  It *must* not store, either: the same fitter holds
  \f$v_{xc}\f$'s fit and the Fock build interleaves with the energy path, so one buffer for both would
  clobber exactly the way `itsXCMix` once did.  The "expansion type" sketched above is therefore NOT needed
  for this, and is parked unless something else asks for it.

  ⇒ The real payoff is not the saved map lookup: the \f$\rho\f$ FORWARD used to go density→basis while the
  \f$H_{xc}\f$ ADJOINT went fitter→basis — two callers reaching one tensor, agreeing by convention.  Both
  now come off ONE held handle, so the adjoint pairing is a class invariant, which is the argument
  `XC_Quadrature`'s own header makes for its two faces one level up.
  ⚠ **Lifetime is now load-bearing**: those closures capture the producing basis, so a STORED copy dangles
  if the basis dies first.  Safe by construction — a fitter co-owns its fit basis by `shared_ptr` — but
  that is an invariant to keep, not an accident of ordering.

  **✅ AND INCREMENT 12 FINISHED THE FACE (user, 2026-08-25).**  ⚠ **NAMING NOTE FOR READERS OF OLDER
  RECORDS: `tDM_CD::DM_RhoAtPoints` is now `tDM_CD::ProjectOnto`.**  Every mention of the old name
  elsewhere in this file and in `doc/OpenWork.md` is a historical record and was true when written.
  - **The rename.**  It returns COEFFICIENTS over the fit basis, not values at points, so "AtPoints" was
    the last of the point vocabulary the 2026-08-23 pass took off everything else — and it made the method
    read like a competing entry point to `operator()(rvec3vec_t)` when it is the identity-carrying version
    of it.  `ProjectOnto(p)` is true of every implementation and pairs with `ScalarProjector::Project`.
    ⛔ `EvalFromOrbitalValues` was considered and rejected on two counts: `Values` is the exact word
    increment 3 removed from the interfaces, and it names the mechanism of ONE branch — the base default
    used no orbital values at all.
  - **The dead default went, and it was measured dead first**: 0 calls across all 760 tests, with a control
    probe on the `IrrepCD` override firing to prove the instrument worked.  Structural, not accidental —
    this face is `tDM_CD`, so arriving here means HAVING a matrix, and a matrix-free density is not a
    `tDM_CD` and never enters (its caller asks `ScalarProjector::Project` directly).  The field route was a
    default on the WRONG FACE.
  - ⇒ **And making it pure immediately earned its keep**: `tPolarized_CDImp` was leaning on that default
    and became abstract, i.e. a whole density family had no stated answer to "project yourself".  It has
    one now (\f$\rho_\uparrow+\rho_\downarrow\f$, the projection being linear).  Nothing called it, so
    it was not a bug — but it was a class that would have silently taken the pointwise field route the
    first time anything did.

  ⇒ **AND THE QUESTION THAT SETTLED THE WHOLE SHAPE (user, 2026-08-25), worth keeping because it is
  general:** why can an `operator()(rvec3vec_t)` override not just BE the fast path?  Because that face
  receives COORDINATES, and \f$\chi_i(r_g)\f$ is not reachable from them — so an override could only
  RECOMPUTE the table, which costs what the pointwise sweep costs.  **Measured** at 4000 points × n=16:
  `op(points)` 77 ms, tensor cold (table build + GEMM) 81 ms, tensor **warm 0.4 ms**.  The ~200–500× is the
  CACHE, not the contraction — materialising Φ is not cheaper than evaluating it, only cheaper the second
  time.  And the override could not cache either: a density is a FRESH OBJECT every SCF iteration
  (`TOrbitalsImp::GetChargeDensity` news one) and is asked for ρ once per iteration, so a density-side
  table has a zero hit rate.  A bare point list is universal precisely BECAUSE it carries no identity, and
  uncacheable for the same reason; passing the projector is how the point set gets one.

- **R1.0e ✅ THE FILE SPLIT IS DONE 2026-09-08; THE SCOPE QUESTION IT EXPOSED IS THE OPEN PART.**
  (Original: USER, 2026-08-23, *"the enormous PWTerms TU is going to need a massive refactoring cleanup
  eventually"*, restated 2026-09-08: *"src/Hamiltonian/Internal/PWTerms.C is huge, again doing too many
  things"*.)

  **What landed:** the 1213-line `Internal/Imp/PWTerms.C` is now FIVE implementation units of the same
  module — the interface-plus-many-Imp-units shape `Internal/Terms.C` has always had:

  | unit | what it holds | lines |
  |---|---|---|
  | `Imp/PWTerms_PP.C` | `Ven_PP_Short` / `_Long` / `_NonLocal` + the G=0 alignment | 212 |
  | `Imp/PWTerms_Hartree.C` | `Vee_Hartree` | 148 |
  | `Imp/PWTerms_XC.C` | `Vxc_Quadrature`, `Vxc_QuadraturePol`, `Vcorr_QuadraturePol`, `MakeVxcTerms` | 173 |
  | `Imp/XCQuadrature_Pair.C` | the PAIR strategy | 319 |
  | `Imp/XCQuadrature_Singles.C` | the SINGLES strategy | 477 |

  Helpers used by more than one unit (`NarrowExact`, `SampledField`) moved to the interface's
  **non-exported** section — module linkage is exactly their scope, and duplicating them per unit would
  have been an ODR trap dressed as tidiness.  Verified line-for-line (6 scaffolding lines differ, no code);
  827/827.

  ★★★ **AND THE SPLIT MADE THE REAL PROBLEM VISIBLE AS A FILE BOUNDARY.**  The user's definition of this
  library (2026-09-08) is the measuring stick:

  > *"At a very high level Hamiltonian is just: charge density in, use orbital basis and fitted functions
  > to evaluate all integrals, spit out \f$H_{ij}(\rho)\f$ and \f$E(\rho)\f$ for each term."*

  By that definition a term is (physics) + (ask the basis for an integral) + (contract with ρ).  **The
  two `XCQuadrature_*` units — 796 of the 1329 implementation lines, 60% — are none of those.**  They are a
  SAMPLING ENGINE, and everything they own is (integration grid) × (fit basis) business, which is
  `doc/Pins.md` pin 2's axis pair, not Hamiltonian business:

  | what `XC_Quadrature` owns today | why it is not Hamiltonian work | where it belongs |
  |---|---|---|
  | ρ sampling + per-density-serial caches (scalar, the {↑,↓} pair, and a separate DM-mix buffer) | a caching policy over a grid, not a term | `qcFitting` |
  | RAW-vs-BALL **route latching** (`LatchRoute`, `SampleOne`, `itsRhoIsRaw`) | which fit route to take is a FITTING decision | `qcFitting` |
  | the **Φ table** contraction and the projector it goes through | χ(r) caching — basis business | `qcBasisSet` / `qcFitting` |
  | **orbit star-averaging** of ρ and the (ρ,m) pair with Shubnikov spin tags (`Symmetrize`, `SymmetrizeSpin`) | crystal symmetry machinery | `qcSymmetry` / `qcFitting` |
  | **site-partitioned moments** + their console/report emission (`SiteMoments`, `PartitionedMoments`, `EmitSiteMoments`) | an OBSERVABLE and its reporting | wherever site observables live — not in a term |
  | **raster geometry** (`Raster()`, voxel counts, the uniform quadrature rule) | grid management | `qcBasisSet` |
  | `bool& ReportGridCharge()` | process-wide MUTABLE state in a library | `theRunPolicy()`, which already carries every other run-scoped switch |

  ✅ **THE MOVE IS LEGAL — CHECKED, NOT ASSUMED (2026-09-08):** nothing in `qcFitting`, `qcBasisSet`,
  `qcMesh`, `qcSymmetry` or `qcChargeDensity` imports `qchem.Hamiltonian.*`, so there is no cycle; and
  `qcFitting` already links `qcBasisSet qcSymmetry qcStructure qcMesh`, i.e. everything the engine touches.
  `qchem.Hamiltonian.Types` — the only Hamiltonian-side thing the engine names — is a pure typedef module
  over `BasisSet::Orbital_1E_IBS<T>` with no Hamiltonian dependency of its own.
  ★ `src/Fitting/Imp/FunctionFitter.C:67` already says its capability question *"is the same question
  `MakeXCQuadrature` asks"* — the duplication was noticed from the other side a year before this.

  ▶ **THE INCREMENT, when it is scheduled:** promote the two `XCQuadrature_*` units + the `XC_Quadrature`
  hierarchy out of `Internal/PWTerms.C` into their own module in `qcFitting` (`qchem.Fitting.XCQuadrature`),
  leaving `PWTerms_XC.C`'s three term classes — which are genuinely thin, and genuinely "functional in,
  \f$H_{ij}\f$ and \f$E\f$ out" — behind.  ⚠ **Do it AFTER `LatticeSum1E`'s ISP split** (item 5 on the
  user's Stage-B list): both refactors touch the collocate/integrate-back seam from opposite sides, and
  landing them together makes neither reviewable.

  ⚠ **STILL OPEN from the original item:** the atom block still derives `Evaluatable_IBS`, and its `op(r)`
  is still the FAKE RADIAL — the promise kept in form and broken in substance, contained (not cured) by
  `ImplicitAngular_IBS`.  That is the remaining scope of step (1): convert `PP_Local::CalculateMatrix` and
  `PP_NonLocal`'s explicit-angular branch to ask the basis for the integral, then the atom density/orbital
  paths, after which the derivation drops and nothing evaluates an atom block from outside.  The relocation
  makes that a LOCAL change per call site instead of a tree-wide one, because the face no longer forces the
  promise on everybody.

- **R1.0f ⚠ THE GUI STILL NEEDS \f$v_{xc}(r)\f$ AND \f$\rho_{DM}-\rho_{fit}\f$, AND THE δ ROUTE HAS NO
  WAY TO GIVE THEM — USER, 2026-08-24.**  Recorded when `FieldEvaluator::EvalField` was deleted (below):
  nothing in the tree evaluates a fitted field any more, but a plotting front end will want both.
  **The two are NOT the same problem, and only one of them is hard:**
  - **\f$v_{xc}(r)\f$ is not a fit-basis question at all.**  It is a POINTWISE functional of the density,
    \f$v_{xc}=f(\rho(r))\f$, and \f$\rho\f$ is evaluatable on any grid the GUI likes.  So the honest route
    is *evaluate \f$\rho\f$ there, apply the functional* — exact, representation-independent, and cheaper
    than interpolating a fit.  Going through a fit expansion would be a worse answer to an easier question.
  - **\f$\rho_{DM}-\rho_{fit}\f$ is representation-dependent, and for \f$\delta\f$ it is IDENTICALLY
    ZERO.**  The δ fit is the identity (\f$c_g=\rho(r_g)\f$), so the residual vanishes at every point where
    the fit is defined and is undefined everywhere else.  The diagnostic measures *how badly does my fit
    basis represent this*, and a δ basis has no fit error by construction — it has QUADRATURE error instead.
    So the δ-side quantity a user actually wants is a DIFFERENT pair: \f$\rho_{DM}\f$ against the
    band-limited \f$\tilde\rho\f$ the Hartree term integrates, which is where the representation error
    actually lives.  Name it for that; do not spell it \f$\rho-\rho_{fit}\f$.
  - ⇒ **If a genuine "evaluate my δ expansion at an arbitrary \a r" is ever needed, it is an INTERPOLATION
    problem on a scattered atom-centred mesh** — expensive, approximate, and owed a name that says so.  It
    must not come back disguised as `op(r)`, which is exactly the trap this whole arc removed.

- **R1.0c ⚠ COVERAGE HOLE: the COMPLEX factored-\f$\rho\f$ contraction is exercised by NO enabled test.**
  **MEASURED 2026-08-23** (instrumented `DeltaFit_IBS::ForwardFactoredT`, counted over the whole
  integration suite): the low-rank route fires **exactly once**, and it is the **`double`** (real TRIM)
  instantiation.  `ForwardFactoredT<dcmplx>` never runs.
  *Not a regression* — the same body was previously inline in `FactoredRho<PeriodicIrrepCD<dcmplx>>`, which
  no ENABLED test instantiates through (the MnO recipes that would are `DISABLED`) — but it is live
  production code on the default route (`QCHEM_DM_LOWRANK` unset ⇒ pivoted Cholesky is ON), so a green
  758/758 must not be read as covering it.
  **Why it matters more than a normal gap:** a wrong factored \f$\rho\f$ is not loud.  It is
  \f$\sum_m|[\Phi L]_{gm}|^2\f$, non-negative by construction, so a mistake shows up as a plausible-looking
  density and a mistuned SCF, not a crash — exactly the failure mode the `FactoredRho` staleness assert was
  added for.
  **Cheapest fix: a UNIT gate, not an SCF one.** For one Bloch block, build \f$D=LL^\dagger\f$ from a
  random thin \f$L\f$ and assert `applyRawFactored(L)` == `applyRaw(D)` to ~1e-12, on both scalars.  That
  pins the identity the whole route rests on without needing a converging magnetic cell, and it would also
  cover the `double` path with something sharper than "MnO happens to run".
  *(A second, weaker option is to enable one small polarized periodic SCF; that costs suite time and still
  only covers whichever scalar that cell happens to use.)*

- **R1.0d ⛔ DEFECT, FOUND 2026-08-23: an IMPOSED-symmetry Becke mesh SILENTLY LOSES ITS SITE BLOCKS, so
  every per-site integrated observable vanishes on exactly the runs that want it.**
  `UnitCell::CreateIntegrationMesh(mp, ops)` (Structure/Imp/UnitCell.C) builds the site-adapted Becke mesh
  via `MakePeriodicBeckeMesh`, which DOES record one block per atom (`MeshBuilder::BeginSite`).  The
  orbit-consistency filter that follows — the pass that drops eps-borderline points whose orbit partners
  were tail-dropped — then rebuilds the whole mesh into a **fresh `qcMesh::MeshBuilder` and never calls
  `BeginSite`**.  The result has `NSites()==0`.
  **MEASURED, three runs:**

  | run | ops | atoms | `NSites()` | moments |
  |---|---|---|---|---|
  | `GPW_SCF.PolarizedRunKeepsItsSpin` | free | 1 | **1** | `[site moments] 0:5.00462  net=5.00462` |
  | `GPW_SCF.O2TripletInBoxMatchesFinite` | 16 (Shubnikov) | 2 | **0** | silent |
  | `GPW_SCF.ImposedShubnikovHoldsAFMThroughSCF_Mn2Box` | 96 (Shubnikov) | 2 | **0** | silent |

  **Why it stayed invisible:** the consumer contract is *"EMPTY when the mesh has no site blocks — ask, do
  not assume"*, which is right for a uniform grid and indistinguishable from this.  So an imposed run
  reports no moments and looks like it correctly had none.  The MnO AFM campaign is imposed by
  construction, i.e. the atomic moments the *"integrated observables, not point probes"* rule exists to
  produce are precisely the ones that are missing.
  ⚠ **NOT a simple `BeginSite()` insertion.**  The filter iterates ORBITS and appends their members, so its
  output is orbit-major; site blocks require contiguous per-site runs.  The fix is to compute a KEEP mask
  from the orbit test and then re-emit in the ORIGINAL mesh order (site-major), calling `BeginSite` at each
  original boundary — which also CHANGES THE POINT ORDER relative to today, hence the summation order of
  every integral over that mesh.  Algebraically identical, last-bits different: it needs its own
  measurement against the imposed-run pins, not a drive-by.

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
- **R1.9 ✅ DONE `3882938e`. Molecular `BasisSetID()` streamed its SEPARATORS as hex addresses.**  **→ doc/CleanupHistory.md**
### R2 — mechanical hygiene

- **R2.21 ✅ FOUND AND FIXED `9da2e825` (2026-08-24). `blaze::conj` IS A NO-OP ON A COMPLEX SCALAR** — it silently broke the Cholesky factor; use `blazem::conjs`.  **→ doc/CleanupHistory.md**
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
- **R2.18 ✅ NAMES DONE `86c5b24d`.  The ENCAPSULATION half is DELIBERATELY LEFT OPEN — user ruling, kept here because it is the open part:**  **→ doc/CleanupHistory.md**
    - **ENCAPSULATION (public vs protected `Make`): low priority, DELIBERATELY LEFT OPEN.**  *"All the
      MakeXXX() functions were originally protected.  For DFT the 3C versions (MakeOverlap3C,
      MakeRepulsion3C) still are.  It seemed like these were purely internal functions ... but that turned
      out to be incorrect in some cases.  Anyway I have no strong policy on this (encapsulation level)
      right now.  Maybe the right policy will emerge as we refactor.  My intuition says that it is a low
      priority decision."*
- **R2.19 ✅ DONE `86c5b24d`. `FittedVxcPol` copied a matrix its child already owned.**  **→ doc/CleanupHistory.md**
- **R2.20 ✅ DONE `7c80e71e` (2026-08-17).**  The four oracle helpers moved into an
  `export namespace qchem` block of `qchem.PeriodicTable` (the user's suggested home); TestUtils.C had
  nothing test-only left and is DELETED with both FILE_SET wirings — scfrun imports no test module.
  721/721 green.  **→ doc/CleanupHistory.md**

- **R2.21 ✅ DONE 2026-08-17 (concurrent-cleanup session), BOTH halves** — the occupation state stores each block's reference under ITS OWN scalar; `RealBlockFillView` deleted.  ONE remainder, kept here because it is the open part:  **→ doc/CleanupHistory.md**
    - **Still open (smaller now):** a `route` ctor argument would let a caller FORCE ball (an A/B instrument
      — today only the `GPW_XCROUTE` env var reports the route, it cannot select it).  That needs a
      capability question on the neutral `Band_FT_IBS` face so the factory can ask without a concrete cast.
      Worth doing when someone actually wants the A/B.
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
- **R2.15 ✅ COMPLETE — the Lebedev DEFAULT FLIP LANDED 2026-08-17, degree-gated** (`nAngular` → degree-typed angular interface; at degree 29 Lebedev delivers 302 directions against GaussLegendre's 450).  **→ doc/CleanupHistory.md**
- **R2.15 (original text) `nAngular` → degree-typed angular interface.**  `nAngular` is a COUNT for Lebedev but a
  DEGREE for GL/EM (and the imposed site-adapted builder consumes it as the degree) — the dual
  semantics BLOCKS flipping the free-run Becke default to the measured-equal Leb-302 (67% of
  GL-29's directions).  Fix: `angularDegree` + per-scheme count resolution; the default flip rides
  along.  (The warn/auto-resolve ergonomics half of this item stays in D5.)

## VERIFY

### V1 — interface-design questions

- **V1.1 ✅ DONE `d49db261`+`4b37221d`+`2bfb83b4` (2026-08-16) — the merge landed; the answer was YES.**
  `ERI3`+`G_ERI3` → **`Projector3<T>`** (one struct, realizations inside — user ruling); the
  `MakeOverlap(f(G))` bridge was already production-dead and died with its using-decl cascade;
  `Band_FT_IBS` deleted — the lattice lineage IS `Orbital_DFT_IBS<dcmplx,dcmplx>`, and
  `Orbital_DFT_IBS<double,dcmplx>` (real TRIM block, complex fit basis) is now a live spelling.
  The metric worry needed no new machinery — (i)'s two axes had already discharged it.
  **→ doc/CleanupHistory.md** (full record + the three commit summaries).
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
- **V1.3 ✅ MECHANISM DONE `72fecf8d`** (both ε-adapters deleted, via `GetEMatrix`).  ⚠ **STILL OPEN: the QUADRATURE-TERM face** — the second list in the original item (`FittedEpsXc`/`FittedVxc` simplification).  **→ doc/CleanupHistory.md**
- **V1.4 ✅ DONE `80fc2ae8`. `DM_RhoAtPoints` Phi key → Irrep (USER RULING 2026-08-05).**.  **→ doc/CleanupHistory.md**
- **V1.5 ✅ DONE `f18a6ee9`+`9ebaebdb` (2026-08-16) — FOUR faces, not three, and a reporting redesign fell
  out first.**  The §K blocker had dissolved piecemeal (grid one-owner landed with #7; V1.1 removed the last
  orbital-flavored method), so the split was executable.  `EmitGridReport` did NOT move onto a face — user
  ruling: providers self-report; every created grid announces at construction, ROLE-labeled by its factory
  (`report::EmitAt` made idempotent so dedup is run-scoped in the REPORT, killing the `static const void*`
  latch in PWTerms and the raw-`cout`-beside-the-report bug).  Then:
  `G_FieldEvaluator` (evaluate) / `G_Quadrature` (FFT engine, via `GriddedScalarFitter::Grid()`) /
  `G_StructureFactor` (seed) / `G_SpectralFilter` (mixer).  **→ doc/CleanupHistory.md.**
  Still open from the old bullet, now standalone: XC_GridEngine's `Lattice_3D::Fold` + dcmplx dependency
  bars molecular reuse of the quadrature engine.
- **V1.6 ✅ DONE `2d0f6982`. `tDM_CD::Accumulate*` — face split, NOT pure-virtual.**  Now `tHF_System_CD` +
  `tHF_Pair_CD`, real path only via `conditional_t`; the complex leaf declares NOTHING (a CRTP mixin, after
  the user pressed that empty bodies are still the interface failing to segregate).  The NDEBUG hazard is
  closed: `Vee`/`Vxc` THROW where they used to build a zeroed J in silence.  **→ doc/CleanupHistory.md**
  *(original)*  Verified override matrix: every
  concrete family relies on the default for exactly 2 of the 4 (IrrepCD lacks `*All`;
  Composite/Polarized lack `*Both`) — pure-virtual just forces 6 new asserting stubs.  The `*Both`
  pair is an internal Composite↔leaf collaboration protocol (only called from
  Imp/CompositeCD.C:48,59) — no business on the public face.  Shape: a whole-system face and a
  pair-partner leaf face (the `tSpinResolved_CD` cross-cast idiom).  NDEBUG hazard on record:
  these void assert-only bodies are silent NO-OPS in Release — a bare IrrepCD through
  `Vee::AccumulateAll` yields a zeroed J and a silently wrong Fock.
- **V1.7 ✅ DONE `2d0f6982`. All NINE denials gone** — three periodic-only CRTP mixins; the mechanism
  (`FourierDensityBase<T>`) was already right there, the families just re-declared outside it.
  **→ doc/CleanupHistory.md**  *(original)*  **The periodic trio (`GetFourierDensity`/`GetRhoOnGrid`/`GetRepulsion3C`) — 9 asserting
  stubs, the largest LSP block in qcChargeDensity.**  Re-declared + NA-asserted on
  Polarized/Composite/IrrepCD for BOTH T (Imp/ChargeDensity.C:97,120,138; Imp/CompositeCD.C:197,
  226,249; Imp/IrrepCD.C:267,286,303).  The correct mechanism ALREADY EXISTS in the same file —
  `FourierDensityBase<T>` gives dcmplx the capability and double an empty base — but the derived
  classes re-declare outside it and assert.  Fix: declarations live only on the dcmplx side
  (if-constexpr-guarded definitions or a `tPeriodic_CD` mixin).
- **V1.8 ✅ DONE `2d0f6982`. The cast EVAPORATED with V1.6**, exactly as the user predicted — and because the
  face is operation-named (`CompleteDirectPair`), not a block accessor.  **→ doc/CleanupHistory.md**
  *(original)*  **`IrrepCD`↔`IrrepCD` concrete same-class casts in the hot path** (Imp/IrrepCD.C:84,98,
  218,227: `Accumulate*Both`/`MixIn`/`GetChangeFrom` take abstract `tDM_CD&` and narrow to the
  concrete leaf to touch `itsDensityMatrix`; the in-file comment names "the IrrepCD↔IrrepCD
  idiom").  Abstract→concrete, the pattern the project rule forbids; also makes MixIn
  unimplementable for any future leaf.  Wants a double-dispatch primitive or an abstract
  density-block face.  (Design with V1.6 — same seam.)
- **V1.9 ✅ DONE `38a1ebd6`. `Structure`→concrete-`UnitCell` down-casts in 4 libraries**.  **→ doc/CleanupHistory.md**
- **V1.10 ✅ DONE `2d0f6982`. Both casts gone.**  The SALC one dissolved into V1.31's `WholeSystemFock_IBS`
  face — its three primitives ARE the steps the cast open-coded; the DHF one became
  `Orbital_RKB_Pair::MakeDirectAgainstL`.  **→ doc/CleanupHistory.md**  *(original)*  **Two abstract→CONCRETE basis casts in src/** — Imp/SymmetryAdapted_IBS.C:109,118
  (Orbital_HF_IBS* → concrete SymmetryAdapted_IBS, solely to reach `itsO`) and
  Internal/Imp/Orbital_DHF_IBS.C:89,109 (Orbital_ERI4_IBS& → Orbital_RKB_HF_IBS_Imp&).  Both are
  "give me your private state" reaches — promote the needed answer to an abstract question on the
  face.  (Unit-test exemption does not apply; these are src/.)
  *(NOT in this list, and deliberately so: `Orbital_ERI4_IBS::Substrate` added by R1.7 is an
  abstract→ABSTRACT cross-cast — the sanctioned direction — and it THROWS naming both bases.)*
- **V1.10b ✅ DONE (see LANDED). Mixer.  **→ doc/CleanupHistory.md**
- **V1.11 ✅ DONE `43bbebad`+`0c818835`+`841eadf2`+`092d1da8`+`2398dd07` (2026-08-17) — the occupation
  seam, five bit-identical increments.**  `OccupationPolicy<T>` in qcElConfig decides every fill
  (`DecideBlockFill` → the two-axis `BlockFill` spec; `HeldOccupationPolicy` is the direct minimiser's
  sibling); the WF carries NO occupation state; ONE `TOrbitals::Fill` replaced the five `TakeElectrons*`
  virtuals; the EC mode bools became the `ReservoirPartition`; the D11 seed-fill hazard closed
  structurally.  **DAG lesson**: qcElConfig→qcOrbitals is a linker cycle — the `OrbitalView<T>` DIP face
  (owned below, implemented above) is the CLAUDE.md inversion example verbatim.
  **→ doc/CleanupHistory.md** (full record).  **With this, ALL SIX doc/RealComplexPlan.md §7
  prerequisites are DONE.**  *Residue (user catch, 2026-08-17): the landed policy is still a
  mode-flag-configured CONCRETE, not D1's abstract interface — the Policy/State split that finishes it is
  filed as **R2.21** (liked, deferred).*
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
  **POST-MERGE (checked 2026-08-17: its three faces are src/WaveFunction + the IrrepBasisSet face +
  SCFIterator — all in the real-TRIM working set).**
- **V1.15 `tBasisSet<T>::Create*FitBasisSet` defaults** — the generic body hard-codes
  `Orbital_DFT_IBS<double>` regardless of T (only the explicit dcmplx specializations save it) and
  derefs an unguarded iterator (a 1E/HF-only basis ⇒ null ⇒ UB in Release).  Also
  `CreateXCQuadrature`'s default body is byte-identical in two places (Imp/BasisSet.C:42 ==
  Band_FT_IBS.C:53).  Decide: pure, or a documented "this basis does not fit" contract; hoist the
  shared default.  (Interacts with the V1.1 CreateXCQuadrature hoist.)
- **V1.16 ✅ DONE `de0292cb`. `ProjectedDensity_AO::GetRepulsion3C` asserting default** — the metric is now
  two refinement faces (`CoulombMetric_ProjectedDensity` / `OverlapMetric_ProjectedDensity`); the base keeps
  only what the FITTER needs.  **The compiler found more than the item predicted:** `tPolarized_CD` and
  `tComposite_CD` were cross-casting to the PLAIN AO face and then calling the Coulomb-only
  `GetRepulsion3C` — the implicit pairing in action, since any `ProjectedDensity_AO` satisfied that cast and
  the poison fired later.  Both now ask for the capability they use.  **→ doc/CleanupHistory.md**
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
- **V1.19 ✅ VISITOR + THROWS DONE 2026-08-17; bit-identical, 734/734.**  ⚠ **ONE DELIBERATE REMAINDER**: the seed's flip-group sub-cell duplication — removing it needs a per-SITE form-factor overload on the basis face, which the item itself weighs against the pseudo-wall pin.  That block is the seed's ONE remaining concrete-`Atom` consumer.  **→ doc/CleanupHistory.md**
- **V1.20 `SymmetryAdapted_IBS` → Internal: decide the library-family rule first.**  Verified: all
  consumers inside the src/BasisSet/ tree BUT in a different CMake target (qcMolecule_BS).  Moving
  to `.Internal.` creates a cross-LIBRARY internal import under the CLAUDE.md rule.  Precedent the
  rule is already bent: `qchem.BasisSet.Internal.GMap` imported by src/Fitting + src/ChargeDensity.
  Decide: is the qcBasisSet* family "one library" for Internal purposes?
- **V1.21 `BandStructure.C`: promote or demote.**  Confirmed test-only (sole import =
  tests/BandStructureUT.C:13); worse, tests/PlaneWaveUT.C:59 defines its OWN local `SolveBands`
  instead of importing.  Either promote (band plots are on the viz roadmap) or demote into the
  test tree; either way kill the duplicate.  **POST-MERGE (checked 2026-08-17: the file is
  src/BasisSet/Lattice_3D/ — the real-TRIM working set).**
- **V1.22 `MakePeriodicBeckeMesh` ε-tail drops vs orbit consistency (W2c find).**  The builder's
  borderline drop decisions (`<eps` screens + `w>0` keep) are per-point and bit-sensitive, so the
  site-adapted caller post-filters orbit-incomplete points (`CreateSiteAdaptedBeckeMesh`).
  Cleaner: make the drop decision ONCE per representative (angular dir × radial shell) and apply
  it to the whole atom orbit inside the builder — removes the second fold pass + the filter.
  **2026-08-17 (concurrent-cleanup session): DEFERRED deliberately** — bit-sensitive on the imposed
  mesh path the running MnO real-TRIM campaign depends on; wrong week to move it.  **New evidence
  while measuring R2.15:** the drop asymmetry is not only an imposed-run problem — the MnO
  seed-mirror gate's orphan check caught the FREE builder keeping a point whose AFM translation
  partner was tail-dropped (w·ρ=0.04, degree-11 Lebedev ⟨111⟩ into a neighbour core).  A
  per-representative drop rule would have made that configuration impossible by construction, so
  the item's case is STRONGER than when filed.
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
- **V1.25 ✅ DONE `e52c00eb`. NOT "minor trims" — there was a LIVE LEAK.**  `GetChargeDensity()` now returns
  `unique_ptr` along the whole chain, so `AtomCalculation::TotalCharge()` stops leaking a composite per call.
  **Three more owning-raw-pointer sites the item never mentioned** turned up on the way: `tPolarized_CDImp`
  had `tSpinDensity`'s exact double-delete shape, and `tComposite_CD::Insert` advertised a raw pointer while
  wrapping it in a `unique_ptr`.  **→ doc/CleanupHistory.md**
  *(original analysis follows)*
  - **`GetChargeDensity()` is a `Get*` that ALLOCATES.**  `tCompositeWF::GetChargeDensity(Spin)` does
    `new tComposite_CD<T>(...)` and inserts every irrep block, on EVERY call; `TOrbitalsImp`, `tIrrepWF`,
    `tPolarizedWF` and `tUnPolarizedWF` all forward or do the same.  The name says accessor; the body is a
    factory.
  - **`AtomCalculation::TotalCharge()` LEAKS it** (Imp/AtomCalculation.C:222): builds a whole composite
    density, reads one number off it, drops the pointer.  `SolidCalculation` happens to get it right
    (`itsImp->cd.reset(cd)`), and the tests mostly remember `delete cd` — so the contract is honoured by
    vigilance, not by the type.  **Cross-ref the RAM question**: a leak per call is exactly the class of
    thing to rule out before profiling anything.
  - **The item's OTHER claim was backwards.**  It said `tSpinDensity` holds "two raw NON-OWNING `tDM_CD*`
    with unmanaged lifetime".  It DOES own them — its dtor deletes both — but declares no copy ctor or
    assignment, so the implicit copy is a DOUBLE-DELETE (rule of three).  And it is only correct today
    because its caller feeds it two freshly-`new`ed densities, i.e. it depends on the accessor's hidden
    factory behaviour.  Three fragilities holding each other up.
  - **Fix (all at once, user 2026-08-12):** return `std::unique_ptr<tDM_CD<T>>` from the accessor chain, so
    ownership is in the TYPE and the compiler finds every caller; `tSpinDensity` holds `unique_ptr`s (copy
    then implicitly deleted, dtor disappears).  Matches CLAUDE.md's rule that a raw `new` should go into a
    smart pointer within a few lines.
  - Original text: non-const `Polarized_CD::GetChargeDensity(Spin)` overload has no external consumer
    (removable); `tSpinDensity` holds two raw `tDM_CD*`.

- **V1.26 ✅ COMPLETE (reconciled 2026-08-17)** — every deliverable landed across the Uniform-vs-Becke selector work.  Three things STILL OPEN after that landing, kept here because they are the open part:  **→ doc/CleanupHistory.md**
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

- **V1.28 ✅ RESOLVED BY THE SHUBNIKOV CAMPAIGN S1–S4** (reconciled 2026-08-17).  **→ doc/CleanupHistory.md**
- **V1.29 ✅ RECONCILED 2026-08-17 — the question is ANSWERED and the dependency DISCHARGED.**  ⚠ What remains UNBUILT is the Hessian-based magnetic-order DISCOVERY loop (Davidson on the stability matrix, swept per non-symmetric irrep block) — deliberately not built for MnO, whose AFM-II order is known and imposed; it earns its keep on materials whose order is unknown.  Depends on V1.28.  **→ doc/CleanupHistory.md**
- **V1.30 ✅ APPEARS ALREADY FIXED — item was STALE** (verified against the tree 2026-08-10).  **→ doc/CleanupHistory.md**
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
  needs a suite sweep since every pinned GPW anchor re-seeds.  **POST-MERGE + the bit-moving batch
  (checked 2026-08-17: both defaults live in TRIM-owned files — GPW_SCF_UT.C:314 and
  SolidCalculation.C:104 — and "every anchor re-seeds" is the V1.22/§K class).**
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
  Raising it moves grids, hence anchors: measure first (D8), same instrument as V2.4.  **POST-MERGE +
  the bit-moving batch (checked 2026-08-17: lives in the GPW evaluator = TRIM working set, and it is
  anchor-moving by its own last sentence).**

- **V2.6 ✅ CLOSED (reconciled 2026-08-17)** — the Becke recipe ladder is fully banked (nRadial=40 right; the angular flip to 17 REFUTED by Al FCC).  Records FOUR refuted guesses; read it before changing either default.  **→ doc/CleanupHistory.md**
- **V2.6a ⛔ ATTEMPTED AND REJECTED 2026-08-07 — flip `angularDegree` 29 → 17.**  Made the one-line change,.  **→ doc/CleanupHistory.md**
- **V2.7 ✅ DONE 2026-08-17.**  `RadialResolutionRatio(mp, alphaMax)` + the `kRadialRatioFloor=3.0` warning from `ResolveXCMesh`.  ⚠ The floor stays INSULATOR-FITTED per its own doc note: re-calibrate on {Si, Mn-atom, Al, MnO} before promoting it beyond a diagnostic.  **→ doc/CleanupHistory.md**
### V3 — repro / campaign bugs (Spin-SAD, 2026-08-04)

- **V3.1 ✅ CLOSED 2026-08-17 — NO LONGER REPRODUCES (dissolved in the interim; V1.11 the likely
  cure).**  Regression anchors added: `ValenceBasisGen.SodiumSeedDensitySpinResolved` +
  `Slater_Low/A_HF_P.Energy/Z1` (UHF H, exact −0.5).  **→ doc/CleanupHistory.md**
- **V3.2 ✅ CLOSED 2026-08-17 — NO LONGER REPRODUCES (same verification pass).**  Anchors:
  `ValenceBasisGen.SodiumSeedDensity{UnpolarizedWithPolarizationShell,SpinResolvedWithPolarizationShell}`.
  **→ doc/CleanupHistory.md**

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

- **D13 Single-species `HGH_LocalPotential` / `HGH_SeparablePotential` silently IGNORE their `int Z`
  argument** (found 2026-08-14 — this exact foot-gun manufactured the retracted "V_long sharp-field
  defect": the `DiffuseDVlongSharpFieldOracle` gate passed `gthMn.local` over an Mn+O cell, so
  production integrated a Mn-q7 field at the O site while the oracle used the true O q6, and Mn's
  empty C list made β=0 so the custom-G-ball path was never even entered).  The per-species classes
  take `Z` only to satisfy the `LocalPotential` face; a caller sweeping `FormFactorLong(a->itsZ,…)`
  over a multi-species Structure with a single-species object gets one species' physics at every
  site, with no diagnostic.  Production is safe (`BuildMultiSpeciesLocal` wraps the router), but
  hand-written gates/probes are not.  Candidate fixes, in the compile-time-over-runtime spirit
  ([[feedback_compile_time_over_runtime]]): store the species Z the object was built for (GetGTH
  knows it) and `assert(Z==itsZ)` in every Z-taking method; or drop `Z` from the single-species
  concrete API entirely and let ONLY `MultiSpecies_*` model the Z-dispatching face.

- **V1.32 Rename the FINITE density leaf `IrrepCD<T>` → non-template `FiniteIrrepCD` (user, 2026-08-17,
  out of the RealComplexPlan 3c-2b split).**  After lineage-as-class, the finite leaf has exactly ONE
  instantiation — `<double>` — and the factory's `if constexpr` guard makes a finite-complex density
  UNREPRESENTABLE (it throws before instantiating).  So the template parameter is vestigial and the
  name no longer says the load-bearing thing: *Finite* (molecules/atoms) is the identity, not the
  scalar.  De-templating also lets its conditional bases (`ProjectedDensityBase<T>`,
  `IrrepHF_PairBase<T,...>`) collapse to the plain double faces.  SMALL: touchers are
  `IrrepCD_Factory` (the double branch), the explicit instantiation + `template <>` member
  definitions in Internal/Imp/IrrepCD.C, the `IrrepCD_HFPair<IrrepCD<double>>` CRTP/friend spellings,
  and tests/Version.C.  Consider the same question for `PeriodicIrrepCD<T>` and DECLINE it there:
  that leaf genuinely has both scalars (the real TRIM block vs general k), so its T is load-bearing.

- **The `OverlapMatrix` static-field integrate-back screens with `CollocMemo::D` — safe only by build
  ORDERING (found 2026-08-18, out of the cross-run pollution hunt).**  `GPW_Evaluator::OverlapMatrix`
  passes `itsCollocMemo->D` (the last collocated density) as the `IntegratePotential` density screen,
  and `MakeLocalPP` routes through it — a STATIC field, whose matrix must not depend on any iteration's
  D.  Today this never bites: each block's static terms are built before that block's first Hartree
  collocation (`tHamiltonianImp::GetMatrix` sums statics before dynamics), so the memo is still invalid
  and the sweep runs unscreened.  But the guarantee is implicit: any flow that rebuilds a static PP
  matrix AFTER a collocation on the same evaluator (a new Structure::ID mid-process, a diagnostic
  probe, a future term reordering) would silently D-sparsify a static matrix — and the I2n cache would
  then serve that sparsified matrix process-wide.  Candidate fix: a `screen` parameter on
  `OverlapMatrix` so the STATIC callers (`MakeLocalPP*`) pass none explicitly, keeping the D-screen an
  opt-in of the per-iteration KS-field path only.  (The cross-run sibling of this hazard — the `Dscr`
  union riding the DBCache'd `Projector3` closures — was fixed 2026-08-18 by instance-scoping the GPW
  3C tensors in `tGPW_IBS`; see `GPW_SCF.CrossRunFirstRunAnomalyProbe`.)

- **V1.34 THE FITTER'S CONTRACTION FACE IS TEMPLATED BUT ONLY HALF-REALISED — an ISP hole with a
  `bad_cast` for a diagnostic** (found 2026-08-28, while decoupling the XC grid from polarization).
  `Fitting::FitContraction<U,TFit>` is templated on the BLOCK scalar and declares one method,
  `hmat_t<U> Overlap(const BasisSet::Orbital_DFT_IBS<U,TFit>&)`.  The ortho scalar fitter implements
  **only** the `<dcmplx,dcmplx>` face.  So `XC_PairQuadrature`'s ball-fit branch, asked for a REAL TRIM
  block, gets a `std::bad_cast` — and carried a hand-written throw saying *"a real TRIM block on the
  legacy ball-fit XC route is not wired"* with the reason recorded as **"nothing real-block reaches
  it"**, which was an observation about reachability, not about the interface.  It stopped being true
  the moment a polarized run could take that route.

  ⛔ **WHY THIS IS AN OOD ITEM AND NOT A MISSING OVERLOAD** (user, 2026-08-28: *"it sounds like we have
  an OOD design problem in the fitter interfaces"*).  A templated interface that exists for two scalars
  and is realised for one is a face that LIES: every caller must either know which instantiation is real
  — knowledge the abstraction exists to remove — or discover it as a `bad_cast` at runtime, in Release,
  mid-SCF.  That is the same failure the project already ruled on elsewhere: *give capabilities only to
  types that have them* (\c feedback_compile_time_over_runtime), and the house style of asking a face
  what it CAN do rather than what it IS.  The current shape supports neither: you cannot ask, and you
  cannot fail to compile.

  ⇒ **The ruling to take, not the patch.**  Three shapes, and picking between them is the item:
    1. **Realise the `<double,dcmplx>` face** — smallest, but it leaves the next unrealised combination
       to be discovered the same way.
    2. **Make the capability ASKABLE** — a `bool CanContract<U>()`-shaped question, or split the face so
       a fitter advertises exactly the scalars it serves.  Matches the `applyRaw`/`applyRawAdjoint`
       pattern the XC pair route already uses (an empty `std::function` IS the capability answer), and
       it is how `MakeXCQuadrature` already decides between its two strategies.
    3. **Make it a compile-time error** — the fitter only exposes the instantiations it defines, so a
       caller that needs another one does not link.  Strongest, and the project's stated preference
       (build-failure over runtime crash); needs the caller side to be scalar-generic in a way it may
       not be.
  ⚠ NOT urgent: the hole is currently unreachable again (the XC adjoint follows the lineage's
  capability, so the ball fit only ever sees complex blocks).  It is on this list because the NEXT
  route change will rediscover it, and because a half-realised templated face is a design defect
  whether or not anything is standing on it today.
  ★ RELATED: [`project_functionfitter_isp_split`] already split `FunctionFitter` into Scalar/Density
  faces on exactly this kind of argument, so this is the same axis, one level down.

### V1.33 — THE BasisSet TAXONOMY IS THE WRONG AXIS (user, 2026-08-20)

`src/BasisSet/{Atom, Molecule, Lattice_3D}` classifies by PHYSICAL SYSTEM, but what the directories
actually contain is classified by BASIS KIND.  The user's proposed axis:

> `Radial` / `Polarized{Cartesian|Spherical}` / `<LocalPeriodic? = GPW>` / `PW`

**The evidence that the current axis has already failed**, found while designing the Φ Bloch point-sum
seam (`doc/OpenWork.md` Step 3):

- `qchem.UnitCell` is imported at **five** sites INSIDE `BasisSet/Molecule/` — `PG_Cart/BasisSet.C`,
  `PG_Cart/Imp/IrrepBasisSet.C`, `Evaluators/PG_Cart_MnD/Evaluator.C`, `PG_Spherical/Imp/LatticeView.C`
  and `LatticeSum1E.C` itself.  A *molecule* has no unit cell.
- The face `Molecule::LatticeSum1E` names a LATTICE inside the MOLECULE namespace, and defines
  `cellphase_t` there.  Any new periodic capability (e.g. a Bloch point-value face) deepens it.
- `PG_Cart::IrrepBasisSet::operator()`'s own comment reads *"the PERIODIC caller (GPW_Evaluator::Eval)"*.

**Root cause (user):** *"Gaussian basis functions/sets are simply not a Molecule specific concept."*
`Molecule` here has come to mean MULTI-CENTRE (as against `Atom` = single-centre), which stopped being
true the moment GPW consumed the same basis periodically.

**SCOPE PRECISION:** the smell is `BasisSet::Molecule` ONLY.  `Symmetry::Molecule` (point groups, against
`Symmetry::Lattice_3D` space groups) is CORRECTLY named — a blanket rename would destroy a real
distinction.

**Two separable increments, both DEFERRED to their own session (user, 2026-08-20):**
1. **Move the periodic capability faces to the system-neutral level** `qchem::BasisSet::`.  Measured as
   dependency-FREE: `qcBasisSet` already links `qcStructure` (where `UnitCell` lives), already imports
   `qchem.Symmetry.Lattice_3D.Fold` (`Internal/GMap.C:19`), and already hosts `Band_DFT_IBS.C` — a
   periodic capability face at that level.  `qcMolecule_BS` depends on `qcBasisSet`, so moving a face UP
   is the existing dependency direction (no cycle).  Blast radius: 9 module importers, 12 files, 68
   textual uses, **no `pybind/` impact**.
2. **Re-cut the taxonomy itself** onto the basis-kind axis above.  Much larger: directories, namespaces,
   every `qchem.BasisSet.Molecule.*` MODULE NAME (see the module-rename dyndep hazard), the
   `qcMolecule_BS` target, the `.vscode` test globs — **and it breaks `pybind/qchem_bridge.cpp`** (2
   references), which is binding-owned: FLAG it, never fix it lib-side (CLAUDE.md).

**RULING (user, 2026-08-20) on the near-term cost of NOT doing this:** `LatticeSum1E` *"has evolved from a
simple three member lattice version of Make{Overlap,Kinetic,Nuclear} into a bit of a monster class … it is
already a mess, making it slightly incrementally messier at this point is not a big concern."*  So the Φ
point-value face is added to `LatticeSum1E` IN PLACE, and the ISP split of that class is deferred here
along with the taxonomy.  **`LatticeSum1E` therefore also wants an ISP review in its own right** — it now
carries collocation, integrate-back and grid machinery that are not one-electron integrals, so even its
NAME is stale.

## `Vxc_QuadraturePol` is dead code (2026-09-04)

`MakeVxcTerms` used to return the exchange/correlation PAIR for a polarized run; since the
`CompositeExFunctional` change it returns ONE term (`Vcorr_QuadraturePol` carrying the composite), so
`Vxc_QuadraturePol` — the exchange-only polarized term — has no constructor call anywhere.  Its
`SiteMoments` already moved to the surviving term.

⇒ Delete it, and while there consider RENAMING `Vcorr_QuadraturePol`: it is no longer "the correlation
half of a pair", it is the spin-native XC term (`Vxc_SpinNative` or similar).  Not done inline because the
rename touches three files and the measurement work wanted a small diff.

## R2.22 ✅ DONE `0210cfb9` (2026-09-06) — MnO annealed **83.2 s → 67.6 s (1.23×)**, `Etot` bit-identical

**Landed exactly as scoped** (step 1 only; the analysis that produced that scope is kept below).  One rule
now covers the iterator's three collaborators — *what it is HANDED it does not delete, what it MAKES it
holds in a `unique_ptr`*: the Hamiltonian is non-owning and the facade keeps ONE for the whole schedule;
the accelerator is non-owning and still replaced per stage by its owner (stale Pulay/DIIS, and the type
changes); the wave function became a `unique_ptr`; `~tSCFIterator` is `= default`.  Applied in all three
facades — the molecular pair keep rebuilding their (cheap) Hamiltonian per `Converge()`, they just own it
now.  The three misleading comments are gone.
- **Acceptance met on both halves:** `Etot=-61.40297529` (bit-identical to the banked value), and
  `setup: hamiltonian ctor` + `setup: becke mesh build` run **once** — the only `[x2]` buckets left are
  the per-stage SCF residues, which should be per stage.  Physics unchanged (m_stag 0.6667, 14+17 iters).
  814/814 green.
- **A pre-existing LEAK fell out of it:** on `RunGpw`'s non-annealed path the Hamiltonian had no owner at
  all once the iterator stopped deleting it — invisible while the iterator was silently cleaning up after
  a caller that never owned anything.
- ⚠ The Benchmark.md MnO row is now stale (83.2 s). **Deliberately NOT re-banked** — stale rows are
  re-banked ONCE at the end of Phase 1, not per increment (`8fbb6132`).

*(the item as filed, and the review that scoped it, follow)*

## R2.22 (original) — `tSCFIterator` DELETES A HAMILTONIAN IT DID NOT CREATE, and it costs 19% of an annealed run (2026-09-06)

`tSCFIterator<T>::~tSCFIterator()` does `delete itsHamiltonian;` on a raw `ham_t*` handed in by the
caller.  CLAUDE.md: *"Raw `new` ops are fine if the pointer quickly goes into a `unique_ptr` or
`shared_ptr` … As a result `delete` should be rare or non-existent."*  This is the opposite — the
Hamiltonian is constructed by the composition root (`SolidCalculation`, or the GPW test harness), handed
over as a bare pointer, and destroyed by an object that is one of its *users*.  `SolidCalculation::Imp`
even documents the smell in a comment: `ham = nullptr;  // owned by the iterator once handed over`.

**WHAT IT COSTS, MEASURED (`doc/ParallelAndOraclePlan.md` 1.1(b)).**  Because the previous stage's
iterator deletes the Hamiltonian when it dies, `SolidCalculation::BuildStage` MUST build a fresh one for
every anneal stage — so an N-stage schedule constructs N Hamiltonians.  The MnO recipe has two, and the
build is **15.5 s** each: **19% of an 83 s threaded run, for an object that is a pure function of
(structure, basis, species, functional, xcMesh, vxcFit)** — none of which change between stages.

⚠ **THE PHYSICS REASON DOES NOT APPLY TO THE HAMILTONIAN.**  The test harness's comment reads *"Fresh
Hamiltonian + accelerator per stage (the iterator OWNS + deletes them; a kT change must not carry stale
DIIS history across the re-seed)"* — and stale DIIS history is a property of the **accelerator**, which
genuinely must be fresh each stage.  The Hamiltonian is there because of the `delete`, nothing else.

⇒ **The fix**: `std::shared_ptr<ham_t>` on the iterator (or a non-owning raw pointer with the facade
owning it), and the `delete` goes.  Call sites: `SolidCalculation` ×3, `GPW_SCF_UT.C` ×5-ish, plus the
molecular `Calculation`/`AtomCalculation` path — all compiler-enforced, so the change is mechanical.
⚠ Before sharing one across stages, confirm the terms carry no per-stage state that must be dropped: the
memo keys are logical density SERIALS (not pointers) since the Dynamic_HT fix, which suggests continuing
the sequence across a stage boundary is fine, but that is an argument, not a measurement.
★ Acceptance: `Etot` bit-identical on the MnO annealed row, and the ledger showing
`setup: hamiltonian ctor [x1]` instead of `[x2]`.

---

### Scope correction (2026-09-06, on review): it is THREE deletes, THREE facades, and the Hamiltonian is a SYMPTOM

**(a) `~tSCFIterator` deletes three raw pointers, and one rule covers all of them.**

| member | origin | today | verdict |
|---|---|---|---|
| `itsHamiltonian` | **given** by the caller | `delete` | ❌ the item above |
| `itsAccelerator` | **given** by the caller | `delete` | ❌ the same defect, unstated |
| `itsWaveFunction` | **created** in the ctor (`WaveFunction::Factory`) | `delete` | ✓ right to free, wrong to be raw |

`itsMixer`, `itsOccPolicy` and `itsOccState` are already `unique_ptr`/by-value, so the class is half
modernised and the three survivors ARE the three deletes.  One rule finishes it — **what you are handed you
do not delete; what you make you hold in a `unique_ptr`** — and `~tSCFIterator` becomes `= default`.
Do the accelerator WITH the Hamiltonian: "a stage needs a fresh accelerator" is a lifetime policy for the
composition root, not a reason for a USER to own it.

**Two comments die with it, both currently false or apologetic:**
- `src/SCFIterator/Imp/SCFIterator.C`, three lines above the dtor: *"Recall that the wavefunction is not
  owned buy this."* — the WF is the one of the three that IS owned.  A comment a future reader would trust
  (the R1.9 lesson), plus the typo.
- `src/WaveFunction/Internal/Imp/CompositeWF.C:224`: `// delete itsAccelerator; NO!!!! SCFiterator deletes
  the accelerator.`  A comment that exists only to document the awkward ownership; it evaporates.

**(b) THREE facades rebuild, not one.**  `Calculation::Converge` (Imp/Calculation.C:195) and
`AtomCalculation::Converge` (Imp/AtomCalculation.C:159) both `delete itsScf` — *"releases the previous
Hamiltonian + accelerator"* — and rebuild.  So a second `Converge()` on a MOLECULE or an ATOM rebuilds its
Hamiltonian too.  Same fix, same rule; the call-site list is SolidCalculation ×3 **+ both molecular facades**.

### WHAT IS A "STAGE BOUNDARY"? (user asked, 2026-09-06) — ONE variable, and it is kT

A stage is `SCFStage = {SCFParams params, SCFAccelerators::Type accelerator}` (SolidCalculation.C:332), so
the TYPE says "any of ~15 SCF parameters may change per stage".  The tree says otherwise.  There is exactly
**one** schedule builder anywhere (`GPW_SCF_UT.C:4580`, the MnO harness — `SolidCalculation::Converge(
vector<SCFStage>)` is public API with no production or CLI caller), every stage starts from one shared
`base`, and precisely three fields are then touched:

| varies per stage | source | independent? |
|---|---|---|
| `SmearingkT` | `MNO_ANNEAL` | ✅ **the one real variable — this IS the anneal** |
| accelerator type | `MNO_ACC` | ❌ a RESPONSE to kT (banked recipe: `MNO_ANNEAL="5e-3,0"` with `MNO_ACC="Ladder,GDM"` — Ladder while smeared, GDM once cold) |
| `MOMSmearPenalty` (Λ) | `MNO_ANNEAL_PENALTY` | ❌ only MEANS anything when kT>0 — a companion to kT |
| `StopOnAccelExhausted` | `(i+1<kTs.size())` | ❌ pure schedule POSITION ("am I the last stage?") |
| every other `SCFParams` field | copied from `base` | — constant |

⇒ **a boundary is: same everything, new kT** (plus one companion and one positional flag).  `SCFStage` is
over-general for its single caller, and `MNO_ANNEAL`/`MNO_ACC` are two parallel env lists that must be the
same length (enforced by an `assert`) — the classic shape of one thing modelled as two.  **PARKED, not
filed as work**: per the cost ruling below this is a test-harness ergonomics wart, not a library defect.

### Is "smeared ⇒ Ladder, cold ⇒ GDM" a LAW? — no, at THREE levels (user, 2026-09-06)

**(1) OUR GDM is kT=0 by IMPLEMENTATION, and the restriction is structural.**  A direct-min leg fills
under `HeldOccupationPolicy`, whose `SmearingkT()` returns 0 **by override, not by configuration** — so
\f$-TS\f$ is identically zero and **a held leg minimises \f$E\f$, never \f$A=E-TS\f$**, which at kT>0 is
the wrong functional.  Measured, not theoretical: `E(t=0)` under the held fill sits **+14.5 Ha** off the
previous smeared iteration on MnO — enough that the line search rejected every `t` until the reference was
re-taken under the same held fill (the "convention shift", `GPW_GDMTRACE`, SCFIterator.C ~L500).

**(2) GDM IS NOT kT=0 BY DEFINITION.**  Geometric/geodesic direct minimisation is formulated on the
Grassmann/Stiefel manifold of occupied SUBSPACES — \f$D=CC^\dagger\f$, idempotent — which is an
integer-occupation object *by construction*.  That is the PARAMETERISATION, not the method.  The
finite-\f$T\f$ generalisation is standard: promote the occupations to variational parameters and minimise
the FREE ENERGY \f$A=E-TS\f$ over orbitals AND occupations (ensemble-DFT direct minimisation; Marzari,
Vanderbilt & Payne 1997 is the canonical reference — their cold smearing exists partly to keep the entropy
term well-behaved, and the occupation-space PRECONDITIONER is the hard part).

**(3) CP2K's OT IS IN EXACTLY OUR POSITION — scaffolded and switched off.**  Checked in their source
(`/home/janr/Code/cp2k`, 2026-09-06), because the folklore "OT cannot smear" is not what the code says:
- `qs_scf_loop_utils.F:274` passes `scf_control%smear` INTO the OT branch; `ot_scf_mini(mos, mo_derivs,
  smear, …)` takes it and `set_mo_occupation(mo_set, smear=smear)` runs after the minimisation.
- The variational occupation axis exists: `&OT`'s `ENERGIES` keyword sets `settings%do_ener`, and under it
  `qs_ot_scf.F:196` builds *"the derivative of the free energy with respect to the evals"* — `rot_mat_u`
  (rotate the subspace) beside `ener_gx = ∂A/∂ε` (vary the occupations), exactly the two-axis shape (2)
  describes.  A sibling keyword `OCCUPATION_PRECONDITIONER` sits next to it.
- **But `qs_ot_types.F:906` reads `! not yet fully implemented` / `CPASSERT(.NOT. settings%do_ener)`.**
  So OT ships as a FIXED-occupation minimiser with the finite-\f$T\f$ axis designed in and asserted off.
⚠ What their OT+`&SMEAR` combination does numerically is NOT established here — only what the source
wires.  (Our own CP2K decks avoid `&OT` anyway; see doc/Benchmark.md §5a.)

⇒ **DESIGN CONSEQUENCE — do not bake "direct minimiser ⇒ integer occupation" in anywhere.**  The thing to
keep general is not which accelerator runs at which kT (ruled below: it stays an explicit stage field).  It
is the CONFLATION inside `HeldOccupationPolicy`, which currently welds together two independent facts:
  (i) **hold the occupied block** — do not re-decide WHICH states are occupied mid-line-search (real, and
      the reason the item exists: a re-ranked fill makes \f$E(t)\f$ discontinuous in \f$t\f$); and
  (ii) **do not smear** — `SmearingkT()≡0`, hence \f$-TS\equiv0\f$.
A realistic finite-\f$T\f$ GDM/OT needs (i) WITHOUT (ii).  R2.21's two-axis policy already accommodates it
— occupancy {Integer, **Fermi**} × ranking {Bare, MOM}, plus the `HoldsStoredBlocks` bit — so the missing
concrete is simply **Fermi-occupancy-with-held-block**, which `HeldOccupationPolicy`'s own header already
names (*"a coupled leg that holds the block but keeps kT is a new sibling with a Fermi occupancy — a new
object, not a new bool"*).  **Nothing in R2.22 may narrow that.**

★ **AND THE PLACE kT=0 IS ACTUALLY HARDCODED IS NOT THE ACCELERATOR — IT IS THE OBJECTIVE ASSEMBLY.**
The direct-min leg builds its objective as
`E(t) = itsHamiltonian->GetTotalEnergy(cd_t) + itsOccPolicy->EntropyTerm()` (SCFIterator.C ~L512/L490):
the energy under the HELD (integer) fill, plus the \f$-TS\f$ of the last REAL fill.  At kT=0 that term is
zero and the expression is exact.  Under a held leg at kT>0 it is a CONSTANT offset across \f$t\f$, so it
cancels in the comparisons and the search still "works" — but the trial's OWN entropy never varies, so it
is not a free-energy minimisation.  **That line is what must change when the Fermi-held policy lands**;
record it here so the next reader does not conclude the accelerator was the blocker.

★ **RULING (unchanged by the above): keep the accelerator an EXPLICIT stage field.**  Deriving it from kT
would bake TODAY's accelerator inventory into the schedule API and be wrong the day the Fermi-held leg
exists.  (This reverses the tentative "collapse `SCFStage` to `{kT, Λ}`" floated in the review.)

### ⇒ PRIORITY (user ruling, 2026-09-06): sort by CONSTRUCTION COST, not by ownership purity

> *"A good reason to rebuild is Pulay history is no longer valid.  Rebuilding objects that are dirt cheap
> to construct is a small (ignorable) problem.  An object cheaply reconstructed for a bad reason like
> *just* changing the stored kT value, maybe just document and leave it for now.  We want focus on
> rebuilding objects that are expensive to construct (Hamiltonian)."*

| rebuilt per stage | cost to construct | reason it is rebuilt | verdict |
|---|---|---|---|
| **Hamiltonian** | **15.5 s/call** — the largest setup bucket in the run | ONLY because the iterator deletes it | ⛔ **FIX. This is the item.** 15.5 s × (N−1) = **19% of an 83 s run** |
| **SCFAccelerator** | cheap (a Factory call) | **GOOD reason** — stale Pulay/DIIS history across a re-seed, AND the type genuinely changes (Ladder→GDM) | ✅ **correct as is; keep rebuilding** |
| WaveFunction + iterator | cheap: all three residue buckets total **0.36 s of 121 s**, and the ledger is a partition (unbucketed 0.03 s).  The stage's first Fock is work that must happen for a new kT anyway | a BAD reason — the accelerator is baked into the WF's per-irrep children (`CompositeWF.C:213`), so a new accelerator forces a new WF forces a new iterator | 📝 **DOCUMENT AND LEAVE** (user).  Cheap object, bad reason: not worth unpicking |

⇒ **R2.22 IS STEP 1 ONLY.**  The earlier two-step sequencing is superseded: "persist the iterator" is
DEMOTED out of this item — it buys 0.36 s and costs a `ResetAccelerators(acc&)` walk plus a density-mixer
history audit.  Revisit only if one of `BuildStage`'s two carries (`AdoptMOMReference(prev…)`,
`AttachProbes()`) is itself found to cause a defect; both are cheap and both currently work.

**So the whole item is:**
1. The facade owns the Hamiltonian; the iterator holds it non-owning; **one Hamiltonian per run**.
2. Riding along for free in the same hunk (compiler-enforced, no behaviour change): the accelerator also
   becomes non-owning (the facade keeps replacing it per stage — that is correct), the WaveFunction becomes
   a `unique_ptr`, `~tSCFIterator` becomes `= default`, and the two misleading comments die.
3. Same fix in all THREE facades (see (b) above).

★ Acceptance: `Etot` bit-identical on the MnO annealed row, and the ledger showing
`setup: hamiltonian ctor [x1]` instead of `[x2]`.
⚠ Still owed before sharing one across stages: confirm no term holds per-stage state.  Memo keys are
density SERIALS (not pointers) since the Dynamic_HT fix, so continuing the sequence across a boundary
should be fine — an argument, not a measurement.

*(Doc correction found while reviewing: R2.18's sub-note "ONE STALE COMMENT LEFT BEHIND ON PURPOSE —
`SCFIterator.C:186` still says Vxc::CalcMatrix … sweep it when that list is released" is DONE BY DRIFT.
`grep -rn "Vxc::CalcMatrix" src/` finds nothing tree-wide.)*

## Setting `MeshParams::cellKind=Becke` alone gives a DEGREE-5 mesh (2026-09-07)

`qcMesh::BeckeXCParams()` is the recipe: nRadial 40, mhlAlpha 2, angularDegree **29**, and it honours
`GPW_BECKE_NR/ALPHA/L` for any argument passed `<0`.  `MeshParams`' own member defaults are **30, 1, 5**.
So a caller that writes `mp.cellKind = UnitCellKind::Becke` — the obvious thing, and what it looks like the
type invites — gets a **degree-5 angular mesh** and no warning.

MEASURED COST OF THE TRAP: it produced a 40 mHa imposed-vs-free discrepancy that read exactly like a
symmetry bug and consumed a diagnosis before the recipe was checked (`doc/SymmetryUpgradePlan.md`,
"SETTLED 2026-09-07").  At the real recipe the same comparison gives 0.046 mHa.

⇒ Options, cheapest first: (a) make `MeshParams`' Becke-relevant defaults MATCH `BeckeXCParams`' so the
two agree however the struct is reached; (b) have whoever resolves a Becke mesh reject an unresolved
recipe rather than silently integrating at degree 5; (c) make `cellKind` unsettable on its own, so asking
for Becke means asking for the recipe.  (c) is the compile-time answer and matches the project's
preference for build failure over a plausible wrong number — a degree-5 XC mesh is exactly a plausible
wrong number.

