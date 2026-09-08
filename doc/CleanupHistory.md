# Cleanup history — the closed record

Split out of `doc/CleanupCandidates.md` on 2026-08-09, when the working doc passed 2000 lines.

NOTHING WAS TRIMMED.  Closed items are moved here in FULL and leave a one-line stub in the worklist,
so cross-references ("see R2.8", "the V1.10b ruling") still resolve and the reader can see that an item
happened without reading why.

Keeping the full text is deliberate, and this session supplied the evidence for it:
  - V1.24(ii) described a GDM fallback bug on 2026-08-03 that the MnO campaign independently re-derived
    six days later, using this doc's own two listed reproducers -- and the entry ALSO recorded that its
    own prescribed fix ("hold position") would have been wrong.
  - V2.6 records FOUR refuted guesses.  Each one is the obvious assumption, and each would have shipped a
    wrong default if the refutation had been trimmed as "already decided".
  - R2.8, R2.2 and V1.27 are the same finding in three libraries (an accidental type correlation standing
    in for the real discriminator).  That pattern is only visible because the earlier two were kept.
A record of what was TRIED AND REJECTED is worth more than a record of what landed, and it is exactly what
gets lost first when a doc is trimmed for length.

---

## LANDED 2026-08-17 — R2.20 `7c80e71e` (concurrent-cleanup session): the oracle helpers out of the test module

`RelativeError` / `RelativeHF/DFT/DHFError` moved from `IntegrationTests/TestUtils.C` (module
`qchem.Unittests.TestUtils`) into an `export namespace qchem` block of `qchem.PeriodicTable` — beside the
Saito/NIST/Dirac tables they wrap, per the user's own suggestion ("could live in PeriodicTable").  TestUtils
had nothing genuinely test-only left, so the module is deleted along with both per-target FILE_SET wirings
(ITMain, scfrun); nine importers swapped.  Two details worth keeping: the helpers carry `<cmath>`+
`<iostream>` in PeriodicTable's own global module fragment (the R1.9 stream-literal lesson), and they are
std includes rather than a `qchem.Math` import so qcCommon keeps its math-free pin.  Note the module's other
exports are GLOBAL-namespace (pre-`qchem::` legacy) — the helpers went into the proper namespace, making
PeriodicTable a mixed-namespace interface until the legacy exports migrate.

*(original item, verbatim)*

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

## CLOSED 2026-08-17 — V3.1 + V3.2 NO LONGER REPRODUCE (concurrent-cleanup session, branch concurrent-cleanup)

Both were filed 2026-08-04 off the Spin-SAD campaign and assigned to the concurrent-cleanup session on
2026-08-17.  Verification (not repair): every reproduction attempt PASSES on origin/main @ `0dffe7d3`,
including configurations HARDER than the filed ones.  What was run:

- **V3.1 as filed (valgen `--spin` surface):** `GenerateSeedDensity` for Na q1 (Zion=1, electrons=1,
  s-window `EvenTemperedWindow(5, 0.03, 2.0)`, `spinResolved=true`) — the polarized pseudo-atom SCF with
  nUp=1, nDown=0.  Converges; charge=1, moment=1 both land.  (The filed one-liner
  `AtomCalculation(11, 0, {...})` is internally inconsistent — charge=0 gives 11 electrons, not the
  described doublet; the faithful facade call is `AtomCalculation(11, 10, {.pseudopotential=true,
  .valence=1, .pol=Polarized})`, which also passes, on both the per-l window AND the Medium accuracy
  pool basis.)
- **V3.1 beyond the filing:** all-electron `AtomCalculation(1, 0, {.type=Slater, .pol=Polarized})` — UHF
  hydrogen, the purest empty-minority-channel case — converges to −0.5 exactly (0.78 ppm vs Saito).
- **V3.2 as filed:** the same Na recipe carrying the unoccupied p window (`{1, EvenTemperedWindow(2,
  0.09, 0.3)}`), both unpolarized and `spinResolved=true`.  Both pass.

**No fix was written here; the defect dissolved between 2026-08-04 and 2026-08-17.**  The most plausible
cure is **V1.11** (the occupation seam, `43bbebad`..`2398dd07`, landed 2026-08-17 — five increments
restructuring exactly the WF/occupation state these failures lived in), with the DIIS `SVTolRel`
relative-prune (runs 31–34) the runner-up candidate for the "Invalid setup of symmetric matrix" symptom.
Per the R2.9(ii) lesson (re-derive the defect before picking from an item's menu), no archaeology was
spent pinning the exact commit: the class is now PINNED BY TESTS instead —

- `ValenceBasisGen.SodiumSeedDensitySpinResolved` (V3.1: nUp=1/nDown=0 through the valgen --spin path)
- `ValenceBasisGen.SodiumSeedDensityUnpolarizedWithPolarizationShell` +
  `SodiumSeedDensitySpinResolvedWithPolarizationShell` (V3.2: the riding-along unoccupied p window)
- `Slater_Low/A_HF_P.Energy/Z1` (the fully-polarized one-electron atom, all-electron; UHF H is EXACT at
  −0.5, so any future regression there is structural, not numerical — Z=1 added to the Slater_Low span)

If either symptom returns, the OLD prescription still stands in the original text below: fix the
empty-channel atom SCF rather than teaching valgen to synthesize the exact pair, and check
`cSCFAcceleratorDIIS::GetNProj` on an empty EC (the filed sibling failure).

*(original items, verbatim)*

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

## LANDED 2026-08-16 — V1.1 COMPLETE (the basis-face merge), 717/717 green

Three commits, each a full green sweep; user rulings of 2026-08-16 folded in: **the MnO session now
WAITS on this work** (real TRIM irreps at (0,0,0), (½,½,½)…), so the old DO-NOT-TOUCH on
`src/BasisSet/Lattice_3D/` was explicitly lifted ("DO-NOT-TOUCH is now please-TOUCH"); **one struct with
realizations inside** was chosen over an abstract tensor interface; **`Projector3`** chosen over `Map3`
(a Map is a keyed container in this code region — ΔG_Map — and this is an operator; "ERI" rejected as
inaccurate since the overlap-metric tensor has no repulsion integral in it).

- `d49db261` — **V1.1(ii) closed: `ERI3` + `G_ERI3` → `Projector3<T>`.**  One 3-centre tensor type with
  the realization a property of the PRODUCING BASIS, not the scalar type (the R2.13/V1.27 lesson made
  structural): `dense` (molecular vector-of-`smat_t`, RAM-hungry), `columns`+`kernel` (plane-wave delta
  support), `apply*` closures (GPW matrix-free).  New lean module `qchem.BasisSet.Internal.Projector3`
  (also home of `IVec3Less`/`ΔG_Map`; `GMap` keeps only the space-group symmetrization utilities and
  re-exports).  `Contract`/`ContractAdjoint` replace `ContractG_ERI3`/`ContractAdjointG_ERI3`.  The
  DB_Cache's two overload-resolved I3C `Get`s merged into one slot; the dcmplx cache's always-empty
  `ERI3` map and the separate `G_ERI3` map died together.  **Load-bearing detail:** the cached accessors
  now key `theCache<TFit>()`, not `theCache<T>()` — identical for every T==TFit instantiation, required
  for `<double,dcmplx>` TRIM blocks.
- `4b37221d` — **V1.1(iii) first half: the `MakeOverlap(f(G))` field bridge retired.**  It was already
  production-dead — Hartree/XC assemble ⟨i|f|j⟩ via `ContractAdjoint` on the SAME tensor the density
  contracted (fit-grid-consistent, doc/GPWPlan §0e step 2); only unit tests still called it, and they now
  use the evaluators' own `OverlapMatrix(f)`.  The load-bearing `using Integrals_Overlap::MakeOverlap`
  un-hiding decl AND the concrete classes' two-base overload-set merges (PlaneWave_IBS/GPW_IBS) became
  noise and were deleted — exactly the cascade the item predicted.
- `2bfb83b4` — **The merge itself: `Band_FT_IBS` deleted; the lattice lineage IS
  `Orbital_DFT_IBS<dcmplx,dcmplx>`.**  `CreateXCQuadrature` hoisted with its neutral default (identical
  to `tBasisSet<T>`'s whole-set default — (iii) second half, which paid only here, as reclassified);
  `EPW_Orbital_DFT_IBS` re-parented; every consumer re-pointed (PWTerms, Seed, DensityMixer,
  IrrepCD_Fourier, `tBasisSet<dcmplx>`'s three `Iterate<>` delegates, OrthoFunctionFitter, tests).
  qcBasisSet now has ONE structure-neutral DFT face; `Orbital_DFT_IBS<double,dcmplx>` is a live spelling.
- **The metric worry ("molecules Dunlap-fit, solids don't") needed NO new machinery** — it was already
  discharged by (i)'s two-axis templating: the metric choice is expressed by the two accessor pairs
  (Repulsion3C=Coulomb, Overlap3C=overlap), the fit-basis argument type, and solve-vs-project behind
  `isOrtho()` in the fitter.  V1.1b (the Dunlap energy expression) is untouched and stays open.
- **What made the answer YES, found on contact:** every consumer outside qcBasisSet touched BOTH tensors
  through exactly two operations (forward D-contraction, adjoint fit-coefficient expansion) — nobody read
  raw data except test oracles.  And the "last separating member" had already gone production-dead
  without anyone recording it.  *Re-derive the defect before picking from an item's menu*, again.

*(original item text follows, moved in full from the worklist)*

- **V1.1 ⚠️ (i) ✅ DONE `57ca229e`; (ii) HALF done; (iii) RECLASSIFIED — one question left.**
  - **(i) DONE, and "trivial" was wrong in the way that mattered** (user, asked directly: *"I might have
    been wrong!!"*).  The one-line version — make the fit args `FIT_*_ABS<T>` — asserts *orbital T == fit
    T*, and that is precisely backwards for the case the type plan turns on.  The fit basis is a property
    of the RUN, not of a block: every call site builds it once from the whole `tBasisSet`, and all irrep/k
    blocks share it.  So when `BlochQN::IsReal()` makes TRIM k-blocks REAL, those real blocks still share
    the run's COMPLEX G-space fit basis.  ⇒ **two independent axes**, now
    `template <class T, class TFit = T>`.  `Orbital_DFT_IBS<double,dcmplx>` had NO SPELLING before.
    Purely additive (the default keeps every existing instantiation), so none of the twelve
    implementers/consumers changed.  Three of four combos are meaningful (user): `<double,double>`,
    `<dcmplx,dcmplx>`, `<double,dcmplx>`; `<dcmplx,double>` is never instantiated.
    - It also retired a self-CONTRADICTION that probably caused the split in the first place:
      `CreateCDFitBasisSet` returned `FIT_CD_ABS<T>*` while `Repulsion3C` accepted only
      `FIT_CD_ABS<double>`, so `Orbital_DFT_IBS<dcmplx>` could not consume the fit basis it created.
    - `Band_FT_IBS`'s own comment already knew — *"never assuming orbital==fit"* — it just had no second
      parameter to say it with.  **A comment that states an invariant the types cannot express is a
      standing invitation to look for the missing parameter.**
  - **(ii) HALF: the tensor element type follows `TFit`.**  A no-op today (every instantiation has
    `T==TFit`) kept for its DOCUMENTATION value (user): it says WHICH AXIS decides, which was the
    ambiguous part.  Correct for all three meaningful combos.  **Still open:** `ERI3` vs `G_ERI3`.
  - **(iii) RECLASSIFIED — not "extras", and not separable.**  Both parts are consequences of the merge:
    the `using Integrals_Overlap<dcmplx>::MakeOverlap` decl is LOAD-BEARING (it un-hides the no-arg form
    against the Fourier `MakeOverlap(f(G))`, so it is only noise once that member moves — and that member
    is the last thing separating the classes); and hoisting `CreateXCQuadrature` pays nothing alone,
    because `Band_FT_IBS` derives from `Orbital_1E_IBS<dcmplx>` NOT `Orbital_DFT_IBS`, and molecules
    already get the identical neutral default from `tBasisSet<double>`.
  - **⇒ ONE question remains: can `Band_FT_IBS` be `Orbital_DFT_IBS<dcmplx,dcmplx>`?**  Plan-level; owner
    is doc/RealComplexPlan.md.
  *(original text follows)*
  **`Orbital_DFT_IBS` ⇄ `Band_FT_IBS` merge.**  User (2026-08-05): Orbital_DFT_IBS simply
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

## LANDED 2026-08-17 — V1.11 COMPLETE (the occupation seam), five increments, 717/717 green each

The LAST doc/RealComplexPlan.md §7 prerequisite.  Design RULED by the user 2026-08-17 (recorded in
SCFStrategyPlan §5b): **policy-owns-state**; **abstract `OccupationPolicy` homed in qcElectronConfiguration**
(D6) with concretes ASSEMBLED from two axes (occupancy {Integer, Fermi, Held} × ranking {bare, MOM}) rather
than multiplied; **the EC class network becomes DATA** (counts + the reservoir partition) while the LIBRARY
keeps its name and hosts the policy; **observer/DIP if the policy needs to cross the DAG**; **MnO real-TRIM
work waits on this seam**.

- `43bbebad` — **inc 1: `FillResult{electronsLeft, minusTS, DPrime}`** kills the `ds_t` tuple whose double
  slot meant leftover electrons for integer fills and −TS for Fermi ones (its own doc had to warn).
- `0c818835` — **inc 2: `ReservoirPartition{spansSpatial, spansSpin, ranksIntegerFill}`** replaces the three
  EC mode bools (8 encodable states, ~4 meaningful; aufbau∧globalFermi now UNREPRESENTABLE);
  `tCompositeWF::FillOrbitals` becomes ONE loop over the reservoir grouping — generalising the map
  `FillOrbitalsSharedFermi` already built internally.  All four fill modes are one algorithm over different
  partitions; the held-fill + kT=0-staircase degrades are one documented flag.
- `841eadf2` — **inc 3: `OccupationPolicy<T>`** — the WF face sheds `SetMOM`/`SetSmearing`/`GetEntropyTerm`/
  `AdoptMOMReference`/`ReleaseMOMReference`; MOM references/counters and the −TS aggregate are policy state;
  the iterator owns the slot (built like the mixer) and its PUBLIC `AdoptMOMReference` face is unchanged.
  The seed fill runs on the policy's EXPLICIT default (prescribed integer, kT=0, no MOM) — closing the
  MEASURED D11 hazard (config arrived only at Iterate while Init filled in the ctor; metals lost charge:
  Al Σw·n=2.25 vs 3).  **DAG lesson: the naive qcElConfig→qcOrbitals link is a linker-rejected CYCLE**
  (qcOrbitals sits above via qcChargeDensity) — resolved with the CLAUDE.md DIP example verbatim: an
  `OrbitalView<T>` face OWNED by qcElConfig, implemented by `TOrbitals` from above.
- `092d1da8` — **inc 4: ONE `TOrbitals::Fill(const BlockFill&)`** replaces the five `TakeElectrons*`
  virtuals (each policy had added one — the OCP violation the item named).  `BlockFill` (T-free, beside
  `OrbitalView`) IS the two-axis product as data; `OccupationPolicy::DecideBlockFill` PRODUCES it (the
  4-way branch + the measured Λ-calibration lore moved off `tIrrepWF`); the five bodies became private
  realizations.
- `2398dd07` — **inc 5: `HeldOccupationPolicy`** replaces the `holdBlock` bool — the §5-flagged
  {occupation × direct-min} cell named.  Wraps the run policy (shared IMOM clocks, shared −TS aggregate);
  stored-order spec, kT structurally 0, capture no-op, `HoldsStoredBlocks()` degrades the reservoir loop.
  The consulted methods went virtual; OT+smearing later is a sibling that holds the block but keeps kT —
  a new policy, not a new bool.  **CORRECTION (user, 2026-08-17, same day): this does NOT make the policy
  the abstract interface D1 ruled** — the base is still a concrete whose behaviour `Configure`'s flags
  select per fill; only `Held` is a true derived policy.  The finishing move (the Policy/State split:
  persistent `OccupationState` + factory-assembled abstract concretes, `Configure` dies) is RULED LIKED
  and deferred — filed as **R2.21**.

Every increment was a full 717/717 sweep with the smeared/metal anchors (AlFCCAnnealedMetal, NaF traces,
MnO shared-spin) in the enabled set.  The refactor never changed WHICH orbitals fill — only who decides.

Still at call sites by design: the shared-μ metal fill's spec (the reservoir-driver migration — a future
session, not a V1.11 obligation).

*(original item text follows, moved in full from the worklist)*

- **V1.11 Occupation seam: two-phase SCF-WF construction + the `TakeElectrons*` family.**
  (i) `SCFWaveFunction::Init/SetMOM/SetSmearing/AdoptMOMReference/ReleaseMOMReference`
  (SCFWaveFunction.C:39-76): run-config known at construction time delivered by post-ctor setters;
  every concrete WF defends against the un-configured state.  (ii) Five `TakeElectrons*` virtuals
  on `TOrbitals<T>` (Orbitals.C:128-152) with a dual-meaning tuple slot (the doc comment itself
  warns "here the double is −TS, NOT leftover electrons"); every new occupation policy so far added
  a virtual to the abstract face (OCP).  One `Fill(const OccupationPolicy&)` returning a named
  struct; the policy object also feeds the WF ctor.  Cross-ref: this IS SCFStrategyPlan's
  "occupation seam" — design it there.

## LANDED 2026-08-16 — V1.5 (the `G_FieldEvaluator` ISP split) + the grid-reporting redesign, 717/717 green

Two commits; the reporting redesign came FIRST because the split's one open design question ("does
`EmitGridReport` belong on the quadrature face?") answered NO and pulled a thread:

- `f18a6ee9` — **Grid reporting: providers self-report, role-labeled, latch-free.**  User rulings:
  *"each object with something to report should decide for itself when and what to report — external calls
  to activate reporting should be the exception"*; *"if grids are created the user wants to know … if every
  new grid reports, we at least have a chance of pruning"* (a created-but-unused grid announcing is a
  FEATURE); the raw-`cout` line beside the report entry *"sounds like a bug — send everything to
  CurrentReport"*; and the PWTerms `static const void*` latch was called out (*"void* has no place in
  modern c++"* — it is CLAUDE.md law already).  Mechanism: `report::EmitAt` is now IDEMPOTENT (identical
  value at the same path = no write, no re-render), so providers announce UNCONDITIONALLY and dedup lives
  in the report, RUN-scoped — the two function-local-static latches (PWTerms' address-keyed `void*`,
  UnitCell's `lastID` string) leaked their dedup across runs and are gone.  The grid's ROLE
  ("xcQuadrature"/"densityFit") is a construction-time fact only its factory knows (R2.16), stamped on
  `PlaneWaveFit_IBS` at creation; the CD grid now announces too — which is exactly what §K's grid
  divergence will need the report to show.
- `9ebaebdb` — **V1.5 proper: one 10-method union → four client-named faces** (the V1.27 lesson):
  `G_FieldEvaluator` (evaluate a fitted map: `EvalField`/`EvalFieldGradient` — the ortho fitters' op(r)),
  `G_Quadrature` (`GridPoints`/`RhoOnGrid`/`ForwardFFT`/`GridCoeff`/`FieldCoeffs`/`Integral` — exposed
  through `GriddedScalarFitter::Grid()`, one owner as #7 left it), `G_StructureFactor`
  (`MakeFourierDensity` — the seed's single ask), `G_SpectralFilter` (`ApplySpectralFilter` — the mixer's
  single ask).  `PW_Grid_Evaluator` implements all four; every consumer casts to exactly its face;
  `Factory(cFIT_SF_ABS)` checks both required faces at the construction seam.  Cut to §K's end state: when
  the fit grids densify, what changes is which object implements `G_Quadrature`, not who asks.
- **Why the §K pin no longer blocked:** the one-owner rule landed with review-finding #7
  (`GriddedScalarFitter::Grid()`), and the V1.1 collapse had already removed the last orbital-assuming
  method from the union — the face was a pure density/potential grid engine before the split touched it.
- **Session status audit of FittingCleanupPlan** (what prompted this): H essentially landed via V1.26/V2.4
  (`MeshParams.eCut` + XCPolicy auto-sizing); I.2 landed (`GridCutoffFactor`→`relCutoff` threaded and
  consumed); §K's `ProjectedScalar_G`→`cvec_t` sub-bullet is OVERTAKEN (the class no longer exists); K's
  core (the densification) is now UNBLOCKED; I.1's residual is the `GetEpsXc()=0.75*GetVxc()` base default
  (exact for Dirac exchange, a silent-wrong inherited default for GGA); C stays dead last.

## LANDED on branch `solid-cleanup` (qchem1, 2026-08-05) — 665/665 ctest green

- `06e23f5d` — **R1.1, R1.2, R1.6, R1.8, R2.1.**  The `te.Exc=0.0` clobber deleted; `=`→`+=` sweep
  (Kinetic, DiracKinetic, RestMass, Ven, PW_Pseudo Een+E_alphaZ, PW_Kinetic, Vee); FittedVee
  re-expressed through locals (its `te.Eee` is a Dunlap combination of its OWN two pieces, so a
  bare `+=` would have read back another term's accumulation); Vee's two dead zero-stores deleted;
  `DM_ContractBlocks` pure virtual + its false "periodic path asserts out" doc fixed; two `Write()`
  pointer-streams dereferenced; FittedVee's missing cast assert added.
  **REFINEMENT TO R1.2 (found while doing it):** `GridChargeLost` must KEEP `=` — Delta_XC_Pol and
  Delta_VcorrPol both write the same value, so `+=` would double a health diagnostic.  A blanket
  sweep breaks it.  That it needs the OPPOSITE rule from every neighbouring field is independent
  evidence for V1.12 (get the diagnostic out of the energy value object).
- `80fc2ae8` — **V1.4** (Phi key → `map<Irrep,...>`, both sides on `Spin::None`).  Verified before
  landing: `Irrep::operator<` compares `SequenceIndex()` by VALUE (not handle identity) and
  `BlochQN::SequenceIndex()` uniquely encodes the k index, so periodic k-blocks keep separate
  tables instead of collapsing into one — the one way this change could have silently mixed one
  block's basis table into another's density.  All GPW/PW anchors unmoved.

- **R2.2** — `Kinetic` + `PW_Kinetic` collapsed to `Kinetic<T>` (new module
  `qchem.Hamiltonian.Internal.Kinetic`, inline, following the `IonIon<T>` recipe); `Internal/Imp/
  Kinetic.C` deleted, `PW_Kinetic` deleted from PWTerms.  Confirmed while doing it that PW_Kinetic's
  stated justification was FALSE: `Integrals_Kinetic<dcmplx>` IS instantiated, and the file that
  instantiates it says in so many words that the periodic dcmplx bases use the cached accessors.  So
  the periodic path now gets the `<p^2>` CACHE for free (it was calling the uncached `MakeKinetic()`
  every build).  Same value either way — the cache key carries k/Ecut/nG — and all PW/GPW anchors are
  unmoved.
  **GOTCHA WORTH KNOWING (now in CLAUDE.md):** moving the body into a template broke the build on
  `0.5*bs->Kinetic()` — *"operator* neither visible in the template definition nor found by ADL"*.
  A module that DEFINES a template using Blaze operators must `import qchem.Blaze` ITSELF; ADL at
  instantiation consults the template's DEFINITION context, not the instantiating TU's imports.
  Verified empirically: the import alone is sufficient and NO `<blaze/Math.h>` include is needed
  (both variants were built and compared).  Non-template code in an `Imp/` TU never hits this, so it
  will recur on every future T-templating collapse (`DiracKinetic`, `RestMass` if they go periodic).

- **R2.4 + V1.9 + R1.3** — the `Structure`→`UnitCell` cast sites are GONE, replaced by the user's
  free pry-out helpers (the `Symmetry::Atom::Getl` pattern), added to `qchem.ReciprocalLattice`
  (which already re-exports `qchem.UnitCell`, and sits in qcStructure below every consumer):
  `isPeriodicCell` (non-throwing probe), `GetUnitCell`, `GetReciprocalCell`, `GetReciprocalLattice`
  (the `Get` forms throw `std::bad_cast` via the reference cast).  All four sites converted:
  `SeedCD`'s anon `ReciprocalOf` DELETED (it was this helper, privately), SeedCD's flip-group cast,
  `MakeDensityMixer`, and `SCFIterator`.
  - **Design point that fell out:** the mixer factory needs a GRACEFUL fallback (it degrades to
    linear D-mixing and warns), so a purely throwing pry-out could not serve it — hence the
    `isPeriodicCell` probe alongside the throwing accessors.  Worth keeping in mind for the other
    capability faces: "can you?" and "give me it" are two different questions.
  - **R1.3 (the slicing copy) is FIXED as a side effect, not stopgapped**: SCFIterator now does
    `st->Clone()`, which is polymorphic (no slice of a derived cell) and returns the
    `shared_ptr<Structure>` the member already held.  V1.10b can still delete the member entirely.
  - Small ISP win: `MakeDensityMixer` was handing the concrete `UnitCell*` to
    `CreateVxcFitBasisSet`, which takes a `Structure*` — it now passes the neutral pointer.
  - Stale-comment batch: `Band_DFT_IBS.C`'s header no longer claims `PlaneWave_IBS` implements it
    (it documents the D1 re-decision instead), `PlaneWaveDFTUT`'s comment now names the real cast
    (`Integrals_Pseudo<dcmplx>`), and FOUR unused imports left SCFIterator (FourierDensity,
    FourierMixCD, Band_FT_IBS, ReciprocalLattice — each had exactly one occurrence in the file:
    its own import line).

- **V1.10b** — mixer creation moved OFF the structure-neutral iterator.  `MakeDensityMixer` (one factory
  running a three-way capability probe with a `cerr` fallback to linear D-mixing) is split into
  `MakeLinearMixer<T>` and `MakePeriodicMixer`; the choice is now a protected virtual
  `tSCFIterator<T>::CreateMixer` (base = linear) overridden by `SolidSCFIterator`.  The class that KNOWS
  it is periodic chooses, so the periodic factory takes Band_FT_IBS + UnitCell + FourierDensity as
  PRECONDITIONS and THROWS — user ruling 2026-08-06: "silent problems always end up consuming a lot of
  time... in the R&D phase failing loudly is encouraged" (a whole multi-hour run could previously
  converge on the wrong mixer behind one `cerr` line).  Everything the override needs is passed in, so
  iterator state stays private.  Verified first that all six Kerker/Pulay call sites are genuine GPW runs
  and that every `tSCFIterator<dcmplx>` in the tree IS a `SolidSCFIterator` (the lone `cSCFIterator` is a
  base-pointer holder), so the override is always reached.
  - **The `isPeriodicCell` prediction was HALF right (worth recording honestly).**  The probe did NOT
    disappear; it survives in exactly two places, and both are now CONTRACT CHECKS rather than
    behavioural branches: the `MakePeriodicMixer` precondition throw, and an `assert` in the base ctor.
    The base ctor keeps the cell snapshot because the `Structure` dangles by `Iterate` time and
    `SolidSCFIterator` inherits its constructors — but the runtime branch there became an assert, since
    the `if constexpr (is_same_v<T,dcmplx>)` above it already answers periodic-vs-molecular at compile
    time.  That doubled question (compile-time gate + runtime re-check) was the actual smell.
  - **Generalizable rule from this:** a periodic-vs-molecular `if` is a decision at the wrong altitude;
    the same probe reads fine as a *precondition*.  "Can you?" and "give me it" stay separate questions
    (V1.9), but "which am I?" should be answered by the type, not re-asked at run time.

- `72fecf8d` — **R1.4 + R1.5 + V1.3 (mechanism)** — three items, one session (2026-08-07); 667/667 ctest green.
  - **R1.4** the two silent-zero `Gradient()` overrides (FourierMixCD, `IrrepCD<dcmplx>`) now THROW.
    The plan said "assert"; asserts are compiled OUT of the build we test (`build/Release` is `-DNDEBUG`
    unless `QCHEM_RELCHECKED=ON`), so an asserting stub would still have returned the silent zero under
    `ctest`.  Throwing is loud in both configurations — the V1.10b ruling applied one level down.
  - **R1.5** `tChargeDensity::EvalBatch` DELETED; the fast batch is now an override of the inherited
    `ScalarFunction::operator()(rvec3vec_t)`.  The fork turned out to be LATENT (the one density-batching
    caller happened to use the `EvalBatch` spelling), not the live perf trap the item claimed.
  - **V1.3** both ε-adapter classes (`FittedEpsXc`, `FittedEpsCPol`) DELETED — but by naming
    `tDynamic_CC`'s method `GetEMatrix`, not by the planned `DM_ContractBlocks` reuse, which turned out
    to be blocked (a `tDynamic_HT` never receives `wholeBasis`, so it cannot enumerate the irrep blocks a
    block map needs — only `tDynamic_HF_HT` can).  The user's own V/E face IS the mechanism: one spelling
    was the entire collision.  Details + the blocked-route writeup under V1.3.
  - **Generalizable rule from R1.4:** "interim: assert" is only a real diagnostic where asserts are LIVE.
    Check the build's NDEBUG status before choosing assert-vs-throw for a *wrong-value* (as opposed to
    crash-soon) failure mode; the same NDEBUG hazard is already on record under V1.6.

**Process note for the next session:** build **`allTests`**, not just `ITMain`.  A first `ctest` in
qchem1 reported 160/590 failures that were ENTIRELY stale per-library `UT*` binaries (undefined
symbols + segfaults from the tree jumping ~30 commits while only ITMain was relinked) — zero real
failures.  After `ninja allTests` the same tree is 665/665.

**R1.3 is now SUPERSEDED-IF-V1.10b-LANDS and was deliberately NOT done** — the mixer-layering fix
deletes `itsKerkerCell` outright, so the one-line `Clone()` stopgap is wasted work if V1.10b lands
in the same session.


- **R1.1 ✅ DONE `06e23f5d`. `FittedVxcPol::GetEnergy` clobbers `te.Exc`** — `te.Exc = 0.0;` before delegating
  (Imp/FittedVxcPol.C:83).  Correct today only because exchange precedes correlation in
  Ham_DFTcorr_P's term list; a silent order-dependent zeroing of any prior XC contributor.
  Fix: remove the zeroing (struct is zero-initialized), let the children `+=`.

- **R1.2 ✅ DONE `06e23f5d` (with one CORRECTION, below). `=` vs `+=` on `EnergyBreakdown`.**  Assigners:
  Kinetic, DiracKinetic, RestMass, Ven, Vee, IonIon, PW_Kinetic, PW_Pseudo (`te.Een=`);
  accumulators: PP_Local/NonLocal, Vxc, FittedVxc, PW_Hartree (`te.Een+=`), Delta_*.  PW_Pseudo(=)
  and PW_Hartree(+=) share `Een` in one Hamiltonian, correct only because the static list is
  iterated before the dynamic list (Imp/HamiltonianImp.C:79-81).  Sweep to `+=` everywhere.

- **R1.3 ✅ DONE `38a1ebd6` — fixed via `Clone()`, not stopgapped. `UnitCell` SLICING copy** — SCFIterator.C:163
  `make_shared<const UnitCell>(*cell)` slices any derived cell (e.g. FCCUnitCell) to a plain
  UnitCell; `Structure::Clone()` exists.  (The cast itself → V1.9.)
  **What the Kerker code does with the cell (user Q, answered 2026-08-05):** the snapshot is only
  passed to `MakeDensityMixer` (SCFIterator.C:244), which uses it for exactly two things
  (DensityMixer.C:327-328): (1) as the plain `Structure*` argument to `CreateVxcFitBasisSet` (no
  concrete-UnitCell need at all), and (2) `MakeReciprocalCell()` → the `ReciprocalLattice` for the
  q²/(q²+q₀²) Kerker filter.  So the entire concrete requirement collapses to "give me your
  reciprocal lattice" — exactly the V1.9 periodic-geometry face.  Slicing severity: BENIGN today
  (FCCUnitCell's ctor just forwards a special cell matrix — no extra state/behavior, so the sliced
  copy keeps itsA intact); latent, becomes live the moment UnitCell grows derived state.  The deep
  copy exists only because `Lattice_3D::GetStructure` returns a temporary that dangles by Iterate
  time (the in-code comment says so) — `Clone()` fixes the slice now; the lifetime fix upstream
  would remove the snapshot entirely.  **SUPERSEDED-IF-V1.10b-LANDS**: the mixer-layering fix
  deletes `itsKerkerCell` outright.  Do the one-line `Clone()` only as a stopgap if V1.10b is not
  in the same session.

- **R1.4 ✅ DONE `72fecf8d` (THROW, not assert — deviation explained). Silent zero `Gradient()` overrides** — FourierMixCD.C:75 and IrrepCD<dcmplx>
  (Internal/Imp/IrrepCD.C:333) return `rvec3_t(0,0,0)` with no assert; a future GGA/plotting
  consumer through the neutral ScalarFunction face gets silently wrong ∇ρ.  ~~Interim: assert.~~
  Real fix (a `DifferentiableField` capability face) → V1.6-adjacent design.
  **Both now `throw std::logic_error` naming the missing implementation.**  The planned `assert` would
  have been a NO-OP in the build we actually test: `build/Release` is `-DNDEBUG` (`QCHEM_RELCHECKED=OFF`,
  CMakeLists.txt:77-93), so an asserting stub still returns the silent zero under `ctest`.  A throw is loud
  in BOTH configurations, and matches the V1.10b "fail loudly in the R&D phase" ruling.  Verified no live
  caller: the density `Gradient` consumers are all `<double>` (SlaterExchange, the fitters, the composite/
  polarized/spin forwarders) — nothing on the periodic path asks for ∇ρ, and 665/665 stays green.

- **R1.5 ✅ DONE `72fecf8d`. `tChargeDensity::EvalBatch` duplicates `ScalarFunction::operator()(rvec3vec_t)`.**
  Identical signature after alias expansion.  `EvalBatch` DELETED (the inherited `ScalarFunction` batch
  op carries the same pointwise-loop default); FourierMixCD's fast factorized-phase path now overrides
  `operator()(rvec3vec_t)`; the four call sites in `XC_GridEngine` (PWTerms.C:375,403,404,411) spell it
  `(*cd)(mesh->Points())`.
  **CORRECTION to the claim above (verified while doing it): the fork was LATENT, not live.**
  OrthoFunctionFitter.C:101 batches `ps.GetScalarFunction()` — a `ProjectedScalar_R`'s field (v_xc, ε_xc),
  never a density — so no FourierMixCD reaches it today.  The ONLY site that batched a density was
  `XC_GridEngine`, and it happened to use the `EvalBatch` spelling, i.e. it got the fast path BY LUCK OF
  SPELLING.  The trap was one neutral-face caller away, not already sprung.

- **R1.6 ✅ DONE `06e23f5d`. `Write()` streams raw POINTERS (hex addresses)** — Imp/FittedVxc.C:111 (`os << itsLDAVxc`)
  and Imp/FittedVxcPol.C:93 (`os << itsUpVxc << itsDownVxc`); `qchem::op<<` binds only
  `const Streamable&`, so these hit `void const*`.  Need `*`.

- **R1.7 ✅ DONE `26af31b6` (2026-08-10). `SymmetryAdapted_IBS::MakeDirect/MakeExchange` return empty
  `ERI4{}` silently**
  (SymmetryAdapted_IBS.C:80-81, comment "never called") — failed QUIETLY: a future generic-ERI4
  caller got a zero Fock contribution.  ~~Interim: assert.~~
  **USER RULING (2026-08-05): SymmetryAdapted_IBS should not be inheriting the (ERI4 part of the)
  `Orbital_HF_IBS<double>` interface.**  The interface's own doc comment already pointed there:
  "Client code rarely need ERIs directly, they only need the contraction over a density matrix."
  Split `Orbital_HF_IBS` into (a) the CONTRACTION face (`Accumulate*` — what the HF terms consume)
  and (b) the ERI4-SUBSTRATE face (`MakeDirect/MakeExchange` + `Direct/Exchange` caches — the
  implementation detail behind the default contraction path).  Concrete AO bases implement both; the
  SALC decorator implements only (a); terms consume only (a); the dummy `ERI4{}` bodies die.

  **What landed — the split exactly as ruled, plus one thing the item did not predict:**
  - `qchem.BasisSet.Orbital_HF_IBS` is now the CONTRACTION face and nothing else: four pure-virtual
    `Accumulate{Direct,Exchange}[Both]`.  No ERI4 anywhere in it.
  - `qchem.BasisSet.Internal.Orbital_ERI4_IBS` (NEW) is the substrate: `MakeDirect`/`MakeExchange`
    pure virtual, the cached `Direct`/`Exchange`, and the ONE implementation of the four
    `Accumulate*` in terms of them.  `src/BasisSet/Imp/Orbital_HF_IBS.C` moved verbatim (bar the
    class name and the cross-cast) to `src/BasisSet/Internal/Imp/Orbital_ERI4_IBS.C`.
  - Concrete AO bases implement both faces.  The two evaluator-templated mixins that BUILD the
    blocks — `Atom::Orbital_HF_IBS<E>` and `Molecule::Orbital_HF_IBS<E>` — were renamed to
    `Orbital_ERI4_IBS<E>`, so a mixin's name says which face it supplies; likewise the RKB/DHF
    `Orbital_RKB_HF_IBS_Imp`.  `SymmetryAdapted_IBS` derives from the contraction face ONLY, and the
    two `{return ERI4();}` bodies are deleted — the question can no longer be asked of it.
  - **The unpredicted win: `Internal.` is the honest address for the substrate.**  Grep says nothing
    outside qcBasisSet ever names an `ERI4`; qcHamiltonian and qcChargeDensity both reach the basis
    as `rohfbs_t = Orbital_HF_IBS<double>` and only ever call `Accumulate*`.  So the substrate went
    into an `Internal.` module and `Orbital_HF_IBS` stopped re-exporting
    `qchem.BasisSet.Internal.ERI4` — two libraries that had the 4-index type in scope for no reason
    no longer do.  That is the DIP half of the item, and it came free with the ISP half.  (Unit
    tests that pin the cache and the bra-ket symmetry import the Internal module explicitly and
    iterate on the new `Real_ERI4_OIBS` typedef — the sanctioned test cheat.)
  - **The one new cast, and why it is the sanctioned kind.**  `Accumulate*` take the partner as the
    CONTRACTION face and C++ has no contravariant parameters, so the substrate implementation
    cross-casts it back through `Orbital_ERI4_IBS::Substrate()` — abstract→abstract, the direction
    the project rule allows — and it THROWS naming both bases rather than asserting (R1.4's ruling:
    `build/Release` is `-DNDEBUG`, so an assert is a no-op in the configuration we actually test).
    The condition it reports is precisely the one the empty `ERI4{}` used to hide: an ERI4 basis
    paired with a basis that has no `(ab|cd)` block spanning the two.
  - **Anchors:** whole suite green; the load-bearing ones are `Cache4Tests.*` (canonical-only cache +
    the `J(a,b)=J(b,a)^T` check against the UNCACHED `MakeDirect`), `PGSymmetry.
    decorator_coulomb_matches_AO_slice` (the SALC path that owns this item), `M_MEvaluator`/
    `M_LibCint` `matrix_3C_4C_match_scalar`, and `M_Calculation.WaterSymmetryLibCint`.
  - Doc references updated in `doc/ERI4Rework.md` (§2 substrate address, §5.4 SALC caveat).

- **V1.31 ✅ DONE `627a4ff9` (2026-08-10).  `SymFockCache` deleted; the SALC path builds ONE whole-AO Fock
  and slices it.**  The item was filed twice wrong before it was right, and both wrong versions are kept in
  doc/CleanupCandidates.md because the corrections are where the value is:
  1. Filed as "caching has escaped `DB_Cache_RAM`" — WRONG on both halves (the two memos hold contracted
     matrices, not J/K tables; and `DB_Cache_RAM` is not immortal, it LRU-evicts ERI4s).  User caught it
     from the description alone.
  2. Re-filed with the user's ruling "version counter over elementwise D compare" — right in general,
     REFUTED at this site: `tPolarized_CD::Version()` forwards to its Up child, so Up and Down would be
     indistinguishable inside one sweep and the Down pass would be served the Up channel's \f$J_{AO}\f$.
     The elementwise compare was load-bearing, not lazy.
  3. The user then asked for a flow diagram, on the hunch that *"this whole thing is just designed wrong.
     We are somehow caching the wrong thing in the wrong place."*  **That was the correct diagnosis**, and
     the chain showed it: the composite's canonical-PAIR loop exists for bases with per-irrep-pair ERI4
     blocks, which a SALC basis does not have at all (R1.7).  Driven through that loop, the decorator
     rebuilt the SAME whole-molecule AO Fock once per irrep; `SymFockCache` existed only to stop it being
     once per PAIR.  So the memo was a prop holding up a loop that should not have been running.

  **What landed — remove the loop, not the staleness test.**
  - New CAPABILITY FACE `BasisSet::WholeSystemFock_IBS<T>` (in `Orbital_HF_IBS.C`, beside the contraction
    face so qcChargeDensity imports nothing new): `AODimension`, `AddAODensity` (\f$D_{AO}\mathrel{+}=ODO^T\f$),
    `MakeAOFock` (the ONE build), `SliceAOFock` (\f$F_{ab}\mathrel{+}=O^TF_{AO}O\f$).  Four pure virtuals,
    no stubs, no default bodies — a basis either has the face or does not, the `tSpinResolved_CD` idiom.
  - `tComposite_CD::Accumulate{Direct,Exchange}All` probe ONCE at the top of the sweep and, when the leaves'
    basis has the face, run: sum the AO densities → ONE build → N slices.  Absent the face, the canonical-pair
    loop is untouched, so every ERI4 basis is bit-for-bit unaffected.
  - `tDM_CD` gains the probe `WholeSystemFock()` and `AddAODensity()`; only the leaf can fold its own D, and
    the SLICE needs no density, so the composite drives that through the basis face directly.  Two virtuals,
    not four.
  - **The exactness is an identity, not an approximation:** J and K are LINEAR in D, so
    \f$\sum_\Gamma F_{AO}(O_\Gamma D_\Gamma O_\Gamma^T) = F_{AO}(\sum_\Gamma O_\Gamma D_\Gamma
    O_\Gamma^T)\f$.  N whole-AO builds per sweep become one.
  - `SymFockCache` is GONE, and with it the elementwise density compare AND its incomplete key (it omitted
    `Ocd`).  **Both defects died of the restructure rather than being fixed** — which is the sign the
    diagnosis was right.
  - **VERIFIED LIVE, not merely green.**  Passing tests prove nothing if the old path still runs, so the
    per-pair route was temporarily booby-trapped to throw and `M_Sym` (9 cases: HF/DFT × polarized/
    unpolarized × symmetry on/off) plus `M_HF_U` were re-run: all 15 still passed, so the SCF never reaches
    the pair route for a SALC basis.  Probe then reverted.  The per-pair methods remain as the direct-call
    route the `PGSymmetry` unit tests exercise.

- **V1.6 / V1.7 / V1.8 / V1.10 ✅ DONE `2d0f6982` (2026-08-16).  Four of the six RealComplexPlan
  prerequisites — one defect in four places: a face declaring a capability half its hierarchy lacks, so
  the other half must write denials.**
  - **V1.6/V1.8, exact exchange.**  `Vee`/`Vxc` are added ONLY by the molecular HF Hamiltonians (the
    periodic `Ham_PW_DFT` adds `Vee_Hartree`, never exact exchange), so the four `Accumulate*` were
    real-only sitting on the general T-templated density face.  Now `tHF_System_CD` (spans every block)
    and `tHF_Pair_CD` (the composite↔leaf pair protocol — never public-face business, its only callers
    were two loops in Imp/CompositeCD.C), inherited through `conditional_t`.  **T-typed, not erased**, per
    the plan: an impossible pairing must fail to COMPILE, not throw.
  - **The user's two corrections did the real work, and both generalise:**
    1. My stated reason for preferring a face over double dispatch (mixed-T pairs multiplying visitor
       arms) was WRONG — mixed-T pairs never arise, every pairwise op pairs a block with its same-irrep
       counterpart, and the site is a `template<> IrrepCD<double>` specialization anyway.  The conclusion
       survived on OCP + data-hiding alone, which is the better argument since it does not depend on the
       T question.  **Check whether your reason and your conclusion are actually connected.**
    2. Make the face OPERATION-named, not a block getter.  `CompleteDirectPair` = "finish this
       contraction with your block".  An abstract block ACCESSOR would have satisfied the compiler while
       reintroducing the `GetDensityMatrix()` that CLAUDE.md cites this very class for NOT having.  That
       naming is also what made V1.8's cast EVAPORATE rather than move — as the user predicted when
       sequencing V1.6 first.
    3. Then: the complex leaf still DECLARED the four and so defined four empty bodies.  *"It sounds like
       we segregate the interface ... so that the complex CD is not forced to make fake functions."*
       Correct — now a real-path-only CRTP mixin (`IrrepCD_HFPair`), CRTP + friendship so the
       implementation reaches the block's own D/basis without any of it becoming public.  The complex leaf
       declares nothing, needs no vtable slots, defines nothing.
  - **The `-DNDEBUG` hazard is closed.**  Those `void` assert-only bodies were silent NO-OPS in the build
    we ship: a bare leaf reaching `Vee::AccumulateAll` yielded a ZEROED J and a wrong Fock, with no
    diagnostic.  Both terms now cross-cast to the system face and THROW.
  - **V1.7, the periodic trio — the largest LSP block here, nine denials, now zero.**  The mechanism was
    already present and correct: `FourierDensity` declared the trio pure-virtual and `FourierDensityBase<T>`
    already handed it to dcmplx alone.  The three families simply RE-declared it in their own bodies for
    both T, so the finite instantiation answered with `assert(false)` in an if-constexpr dead branch.  Three
    periodic-only CRTP mixins; consumers already cross-cast to the `FourierDensity` face, so nothing outside
    changed.  **When a denial appears, look for the conditional base that already exists.**
  - **V1.10, both abstract→concrete basis casts.**  The `SymmetryAdapted_IBS` one DISSOLVED into a face
    added three items earlier: V1.31's `AddAODensity`/`MakeAOFock`/`SliceAOFock` are literally the three
    steps the cast was open-coding, so the partner folds its OWN block up and nobody hands out SALC columns.
    The DHF one became `Orbital_RKB_Pair::MakeDirectAgainstL` — the partner builds against the caller's
    large component instead of surrendering `itsRKBL`.  **A cast that reaches for private state is usually
    an operation that already exists somewhere as a face.**
  - 716/716 green.  An earlier sweep showed 22 failures and was NOT a regression — it ran while ninja
    relinked underneath it.  Re-verified on a sweep started after the final build with nothing touched.

- **R2.18 ✅ NAMES DONE `86c5b24d` (2026-08-10); encapsulation half deliberately left open.  The `Make`/`Get`
  pair in qcHamiltonian.**  USER: *"GetMatrix goes through caching ... if there is no cache it calls
  MakeMatrix() which does return by value.  So if you need a return by value override it should be the
  MakeXXX() call."*  qcBasisSet follows the pair without exception (11 `Make` verbs); qcHamiltonian spelled
  the same role `CalculateMatrix` on the static base and `CalcMatrix` on both dynamic ones — two names for
  one role in one file, neither the project verb.  All renamed to `MakeMatrix` (plus
  `CalculateMatrixRadial`→`MakeMatrixRadial`); the compiler found every override.
  - **The visibility half was filed on a FALSE premise, and the user corrected it.**  The item said
    qcBasisSet makes `Make` public and qcHamiltonian makes it protected, so qcHamiltonian is the outlier.
    **Backwards:** *"All the MakeXXX() functions were originally protected.  For DFT the 3C versions
    (MakeOverlap3C, MakeRepulsion3C) still are.  It seemed like these were purely internal functions ... but
    that turned out to be incorrect in some cases."*  Protected is the ORIGINAL state; the public ones are
    DRIFT, each from a case where "purely internal" proved wrong.  qcHamiltonian is not an outlier, it is
    un-drifted.  **Standardising on qcBasisSet would have standardised on the drift** — which is what the
    item, as filed, would have led someone to do.
  - **Ruled LOW priority and left open** (*"I have no strong policy on this right now.  Maybe the right
    policy will emerge as we refactor.  My intuition says it is a low priority decision."*), with the
    standing practice that a `Make` going public must record WHY at the declaration — see the convention
    box at the top of doc/CleanupCandidates.md.  Cheap to defer, expensive to guess.
  - One comment at `SCFIterator.C:186` was initially left stale (MnO DO-NOT-TOUCH list); fixed in the same
    commit once the user confirmed the campaign was paused.

- **R2.19 ✅ DONE `86c5b24d` (2026-08-10).  `FittedVxcPol` copied a matrix its child already owned.**
  Found by following the user's Get/Make remark into the code.  `FittedVxcPol::CalcMatrix` was a PURE
  FORWARDER — both branches ended in `(s==Spin::Up ? itsUpVxc : itsDownVxc)->GetMatrix(...)`, returning BY
  VALUE a matrix the child's own per-Irrep cache already held stably, which the `tDynamic_HT_Imp_NoCache`
  base then stored in scratch purely to have something to return a reference to.  One full matrix copy per
  call, per spin, per irrep, per SCF iteration, to satisfy a signature.  Its direct HF analogue `VxcPol`
  (Imp/VxcPol.C:39-51) had always done it right: override `GetMatrix`, return the child's reference.
  - **Fix:** `FittedVxcPol` now overrides `GetMatrix`, forwards the child's reference, and **no longer
    derives from `tDynamic_HT_Imp_NoCache` — it has no `MakeMatrix` at all.**
  - **The generalisable test, and the reason this item existed:** *if your `MakeXxx` does not COMPUTE
    anything, you are not an implementer, you are a forwarder* — override `Get` and hand back what you
    forward to.  The old code had to invent a `Make` (and a scratch slot for its result) purely to satisfy
    a base it should not have had.  The Get/Make split is what makes that visible; before it, "returns a
    reference to shared scratch" looked like a lifetime question rather than a wrong-base-class question.
  - **Verified NOT applicable to its sibling:** `FittedVcorrPol` KEEPS `NoCache` — it genuinely recomputes
    (one `itsVcFitter` shared across both channels, refit per spin), so R2.9(ii)'s Irrep-keyed scratch still
    earns its place and is now that one class's private business.

- **R2.9 ✅ DONE `268473b9` (2026-08-10).  Small Hamiltonian hardening — all three sub-items.**
  Original text: (i) `XC_GridEngine` constness laundering — Rho/RhoPol/Matrix/Phi non-const, called from
  const term methods via a non-const `shared_ptr`; every other cache in the module is `mutable`+const-method
  — align it (and note: no cross-invalidation between its two rho caches, scalar + Up/Dn, which can be live
  simultaneously).  (ii) `tDynamic_HT_Imp_NoCache::GetMatrix` returns a reference to shared scratch
  (HamiltonianTerm.C:104-107) — safe today only by immediate consumption; document or return by value.
  (iii) `Dynamic_HF_HT_Imp::itsWholeBasis` latches on first use (Imp/HF_HT.C:31) — assert on change.

  Done in the order (iii), (ii), (i) — by VERIFICATION cost, not importance: (iii) and (ii) are checked by
  the fast molecular tests, (i) needs a GPW sweep.

  - **(iii) THROW, not the assert the item asked for.**  Same deviation as R1.4, same reason: `build/Release`
    is `-DNDEBUG`, so an assert would be a no-op in the configuration we actually test — and the null-basis
    precondition three lines above, on the SAME parameter, already throws.  One function should not check two
    preconditions on one argument at two different volumes.  What a change would mean is a wrong answer, not
    a crash: `itsJKs` is keyed by BasisSetID and built by walking the latched composite, so a second basis
    would be served the first one's contraction — a plausible Fock from the wrong cross-irrep view.
    - **Follow-up worth filing separately (R2.16 shape, NOT done here):** the latch exists because
      `GetEnergy` has no basis argument, i.e. a run-stable fact is discovered lazily instead of being
      supplied at construction — exactly the pattern R2.16 rules against.  The precedent already exists:
      `Ham_DFT_U`/`Ham_DFTcorr_U` ctors take `const rbs_t* bs`; `Ham_HF_U`/`Ham_HF_P` do not.  Making it a
      ctor parameter is the real fix and is mechanical for the DFT family, but it changes two public
      Hamiltonian ctors and their call sites — more than "small hardening" scopes, hence left.
  - **(ii) NEITHER of the two options the item offered.**  Both were considered and both are wrong:
    - *Return by value* is not available: `GetMatrix` OVERRIDES `tDynamic_HT<T>::GetMatrix`, which returns
      `const hmat_t<T>&`, and C++ has no non-covariant return-type change on an override.  Changing the BASE
      to by-value would cost the CACHING sibling a full matrix copy per term per irrep per iteration, buying
      nothing — its whole point is to avoid recomputation, not copies.
      **USER, 2026-08-10 — the deeper reason, which reframes the whole option:** *"GetMatrix goes through
      caching ... if there is no cache it calls MakeMatrix() which does return by value.  So if you need a
      return by value override it should be the MakeXXX() call."*  Right: `Get`/`Make` is a PAIR, `Get`
      returns a reference and `Make` returns by value, and a by-value need is served by asking the other
      half — never by changing `Get`'s return type.  So "return by value" was not merely impractical here,
      it was reaching for the wrong half of an existing convention.  Two items came out of that remark:
      **R2.18** (qcHamiltonian spells the `Make` verb as `CalculateMatrix`/`CalcMatrix`, two names for one
      role, and makes it `protected` where qcBasisSet makes it public — so the convention's escape hatch is
      not actually open to a term's client) and **R2.19** (`FittedVxcPol` is a pure FORWARDER that copies a
      matrix its child's cache already owns; its HF twin `VxcPol` already returns the child's reference).
    - *Document it* understates the defect.  The real problem is not "a reference to shared scratch"; it is
      that **two implementations of one interface made DIFFERENT lifetime promises, invisibly.**  The caching
      sibling keys on `Irrep`, which folds in Spin via `GetIrrep(s)`, so its Up and Down blocks sit in
      separate slots; the NoCache one had a single `itsMat`, so they aliased.  A caller holding
      `tDynamic_HT<T>*` cannot tell which it has, and the natural polarized idiom — bind `up`, then bind
      `dn` — is correct against one and silently wrong against the other.  `FittedVxcPol`, ONE object
      serving both spin channels, is exactly where it would have been sprung.
    - **So: key the scratch by `Irrep` (reusing `tHT_Common`), and ALWAYS assign — never look up.**  Storage
      is shared with the caching sibling; the CACHING is not.  Both siblings now promise the same thing:
      valid until the next call for the same Irrep.  Cost is one matrix per Irrep, the bound the sibling
      already carried.  Verified latent, not live: every consumer in the tree is an immediate
      `H += t->GetMatrix(...)`; nothing stores the reference.
  - **(i) const + `mutable`, and the `shared_ptr` made `const` too.**  All four accessors (`Rho`, `RhoPol`,
    `Matrix`, `Phi`) are now const and every cache they touch is `mutable` — the idiom `tHT_Common::itsCache`,
    `tDynamic_HT_Imp::itsCacheVersion` and `Dynamic_HF_HT_Imp::itsJKs` already use.  `itsMesh`/`itsFold` are
    deliberately NOT mutable: construction-time, must not move.  The three terms' `engine_t` became
    `shared_ptr<const XC_GridEngine>`, which is what actually retires the laundering — previously the
    constness was defeated by holding a non-const pointer, so the const-ness of the term methods said nothing.
  - **The two-rho-cache hazard the item flagged is now PINNED, not just noted.**  `itsRho` and the
    `itsRhoUp/Dn` pair guard only their own serial and never invalidate each other, so an engine driven on
    both routes for different densities would report "fresh" while holding one stale raster.  It is
    unreachable today, and CHECKED rather than assumed: there is exactly ONE construction site
    (Imp/Hamiltonians.C:240) and its `if (polarized)` branch adds either two Pol terms (RhoPol only) or two
    scalar terms (Rho only) to a single shared engine — R2.16's own good pattern, where absence of a
    capability means the term is not in the list.  So each accessor now ASSERTS the other route was never
    used, and the members carry a `\warning` that any mixed/GGA route must add real cross-invalidation
    first.  Assert (not throw) is deliberate, and is the OPPOSITE call from (iii) in the same commit: (iii)
    guards a condition a caller could reach, so it must be loud in Release; this one is excluded by
    construction, so it is an invariant pinned for a future change — and a throw would cost a branch in the
    per-iteration rho path to defend against something no caller can do.

- **R1.8 ✅ DONE `06e23f5d`. `FittedVee` casts `bs` and dereferences with NO assert** (Imp/FittedVee.C:41-42) — the
  sibling sites at least assert.  (The odftbs_t casts themselves are sanctioned abstract→abstract.)


- **R2.1 ✅ DONE `06e23f5d`. `tDM_CD::DM_ContractBlocks` → pure virtual.**  The asserting default is DEAD — all three
  concrete families override it (IrrepCD both T, Composite, Polarized); pure-virtual is free today.
  Also fix its false doxygen ("the periodic path asserts out" — the periodic path CALLS it,
  PWTerms.C:158).

- **R2.2 ✅ DONE `48e25b74`. Collapse `Kinetic` + `PW_Kinetic` → `Kinetic<T>`.**  Both are 0.5×(kinetic matrix);
  PW_Kinetic's stated justification (PWTerms.C:59-60, "symmetric cache bypassed for complex") is
  FALSE — `Integrals_Kinetic<dcmplx>` is instantiated (Imp/Orbital_1E_IBS.C:30).  Follow the
  `IonIon<T>` collapse recipe (Internal/IonIon.C documents it); PW_Kinetic dies.

- **R2.3 ⛔ WITHDRAWN — NOT a free dedup; re-filed as part of V1.1 (verified 2026-08-07).**
  ~~Dedup the literal 4-line Imp mirror — Imp/Band_FT_IBS.C:15-25 == Imp/Orbital_DFT_IBS.C:10-20 (only
  `theCache<T>` differs).  Mechanically removable independent of the V1.1 merge.~~
  The claim "only `theCache<T>` differs" is FALSE.  Read side by side, the two bodies differ in THREE ways:
  (i) `theCache<dcmplx>` vs `theCache<T>`; (ii) the RETURN type — `const G_ERI3&` vs `const ERI3<T>&`;
  (iii) the ARGUMENT types — `cFIT_SF_ABS`/`cFIT_CD_ABS` vs `rFIT_SF_ABS`/`rFIT_CD_ABS`.  (ii) and (iii)
  are precisely V1.1's two listed blockers (the `ERI3<T>` vs `G_ERI3` return-type question and the
  hard-wired-REAL fit-face arguments).  So there is no dedup to do here that is not the V1.1 merge itself —
  a shared body cannot be written until those two are settled.  **Do it inside V1.1, not before it.**

- **R2.4 ✅ DONE `38a1ebd6`. Stale-comment/import batch**: Band_DFT_IBS.C header claims PlaneWave_IBS implements it
  (false); PlaneWaveDFTUT.C:1079 claims PW_Pseudo routes through Band_DFT_IBS (it casts
  `Integrals_Pseudo<dcmplx>`, PWTerms.C:46); SCFIterator.C:30-31 imports of
  FourierDensity/FourierMixCD are stale (not consumers); PWTerms.C:59-60 false cache comment dies
  with R2.2.

- **R2.6 ✅ DONE 2026-08-07. The `LDAVxc` bundle** — a "Hamiltonian term" whose `CalcMatrix`/`GetEnergy` call
  `exit(-1)`; its only real job is `GetScalarFunction()` (the fit callback).  `FittedVxc` holds it
  as a raw OWNING pointer to the CONCRETE class (Terms.C:291, DIP violation);
  `tFittablePotential` (HamiltonianTerm.C:113-119) exists solely for it; `FittedVxc::
  UseChargeDensity` is dead (overrides nothing, no caller — `newCD()` already triggers the refit).
  One fix kills all four: collapse LDAVxc to a plain `Fitting::ProjectedScalar_R` adapter, hold as
  `unique_ptr<ProjectedScalar_R>`, delete `tFittablePotential`.
  **All four verified before the cut, and all four are gone.**  The module
  `qchem.Hamiltonian.Internal.LDAVxc` is DELETED outright (both TUs + the two CMake entries); LDAVxc had
  exactly one user in the whole tree, `FittedVxc`.
  - **Went one step further than the plan, for free.**  Rather than an adapter that hands the fitter the
    `ExFunctional` (what LDAVxc did), the adapter SELF-EVALUATES: `VxcDensity{ex,cd}` returning
    `ex->GetVxc((*cd)(r))`, an exact mirror of the `EpsXcDensity` already sitting beside it (and of
    FittedVcorrPol's `PolVcDensity`/`PolEpsCDensity` pair).  Verified byte-identical first: ALL THREE
    `ExFunctional` subclasses define `operator()(r)` as literally `GetVxc((*itsChargeDensity)(r))`
    (SlaterExchange, VWN_Correlation, Libxc_LDA).  The V and E fields are now visibly the same shape,
    differing in one method call.
  - **Consequence worth acting on: this makes V1.13 a DELETION rather than a redesign.**
    `ExFunctional::InsertChargeDensity` now has ZERO callers tree-wide (it had exactly one, in LDAVxc),
    so `itsChargeDensity` is never set and the whole FIELD face of `ExFunctional`
    (`operator()`, `Gradient`, the member) is dead code.  The fitter reaches the field through the
    adapters instead, which take the density as a ctor argument — the hidden-init landmine V1.13
    describes cannot happen on this path any more.  NOT deleted here: removing the field face means
    `ExFunctional` stops being a `ScalarFunction<double>`, which IS V1.13's value-face/field-face split,
    and it also retires `SlaterExchange::Gradient`.  Do it as V1.13, now cheap and compiler-verifiable.
  - Not touched (still V1.13): `SetPolarized`/`isPolarized` (still no caller, so `isPolarized` is still
    permanently true).

- **R2.7 ✅ DONE 2026-08-07. `FittedCD::Clone()` — delete.**  Pure virtual (FittedCD.C:28) whose SOLE implementation
  asserts false and returns nullptr (Imp/FittedCDImp.C:58-64).  Dead contract clause; restore when
  the polarized-from-unpolarized path exists (real blocker per the assert message: a cloneable
  FunctionFitter).

- **R2.8 ✅ DONE 2026-08-07. `InsertStandardTerms<dcmplx>` = assert(false)** (Imp/HamiltonianImp.C:49-53) — a
  molecular-only convenience on the T-generic base; move it down to the real lineage.
  **RESOLVED BY DELETION, not relocation (USER RULING 2026-08-07: "the whole idea of having
  InsertStandardTerms at all was too clever by half").  `InsertStandardTerms` IS GONE; `rHamiltonianImp` is
  a plain alias again.**  Each of the seven callers now spells its core out in three `Add`s.
  - **The user's diagnosis, and it is the important part: `double` vs `dcmplx` was never the discriminator.**
    The bit that decides the core is BARE vs PSEUDISED nuclei, and it decides TWO things at once — the
    ion charge (Z vs Zion) AND which electron-nuclear term(s) exist.  Decisive evidence already in the tree:
    **`Ham_PP` is `<double>` and never called the helper either** (Imp/Hamiltonians.C:108-111 builds
    `Kinetic<double>` + `IonIon<double>(vloc->ZionFn())` + `PP_Local` [+ `PP_NonLocal`]).  So `double` sits
    on BOTH sides of the split; it is not a molecule-vs-solid decision either.
  - **And the basis lineage is a SECOND, independent axis** — it picks WHICH electron-nuclear
    implementation (`PP_Local`'s mesh quadrature vs `Ven_PP_Short`/`_Long`'s G-space route), not WHETHER the
    nuclei are pseudised.  The old comment "the complex Hamiltonian assembles its terms explicitly" named the
    symptom; the cause is that it is a PSEUDOPOTENTIAL Hamiltonian, exactly like the `<double>` `Ham_PP`.
  - Add the Dirac lineage (DiracKinetic + RestMass, and no ion-ion at all) and the "standard" set was
    standard for 7 of 11 Hamiltonians.  Three explicit `Add`s read better than a name that hides them.
  - *(Kept for the record, since it is a genuine C++ gotcha that shaped the original stub:)*
    `template class tHamiltonianImp<dcmplx>;` is an EXPLICIT instantiation, which instantiates every member
    definition — so while the member lived on the T-generic base, the dcmplx side was *required* to have a
    body, which is why it was an `assert(false)` rather than simply absent.  "Just don't define it for that
    T" is not available under explicit class instantiation.

- **R2.10 ✅ DONE 2026-08-07. `Fit_IBS::SetMesh` → ctor parameter.**  Two-phase construction; the construction-time
  principle was already settled for the XC quadrature — build the mesh first, hand it in.
  `SetMesh` deleted; `Fit_IBS(const Structure&, const MeshParams&)` builds and owns the mesh.  All three
  `EFit_IBS` lineages (Atom / PG_Cart / PG_Spherical) forward to it — `Fit_IBS` is a VIRTUAL base of each,
  so the most-derived `EFit_IBS` initialises it, which is exactly where the creators already have both
  arguments.  The six `CreateCD/VxcFitBasisSet` sites became one-liners.
  **The invariant is now ESTABLISHED rather than re-checked:** the two asserts that guarded the half-built
  state inside `Norm()` and `Overlap(f)` ("SetMesh must run before any numerical integral") are replaced by
  one ctor postcondition.  Every numerical integral this class offers runs over that mesh, so there was
  never a valid state between "constructed" and "has a mesh".

- **R2.11 ✅ DONE 2026-08-07. `DB_Cache_RAM.C`** — a screenful of `-Winconsistent-missing-override` warnings on every
  qcBasisSet build (`Get`/`Register`/`GetCache*`).  Mechanical `override` sweep — 15 members marked.
  Swept the qcHamiltonian ones too (`FittedVxc`/`FittedVcorrPol`, newly inconsistent because V1.3's
  `GetEMatrix` arrived with `override` while its siblings had none).  **`ninja allTests` is now
  warning-free**, which is the real win: a screenful of known-benign warnings is where a NEW one hides.

- **R2.12 ✅ DONE 2026-08-07. `UnmatchedCounts`/fold `tol` defaults** — 1e-8 fractional as a literal in three places
  (GPW `CreateXCQuadrature`, SymmetrizeMesh overloads); name it once.
  Now `qchem::kMeshMatchTol` in SymmetrizeMesh.C, used by all FOUR SpaceGroup overload defaults (the item
  said three; it was four) plus the GPW `CreateXCQuadrature` site.  They all decide the same question —
  "are these two mesh points the same point?" — so they must agree, and a reader should not have to diff
  four default arguments to confirm that they do.

- **R2.13 ✅ DONE 2026-08-07. Becke strings/labels rename in `Delta_*`/`XC_GridEngine`.**  Verified: the classes are
  already behaviorally mesh-neutral — ZERO branching on Becke; only 3 `Write()` strings + 3
  profiler labels (the one factory branch lives in Imp/Hamiltonians.C:201 where policy belongs).
  Rename "becke" → "XC mesh"/"XC quadrature".  (The REAL non-neutrality — the
  `Symmetry::Lattice_3D::Fold` + dcmplx dependency — is V1.5's problem.)
  Six renamed (3 `Write()` strings + 3 profiler-label occurrences, as the item predicted) → "XC-mesh".
  **Deliberately NOT renamed: the `[Becke grid]` console lines** in `Imp/UnitCell.C` and `Imp/GPW_IBS.C`.
  Those belong to the mesh BUILDER, where Becke is genuinely the scheme being built — an accurate name, not
  a leaked assumption.  The rename targets classes that are mesh-NEUTRAL but were labelled Becke; it should
  not strip the name from the one place it is true.

- **V1.4 ✅ DONE `80fc2ae8`. `DM_RhoAtPoints` Phi key → Irrep (USER RULING 2026-08-05).**
  User model: a CompositeCD is a list of IrrepCDs — not BasisSetCDs; one irrep points to one IBS.
  `Irrep` has meaning in the real world outside of code; `BasisSetID` is purely a code construct,
  whose job is the DB cache's CROSS-RUN caching (same irrep, different radial functions) — caching
  grids inside one SCF run is a different concern and should not borrow its key.  Verified
  enablers: `Irrep` is designed as a std::map key (`operator<`, src/Symmetry/Irrep.C:25) and every
  IBS exposes it (`IrrepBasisSet::GetIrrep(const Spin&)`, IrrepBasisSet.C:62).  The spin detail
  is resolved (user, 2026-08-05): **Spin=None is a valid spin state** — key on
  `GetIrrep(Spin::None)`, the spatial irrep.  → **READY**: change `XC_GridEngine::itsPhi` +
  `DM_RhoAtPoints` signature from `map<std::string,...>` to `map<Irrep,...>`.

- **V1.9 ✅ DONE `38a1ebd6`. `Structure`→concrete-`UnitCell` down-casts in 4 libraries**
  (SCFIterator.C:163 [+ R1.3 slicing], DensityMixer.C:319, Imp/SeedCD.C:26,94).  UnitCell.C:38-41
  itself documents the pattern as the thing to avoid.
  **USER FIX SHAPE (2026-08-05): free pry-out HELPERS, exactly like the atom-symmetry precedent.**
  `Symmetry::Atom::Getl(const Symmetry&)` (src/Symmetry/Atom/Spherical.C:68-73) exists because l is
  needed literally *everywhere*; the free helpers downcast the abstract base to the atomic concrete
  and throw `std::bad_cast` on mismatch — "casting and error handling in ONE place" and readable
  client code.  Do the same for Lattice_3D: `GetReciprocalCell(const Structure*)` /
  `GetReciprocalLattice(const Structure*)` (+ others as they show up — cell volume, frac↔cart).
  Semantically identical to the cast, but the 4 sites become one-liners naming what they WANT.
  Verified demand (the R1.3 answer): every consumer of the Kerker cell snapshot wants exactly
  "give me your reciprocal lattice" — `CreateVxcFitBasisSet` takes the plain `Structure*`, and the
  only concrete use is `MakeReciprocalCell()` → `ReciprocalLattice` (DensityMixer.C:327-328).
  Open (minor): helper-in-a-lib placement — the helpers need `UnitCell` visible, so they live
  wherever `Getl`'s analogue would (qcStructure or a Lattice_3D helper module), NOT in the
  consuming libraries.  Note this is the pragmatic sibling of the abstract periodic-geometry
  capability face; the helpers can land FIRST and the face later behind them (the helper body is
  the only thing that changes).

- **V1.10b ✅ DONE (see LANDED). Mixer
  LAYERING: `SCFIterator` should only ever see `tDensityMixer` (USER DESIGN,
  2026-08-05).**  Today the iterator knows about Kerker specifically: it holds a Kerker-only
  `itsKerkerCell` snapshot (SCFIterator.C:182, built by the R1.3 cast at :163) and calls
  `MakeDensityMixer(..., itsBS, itsKerkerCell.get(), itsCD.get())` at :244 — i.e. mixer CREATION
  (and its periodic-basis/cell/seed feasibility probe, DensityMixer.C:319-328, which `cerr`s
  "[Mixer] DISABLED" when the probe fails) lives inside the structure-neutral iterator.  Correct
  placement: build the mixer ABOVE the iterator (CalculateSolid / RunGPW facade level) and inject
  the `tDensityMixer*`, or resolve it in a `SolidSCFIterator` via virtual dispatch.  **This
  SUBSUMES R1.3 and one of V1.9's four casts**: with the mixer arriving pre-built, `itsKerkerCell`
  and its slicing deep-copy disappear entirely (the copy exists only to dodge the dangling
  `Lattice_3D::GetStructure` temporary — a problem the facade doesn't have).  Also relocates the
  DISABLED-fallback decision to a layer that can report it properly.

- **V1.13 ✅ DONE 2026-08-07 — executed as the compiler-verified DELETION R2.6 made possible.  `ExFunctional` — data + setters on an abstract face, with a hidden-init landmine.**
  **What went:** the `ScalarFunction<double>` base (the FIELD face), `operator()`/`Gradient` on all three
  implementations (SlaterExchange, VWN_Correlation, Libxc_LDA), `InsertChargeDensity`,
  `itsChargeDensity`, `SetPolarized`, `isPolarized` — and with them the whole TU
  `Internal/Imp/ExchangeFunctional.C` (it held only the ctor that initialised the two dead members) plus
  five now-dead imports.  `ExFunctional` is now a DATA-FREE value face: `GetVxc`, `GetEpsXc`,
  `GridCutoffFactor`.
  - The claim "the field face is dead" was checked BY THE COMPILER, not by inspection: deleting it built
    clean, so nothing needed it.
  - **Why it was safe to delete rather than reimplement:** every one of the three implementations defined
    `operator()(r)` as literally `GetVxc((*itsChargeDensity)(r))` — the value face composed with a density
    the object should never have owned.  A field is now built where it is USED, by adapters taking
    (functional, density) as CONSTRUCTOR arguments: `VxcDensity`/`EpsXcDensity` (Imp/FittedVxc.C),
    `PolVcDensity`/`PolEpsCDensity` (Imp/FittedVcorrPol.C), `PWVxcField` (Imp/PWTerms.C).  So the density
    arrives as an argument, and the hidden-init landmine cannot recur.
  - `SlaterExchange::Gradient` went too — it was the only real implementation of the field face, unused,
    and the one place that read `isPolarized` (permanently true, so its unpol branch was dead code).
  - Polarization is now expressed where it belongs: the spin-native `SpinCorrelation` face, and
    `SlaterExchange::itsSpin` (which is what `GetVxc` ALREADY branched on — the bool and the flag were two
    different answers to one question).
  *(Original analysis:)*
  Carries `itsChargeDensity*` + `isPolarized` (violates the data-free-interface convention).  Its
  ScalarFunction face dereferences `itsChargeDensity`, set ONLY by `LDAVxc::UseChargeDensity` — on
  the PW/GPW path it is NULL and `op()(r)` would segfault.  `SetPolarized` has NO caller ⇒
  `isPolarized` is permanently true and Gradient's unpol branch is dead; meanwhile GetVxc branches
  on a DIFFERENT flag (SlaterExchange::itsSpin).  Split the value face (GetVxc/GetEpsXc(ρ)) from
  the field face; polarized = a type, not a bool (the `SpinCorrelation` face below it is the
  correct data-free shape).  ~~(Coordinate with R2.6 — LDAVxc is the only setter caller.)~~
  **PROMOTED after R2.6 landed (2026-08-07): this is now a DELETION, not a redesign.**  R2.6 removed the
  ONLY caller of `InsertChargeDensity` (it was `LDAVxc::UseChargeDensity`), so `itsChargeDensity` is never
  set and the entire FIELD face — `ExFunctional::operator()`, `Gradient`, the member — is dead tree-wide.
  The fitter now reaches the field through the `VxcDensity`/`EpsXcDensity` adapters, which take the density
  as a CTOR ARGUMENT, so the hidden init cannot recur.  Executing it = delete the member + the two
  virtuals + `InsertChargeDensity`, drop `ExFunctional`'s `ScalarFunction<double>` base, and retire
  `SlaterExchange::Gradient` (its only reason to exist was that face).  The COMPILER verifies the claim:
  anything still needing the field face fails to build.  `SetPolarized`/`isPolarized` (still callerless,
  hence permanently true) rides along in the same pass.

- **V2.4 ✅ DONE 2026-08-08 — margin validated, selector ARMED.**  Converged-run A/B on both systems the
  selector routes to uniform (`GPW_SCF.DISABLED_GridRouteAB_{SiGamma,AlFCC}`; scored on the converged DENSITY
  per D8, never ΔE_total).  Uniform at its own cutoff matched a fine Becke reference at least as well as the
  PRODUCTION Becke mesh, for 4–15x fewer points — and much better on total energy (Si 30x, Al 14x), because
  B-prod's residual is the angular error V2.6 measured.  `U-sel ≡ U-4x` on both, so the uniform route is
  genuinely converged at the selector's cutoff.  `kUniformMargin=2.0` stands.  `Auto` now ACTS
  (`GPW_XCGRID_NOSELECT=1` is the A/B valve); 683/683 green.
  - **No contradiction with the 2026-08-01 Becke default**, which was justified for diffuse bases and sharp
    cores: those systems (F α_max=40, MnO) are exactly the ones the selector already routes to Becke.  The
    selector reproduces that split automatically instead of applying one answer to both regimes.
  - **One anchor moved, and to a value already in the tree:** `AlFCCMetalIBZExact` −2.1174805 → **−2.1169707**,
    which is verbatim the "uniform-route pair" its own comment had recorded since 2026-08-02.  The re-pin
    switches which of two long-known values is the default; it introduces no new quantity.  IBZ folding stays
    exact (`AlFCCMetalGlobalMu` prints the same value).
  - **⚠️ ARMING EXPOSED A REAL BUG, which is the main find.**  `VxcFit::Auto` paired uniform with the
    PLANE-WAVE fit unconditionally — and `PWFittedVxc` is not spin-native, so the two Si spin-collapse gates
    (`Polarized{,Seed}SingletMatchesUnpolarizedSiGamma`) threw the moment Auto could route a soft system to
    uniform.  The throw message named its own fix: Delta works on EITHER grid, and the code already said so.
    `Auto` now picks Delta whenever PW cannot do the job — Becke grid OR polarized run — so the throw is
    reachable only via an explicit `VxcFit::PlaneWave`.  **This was latent before V2.4, not caused by it:**
    a polarized uniform-grid run was simply unreachable while Auto always chose Becke.
  - **Recorded honestly: the one place Becke still wins is \f$\|\Delta\rho\|_\infty\f$** — 4.7x better at
    the WORST point on Si, which is the core.  A property that samples the core (hyperfine, EFG, core-level
    shifts) should prefer Becke even where the selector says uniform.  The integrated norm is not the whole
    story, and the selector optimises the integrated norm.
  - Method note that changed the answer: the first Si A/B had 3 of 4 arms hit FIT-FLOOR STALL, because
    `imposeSymmetry=false` was copied from the V2.6 ladders — right there (measure the bare angular rule),
    wrong here (free-mesh Si/Γ oscillates in its degenerate manifold, so the densities were not converged and
    the scores read 20–40x larger).  A converged-density metric means nothing until every arm converges.
  *(original item text:)*  **Calibrate `kUniformMargin` and ARM the V1.26 selector.**  The cost model says uniform beats Becke
  15x on Si; the 2026-08-01 measurement says Becke is the right grid there.  Both cannot be right, and the
  margin (a guess at 2.0) is where the discrepancy has to be absorbed.  Instrument: the
  `[XC grid choice] Auto:` line, now emitted by every GPW run.  Measurement per the D8 pin — grid-convergence
  of ρ/property against a fine reference, NEVER ΔE_total.  Then update
  `XCPolicy.AutoStillResolvesToBeckeWhileTheSelectorIsDisarmed`, deliberately, and arm.


- **V2.6a ⛔ ATTEMPTED AND REJECTED 2026-08-07 — flip `angularDegree` 29 → 17.**  Made the one-line change,
  ran the suite, and BACKED IT OUT on the evidence.  Default stays 29.  Three things came back:
  1. **A FOURTH system refutes the three-system recommendation: Al FCC (simple metal).**  Both Al anchors
     moved **6.4e−4** (6x `AlFCCMetalIBZExact`'s tolerance) — `−2.1174805` → `−2.11812`.  Checked the obvious
     suspicion first and it was wrong: the full-mesh sibling `AlFCCMetalGlobalMu` prints the SAME `−2.11812`,
     so IBZ folding stays exact and the site-adapted mesh did not degrade — both routes moved together.  It
     is a genuine XC shift, and Al is simply harder than every system in the ladder set.
     - Ladder confirms it, and Al is the worst case on BOTH axes: `max|dVxc|` at GL-15 = 2.4e−3 (out of
       tolerance), GL-17 = 3.9e−4, GL-29 = 2.6e−4; radial nR=30 = 9.3e−3 (out), nR=40 = 2.6e−4, nR=60 = 9.3e−5.
     - **And it is NON-MONOTONIC**: GL-11 (1.8e−2) is worse than GL-9 (8.3e−3); GL-23 (5.2e−4) worse than
       GL-17 (3.9e−4).  So no clean "degree N suffices" statement survives Al at all — the error rattles
       around 3–5e−4 from 17 up to 29.
     - **Why a metal is the worst case is exactly the mechanism V2.6 identified**: the angular requirement is
       set by the fuzzy-Voronoi PARTITION SURFACE, and a nearly-free-electron valence density is the one that
       puts substantial charge ON that surface.  Si/NaF/Mn all concentrate charge near nuclei.  So the theory
       held; the sample just didn't contain its own worst case.
  2. **`AlFCCDegenerateShellAufbauStalls` flipped QUALITATIVELY: it now converges, and that is a FAILURE.**
     That test asserts an integer-aufbau degenerate 3p shell CANNOT converge Δρ (the density rotates freely
     in the degenerate manifold).  A coarser angular grid has larger orientation-dependent quadrature error,
     which PINS the rotation — so the run "converges" to a state held in place by grid error.  Re-pinning
     that `EXPECT_FALSE` to `EXPECT_TRUE` would have recorded a grid artifact as a physics result.  Same
     mechanism as the `SiPseudoAtomInBoxMatchesFinite` caveat under V1.26, now biting from the other side.
  3. **⚠️ THE METHODOLOGICAL FINDING, and the most transferable thing in this whole item: a frozen-density
     quadrature error UNDERSTATES the self-consistent shift on a METAL.**  Al at degree 17 measures
     `max|dVxc|=3.9e-4`, comfortably inside the gate — yet its SCF total moved **6.4e−4**.  Grid error feeds
     back through the density, the Fermi level and the occupations; on an insulator the two numbers agree
     (Si/NaF anchors did not move), across a Fermi surface they do not.  **A ladder measures the QUADRATURE;
     for a metal that is a LOWER BOUND on what the SCF does with it.**  The ladder method (D8-compliant, and
     still the right instrument) simply cannot see this — it needs a converged-run A/B alongside.
  - **Where this leaves the recommendation:** degree 17 is defensible for INSULATORS and is not safe as a
     global default.  The honest options are (a) keep 29, (b) make the degree adaptive (metal-vs-insulator is
     already known at the facade — `globalFermi`/smearing), or (c) re-measure with converged-run A/Bs rather
     than frozen-density ladders.  (b) is the interesting one and is genuinely new work.
  - **Consequence for V2.4: it is UNBLOCKED, not gated.**  The Becke side stays at 450 directions, so every
     crossover number recorded under V1.26 stands as written — no refit needed.


- **D6 ✅ DONE 2026-08-07. `BeckeXCParams()` lives in the TEST file + `ResolveXCMesh` (test driver)** — the
  de-facto PRODUCTION Becke recipe and the run-policy resolution (Auto grid × imposeSymmetry interplay)
  living in the integration-test harness; both belong with the facade/driver once the policy
  object exists (library beside `MeshParams`; tests read it, not the other way round).
  **Both moved verbatim into `qchem.Mesh` (src/Mesh/Mesh.C), declared immediately after `MeshParams`;
  the anonymous-namespace copies in `IntegrationTests/GPW_SCF_UT.C` are gone and its 15 call sites now
  spell `qcMesh::`.  670/670 ctest green, every GPW anchor unmoved (the moved bodies are logic-identical).**
  - **The reason it did NOT need the policy object first — worth recording, because the item's own
    wording ("once the policy object exists") assumed it would.**  `ResolveXCMesh` consults **NO run
    context** any more: its last context-dependent branch was the BZ-reduced carve-out, retired
    2026-08-02 once the site-adapted invariant mesh was gate-verified.  What was left is a
    context-FREE default (`Auto` → the calibrated recipe), which is exactly the kind of thing that
    can sit beside the enum it resolves.  Same shape as the session's standing heuristic: the thing
    that looked like it needed extra machinery turned out to need less, not more.
  - **Why `qchem.Mesh` and not `qcHamiltonian`:** `UnitCellKind::Auto` is DECLARED in `qcMesh`, and its
    own doc comment said "a POLICY layer that knows the run context resolves it" — i.e. the library
    shipped a value only a test file could interpret.  An enum that offers `Auto` should ship the
    canonical resolution of `Auto`; that comment is now updated to point at `ResolveXCMesh`.
  - **The `GPW_BECKE_*` env instruments moved WITH the recipe, deliberately** — `getenv` in library code
    is already the house idiom for GPW sweep knobs (`GPW_XCROUTE`, `GPW_STREAM_FOLD`, `GPW_RELCUTOFF`,
    `GPW_OMP_THREADS`, …, ~25 sites in `src/`).  Splitting them off would have BROKEN the instrument:
    most runs reach the recipe through `Auto`, so a library resolver returning the un-overridden default
    would make `GPW_BECKE_L`/`_NR`/`_ALPHA`/`_ANG`/`_ROT` no-ops for exactly the runs they exist to sweep.
    Verified live from the new home: `GPW_BECKE_NR=12` moves the `[Becke grid]` line from nR=40 to nR=12.
  - Not changed (and NOT a regression of this move): under `imposeSymmetry` the site-adapted builder still
    silently replaces `mp.angular`, so `GPW_BECKE_ANG` has no visible effect on an imposed run — that is
    D5's complaint, unchanged.
  - **FOLLOW-UP RULING, same day → V1.26: `ResolveXCMesh` is not finished, it is SEEDED.**  The user's answer
    to "who decides Uniform vs Becke" is *the user, at the solid facade* — with `Auto` becoming a COST
    SELECTOR (uniform \f$n^3\f$ vs the Becke point count) rather than the unconditional `Auto`→Becke landed
    here, and a diagnostic when an explicit choice is the losing one.  So this function acquires a real body;
    the MOVE done here is what makes that a one-place library edit instead of a test-file edit.
  - **Deliberately NOT done: resolving `Auto` inside `Ham_PW_DFT::BuildTerms`.**  It is tempting (an
    unresolved `Auto` silently reads as `Uniform` downstream, since every consumer tests `==Becke`), but
    `xcMesh` reaches several Hamiltonian lineages and the resolution point is a facade decision, not a
    term-builder one.  Flagged in the `UnitCellKind` doc instead: resolve where the mesh spec ENTERS the
    Hamiltonian.  Revisit with the `SymmetryPolicy`/facade pass that D5 also waits on.

## DONE (record)

- **BZ creep on the neutral `Symmetry` base (ISP) — DONE 2026-08-04 (same day): the ATOM SHELL
  CONVENTION landed.**  `BlochQN::GetDegeneracy()` = star (ctor converts the k-mesh layer's w_k;
  asserted integer), `GetWeight()` = uniform 1/N_mesh, `StarSize()` + the `EnergyLevel` report
  scaling DELETED (the 8/8 wedge display is now the physical occupation), `Crystal_EC::GetN` =
  star×Nval per block in insulator mode / plain Nval for the global-μ total.  No factory/KBlock/
  call-site changes (BlochFactory keeps taking w_k); free meshes (star=1) bit-identical; g=w·degen
  invariance verified (Al metal); imposed Si diamond IBZ reproduces E to all digits; 656/656.
  The base is down to TWO Bloch-flavored members (GetWeight, MergeAcrossIrreps) — the remaining
  capability-face split is now small enough to fold into any future Symmetry-base touch.
  *(Original analysis, for the record:)*  The base carried three Bloch-flavored defaulted virtuals
  (`GetWeight()` pre-existing; `MergeAcrossIrreps()` + `StarSize()` from the 2026-08-03
  MergeTol/IBZ-table fixes).  Rather than split a capability face, the user's design was adopted:
  a wedge k-block is a SHELL exactly like an atom's l-block (one stored representative,
  degeneracy = the symmetry copies), so **`BlochQN::GetDegeneracy()` = star size** (spatial; spin
  stays layered by `Irrep`) and `StarSize()` is DELETED.  A coordinated convention switch — the
  invariant is w_k × (per-block quantity), so the star factor must leave the weight in the same
  commit: `GetWeight()` becomes plain 1/N_mesh (most of the weight plumbing through `ReduceToIBZ`
  → `KBlock` → `BlochQN` evaporates); `Crystal_EC::GetN` per wedge block becomes star×Nelec (each
  band holds star×2 natively; the level table needs NO report-layer scaling); every w_k consumer
  rebalanced (density sums, `tIrrepWF::GetEntropyTerm`; `FillOrbitalsGlobalFermi`'s g = w·degen is
  INVARIANT — sanity check).  `MergeAcrossIrreps()` stays (unfolded meshes still carry star
  partners as separate blocks).

## Answered questions (from the pre-reorg sections)

- **"XC_GridEngine — I don't know what this is."**  It is the shared XC quadrature engine: holds
  the mesh, the fold, the geometry-fixed Φ tables, and the ρ caches, SHARED via one `shared_ptr`
  across Delta_XC / Delta_XC_Pol / Delta_VcorrPol (Imp/Hamiltonians.C:217) so the Φ tables are
  built once.  "Engine" earns its keep through the sharing — merging it into Delta_XC would
  triple the Φ tables.  Rename candidate: `XC_Quadrature`.  The "Becke route is pure QUADRATURE"
  phrasing = fit basis is (pseudo) delta functions, per the user's note.
- **"PW_XC: what distinguishes it from FittedVxc?"**  Genuinely different mechanics, not just
  PWs: raster ρ by inverse FFT per density serial, the RAW-adjoint vs BALL-fit fork, energy by
  direct grid quadrature (no ε fit, no DM_Contract).  And yes, GPW uses it: there is no Ham_GPW —
  `Ham_PW_DFT` + GPW basis with `UnitCellKind::Uniform` instantiates PW_XC
  (Imp/Hamiltonians.C:246-249).
- **"Band_DFT_IBS — is this no longer used?"**  Dead in code; deliberately kept — see D1.
- **"BandStructure.C only consumers seem to be unit tests."**  Confirmed — see V1.21.

---

# Appendix — V1.2 Orbital_PP_IBS feasibility probe (2026-08-05, read-only)

**Verdict: NO hard blocker.  The inversion works; four real roadblocks (all with mitigations)
and two minor ones.**

What the basis side ACTUALLY consumes from `LocalPotential`/`SeparablePotential` (verified — no
implementation reaches for rloc/C1..C4/h_ij/Zion):
- Local: `FormFactor{,Long,Short}(Z,G2)`, `FormFactorG0Short(Z)`, `Vloc(Z,r)`, and the optional
  `ShortRangeGaussian(Z) -> {c,n,alpha}` cross-cast face.
- Separable: `NumProjectors/AngularMomentum/Coefficient` (Coefficient is already a DIAGONALIZED
  per-projector scalar — no D-matrix reaches the basis), `Projector(Z,p,q)`, `BetaR(Z,p,r)`, and
  the optional `BetaGaussian` face.
- `Zion` is consumed only term-side (IonIon).

Roadblocks, ranked:
1. **The analytic-FT fast path is load-bearing, not optional.**  PW has NO real-space local path
  at all, and GPW's G-space local route is a BOX-INDEPENDENCE requirement (Evaluator.C:836-841),
  not a speed hack.  An opaque-radial-function argument is a non-starter.  Mitigation: the
  neutral field argument is DUAL-SPECTRAL — `ValueR(Z,r,range)` AND `ValueQ(Z,q2,range)` — which
  is what `LocalPotential_Q` already is, renamed (the TYPES are already neutral; only the NAMES
  are PP-flavored).
2. **Three optional capability faces select the fast paths via dynamic_cast**
  (`LocalPotential_Gaussian`, `SeparablePotential_R`, `SeparablePotential_Gaussian`) — re-express
  as neutral optional faces; preserve the `MultiSpecies_*` forwarding.  Bulk of the mechanical work.
3. **The G=0 alignment leaks into the basis** — `FormFactorG0Short` is read basis-side
  (Evaluator.C:1039) to keep the analytic lattice sum consistent with the drop-ΔG=0 convention.
  Neutral face needs a `CellMeanQ(Z,range)` scalar or the two GPW local paths silently disagree.
4. **The β sharpness hint is a parameter leak** — `ShortRangeGaussian(Z)[0].alpha` pulled purely
  to size grid levels (Evaluator.C:901-906, the NaF diving-ghost note).  Neutral face needs an
  explicit `EffectiveExponent(Z,range)` (0 = unknown → old behavior, already the fallback).
5. (minor) **Projector BRACKETS are the wrong primitive for PW** — `MakeSeparablePotential` uses
  the (2l+1)P_l(cosγ) addition theorem to avoid the m-loop and never forms ⟨i|β_lm⟩ vectors.
  Keep the projector primitive at MATRIX level (`MakeProjectorMatrix(st, projectorSet)`, D
  inside); GPW's bracket loop stays an implementation detail.
6. (minor) `Integrals_Pseudo::MakeLocalPotential` (unsplit) is unit-test-only — drop, don't port
  (re-plumb its 7 test call sites).

DAG findings (better than hoped):
- qcPseudopotential imports NOTHING from qcBasisSet (only qcMath/qcCommon/qcStructure) — the
  inversion edge qcPseudopotential→qcBasisSet is cycle-free.  BUT the adapter that converts
  `LocalPotential` → the neutral field argument can live in **qcHamiltonian** (which already
  links both and owns the terms) — then qcPseudopotential stays a pure LEAF with zero new edges.
  Strictly better DAG.
- The face + argument interfaces → qcBasisSet; the pure-data `RadialGaussianTerm{c,n,alpha}`
  struct → qcMath (already hosts Monomial/CartTerm, links nothing).  Do NOT reuse
  `Molecule::LatticeSum1E::GaussianFunction` as the public type (qcMolecule_BS links UP to
  qcBasisSet → cycle).
- After the retarget, `PWTerms.C`'s two casts go to `BasisSet::Orbital_PP_IBS<dcmplx>` and
  qcPseudopotential drops off qcLattice_BS's link line (modulo the two test TUs importing
  GTH_Potentials directly).

Structure-neutrality check: the primitives are NOT inherently periodic — `PP_Local`/`PP_NonLocal`
compute exactly ⟨i|V_species|j⟩ and the D|β⟩⟨β| loop by mesh quadrature, so a molecular basis can
implement the same face later.  Three PERIODIC CONVENTIONS must be parameters, not baked in:
(1) dropping ΔG=0; (2) the cell-mean subtraction (already gated on `isFinite()`); (3) the
long/short RANGE SPLIT itself (exists to route V_long into the Hartree Poisson — a molecular
implementor answers `Full` and ignores it, exactly what `FormFactorShort`'s default-0 already does).

Minimal neutral signature sketch (2 virtuals instead of 4 — Long/Short/Full collapse into a
`FieldRange` argument):

```cpp
// qcMath — pure data, no deps
struct RadialGaussianTerm { double c; int n; double alpha; };   // c r^{l+2n} e^{-alpha r^2}

// qcBasisSet
enum class FieldRange { Full, Long, Short };                    // a range split, NOT a PP concept

class SpeciesRadialField {
  virtual double ValueR   (int Z, double r,  FieldRange) const =0;  // <- Vloc
  virtual double ValueQ   (int Z, double q2, FieldRange) const =0;  // <- FormFactorLong/Short (fast path kept)
  virtual double CellMeanQ(int Z,            FieldRange) const {return 0;}  // <- FormFactorG0*
  virtual double EffectiveExponent(int Z,    FieldRange) const {return 0;}  // <- the beta grid hint
};
class SpeciesRadialField_Gaussian   // optional capability, cross-cast
  { virtual std::vector<Math::RadialGaussianTerm> AsGaussians(int Z, FieldRange) const =0; };

class SpeciesProjectorSet {
  virtual size_t Count (int Z) const =0;
  virtual double Weight(int Z, size_t p) const =0;              // <- Coefficient (scalar; no D-matrix)
  virtual int    L     (int Z, size_t p) const =0;
  virtual double RadialQ(int Z, size_t p, double q) const =0;   // <- Projector(q)
};
class SpeciesProjectorSet_R        { virtual double RadialR(int Z, size_t p, double r) const =0; };
class SpeciesProjectorSet_Gaussian { virtual std::vector<Math::RadialGaussianTerm> AsGaussians(int Z, size_t p) const =0; };

template<class T> class Orbital_PP_IBS {                        // the face, in qcBasisSet
  virtual hmat_t<T> MakeSpeciesFieldMatrix(const Structure*, const SpeciesRadialField&, FieldRange) const =0;
  virtual hmat_t<T> MakeProjectorMatrix   (const Structure*, const SpeciesProjectorSet&)            const =0;
};
```

`HGH_LocalPotential`/`HGH_SeparablePotential` then simply ALSO implement the neutral faces
(rename-level: `FormFactorLong(Z,G2)` → `ValueQ(Z,G2,Long)`).

---

# HARVEST 2026-09-08 — fourteen items closed in place but never moved

These carried their DONE verdict inside `doc/CleanupCandidates.md` and were never harvested, so the
worklist held **1165 lines of closed record** among 34 open items — the reason it had become unreadable.
Moved here verbatim, each leaving the usual one-line stub.  Where an item had a SURVIVING remainder, the
remainder stayed behind in the worklist (it is open work) and only the closed record moved.

---

## R1.9 — harvested 2026-09-08 (was CleanupCandidates.md L1123-1162)

*(verbatim)*

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

---

## R2.21a — harvested 2026-09-08 (was CleanupCandidates.md L1165-1203)

*(verbatim)*

- **R2.21 ✅ FOUND AND FIXED `9da2e825` (2026-08-24). `blaze::conj` IS A NO-OP ON A COMPLEX SCALAR — and it
  had silently corrupted the pivoted-Cholesky factor.**  Found by the R1.0c gate on its first run.

  **THE TRAP.**  `blaze::conj` conjugates a MATRIX or a VECTOR but is the **identity** on a
  `std::complex` scalar.  Measured directly against this submodule (clang 21):

  | expression | result |
  |---|---|
  | `blaze::conj(M)(0,0)` | `(1,-2)` ✅ |
  | `blaze::conj(M(0,0))` | `(1, 2)` ⛔ — conjugating an ELEMENT |
  | `std::conj(z)` | `(1,-2)` ✅ |

  So the natural-looking `blazem::conj(M(i,j))` compiles, reads correctly, and does nothing.

  **WHAT IT BROKE.**  `LowRankFactor` built \f$L=PU^{T}\f$ instead of \f$PU^{\dagger}\f$, hence
  \f$LL^{\dagger}=\bar D=D^{T}\f$, and
  \f[ \rho'-\rho=\sum_{i<j}4\,\mathrm{Im}(D_{ij})\,\mathrm{Im}\!\left(\Phi_{gi}\overline{\Phi_{gj}}\right). \f]
  **Both factors must be nonzero to see it**, and that is the whole story of why it survived: the second
  vanishes for REAL \f$\Phi\f$ — every Γ/TRIM block, which is all the enabled suite ever ran. At
  \f$k=(1/4,0,0)\f$ on a 4×4×4 mesh it is wrong by **~1e-1 relative in \f$\rho\f$**.
  - **It was LIVE production code**: `QCHEM_DM_LOWRANK` defaults ON, so `FactoredRho` is the default for
    every Bloch block.
  - **And it was silent.**  The existing PSD guard compares \f$\mathrm{Tr}(LL^\dagger)\f$ to
    \f$\mathrm{Tr}(D)\f$ through \f$|L|^2\f$ — and a conjugation error does not move a modulus.  A guard
    can only catch what it is sensitive to.

  **FIX.**  `blazem::conjs(x)` — an explicitly OVERLOADED scalar conjugate (`double`→`double`, so the real
  case cannot pick up a `std::complex` return) added to `qchem.Blaze` with the trap documented on it.  Two
  call sites corrected: `LowRankFactor` (the bug) and the `[DM subspace]` diagnostic (a wrong printed
  overlap only).  **A survey of the other seven `blazem::conj` uses in the tree found them all to be
  matrix/vector expressions, which are correct** — so the blast radius was exactly these two.

  ⚠ **THE TRANSFERABLE LESSON, and it is about the GATE, not the bug.**  The pre-existing complex rig could
  not have caught this no matter what it asserted, because every k on a 2×2×2 mesh is TRIM
  (\f$2k\f$ is a reciprocal lattice vector), so its "complex" block carried REAL \f$\Phi\f$ in a complex
  matrix.  **A complex TYPE is not a complex VALUE**: a gate meant to exercise complex arithmetic has to be
  built at a geometry where the numbers are actually complex, and that has to be stated in the test rather
  than assumed from the scalar.  The gate now says so, and names the k it uses and why.

---

## R2.18 — harvested 2026-09-08 (was CleanupCandidates.md L1316-1357)

*(verbatim)*

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
  - **ONE STALE COMMENT LEFT BEHIND ON PURPOSE — ✅ GONE BY DRIFT (verified 2026-09-06: `grep -rn
    "Vxc::CalcMatrix" src/` finds nothing).**  It said `src/SCFIterator/Imp/SCFIterator.C:186` still reads
    "Vxc::CalcMatrix"; it was a comment only, and SCFIterator was on the MnO campaign's DO-NOT-TOUCH list,
    so it was NOT edited — reaching into their working set for a comment is not worth a collision.  The
    list is long released and some later edit swept it; nothing to do.
  - Found while writing up R2.9(ii), where "return by value" was considered as a fix to `GetMatrix` — the
    convention says that was the wrong half of the pair to reach for.

---

## R2.19 — harvested 2026-09-08 (was CleanupCandidates.md L1358-1386)

*(verbatim)*

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

---

## R2.21b — harvested 2026-09-08 (was CleanupCandidates.md L1392-1509)

*(verbatim)*

- **R2.21 ✅ DONE 2026-08-17 (concurrent-cleanup session), BOTH halves — and it was NOT optional after
  all: the state half is what unblocked MOM on a real TRIM block.**  `OccupationState` landed as a scalar-INDEPENDENT
  persistent ledger (per-block MOM references, fill clocks, cross-irrep arming, the −TS aggregate), owned
  by the SCFIterator beside its policy; `OccupationPolicy<T>` is built over it and keeps only the
  decision.  **The load-bearing detail: each block's reference is stored under the BLOCK's own scalar (a
  variant), not the run's.**  That is what a mixed real/complex mesh needs — before it, a real TRIM block
  in a complex run had nowhere to put a `mat_t<double>` reference and its capture THREW
  ("run with forceComplex"), i.e. switching MOM on turned into a mid-SCF failure after the basis, the
  seed and the first Fock were already paid for.  Measured after: Si (3,1,1) mixed mesh with MOM on
  converges REAL and matches its forced-complex twin to **ΔE = 1.8e-15**.  Two things fell out rather
  than being built: `HeldOccupationPolicy` stopped wrapping the run policy (the shared clocks and the −TS
  aggregate are the state's, so only its three DECISIONS remain), and `RealBlockFillView` — the
  forwarding adapter whose capture did the throwing — **deleted**, since a real block now takes a genuine
  `OccupationPolicy<double>` over the shared state.  New pins: `OccupationState.*` (5 unit tests in
  `UTElConfig`).  746/746 green.
  **SHAPE HALF also landed (same session):** `OccupationPolicy<T>` is now genuinely ABSTRACT and the four
  behaviours `Configure` used to select between are four pairs of objects — occupancy
  {`IntegerOccupancy`, `FermiOccupancy(kT)`} × ranking {`BareRanking` (a null OBJECT, not a "MOM off"
  flag), `MOMRanking(startIter,Λ)`} — assembled by `MakeOccupationPolicy<T>(OccupationConfig, state&)`
  once per Iterate, so `kT>0` is answered by WHICH OBJECT EXISTS instead of being re-branched inside every
  fill.  **`Configure` is dead.**  `HeldOccupationPolicy` is Integer×Bare plus `HoldsStoredBlocks` — a
  sibling, not a wrapper.  The two end-of-fill calls `tIrrepWF` made collapsed into one `OnBlockFilled`
  hook (the shared-μ metal path deliberately keeps a bare `CountFill`: its μ was solved on bare ε, so
  capturing a reference from it would snapshot a subspace the ranking never shaped — documented at the
  declaration).  **Configuration is a VALUE** (`OccupationConfig`), which is what lets the mixed-mesh
  cross arm rebuild the run's own policy one scalar over with a single Factory call.  Deviation from the
  item's sketch, with reason: the Factory takes `OccupationConfig`, not `SCFParams` — `SCFParams` lives in
  qcSCFIterator, ABOVE this library in the DAG, so the iterator converts.  MOM-on-mixed-mesh re-verified
  after the reshape (ΔE = 1.8e-15, unchanged); 746/746.
  **A note for the END-TO-END gate, which this session deliberately did NOT add:** the acceptance run
  above lived in a temporary probe in `GPW_SCF_UT.C` (the real-TRIM session's file) and was removed
  before commit.  Worth adding there at the merge: the (3,1,1) mixed mesh with `UseMOM=true`,
  `MOMStartIter=2`, real vs `forceComplex` twins, `EXPECT_NEAR(..., 1e-9)`.
  *(original item follows)*
  **The OccupationPolicy/OccupationState split (USER 2026-08-17: "I really like your Policy/State
  split.  But we don't need it right now").**  V1.11's landed `OccupationPolicy<T>` is NOT the abstract
  interface D1 ruled — it is a CONCRETE class whose behaviour is selected by
  `Configure(useMOM, momStartIter, kT, momPenalty)`: `DecideBlockFill` branches on `kT>0` and `penalty>0`
  PER FILL, i.e. four policies folded into one object, mode-selected (the user's catch).  Only
  `HeldOccupationPolicy` is a genuine derived policy.  The flag shape was forced by a real constraint —
  the object must survive RECONFIGURATION (grid-continuation adopts MOM references after construction but
  before Iterate; annealed runs call Iterate per stage with DIFFERENT kT), so a rebuild-from-SCFParams
  factory would have dropped the references — and the fix is to split what that conflated:
  - **`OccupationState<T>`** — persistent, in the iterator's slot: the per-block MOM references + fill
    counts, the cross-irrep arming, the −TS aggregate.  `AdoptMOMReference`/`ReleaseReferences`/
    `EntropyTerm` talk to THIS; it lives from construction to the last stage.
  - **`OccupationPolicy<T>` goes genuinely abstract** — `DecideBlockFill`, `HoldsStoredBlocks`,
    `SmearingkT` (the reservoir driver's one remaining ask, default 0), and ONE `OnBlockFilled` hook
    (count + capture-if-due; collapses the two calls `tIrrepWF` makes today).  Concretes = the ruled
    two-axis assembly: occupancy {`Integer`, `Fermi(kT)`} composed with a ranking component
    {`Bare` (the null), `MOM(startIter, Λ, state&)`} — the SCFStrategyPlan null-object idiom, not a
    nullable flag.  `Held(state&)` references the state DIRECTLY (cleaner than today's run-policy wrap).
  - **`Configure` dies.**  A `Factory(SCFParams, state&)` assembles concretes at the top of each Iterate —
    the `kT>0` branch runs ONCE at assembly, never per fill.  The seed fill gets an explicit
    `Integer×Bare` at construction (the D11 semantics, unchanged); annealed stage transitions get a fresh
    assembly against the SAME state — today's semantics exactly, with selection at construction instead of
    interrogation.
  - Blast radius: the policy module, the iterator slot, two `tIrrepWF` call sites.  The WF faces
    (`FillOrbitals(pol,…)`) do NOT change.  Bit-identical by the V1.11 discipline; smeared/metal anchors
    guard did-E-move.  NB the delayed-IMOM adaptivity (MOM behaves as bare aufbau until a reference
    exists) is STATE-dependence, not config — it stays inside the MOM ranking component and never
    justified the kT branch.
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

---

## R2.15 — harvested 2026-09-08 (was CleanupCandidates.md L1581-1709)

*(verbatim)*

- **R2.15 ✅ COMPLETE — DEFAULT FLIP LANDED 2026-08-17 (concurrent-cleanup session), DEGREE-GATED.**
  `BeckeXCParams` now defaults the angular scheme to **Lebedev at degree ≥ 29** (302 vs GL's 450
  directions = 67%), GaussLegendre below — the gate exists because the first, ungated sweep FAILED the
  MnO seed-mirror gate's own degree-11 recipe (Leb-50's ⟨111⟩ orbit dives into neighbour Mn cores:
  orphan w·ρ=0.04 vs the 1e-8 eps-tail contract) — the recorded low-degree alignment poison caught
  LIVE.  Degrees 15–23 stay GL until measured.  The flip's measurement (decision 7's precondition, all
  D8-compliant): Si gate passes with Leb-302 (dExc +7.5e-5 vs GL's +1.1e-4, dVxc equal); NaF internal
  convergence equal (1.8e-5 vs 1.5e-5); Mn angularly trivial (V2.6); **Al metal by converged-run A/B
  per the V2.6a lesson: Leb-302 lands 3× closer to its fine reference than GL-29 to its own** (dEtot
  +6.0e-5 vs +1.9e-4).  Consequence accepted: the cheaper Becke side moves the selector crossover —
  Si/sipp re-routes Uniform→Becke (inside the 2× margin; safe direction).  726/726 green including all
  anchors.  `GPW_BECKE_ANG=gl|lebedev` forces either scheme at any degree (the A/B valve).
  *(the interface half, landed 2026-08-07, and the decision list follow)*
  **INTERFACE DONE 2026-08-07 (`nAngular` → `angularDegree`, behaviour-preserving).
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

---

## V1.3 — harvested 2026-09-08 (was CleanupCandidates.md L1789-1857)

*(verbatim)*

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

---

## V1.19 — harvested 2026-09-08 (was CleanupCandidates.md L1958-1977)

*(verbatim)*

- **V1.19 ✅ VISITOR + THROWS DONE 2026-08-17 (concurrent-cleanup session); ONE deliberate remainder.**
  `Structure::ForEachSite(fn(Z,R,spinFlip))` landed beside its precedent `SumFormFactors` — one place
  reads the atom fields, five consumer loops converted (SAD assembly, `IonicSADTargets`,
  `MagneticDecoration`, the SeedCD ctor incl. its separate anyFlip pass), and SeedCD's point-eval loops
  now use a per-atom-parallel scale table instead of re-asking the structure per point (a map lookup per
  atom per mesh point, gone).  The assert(false)+nullptr arms THROW (an assert-only arm was a silent
  core guess under `-DNDEBUG`).  Bit-identical; 734/734.  **REMAINDER (deferred, recorded at the site):**
  the flip-group sub-cell duplication — removing it needs a per-SITE form-factor overload on the basis
  face, which the item itself weighs against the pseudo-wall pin; that block is now the seed's ONE
  remaining concrete-atom consumer.  *(original text follows)*
  **Seed assembly: give `Structure` the question.**  Seed code reads concrete `Atom` public
  fields (Imp/Seed.C:102-104 `a->itsZ`,`a->itsR`; Imp/SeedCD.C:91 `itsSpinFlip`) and clones the
  UnitCell into (unflipped, flipped) groups because `MakeFourierDensity(st, formFactor(Z,g2))` is
  species-keyed — the flip-group sub-cell duplication is the SYMPTOM; the neutral face yielding
  raw atoms is the CAUSE.  A `ForEachSite(fn(Z,R,flip))` / per-atom form-factor visitor (or the
  per-atom-index `G_FieldEvaluator` overload — weigh against the pseudo-wall pin) removes the
  field access, the UnitCell casts, AND the sub-cell duplication.  Also `MakeSeedDensity`'s closed
  enum switch has assert(false)+nullptr arms (Imp/Seed.C:152,158) — in Release an unsupported
  (strategy × T) cell returns null and SCFIterator silently runs a core guess;
  compile-time-over-runtime says make unsupported combos unconstructible or throw.

---

## V1.26 — harvested 2026-09-08 (was CleanupCandidates.md L2054-2253)

*(verbatim)*

- **V1.26 ✅ COMPLETE (reconciled 2026-08-17; analysis kept below).**  Every deliverable landed across
  three sessions: the cost model + `XCMeshSharpness` + `ResolveXCMesh` with the asymmetric diagnostics
  (D6/V1.26 sessions), the Nyquist bridge (`UniformDivisions`/`UniformCutoff`), the ARMED selector with
  its converged-run-validated margin (V2.4, 2026-08-08), the sized uniform verdict, and the radial
  sibling warning (V2.7, 2026-08-17).  Post-flip note: R2.15's Lebedev default cheapened the Becke side
  33%, moving the crossover — Si/sipp now routes Auto→Becke (inside the 2× margin, the safe direction);
  the mechanism tests were updated with the reason recorded.  *(original ruling + analysis follow)*
  **Uniform-vs-Becke: a SMOOTHNESS question that reduces to a COST CROSSOVER — so `Auto` becomes a
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

---

## V1.28 — harvested 2026-09-08 (was CleanupCandidates.md L2364-2468)

*(verbatim)*

- **V1.28 ✅ RESOLVED BY THE SHUBNIKOV CAMPAIGN S1–S4 (reconciled 2026-08-17; analysis kept below — it
  was the design driver).**  The item's predicted hazard never shipped: S1 (`ShubnikovOps` + the
  anti-translation coset, `5ad45c06`) gave the op set exactly the σ the item said `ReciprocalOp` lacked;
  S2 (`d9fe59fe`) landed the signed (ρ,m) projectors + the flip-fixed audit + `MagneticSymmetryDefects`;
  S3 (`f8fd4fa0`) the decoration→siteSpins→factory resolution; and S4's run 38 is THE CONVERGED MnO
  AFM-II **under imposition** at default knobs — the magnetic star-average CLOSED the tie floor instead
  of erasing the order.  The item's own refinements held up: the staggered-vs-not discriminator became
  S4's "legacy op faces = σ=None subgroup ONLY" pin, and "detected grey is sublattice-preserving
  (erasure unreachable)" is the run-37 live catch.  *(original analysis follows)*
  **⚠️ IMPOSING SYMMETRY ON AN AFM STRUCTURE WOULD DESTROY THE AFM ORDER — the density star-average is
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

---

## V1.29 — harvested 2026-09-08 (was CleanupCandidates.md L2469-2520)

*(verbatim)*

- **V1.29 ✅ RECONCILED 2026-08-17 — the question is ANSWERED, the dependency is DISCHARGED, and the
  workflow itself is deliberately NOT built.**  The literature answer stands in the item (SCF stability
  analysis / broken-symmetry DFT); the "depends on V1.28" line is discharged (Shubnikov S1–S4 landed —
  step 6 is expressible today, and S3's decoration→siteSpins resolution IS steps 5–6 for the
  collinear-known case); and MnO followed the item's own advice ("do not over-build it for MnO: impose
  the KNOWN structure").  What remains unbuilt is the Hessian-based DISCOVERY loop for materials of
  UNKNOWN order — a V4-class watch trigger (build it when such a material enters the queue; the
  irrep-block Davidson sketch below is the design).  *(original text follows)*
  **The spontaneous-symmetry-breaking DISCOVERY workflow — it is an established method, and step 4
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

---

## V1.30 — harvested 2026-09-08 (was CleanupCandidates.md L2521-2546)

*(verbatim)*

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

---

## V2.6 — harvested 2026-09-08 (was CleanupCandidates.md L2759-2860)

*(verbatim)*

- **V2.6 ✅ CLOSED (reconciled 2026-08-17) — the measurement is fully banked and every product landed:**
  the four-system ladder verdict (radial 40 RIGHT, angular 29 kept after V2.6a's refuted flip to 17),
  the V2.4 armed selector, the V2.7 radial-adequacy warning, and R2.15's degree-gated Lebedev flip
  (equal degree, 67% of the directions — the cost cut V2.6 wanted, taken on the SCHEME axis the degree
  axis refused to give).  The data tables + four refuted guesses live in the section header above
  `GPW_SCF.DISABLED_BeckeRecipeLadder_*`; the standing rule ("calibrate a grid criterion on a simple
  METAL, or do not ship it as a global default") is now cited by both V2.7 and R2.15.  *(original
  question + analysis follow)*  **Are the Becke recipe's `nRadial`/`angularDegree` defaults
  over-generous?  USER 2026-08-07:**
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

---

## V2.7 — harvested 2026-09-08 (was CleanupCandidates.md L2862-2873)

*(verbatim)*

- **V2.7 ✅ DONE 2026-08-17 (concurrent-cleanup session).**  `RadialResolutionRatio(mp, alphaMax)` in
  qchem.Mesh.XCPolicy — the closed-form \f$r_{peak}/\Delta r(r_{peak})\f$ statistic, MHL-only (+∞ = no
  claim off-MHL) — plus a ONE-level warning (`kRadialRatioFloor=3.0`) spoken from `ResolveXCMesh`
  whenever a Becke mesh is the outcome and sharpness is known.  The closed form reproduces the V2.6
  record's own quoted anchors (Si 3.4@nR=30, Al 3.2@nR=30 — pinned in `XCPolicyTests`), production
  configs all clear the floor (NaF, the tightest, at 3.09) so healthy runs stay quiet, and the measured
  4-50x-inadequate nR=20 class all fire.  A two-level version with a metal info-line (<4.2, Al's first
  adequate rung) was REJECTED — it would print on every production NaF run to warn about a metallicity
  the cost model cannot see; the metal fact rides in the warning text instead.  Floor remains
  insulator-fitted per its doc note: re-calibrate on {Si, Mn-atom, Al, MnO} before promoting it beyond
  a diagnostic.

