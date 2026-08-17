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
