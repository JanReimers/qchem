> Verbatim archive of CleanupCandidates v2 as of 2026-10-01 (RECORD, no open work).  Live successors: `doc/OOD-SOLID-Cleanup.md` and `doc/CleanCode.md`.  Row ids cited as `H3 <id>`; the v1 worklist is `CleanupHistory2.md`.

# Cleanup Candidates — the SOLID/OOD debt worklist (v2, rebuilt 2026-09-19)

**What this file is.**  The OOD debt in the abstract interfaces — faces that lie, casts that should be
types, axes that got fused — with the user's charter at the top and one short row per OPEN item below.
Tooling/build/hygiene debt is `doc/OpenWork.md` §3, not here; features are its §2.

**Where the history went.**  The v1 worklist (2026-08 → 2026-09-16, 3257 lines, ~70 closed rows) is
`doc/Records/CleanupHistory2.md` **verbatim** — every row below cites it by id (`H2 R1.0`), and every
closed row's full argument is there.  The earlier harvest (2026-08-07 → 09-15) is `doc/Records/CleanupHistory.md`.
Three conventions that were living here moved to `CLAUDE.md` Design (Get/Make; say why when loosening
encapsulation; name a capability query for what the client consumes).

**Rules.**  (1) A ✅ goes to `doc/Records/CleanupHistory2.md` the day it is written, with a one-line stub
here — or simply deleted from here, since the verbatim record already carries the row.  (2) Check the tree
before believing a row: on 2026-09-19 three "open" rows (D11, V2.1, the `Vxc_QuadraturePol` dead code) were
already closed by V1.37 and the occupation-policy work.  (3) A durable ruling becomes a `doc/Pins.md` pin;
the argument stays in the history.

Legend: ✅ done · ⛔ refuted/withdrawn · 🔶 design ruled, not built · ⏸ deliberately deferred · unmarked = open.

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

---

## R — READY (correctness-adjacent first, then hygiene)

| id | what is open | next concrete action · record |
|---|---|---|
| **R1.0r** | **ρ is star-averaged under a BIGGER group than the k-mesh has.**  `DetectPointOps` hands the density the full crystal group while `FoldGrid` (KP-0, pin 13) folds only under the mesh-invariant ops; on an anisotropic mesh (Si 2×1×1: 2 ops vs 48) the SCF fixed point is not the one the sampling defines.  Not a count bug (a star-average is norm-preserving) — a JUDGEMENT: VASP/QE symmetrize under the INTERSECTION group | needs a ruling; cheap once decided — expose `FoldGrid`'s mesh-symmetry test as a `MeshSymmetryOps(N, shift, ops)` filter and intersect in `DetectPointOps`.  Γ-only and isotropic meshes unaffected · H2 R1.0r |
| **R1.0 (metric axis)** | **`FIT_SF_Ortho<T>` — separate the metric axis into faces, BOTH sides together** (specced 2026-08-23, not built; was OpenWork item 6).  `OverlapDiagonal()` sits on the metric-NEUTRAL face so `Fit_IBS` invents an answer in the wrong normalisation.  ⚠ ACCEPTANCE: must NOT become `if (dynamic_cast<FIT_SF_NonOrtho*>)…` — a type switch wearing a cast is worse than the bool.  Measured: all eight `isOrtho()` sites are asserts, zero live branches; narrowing a PARAMETER type is the sanctioned substitute; the fitter Factory's representation branch is NOT a metric branch and stays | move `OverlapDiagonal()` to `FIT_SF_Ortho<T>` (δ, PW derive it), `Fit_IBS` loses it outright, `OrthogonalFit(const FIT_SF_Ortho<T>&)`; mirror `FIT_CD_ABS::isOrtho()` in the same increment; no orthonormal marker face · H2 R1.0 "NEXT INCREMENT, SPECCED" |
| **R1.0 (step 2)** | ⏸ **Make atomic `op(r)` HONEST** (real Y_lm combinations) and retire `ImplicitAngular_IBS` — explicitly SCOPED OUT by the user (2026-08-22): step (1) removed the consumers, so the fake sits behind a face nobody calls | only if something needs it; settle first what `op(r)` returns per (i, m) given the spherical solver keeps m-degeneracy in the occupations · H2 R1.0 |
| **R1.0q** | ★★ **Standing target: every `Dynamic_HT` term on the `MatrixForward<T>`/`MatrixAdjoint<T>` pair** (user 2026-09-09).  Static terms are explicitly OUT (no D to push, no field to pull) — the same partition `RefreshForDensity` folds over.  ✅ **FIRST TERM DONE 2026-09-19 (`6f71e1a6`) — `FittedVee`, the molecular Coulomb-metric density fit, and it PROVES R1.0p by construction**: `DenseProjector3Integrator<T>` realises both halves over the basis's analytic \f$\langle ab|c\rangle\f$ tensor with no grid anywhere; `ConstrainedFF` is "actor 2" (one integrator per block, forward vended through the new `Fitting::DensityProjector` — the Coulomb-metric mirror of `ScalarProjector` — adjoint used in `Repulsion`); the J⁻¹ solve moved from the projection to the fitter; `FiniteIrrepCD::GetRepulsion3C` no longer touches `.dense`.  Bit-identical (every M_DFT iteration line), 4 unit gates, 885/885 | **NEXT:** `FittedVxc` (molecular XC) — harder, and INSTRUCTIVE: its forward is ρ at the fit basis's mesh points (a `DenseMatrixIntegrator` over the `Fit_IBS` mesh), the nonlinear functional sits between, the adjoint is the Gaussian overlap fit + `Overlap3C`.  ⚠ It hits `IrrepCD_Core<T>::ProjectOnto`'s hard-coded `Orbital_DFT_IBS<T,dcmplx>` and `ScalarProjector`'s `dcmplx`-only overloads: the TFit axis is fused into the density's projection route — V1.35's fusion one library down, and the reason V1.35 is a two-library campaign.  Do it after V1.34's ruling lands (below), not before.  Then `Vee_Hartree` + the quadrature XC terms (already on the sampler/projector; mostly naming), ✅ **+U BORN on it 2026-09-20 (`6aa8170d`)**: `LowdinProjector<TBlock>` is the pair, the density projects through `ScalarProjector`, the term keeps the adjoint · H2 R1.0q, R1.0p |
| **R1.0b** | **Shared-radial (SP/"L"-shell) support in the Gaussian94 reader** — the user's structural cure for the Cartesian-d s-contaminant rank deficiency (CP2K encodes it natively; MOLOPT shares one exponent set across s, p, d) | ⏸ BLOCKED on the flagged `PG_Cart::IrrepBasisSet` bug (the reader MERGES same-exponent shells across l — hence valgen's "keep exponents disjoint" rule).  Measure, don't assume: sharing exponents buys compactness, whether it buys conditioning is an experiment (`GPW_SCF.MnAtomInBoxDChannel`) · H2 R1.0b |
| **R1.0f** | The GUI will want \f$v_{xc}(r)\f$ and \f$\rho_{DM}-\rho_{fit}\f$; nothing in the tree evaluates a fitted field any more and the δ route cannot.  \f$v_{xc}\f$ is NOT a fit-basis question (evaluate ρ anywhere, apply the functional); for δ the residual is identically zero — the honest δ diagnostic is \f$\rho_{DM}\f$ vs the band-limited \f$\tilde\rho\f$ Hartree integrates | GUI-side, on demand; it must not come back disguised as `op(r)`; a δ expansion at arbitrary r is an INTERPOLATION problem and is owed a name that says so · H2 R1.0f |
| **§K + I.1** (from FittingCleanupPlan) | **§K** fit-{G} densification (`ProjectedScalar_G` sub-bullet overtaken; the densification itself remains) — UNBLOCKED, deferred by user ruling as anchor-moving (sprint S row A3).  **I.1 residual**: `GetEpsXc()=0.75*GetVxc()` base default (`ExchangeFunctional.C`) is exact for Dirac exchange only — a silent-wrong inherited default the day a GGA forgets to override | §K in the sprint window; I.1 with whatever touches the functionals first (GGA, OpenWork §2) · H2 header L60–75 |
| **BM(1)** | **`MeshParams::cellKind=Becke` alone gives a DEGREE-5 mesh** (struct defaults 30/1/5 vs `BeckeXCParams` 40/2/29) — cost a bogus 40 mHa "symmetry bug" (2026-09-07).  ⚠ `GPW_BECKE_L/NR/ALPHA` are consulted only for arguments `<0` | (c) make `cellKind` unsettable on its own so asking for Becke means asking for the recipe — the compile-time answer; (a)/(b) are the cheaper fallbacks · H2 "cellKind=Becke" |
| **R2.14** (remainder) | The two XC-term renames waited on V2.1: `Delta_XC → DeltaFittedVxc`, `PW_XC → PWFittedVxc` ("FittedVxc + WHICH fit basis"; PW is the binding requirement — the `G_FieldEvaluator` quadrature axis — not `isOrtho`) | ✅ OVERTAKEN by V1.37 (one term per operator; neither class exists).  Check the surviving names read as "FittedVxc + fit basis" and close · H2 R2.14 |
| **R2.18** (remainder) | ⏸ Encapsulation of the `Make*` family (public vs protected) — user: *"no strong policy … maybe the right policy will emerge as we refactor"* | deliberately open; the CLAUDE.md rule "say WHY at the declaration when loosening" is how the evidence accrues · H2 R2.18 |
| **R2.21** (remainder) | A `route` ctor argument to FORCE the ball route for an A/B (today `GPW_XCROUTE` only reports); needs a capability question on the neutral `Band_FT_IBS` face | when someone wants the A/B · H2 R2.21 |
| **R2.23** | The atomic Rk cache sizes by incremental `LMax` discovery + evict-on-register (`Cache4::Register`, `Rk::isSupported`); ERI4Rework §6's declare-ceiling / demand-grow fix never started | demand-grow (preferred): the basis declares its worst-case `sym_t` at registration, `Rk` grows monotonically, eviction deleted; atomic layer only; bit-identical; RAM is the only measurable · H2 R2.23 |
| **R2.24** | `qcSymmetry`'s directories are named for the SYSTEM (`Atom/Molecule/Lattice_3D`), not the GROUP (`O3/Point/Lattice`) they mirror (pin 14) | mechanical rename, one commit, when nothing is in flight in `src/Symmetry/` · H2 R2.24 |
| **`Vcorr_QuadraturePol` name** | `Vxc_QuadraturePol` is gone (V1.37); the survivor `Vcorr_QuadraturePol` (`PWTerms.C`, `Hamiltonians.C`) is no longer "the correlation half of a pair" — it is the spin-native XC term | rename (`Vxc_SpinNative` or whatever V1.37's family reads as); three files · H2 "Vxc_QuadraturePol is dead code" |
| **`MinVirial`** | ✅ CLOSED in passing 2026-09-19: the default is `1e30` "DEFAULT OFF" (`SCFParams.C:21`) — the 1e-13 comment/default clash V1.27 recorded is gone | — |

## V — VERIFY: design questions that need a ruling or a measurement

| id | what is open | next concrete action · record |
|---|---|---|
| **V1.35** ★★ | **`Dynamic`/`Static` and the BLOCK SCALAR are orthogonal axes and the term hierarchy has FUSED them**: `Dynamic_HT_RealBlock` is the mixed corner (block `double`, run `dcmplx`) given its own name, so every cross-axis capability is declared twice (R1.0h's two slot hooks).  Cure = the tree's own answer, a two-parameter template `tDynamic_HT<TBlock,TRun>` (as `Orbital_DFT_IBS<U,TFit>`), the `*_RealBlock` faces disappear — a diamond done correctly (CLAUDE.md).  ⛔ It is a TWO-LIBRARY campaign: the same fusion sits in `ChargeDensity::Dynamic_CC_RealBlock` and every density that implements it | **BEFORE IT, in order (ruled 2026-09-19):** (1) **R1.0q's first molecular term** — a term whose matrix comes off a block-scalar projector has a scalar-generic `MakeMatrix`, so its `MakeMatrixR` degenerates into a copy and V1.35 becomes a mechanical collapse; (2) **V1.34's ruling, shape 3** — V1.35 makes real blocks first-class instantiations, so `FitContraction<double,dcmplx>` must exist or V1.35 lands on a `bad_cast`; (3) ✅ **+U written SCALAR-GENERIC from day one — DONE 2026-09-20 (`6aa8170d`)**: `Hubbard_U` is one `template<class TBlock> MakeMatrixT` body, `MakeMatrix`/`MakeMatrixR` one line each, `LowdinProjector<TBlock>` the pair; its `PrepareSlots` is the two-base example the design note below describes.  THEN plan V1.35; do not start it as an increment.  Fingerprint and finish line: the **31 `*R`-suffixed methods** (`GetMatrixR`, `MakeMatrixR`…) across 12 files — they exist only because the corners are types; done = the suffixes are gone.  Design first: a term inheriting both instantiations writes one explicit `PrepareSlots` calling both bases (a compile error until it declares its two caches — a feature) · H2 V1.35, V1.35a |
| **V1.34** 🔶 | **The fitter's contraction face is templated but half-realised**: `FitContraction<U,TFit>` exists for two scalars, is realised for `<dcmplx,dcmplx>` only, and fails as a `bad_cast` mid-SCF in Release.  A face that lies | **RULED 2026-09-19: shape (3), compile-time** — a fitter exposes exactly the instantiations it defines, so a caller needing another one does not LINK; no askable-capability bool, no runtime cast.  V1.35 then makes `<double,dcmplx>` a real instantiation rather than a corner to discover.  Design note from R1.0q's first term: `ConstrainedFF` now carries a `static_assert(T==double)` naming its lineage — the same idea one step earlier · H2 V1.34 |
| **V1.27** (remainder) 🔶 | **`MolecularSCFIterator`/`SolidSCFIterator` are named for the STRUCTURE and discriminate on PP-ness + grid-ness** (user: all four {Molecular,Solid}×{PP,non-PP} will be needed; mixins).  The virial axis collapsed to a derived bool (`IsVirialValid()`, landed); what remains is the column set = a TRACE POLICY chosen from what the run IS (pseudised? variational? gapless?), not from the iterator's type; `MolecularSCFIterator` is an EMPTY subclass | a small `TraceColumns` value the facade supplies; sequence with whatever next touches the iterators (the MD `tSCFIterator` type in OpenWork §2 is a natural moment) · H2 V1.27 |
| **V1.24** (i)+(iii) | **GDM**: (i) `FDMax` reads as a step size but is the ENGAGEMENT gate on ‖[F′,D′]‖ — rename (`EngageBelowFD`) and consider an intensive norm; (iii) SOFT-DIRECTION preconditioning — the 1/(ε_a−ε_i) diagonal Hessian blows up along near-degenerate diffuse modes, and the landed precondition check only makes GDM DECLINE on a wrong occupation.  (ii) the fallback commit is ✅ | with the OT build (OpenWork §2 row OT) — same minimiser family, same seams · H2 V1.24 |
| **V1.22** | **Becke partition drop decisions are per-point and bit-sensitive**, so the site-adapted caller post-filters orbit-incomplete points; the free builder once kept a point whose AFM translation partner was tail-dropped.  Per-REPRESENTATIVE drop rule (angular dir × radial shell, applied to the whole orbit) removes the second fold pass + the filter — and is up to ~23× on the partition cost (98816 points, 4290 orbits) | sprint S row A2 (anchor-moving) · H2 V1.22 |
| **V1.1b** 🔶 | `Eee = 2·EeeFit − EeeFitFit` is DUNLAP-SPECIFIC (robust only under a COULOMB-metric fit); the overlap-metric seed path (`NumericCD`) is kept off the energy expression only by an accident of typing (`tDM_CD*` parameters).  ⚠ GPW's fit is NOT exact (Gaussian products have infinite bandwidth) — it is safe because it never uses the Dunlap form, not because ρ̃=ρ | name the invariant where it lives (the fitter/energy pair) and fold "which metric ⇒ which energy expression" into the metric face (R1.0 metric axis).  Awaiting the user's re-read of the paper · H2 V1.1b |
| **V1.3** (remainder) | the QUADRATURE-TERM face — `FittedEpsXc`/`FittedVxc` simplification (the second list of the original item) | ✅ likely OVERTAKEN by V1.37 (no `FittedEpsXc`/`itsEpsFitter` in the tree); verify the molecular term samples ρ once for both \f$v_{xc}\f$ and \f$\epsilon_{xc}\f$ (the R1.0 ruling) and close · H2 V1.3 |
| **V1.19** (remainder) | ⏸ the seed's flip-group sub-cell duplication — removing it needs a per-SITE form-factor overload on the basis face, which the pseudo-wall pin (CLAUDE.md, Pins 9) weighs against; the seed's ONE remaining concrete-`Atom` consumer | deliberate; revisit only if a second consumer wants per-site form factors · H2 V1.19 |
| **V1.29** (remainder) | the Hessian-based magnetic-ORDER DISCOVERY loop (Davidson on the stability matrix per non-symmetric irrep block) — unbuilt, not needed for MnO whose order is imposed | = the SSB-descent feature (OpenWork §2); on a material whose order is unknown · H2 V1.29 |
| **V1.38** ⏸ | The Point spec in the core + one thin IBS class per (G, engine) — the molecular tier's evaluator-injection sequel; 2–4 sessions, bit-identical.  STASHED by agreement | triggers: a second Gaussian engine needing the lattice role; the NAO family; a fourth hand-written `PG_*::Orbital_IBS`; or the `LatticeSum1E` ISP split landing (which goes FIRST) · H2 V1.38 |
| **V1.40** | **`Hamiltonian::ManifoldSymmetry` is three-quarters group theory living in qcHamiltonian** (user, 2026-09-21, on reading `Hubbard.C`).  What it does with no Hubbard concept in sight: decompose a finite orthogonal matrix rep into irreps by CLUSTERING (character vectors, no table — the molecular tables are abelian-only, and this gives 3-D/2-D irreps for free), build isotypic projectors (scale from Tr A²/Tr A, complex-type-safe), BRANCH a subgroup's irreps against a parent's (\f$\dim=\mathrm{Tr}\,P_kP_p\f$), name a vector by dominant isotypic weight, and rotate a degenerate eigen-cluster onto the projectors.  Only three things are Hubbard's: the word "U slot", `WriteSlots(U, Uirrep)`, and `DudarevInEigenbasis` (already a free function).  `Rep(AoShells, ops)` is qcSymmetry's too — it is `BuildOperationRep` for one site with no permutation.  Two clients are already waiting in qcSymmetry: the site-adapted meshes (non-abelian site irreps) and the I4 lattice SALC over the spherical view (a rep → irreps decomposition is exactly what it lacks).  ⚠ **The one real obstacle: qcSymmetry is LAPACK-free by design** (`SphericalRep.C`, links only qcMath+qcCommon) and the clustering needs a symmetric eigensolver (`blazem::eigen`) | **Not now — let the +U work finish and settle first (increments 3–4 will exercise it).**  Then: split `ManifoldSymmetry` into `Symmetry::IsotypicDecomposition` (rep → irreps, projectors, branching, labelling; the qcSymmetry half) and a thin Hubbard `SlotTable` over it (slot order, U per slot, the printed table); **LAPACK question RULED 2026-09-21 (user): "fair game for any library that needs it … a totally harmless dependency"** — qcSymmetry takes `blazem::eigen` directly; no injection, and the `SphericalRep.C` "LAPACK-free" remark is a note about the past, not a constraint.  Gates move with it (`UTHamiltonian ManifoldSymmetry.*` → `UTSymmetry`); the MnO slot table `{7,14,14}` and the O_h characters are the bit-for-bit finish line · H2 V1.40 |
| **V1.39** | **The SUM OF MATRICES could be the MATRIX OF THE SUM**: terms sharing a quadrature hand back finished matrices and each pays its own gather; a term could contribute its FIELD to a shared quadrature that gathers once.  Expressing that without destroying the term seam is the work | not urgent; live when a third grid-shared term (GGA's gradient terms; NOT +U, a matrix-space term) arrives.  Do not paper over it with a cache · H2 V1.39 |
| **V4.1 / V4.2** (watch triggers) | `CollocMemo` dual duty (replay memo + adjoint D-screen, different lifetimes) — split at a THIRD consumer.  `SolveSPD`/NNLS module-private in `SymmetrizeMesh.C` — promote to qcMath at a SECOND consumer (Blaze has no NNLS) | act when the trigger appears · H2 V4 |

## D — decided elsewhere / structural debt awaiting its campaign

| id | what is open | next concrete action · record |
|---|---|---|
| **D-BASIS-PP** | ⛔ **Nothing checks the BASIS against the PSEUDOPOTENTIAL.**  A run declares its PP variant via `SolidCalcOptions::species` (`{"Li",3}`) but a `.bsd` block is keyed by ELEMENT alone, so a q3 run can silently be handed the q1 block — a basis validated for a one-electron valence describing three, with no diagnostic anywhere.  Latent today only because every element in `valence_lowq_*` happens to be the variant its runs use | found 2026-09-23 while minting the Li bases (`doc/HubbardUPlan.md` §5a.2).  The `.bsd` headers are already rich; add a machine-readable per-element provenance line (invisible to a Gaussian94 reader) and have the factory assert each block's q against the declared valence, THROWing on a mismatch.  Sibling of D-SEED1; goes live the moment `valence_semicore.bsd` exists |
| **D-SEED1** | **The atomic seed library is keyed by (Z, functional) but must distinguish (Z, functional, **q**)** — `FindAtomicEntry` with `Nval<0` returns the FIRST matching entry, and `IonicSADTargets` calls it that way to read the neutral valence count.  So "which entry is the neutral one" is decided by FILE ORDER and stated nowhere: it is correct today only because every ion happens to be appended after its neutral.  ⛔ Li makes it genuinely ambiguous — q1 neutral is 1 electron, q3 neutral is 3, both are "neutral Li" — so the q1-vs-q3 discriminator (`doc/HubbardUPlan.md`) cannot A/B the seed until this is keyed properly | found 2026-09-23 while minting the Li basis.  Either key the lookup on the PP variant the caller already knows (`SolidCalcOptions::species` carries q), or make an ambiguous match THROW rather than pick — `throw` is a marker.  A silent first-match on a physics input is how a wrong number gets banked |
| **D9** | **The ρ̃ mixer FUSES the preconditioner with the extrapolator** (`PulayMixer` owns both the Kerker filter and the history; `PolarizedDensityMixer` is a `tDensityMixer` rather than a filter stage).  End state = `residual → per-channel preconditioner → ONE joint extrapolator` — the VASP/QE/CP2K architecture, and the shape N3 (charge vs spin channels) and the non-collinear 2×2 spin density need | ★ lands WITH N3 / +U (OpenWork §1 step 5.2); makes Broyden a drop-in beside Pulay.  Re-measure run 11's (ρ,m) verdict when it lands (it was refuted at `PulayDepth=0`, i.e. with no extrapolator to shape) · H2 D9 |
| **D12** | **`GPW_Evaluator::Eval` truncates by RADIUS, not by MAGNITUDE** — one global `itsMaxReach` from the most diffuse exponent + a geometric `BuildImages` bound that has already been patched once (oblique MnO cell).  A **pin 1 violation** in an evaluator whose analytic lattice sums already screen per pair on magnitude and REPORT the reach | screen each (function, image) term on its own magnitude, per pair like the analytic path; the radius becomes a reported consequence.  Refuted as the MnO sublattice cause — this is the design debt on its own · H2 D12 |
| **D13** | **Single-species `HGH_LocalPotential`/`HGH_SeparablePotential` silently IGNORE their `int Z`** — the foot-gun that manufactured the retracted "V_long sharp-field defect" (a Mn q7 field integrated at the O site).  Production is safe (`MultiSpecies_*` routers); hand-written probes are not | compile-time-over-runtime: drop `Z` from the single-species concrete API and let ONLY `MultiSpecies_*` model the Z-dispatching face (or store `itsZ` and throw on mismatch) · H2 D13 |
| **D10** | **`DisplayEigen` builds rows and renders in ONE pass, so nothing is assertable** — three real bugs found by reading logs (dropped doubly-empty levels; `setprecision(0)` on 0.996; an absent channel's ε filled from the other, which manufactured MnO's spin-up "hole").  Deeper: rows are paired by `(n, sym)` but `n` indexes a DEGENERATE GROUP that differs between channels once ε↑≠ε↓ | extract the row build as a pure function of the two `EnergyLevels`, pair by energy/character; then each bug is a unit test · H2 D10 |
| **D4** | **The basis carries the run's SYMMETRY POLICY** (`GetReciprocalPointOps`/`GetDetectedReciprocalOps` on `tBasisSet`, `SetSymmetryOps` on the GPW evaluator, `imposeSymmetry` per factory) — 9 sites | collapses to consumers reading ONE `SymmetryPolicy` context when `SymmetryUpgradePlan.md` §3's object is plumbed; the per-factory bool dies · H2 D4 |
| **D3** | `FIT_SF_ABS::SymmetrizeRaster` — a symmetry op on a fit-basis face (SRP); 3 sites | removable once a τ-acting DIRECT-raster fold variant of `FoldGrid` exists: precompute, ctor-inject into `tComposite_CD`, delete the virtual · H2 D3 |
| **D5** | Becke `MeshParams` ergonomics under `imposeSymmetry`: `mp.angular` silently REPLACED by the site-adapted rule; an unachievable L for a low-symmetry site is a bare assert | with the `SymmetryPolicy`/facade pass: auto-resolve, warn on scheme override, a real error on exhaustion (the degree-typed half landed as R2.15) · H2 D5 |
| **D2** | "Why is a `FourierCD` different from a `tChargeDensity`?" — governed by pin 6 (no "Fourier" in abstract faces).  `FourierDensity` is a cross-cast MIXIN; the seeds (`NumericCD`/`SeedCD`/`PolarizedSeedCD`) are a third leg beside `DM_CD`/`FittedCD` | decide with the metric face (R1.0) and V1.1b — the seeds' overlap-metric fit is the same boundary · H2 D2 |
| **D11** | ✅ CLOSED 2026-09-19 (verified): the occupation policy is passed at CONSTRUCTION and the seed fill runs on its explicit default state (`SCFIterator.C:222`) — the ordering defect is gone | — |
| **V2.1** | ✅ OVERTAKEN by V1.37 (one spin-native XC term per operator; `Delta_XC*` gone).  Its DURABLE part is a ruling for the non-collinear future: **represent spin as SU(2)/matrix, not SO(3)/direction** — real diagonals + ONE complex off-diagonal (4 real DOF = (n, m)); the collinear tier is the diagonal case; symmetry acts as \f$U\rho U^\dagger\f$ and COMPOSES; do not let a `Spin` label leak into an XC functional's contract as a permanent parameter | ruling kept here until non-collinear work starts, then a pin · H2 V2.1 |

## Rulings that live in this file's history and are cited from code

`R1.0` (a fit basis provides the integrals to do the fit, and integrals over ITSELF, nothing else; house
names `Overlap`/`Charge`/`Repulsion`; a value table has no business in an interface; nobody holds both
halves of a forward/adjoint pair — build one object, hand each client its half), `R1.0h` (the per-iteration
scope object DECLINED), `R1.0j` (the density↔operator seam is not "XC" and not "quadrature"; `LatchRoute` is
not redundant; the exch/corr one-collocation sharing is worth 4.8 s/iteration on NaF — never push ρ back
into the terms), `R2.5b` (`throw` is a marker), `V1.20` (`.Internal.` = the FAMILY boundary), `V1.37`
(Pol/UnPol are imposed subgroups).  Read them in `doc/Records/CleanupHistory2.md` by id.

## 2026-10-02 — D-BANNER, D-LOWRANK closed
- **D-BANNER ✅**: the solid banner (`SolidCalculation.C` EmitSCFBanner) already named Kerker + Pulay; the remaining defect was the JSON report's `scf.standard.mixer` tag in `Calculation.C` (`Pul` hid Kerker).  Now `Lin`/`Ker`/`Lin+Pul`/`Ker+Pul`.  M_Calculation 12/12.
- **D-LOWRANK ✅**: `IrrepCD.C` low-rank-D doc now qualifies the "rank same from tol 1e-6 to 1e-12" claim as kT=0 only.  (Making pivoted-Cholesky failure loud was already pin 21 territory; not re-opened.)
- **D-LOWRANK (loud half) ✅**: `LowRankFactor` (`IrrepCD.C`) now prints `[DM factor] WARNING` (first 5 per run) when D fails the PSD/trace guard; the exact full-rank route was already the fallback.  D-BANNER's "linear fallback prints no mixer line" was already stale (solid banner says LINEAR D-mixing).  Full `ctest -j8`: 982/983, only the known D-FLAKE timing test (`M_PG_BoxWalk.WhereTheContractionSpendsItsTime`) failed.
- **D-FLAKE ✅**: `M_PG_BoxWalk.WhereTheContractionSpendsItsTime` failed the sweep at 27% vs a 25% bound (quiet box ~16%).  Fix = the no-cry-wolf recipe, written into the test: ratio of back-to-back timings, THREAD CPU time (not wall), MIN of trials (noise is one-sided), retry the failing row with 15 trials before failing.  12 concurrent copies pass, no retry fired.
- **`Vcorr_QuadraturePol` name ✅ (stale row)**: the class is already `Vxc_Quadrature` (V1.37 merged the family); only a stale comment in `Hamiltonians.C` still named the old pair — fixed.
- **D-SKCOND ✅**: a facade run with no report open now prints one `[basis cond]` line (worst min-eig S, worst cond(S), owning k-blocks) via `PrintGpwConditioning`; with a report open the per-irrep rows (`basis.perIrrep`) already carried it.  Full ctest -j8 983/983.
- **D-CUBE0 ✅**: instead of a recurring `GPW_CONTRACT_CUBE=0` sweep, `NR_Evaluator::ContractCubeOverride()` (test hook) + `M_TransitionCollocation.WalkAndContractionRoutesAgree` run the production CollocateDensity/IntegratePotential on both routes in one process (agree to ~1e-14, gate 1e-10).  The env sweep also passed 983/983 on 2026-10-02.
- **D-MAKEWATER ✅**: the five integration tests (M_DFT, M_Sym, M_Response, M_Calculation, M_HF_U) now use `Materials::GetMolecule("H2O")`.  The two BasisSet unit tests keep their own geometry (Å, different).
- **D-STRUCTDATA step 1 ✅ (2026-10-02)**: new `qchem.StructureData` module in qcStructure (`GetMolecule`/`GetCell`/`KindOf`/`CellNames`/`MoleculeNames`/`AtomInBox`/`DimerInBox`; concrete return types, kind mismatch throws); data files moved to `src/Structure/Data/` (macro `STRUCTURE_DATA_PATH`); `Materials` (Calculation) now wraps it and adds only the valence counts; `Materials::GetMolecule` removed, 5 integration tests + the two BasisSet unit tests (`M_LibCint`, `M_MEvaluator`, now on the shared bohr water) use it; 4 new `UTStructure` tests.  ctest 988/988.  Step 2 (PP data → own file, + D-SEED1, D-BASIS-PP) stays in the tracker.
- **D-LONG ✅ (2026-10-02)**: budget 60 s CPU (user).  `qchem_discover_tests` (root CMakeLists) discovers a test whose name ENDS in `_Long` in a second pass labelled `long`: sweep = `ctest -j8 -LE long` (987/987), long tier = `ctest -j3 -L long` (3/3, 377 s wall): `GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long` and `Γ_U_Shub_Pol_Smear_CP2K_Long` (were DISABLED_, now run), `ResponsePolarizability.GPW_Si_k211_U2_FrozenChi_eqFiniteDifferenceLRT_Long` (191 s under -j8).  `scripts/testgrid` treats `_Long` as a suffix.  Pre-existing testgrid violations (+U seed tokens ACBN0/Atomic3p, not axis tokens) are unrelated and still open.
- **D-NAFGRID ✅ (2026-10-02)**: ported to the facade: `GPW_NaF.Γ_Imp_eqColdStart` (coarse Ecut=40 alias-free -> `SolidCalculation::Restart` onto the auto fine grid lands E=-24.43040 == the cold fine run to 1e-4; 12 iterations vs the cold run's 21; coarse stage -24.43245, 17 iters) and `GPW_NaF.Γ_Imp_eqExactResume` (same-grid restart: E reproduced to 6e-9; first iterate bounded at 1e-3 -- it is 1.3e-4 off, an OpenWork §3 row).  ~33 s total (not `_Long`).  The raw-iterator `DISABLED_Γ_GridContinuation` is deleted.
- **D-SEED1 ⛔→partly ✅ (2026-10-02)**: the NEUTRAL ambiguity already threw (9f926e37, 2026-09-23); an EXPLICIT `Nval` was still first-match-wins (Li q3's +2 ion and q1's neutral both hold 1 electron) -- it now throws too, with a fixture-library unit test (`AtomicDensity.AmbiguousPseudopotentialVariantThrowsInsteadOfPickingTheFirst`).  Remainder (key the request on q) is in D-STRUCTDATA step 2.
- **D10 (build/render split) ✅ (2026-10-02)**: `PolarizedEigenRows`/`UnpolarizedEigenRows`/`OccupationCell` in new module `qchem.WaveFunction.EigenRows` (qcWaveFunction); `tCompositeWF::DisplayEigen*` only draws; `EnergyLevel` gained a value ctor for tests; new exe `UTWaveFunction` (7 tests, in the `allTests` DEPENDS list).  Rendered tables verified unchanged (NaF unpolarized byte-identical; polarized water shape unchanged).  Pairing-by-energy stays open in CleanCode.
- **D-TESTER ✅ stale (2026-10-02, user)**: `QchemTester` is long gone (only "migrated off" comments remain) and the symmetry-naming cleanup is done.  The one live fragment, `MolecularSym_EC` (a stop-gap per-irrep occupation EC for the first symmetric-HF test), was DEAD CODE -- no module imports `qchem.ElectronConfiguration.MolecularSym` -- so it was deleted rather than renamed to `FixedIrrepOcc_EC`.  D-RUNDATA renamed D-CLOUD (archive account is the only remainder).
