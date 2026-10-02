# OOD / SOLID debt — open rows only (split from CleanupCandidates v2, 2026-10-01)

**Charter (user).**  The abstract interfaces are what client code sees.  The goal: ~95% of the high-level
faces (charge density, Hamiltonian, orbitals, wavefunction, SCF iterator, accelerators) structure-neutral and
basis-neutral (atoms, molecules, 3D lattices; 1D/2D later).  Structure-specific classes creeping into
qcChargeDensity/qcHamiltonian/qcWaveFunction/qcOrbitals are the debt to pull back. Concrete classes that are structure-specific is a small (or no) problem as these should not be seen by client code.  Abstract interfaces that the client code sees is the priority for cleanliness.  Principles in short:
SRP applies to mixins, not to final evaluators; OCP — change a neutral face only for a reason, and ask what it
means for atoms and molecules; LSP — a virtual defaulting to "not implemented" is a violation; ISP — a
getter needs a client that will DO something with it (ask the owner instead); a setter wants to be a ctor arg;
DIP — callback interface owned by the lower layer (`tStatic_CC`/`tDynamic_CC` are the model).
**Heuristic:** a "special case" that is SIMPLER than the general case is the general case with terms zeroed;
one that needs machinery the general case lacks is real and belongs on the types that have it
(`tSpinResolved_CD` idiom), not on the base with a stub.

**Rules.**  File OOD debt HERE, never fix inline.  A ✅ leaves the same day to `doc/Records/CleanupHistory3.md`
(or later) — delete the row.  Check the tree before believing a row.  A durable ruling becomes a
`doc/Pins.md` pin.  Full argument for every row: `doc/Records/CleanupHistory3.md` (= `H3`) and
`CleanupHistory2.md` (= `H2`, by id).  Legend: 🔶 ruled not built · ⏸ deferred · unmarked open.
Hygiene/tooling rows: `doc/CleanCode.md`.

| id | what is open | next action · record |
|---|---|---|
| **V1.35** ★★ | `Dynamic`/`Static` and the BLOCK SCALAR are orthogonal axes, fused: `Dynamic_HT_RealBlock` is the mixed corner (block double, run dcmplx), so cross-axis capabilities are declared twice.  Cure: two-parameter `tDynamic_HT<TBlock,TRun>`; `*_RealBlock` faces vanish.  TWO-library campaign (`Dynamic_CC_RealBlock` too) | Prereqs: R1.0q first term ✅, +U scalar-generic ✅, **V1.34 ruling (shape 3) still to be realised**.  THEN plan V1.35 (not an increment).  Finish line: the 31 `*R`-suffixed methods (`GetMatrixR`, `MakeMatrixR`…) gone.  A term inheriting both instantiations writes one explicit `PrepareSlots` · H2 V1.35, V1.35a |
| **R1.0q** ★★ | Standing target: every `Dynamic_HT` term on the `MatrixForward<T>`/`MatrixAdjoint<T>` pair (static terms OUT).  Done: `FittedVee` (6f71e1a6), `Hubbard_U` (6aa8170d) | **Next `FittedVxc`** (molecular XC): hits `IrrepCD_Core<T>::ProjectOnto`'s hard-coded `Orbital_DFT_IBS<T,dcmplx>` and `ScalarProjector`'s dcmplx-only overloads (V1.35's fusion one library down); do after V1.34 is realised.  Then `Vee_Hartree` + quadrature XC terms (mostly naming) · H2 R1.0q, R1.0p |
| **V1.34** 🔶 | `FitContraction<U,TFit>` templated but realised for `<dcmplx,dcmplx>` only; other instantiations fail as a `bad_cast` mid-SCF | RULED shape (3), compile-time: a fitter exposes exactly the instantiations it defines, a caller needing another does not LINK; no capability bool, no cast.  `ConstrainedFF` has the `static_assert(T==double)` precedent · H2 V1.34 |
| **R1.0 (metric axis)** | `OverlapDiagonal()` sits on the metric-NEUTRAL face so `Fit_IBS` invents an answer in the wrong normalisation.  ⚠ Must NOT become `dynamic_cast<FIT_SF_NonOrtho*>`; all 8 `isOrtho()` sites are asserts | Move `OverlapDiagonal()` to new `FIT_SF_Ortho<T>` (δ, PW derive), `Fit_IBS` loses it, `OrthogonalFit(const FIT_SF_Ortho<T>&)`; mirror `FIT_CD_ABS::isOrtho()` same increment · H2 R1.0 "NEXT INCREMENT, SPECCED" |
| **R1.0 (step 2)** ⏸ | Atomic `op(r)` is fake (retire `ImplicitAngular_IBS`); SCOPED OUT by user 2026-08-22 | Only if something needs it · H2 R1.0 |
| **R1.0r** | ρ star-averaged under a BIGGER group than the k-mesh has (`DetectPointOps` full crystal group vs `FoldGrid` mesh-invariant ops; Si 2×1×1: 2 vs 48 ops).  VASP/QE use the INTERSECTION | Needs a ruling; then expose `MeshSymmetryOps(N, shift, ops)` and intersect in `DetectPointOps` · H2 R1.0r |
| **V1.27** 🔶 | `MolecularSCFIterator`/`SolidSCFIterator` named for STRUCTURE but discriminate on PP-ness/grid-ness; `MolecularSCFIterator` is an empty subclass.  Remaining: trace-column set should be a policy from what the run IS | Small `TraceColumns` value supplied by the facade; do when next touching the iterators · H2 V1.27 |
| **V1.1b** 🔶 | `Eee = 2·EeeFit − EeeFitFit` is Dunlap/Coulomb-metric-specific; overlap-metric seed path kept off it only by accident of typing | Name the invariant in the fitter/energy pair; fold "metric ⇒ energy expression" into the metric face (R1.0).  Awaiting user re-read of the paper · H2 V1.1b |
| **V1.40** | `Hamiltonian::ManifoldSymmetry` is 3/4 group theory living in qcHamiltonian (irrep clustering, isotypic projectors, branching, labelling); `Rep(AoShells, ops)` is qcSymmetry's too | After +U settles: split into `Symmetry::IsotypicDecomposition` (qcSymmetry; may link `blazem::eigen` — ruled harmless) + thin Hubbard `SlotTable`; gates move `UTHamiltonian ManifoldSymmetry.*` → `UTSymmetry`; MnO slot table `{7,14,14}` + O_h characters are the finish line · H2 V1.40 |
| **V1.39** | Sum of matrices could be matrix of the sum: grid-sharing terms each pay their own gather | Live when a third grid-shared term arrives (GGA gradient terms).  No cache band-aid · H2 V1.39 |
| **V1.38** ⏸ | Point spec in core + one thin IBS class per (G, engine) | Triggers: second Gaussian engine for lattice role; NAO family; fourth hand-written `PG_*::Orbital_IBS`; `LatticeSum1E` ISP split (goes first) · H2 V1.38 |
| **V4.1 / V4.2** | `CollocMemo` dual duty (replay memo + adjoint D-screen); `SolveSPD`/NNLS module-private in `SymmetrizeMesh.C` | Split at a THIRD consumer / promote to qcMath at a SECOND · H2 V4 |
| **BM(1)** | `MeshParams::cellKind=Becke` alone gives a DEGREE-5 mesh (defaults 30/1/5 vs `BeckeXCParams` 40/2/29); cost a bogus 40 mHa "symmetry bug" | Make `cellKind` unsettable alone so asking for Becke means the recipe (compile-time) · H2 "cellKind=Becke" |
| **§I.1** | `GetEpsXc()=0.75*GetVxc()` base default (`ExchangeFunctional.C`) exact for Dirac only — silent-wrong inherited default (LSP) | With whatever touches the functionals first (GGA) · H2 header L60–75 |
| **R2.18** ⏸ | Public vs protected `Make*` — no policy yet | CLAUDE.md "say WHY at the declaration" accrues the evidence · H2 R2.18 |
| **D9** | ρ̃ mixer FUSES preconditioner with extrapolator (`PulayMixer` owns Kerker + history; `PolarizedDensityMixer` is a `tDensityMixer`).  End state: residual → per-channel preconditioner → ONE joint extrapolator | Lands WITH N3 / +U; Broyden drop-in; re-measure run 11's (ρ,m) verdict · H2 D9 |
| **D4** | Basis carries the run's SYMMETRY POLICY (`GetDetectedReciprocalOps` on `tBasisSet`, `SetSymmetryOps` on GPW evaluator, `imposeSymmetry` per factory; 9 sites) | Collapses to one `SymmetryPolicy` context when `Records/SymmetryUpgradePlan.md` §3 is plumbed · H2 D4 |
| **D3** | `FIT_SF_ABS::SymmetrizeRaster` — symmetry op on a fit-basis face (SRP; 3 sites) | Needs a τ-acting direct-raster `FoldGrid` variant; ctor-inject into `tComposite_CD`, delete the virtual · H2 D3 |
| **D2** | Why is `FourierCD` different from `tChargeDensity`? (pin 6); seeds are a third leg beside `DM_CD`/`FittedCD` | Decide with R1.0 metric face and V1.1b · H2 D2 |
| **D13** | Single-species `HGH_LocalPotential`/`HGH_SeparablePotential` silently IGNORE `int Z` (foot-gun behind the retracted "V_long defect") | Drop `Z` from single-species API; only `MultiSpecies_*` models the Z-dispatching face (or store+throw) · H2 D13 |
| **V1.19** ⏸ | Seed's flip-group sub-cell duplication needs per-SITE form-factor overload (weighed against pseudo-wall pin) | Revisit only for a second consumer · H2 V1.19 |
| **D-PW-FACADE** | `SolidCalculation` is built over a Gaussian basis; the 2 PW SCF drivers live in `IntegrationTests/PW/Harness.C` | grow a basis-family axis on the facade, then retire the harness (moved from OpenWork §3) |
| **D-GPWPAIR** | the GPW (orbital family × fit family) pairing is FROZEN in the engine (`GPW_Evaluator` owns its PW density-fit grid; pin 2: pairings are POLICY at the factory) | becomes a defect the day a second fit family wants the same orbital engine |
| **R1.0n/o** | basis-side nullable vendor: `GetRhoOnGrid` signals "none" by empty vector | needs a ruling · H2 R1.0n/o |

## Durable ruling kept here (until non-collinear work starts, then a pin)

**V2.1:** represent spin as SU(2)/matrix, not SO(3)/direction — real diagonals + ONE complex off-diagonal
(4 real DOF = (n, m)); collinear = diagonal case; symmetry acts as UρU† and composes; never let a `Spin`
label leak into an XC functional's contract.  Rulings cited from code (read in `H2`/`H3` by id): `R1.0`
(fit basis provides integrals over ITSELF only; nobody holds both halves of a forward/adjoint pair),
`R1.0h`, `R1.0j` (never push ρ back into the terms), `R2.5b` (`throw` is a marker), `V1.20` (`.Internal.` =
FAMILY boundary), `V1.37` (Pol/UnPol are imposed subgroups).
| **V-CALCNET** | The `{Atom,Molecule,Solid}Calculation` class network is an EMERGING phenomenon, not a designed one.  Open placement questions it raises: where PP data (valence counts, basis provenance, seed variant) is read -- the Calculation level (current plan for `pseudopotentials.json`, `doc/CleanCode.md` D-STRUCTDATA step 2) may be too HIGH; it may be an implementation detail of the Hamiltonian PP terms or of `BasisSet::Orbital_PP` (user, 2026-10-02) | **DEFERRED by user**: do a full design review of the whole Calculation network AFTER the major `doc/OpenWork.md` §2 features land; do not litigate piecemeal.  Until then D-STRUCTDATA step 2 may proceed at the Calculation level knowing it can move |
