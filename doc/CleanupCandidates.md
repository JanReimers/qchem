# Cleanup candidates (Claude's running list for the SOLID/OOD/CleanCode session)

Things noticed in passing while adding features — flagged here instead of fixed inline, so the
refactoring session can batch them.  (User keeps the master list; merge freely.)

- I am mostly concerned about abstract interfaces.  These are what the client code for any module is supposed to see.
The code is an attempt (probably never been done before) to capture one set of high level requirements for any structure: atoms, molecules and 3D lattices (also 1D polymers, 2D graphene).  In other words 95% of the high level abstract interfaces (Charge density, Hamiltonian, Orbitals, Wavefunction, SCF iterator, Accelerators) should be structure neutral.  Even the basis set interfaces (qcBasisSet) are structure neutral.  
- I have so far identified 4 libraries that have a structure neutral face with specific structure specializations:
  1. qcStructure: Structure has derived classes for Atom, Molcule, UnitCell.  Lattice_3D and ReciprocalLattice have no inheritance relation with Structure.
  2. qcSymmetry: Separate folders for Atom/Molecule/Lattice_3D symmetry types.  Spin degrees of freedom are orthogonal to spatial symmetry ... but that only works for non relativistic irreps.  
  3. ElectronCOnfigurations: We need separate classes for non-aufbau (specific # of electrons per irrep) filling.  We also have separate classes Atom, Molecule and Crystal aufbau filling ... can these be combinedid into one generail Aufbau filler?
  4. qcBasisSet->qcAtom_BS,qcMolecule_BS, qc_Lattice_3D_BS
  - With the introduction of lattice calculations some structure specific classes have been creeping into the qcChargeDensity and qcHamiltonian libraries.  Also possibly into qcWaveFunction and qcOrbitals.  If we can find ways of getting these last for back to more generic status that would be a big win.



## SOLID principles 

1. SRP: Single Responsibility Principle
  - This principle is more about concrete classes than it is about interfaces.
  - We have a group of classes called Evaluators that do the low leveer work evaluating all required integrals, op(r)/grad(r) and a couple of other minor things, for a particulkar basis function type.  These have multiple responsibilities but they are mostly built up through a network of mixins in order achieve the end result.
2. OCP: Open Closed Princple
  - We are still in the R&D stages for this so abstract interfaces will be modified.  We just need a good reason to do so.  
  - For example if the extend and structure neutral interface in order get some Lattice_3D feature working we need ask the question: What does this change mean for atoms and molcules?
 
3. LSP: Liskov Substitution Principle
  - virtual functions that default to some sort of "not implemented" behaviour, violate the LSP.
4. ISP: Interface Segregation Principle
  - There is always an urge to add getters and setters.  For setters, always ask: Is this something I can set and construction time?  For getters always ask: What are we going to do with the getter data?  Why not ask the owning class to do that task instead (maintain encapsulation).  
  - I do not claim to that all the abstract interfaces in qchem are truely segregated.  A good example is the irrep basis set interfaces:  They do three things:
    a. Evlautate integrals.
    b. Evlauate op(r) and grad(r)
    c. Expose the irrep symmetry, as an abstract interface pointer.
  - Is it useful to break these up in order to faithfully obey the ISP?  I don't know but my guess in no ... or at least that we have "bigger fish to fry".
  - What we do need is that anytime our abstract interfaces get augmented those new functions should be part of an existing responsibility, not introducing a new one.
5. DIP: Dependency Inversion Princple
  - THis is possibly the most powerful concept, and often eneables adherance to the previous 4 princples.
  - This applies equally well for C++ classes, C++ modules (DAG enforced by compliler) and libraries (DAG enforced by linker).
  - example: Lib A depends on lib B.  B needs a way to send info to A.  Create an abstract interface AI (ha ha)in B, Classes in A derive from B::AI, when instances of the A classes are passed to B (as an AI*) B can then call back into A. 
  - src/ChargeDensity/ChargeDensity.C tStatic_CC, tDynamic_CC are a working example that the Hamiltonian library classes derive from and pass back into the ChargeDensity library.
  - As well as being a dependency inversion, this probably has a "Gang of Four" pattern name.

There are 6 other principles related to package cohesion and coupling.  My experience is that the first 5 are the most commonly misunderstood and misapplied.  


## Fit-basis / mesh seams (from the 2026-08-01/02 symmetry+XC campaign)
- **`FIT_SF_ABS::SymmetrizeRaster`** — (SRP) symmetry op living on a fit-basis interface.  Removable once a
  τ-acting DIRECT-raster fold variant of `FoldGrid` exists (T3 groundwork): precompute the raster fold
  and ctor-inject it into `tComposite_CD` beside `itsPointOps`, then delete the virtual (uniform route's
  voxel-shift moves out of `PlaneWaveFit_IBS`).
- **`Fit_IBS::SetMesh`** — two-phase construction (create, then SetMesh).  Should be a ctor parameter
  (build the mesh first, hand it in), per the construction-time principle settled for the XC quadrature.
- **The basis-as-policy-carrier virtual family** — `GetReciprocalPointOps` / `GetDetectedReciprocalOps`
  on `tBasisSet` + `SetSymmetryOps` on the GPW evaluator: three op-flavored accessors that exist because
  the run's §3 symmetry policy rides the basis.  When the run-level policy object (plan §3 /
  `Lattice_3D`-owned symmetry) is fully plumbed, these can collapse to consumers reading ONE context.
- **`GPWParams.imposeSymmetry` (bool; renamed from `reduceBZ` 2026-08-02, user request — the switch
  imposes the whole space group, not just the k-fold)** — still the §3 assert-switch as a bool on one
  basis family's params; should become the shared `SymmetryPolicy` once PW/APW/LAPW factories read
  `Lattice_3D::GetSpaceGroup()` + a policy, not a per-factory flag.
- **Becke `MeshParams` ergonomics under `imposeSymmetry` (user, 2026-08-02)** — `mp.nAngular` IS
  honored (the exactness degree L on both routes) but `mp.angular` (the GL/Lebedev scheme) is silently
  REPLACED by the site-adapted rule (announced, not warned).  Wanted: `nAngular<=0` → auto-resolve the
  calibrated default; a console warn when the scheme is overridden; a real error (not the bare assert)
  when the requested L is unachievable for a low-symmetry site (C1/Cs seed-pool exhaustion).  Land with
  the `SymmetryPolicy`/facade pass.
  **GREW A BLOCKER ROLE (2026-08-02):** `nAngular` is a COUNT for Lebedev but a DEGREE for GL/EM (and
  the imposed site-adapted builder consumes it as the degree) — this dual semantics is what BLOCKS
  flipping the free-run Becke default to the measured-equal Leb-302 (67% of GL-29's directions; plan
  §6a).  The fix wants a degree-typed interface (`angularDegree` + per-scheme count resolution), then
  the default flip rides along.

## Spin-native solid pipeline (from the 2026-08-04 tier-4b campaign)
- **`Delta_XC` vs `Delta_XC_Pol`+`Delta_VcorrPol` — the unpol pair should become the ζ=0 collapse.**
  Per the spin-native-is-primary bias: today the Becke route has THREE term classes (unpolarized
  Delta_XC×2, polarized Delta_XC_Pol + Delta_VcorrPol) where the polarized pair evaluated at
  ρ↑=ρ↓=ρ/2 IS the unpolarized answer (the singlet-collapse gate proves it to 0.12 mHa).  Collapsing
  to the polarized pair only (unpol = two identical channels) halves the class count at ~2× the XC
  pointwise work for closed shells — the same trade Ham_PP documents.  Decide with a perf measure.
- **`XC_GridEngine` now carries TWO rho caches** (scalar `itsRho` for the unpol route, the `{Up,Dn}`
  pair for RhoPol) — falls out for free if the previous item lands (pair only; total = up+dn).
- **Polarized PLANE-WAVE Vxc fit route** — `Ham_PW_DFT` polarized currently THROWS for
  `VxcFit::PlaneWave`: per-channel PW_XC needs per-spin rho-grid caches (PW_XC's `itsRhoGrid` is
  keyed on `cd->Version()` alone, which a polarized density aliases across channels — the same trap
  the engine's RhoPol pair-cache fixes).  Design note in the throw message.
- **`CollocMemo` dual duty** — the replay memo (exact-D level densities) and the adjoint D-screen
  now live in one struct with different lifetimes (last-D vs union-Dscr).  Fine today; if a third
  consumer appears, split the two concerns.
- **GPW default seed policy** — `GpwOptions.seed` defaults to `Uniform`, which the Na-doublet
  campaign showed has a STABLE wrong basin for electron-sparse systems (the lone-electron doublet
  converged 72 mHa above the true minimum with every health metric green).  The molecular facade
  already defaults DFT runs to SAD.  Candidate: default GPW to `IonicSAD` (or SAD-family), with
  Uniform as the explicit opt-in — needs a suite sweep since every pinned GPW anchor would re-seed.

## Spin-SAD seed campaign (2026-08-04, SCFSeedingPlan §10 increments A/B)
- **Polarized ATOMIC solver fails on an empty minority channel** — `AtomCalculation(11, 0,
  {.pseudopotential=true, .valence=1, .pol=Polarized})` (Na q1 doublet: nUp=1, nDown=0) dies with
  "Invalid setup of symmetric matrix" before the first Fock.  Same failure class as the tier-4b
  driver finding (`cSCFAcceleratorDIIS::GetNProj` segfault on an empty EC).  Worked around for the
  Na library entry (1 valence electron ⇒ the pair is EXACT without an SCF: rho_up=rho, rho_dn=0,
  hand-constructed in atomic_valence_densities.json), but any future fully-polarized one-electron
  species hits it again — the valgen `--spin` path should either fix the empty-channel atom SCF or
  synthesize the exact pair itself when nDown==0.
- **Unoccupied-l shells break the atom SCF** — the same Na valgen run ALSO failed when the recipe
  carried the (unoccupied) p window from valence_lowq_sr (`--shell 1:2:0.09:0.3`): "Invalid setup
  of symmetric matrix" even UNPOLARIZED-adjacent.  GenerateValenceBasis documents "higher-l
  polarization shells ride along un-validated (as intended)" — but at least the polarized
  GenerateSeedDensity path chokes on them.  Reproduce + fix or document the restriction.
- **`SeedCD` flip-group sub-cells duplicate the Structure** — the per-channel AFM assembly clones
  the UnitCell into (unflipped, flipped) groups because `MakeFourierDensity(st, formFactor(Z,g2))`
  is species-keyed.  If a third per-atom attribute ever needs G-assembly (site charges, isotopes),
  consider a per-atom-index form-factor overload on `G_FieldEvaluator` instead (basis-interface
  change — weigh against the pseudo-wall pin).


## qcBasisSet* library family
- We have a number of Orbital_*_IBS interfaces that mostly define what types of integrals are required for 1E/DFT/HF/DHF calculation.  I think would help to define an Orbital_PP_IBS interface defining the types of integrals required specifically for PP calculations.  This might be a superset of the Orbital_DFT_IBS interface.
- Band_DFT_IBS is this no longer used?
- For lattice we created src/BasisSet/Band_FT_IBS.C which does almost the same thing as the generic src/BasisSet/Orbital_DFT_IBS.C that works for atoms and molecules.  If we are able to merge Orbital_DFT_IBS and Band_FT_IBS then that proves that there was no need to create a new separate interface just for lattice basis sets.
- SymmetryAdapted_IBS this looks like it could be moved into the Internal folder.
- G_FieldEvaluator this is tightly coupled to 4 consumers DensityMixing, SeedCD, Fitting, and PWVxcField
  - Do all these consumers use the full class, or do they just need a small interface? In other words is the there DIP improvement available?
- src/BasisSet/Lattice_3D/BandStructure.C only consumers seem to be unit tests.
- BasisSet::GetReciprocalPointOps Only consumer is tCompositeWF<T>::GetChargeDensity(Spin s). This high high level structure neutral code calling a lattice specific member function in structure neutral BasisSet. 
  - What is tCompositeWF trying to do, and how would transfer to molcules and atoms?
- BasisSet::GetDetectedReciprocalOps() const {return {};} only consumer is a unit test.


## qcChargeDensity library
OOD: Charge density has become a complex forest of classes.  
### Aspects:
- Two primary varients:
  1. DM_CD which is capbale of 2ERIs.  THis is the first stop after any SCF iteration, even if we never ask for HF-ERIs.
  2. FittedCDs for DFT, Not HF-ERI capable.  Fit basis = Gaussian/PW/Delta are all potential fit basis sets.  PW and delta are not conventionaly thought of as fit basis sets (ortho->no metric, trivial blah blah).  But we can capture common behaviour in software if we shirk convention in this case.
  - Both are ScalarFunction<T> for op(r)/grad(r).
  - Other than ScalarFunction<T> DM_CD and FittedCD do not share a common "charge density" base.
  - Currently FittedCDs are very transient: Create,DoFit,GetRepusion(orbtial basis set).
- Polarized/Unpolarized: This should be orthogonal to DM/Fitted aspects.
- Composite (multiple irreps).  Is this really different than Polarized (irreps can hold spin)?
- Version traking
- Mixing algos.
- SPin density (for plotting)
- Seeding: NumericalCD and SeedCD: No density matrix.

### Issues
- Why is a FourierCD different than a tChareDensity?
  - At a high level Fourier/PW is an implementation detail.
  - Consumers: src/Hamiltonian/Internal/Imp/PWTerms.C, src/SCFIterator/Imp/SCFIterator.C
  - Granted that for efficiency return types for FourierCD may require different data structures than Gaussians. src/BasisSet/Internal/GMap.C G_ERI3 is a multi purpose data structure that attempts to work around this problem.  It is not a perfect solution but it allows us to harmonize interfaces and move in SOLID-clean direction.  G_ERI3 ands up being a spec for whats needed in terms of data structure flexibility.  In the end high level code should not have to worry about G or no G, PW, Fourier.
- tChargeDensity::EvalBatch just duplicates ScalareFunction::operator()(const rvec3vec_t& rs) const
- LSP violation: tDM_CD::DM_ContractBlocks
  { assert(false && "DM_ContractBlocks: only a finite (Gaussian) density matrix implements this"); return 0.0; }
- tDM_CD::DM_RhoAtPoints(const rvec3vec_t& r, const std::map<std::string,mat_t<T>>& /*Phi*/)
  Phi uses BasisSetID as map key.  Why not irrep?
- tDM_CD::Accumulate*  4 functions ... all LSP violations.  They should be pure virtual.

## Hamiltonian

- I tried to use some conventions for naming terms:
  - Enn/Vnn  nuclear-nuclear repulsion
  - Een/Ven  electron-nuclear attraction
  - Eee/Vee  electron-electron repulsion
  - Eex/Vex  exchange
  - Ecorr/Vcorr correlation
  - Exc/Vxc  exchange-correlation
- For solids these don't work perfectly because at G=0 we can only define certain combinations of Enn+Een+Eee as finite.  Also PPs are a combination of Ven+Vee + KB projectors. I use Ven for PPs, where n is understood to be a shielded nucleus.
- Current names in src/Hamiltonian/Internal/Terms.C
  - PP_Local: Can we call this Ven_Local?
  - PP_NonLocal: Can we call this Ven_NonLocal?
  - FittedEpsXc: Can we call this FittedExc:
    -Only one consumer: FittedVxc can we merge FittedEpsXc into FittedVxc?
  
- Current names in src/Hamiltonian/Internal/PWTerms.C
  - PW_Pseudo  I don't like PW prefix. Can we call this Ven-pp ?
  - PW_Kinetic this class does not do any thing different from the generaic Kinetic term. bs->MakeKinetic(); knows how to do the PW kinetic integrals.  
  - PW_Hartree I don't like PW prefix. Can we call this Vee-pp ?
  - PW_XC: What is distinguising this from FittedVxc?  Is it PWs?  I think GPW also uses this class.
  - XC_GridEngine I don't know what this is.  Engine in the names suggest it is going to do some work ... but that does not tell me anything useful.  "The Becke route is a pure QUADRATURE " this is conflating unit dell grid layout with fit basis for Vxc.  "pure QUADRATURE" = fit basis is just (pseudo) delta functions. 
  - Delta_XC: This is a FittedVxc using delta functions as the fit basis. Mayde call it DeltaFittedVxc.
    - I understand that it makes no sense to try an pump a delta function basis set through the generic FittedVxc plumbing, it needs special implementation for efficienciy.
    - BUT maybe special handling can be moved into the Fitting library.  This will require deeper analysis.
    -Why do we need separate "Engine" class.  Can we merge DeltaFittedVxc and XC_GridEngine?
    - All refernces to the Becke grid/Mesh should be eleminated from both of these classes ... they should be grid/Mesh neutral.
  - Delta_XC_Pol:  Can we make a generic pol Vxc that holds to abstract base FittedVxc* and that covers Molecules/PW/GPW/DeltaFittedVxc etc.
  - Delta_VcorrPol: Again generic?
  

## Older / unrelated spots hit while working nearby
- **`DB_Cache_RAM.C`** — a screenful of `-Winconsistent-missing-override` warnings on every qcBasisSet
  build (`Get`/`Register`/`GetCache*`).  Mechanical `override` sweep.
- **`BeckeXCParams()` lives in the TEST file** (GPW_SCF_UT.C anon namespace) but is the de-facto
  PRODUCTION Becke recipe (`ResolveXCMesh` default).  The calibrated recipe belongs in the library
  (e.g. beside `MeshParams`), tests reading it — not the other way round.
- **`ResolveXCMesh` (test driver)** — run-policy resolution (Auto grid × reduceBZ interplay) living in
  the integration-test harness; belongs with the facade/driver once the policy object exists.
- **`UnmatchedCounts`/fold `tol` defaults** — 1e-8 fractional appears as a literal in three places
  (GPW `CreateXCQuadrature`, SymmetrizeMesh overloads); name it once.
- **`Symmetry::Lattice_3D::DirectOf`** — currently unused after the CreateXCQuadrature move (GPW uses
  its native direct ops).  Keep (documents the U=Wᵀ convention; T3 will want it) or fold its doc into
  `ReciprocalOp` and drop — decide at the refactor session.
- **`MakePeriodicBeckeMesh` ε-tail drops vs orbit consistency (W2c find)** — the builder's
  borderline drop decisions (`<eps` screens + `w>0` keep) are per-point and bit-sensitive, so the
  site-adapted caller post-filters orbit-incomplete points (`CreateSiteAdaptedBeckeMesh`).  Cleaner
  long-term: make the drop decision ONCE per representative (angular dir × radial shell) and apply
  it to the whole atom orbit inside the builder — removes the second fold pass + the filter.
- **`SolveSPD`/NNLS in `SymmetrizeMesh.C`** — module-private dense Cholesky + Lawson-Hanson NNLS
  hand-rolled beside the mesh code; if another consumer of small dense LS/NNLS appears, promote to
  qcMath (Blaze has no NNLS).
- **`GDMParams::FDMax` naming + the fallback commit (2026-08-03 imposed×GDM investigation)** —
  (i) `FDMax` reads as a step size but is the ENGAGEMENT gate on ‖[F′,D′]‖ (the geodesic step size
  is the quadratic-model `itsStdef` capped by `Trust` radians); rename to something like
  `EngageBelowFD` and consider making the norm intensive (per-element or per-electron) so one value
  transfers across basis sizes.  (ii) `DirectMinStep`'s 12-backtrack line search COMMITS a tiny
  non-descent step on fallback (SCFIterator.C ~L424) — the measured uphill leak on imposed NaF-SR
  (+23–56 mHa over 100 iterations along projector-curved diffuse directions).  Fix: hold position on
  fallback (or accept `best` only within a noise floor); pair with soft-direction preconditioning —
  the 1/(ε_a−ε_i) diagonal Hessian blows up the step exactly along the near-degenerate diffuse
  modes.  Reproducers: DISABLED_ImposedGDMProbe_SiDiamondIBZ (healthy), DISABLED_NaFImposedGDMSmearProbe
  (pathological, NAFGDM_* knobs); GPW_GDMTRACE=1 shows DESCENT/FALLBACK per step.
- **BZ creep on the neutral `Symmetry` base (ISP) — DONE 2026-08-04 (same day): the ATOM SHELL
  CONVENTION landed.**  `BlochQN::GetDegeneracy()` = star (ctor converts the k-mesh layer's w_k;
  asserted integer), `GetWeight()` = uniform 1/N_mesh, `StarSize()` + the `EnergyLevel` report
  scaling DELETED (the 8/8 wedge display is now the physical occupation), `Crystal_EC::GetN` =
  star×Nval per block in insulator mode / plain Nval for the global-μ total.  No factory/KBlock/
  call-site changes (BlochFactory keeps taking w_k); free meshes (star=1) bit-identical; g=w·degen
  invariance verified (Al metal); imposed Si diamond IBZ reproduces E to all digits; 656/656.
  The base is down to TWO Bloch-flavored members (GetWeight, MergeAcrossIrreps) — the remaining
  capability-face split is now small enough to fold into any future Symmetry-base touch.
  *(Original analysis, for the record:)*  The base now carries three Bloch-flavored defaulted virtuals (`GetWeight()`
  pre-existing; `MergeAcrossIrreps()` + `StarSize()` from the 2026-08-03 MergeTol/IBZ-table fixes).
  Rather than split a capability face, adopt the user's design: a wedge k-block is a SHELL exactly
  like an atom's l-block (one stored representative, degeneracy = the symmetry copies), so
  **`BlochQN::GetDegeneracy()` = star size** (spatial; spin stays layered by `Irrep`) and
  `StarSize()` is DELETED.  This is a coordinated convention switch — the invariant is
  w_k × (per-block quantity), so the star factor must leave the weight in the same commit:
  `GetWeight()` becomes plain 1/N_mesh (derivable from `N` — most of the weight plumbing through
  `ReduceToIBZ` → `KBlock` → `BlochQN` evaporates); `Crystal_EC::GetN` per wedge block becomes
  star×Nelec (all star copies fill through the one block — each band then holds star×2 natively and
  the level table needs NO report-layer scaling: revert the `EnergyLevel` ctor scaling);
  rebalance every w_k consumer (density sums, `tIrrepWF::GetEntropyTerm`;
  `FillOrbitalsGlobalFermi`'s g = w·degen is INVARIANT under the switch — sanity check).
  `MergeAcrossIrreps()` stays (unfolded meshes still carry star partners as separate blocks).
  Re-verify every IBZ energy anchor.  Dedicated refactor session — multi-seam.
