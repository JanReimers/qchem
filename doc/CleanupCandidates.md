# Cleanup candidates (Claude's running list for the SOLID/OOD/CleanCode session)

Things noticed in passing while adding features — flagged here instead of fixed inline, so the
refactoring session can batch them.  (User keeps the master list; merge freely.)
## SOLID principles 

- I am mostly concerned about abstract interfaces.  These are what the client code for any module is supposed to see.
The code is an attempt (probably never been done before) to capture one set of high level requirements for any structure: atoms, molecules and 3D lattices (also 1D polymers, 2D graphene).  In other words 95% of the high level abstract interfaces (Charge density, Hamiltonian, Orbitals, Wavefunction, SCF iterator, Accelerators) should be structure neutral.  Even the basis set interfaces (qcBasisSet) are structure neutral.  


1. SRP: Single Responsibility Principle
  - I do not claim to that all the abstract interfaces in qchem are truely single responsibility.  A good example is the irrep basis set interfaces:  They do three things:
    a. Evlautate integrals.
    b. Evlauate op(r) and grad(r)
    c. Expose the irrep symmetry, as an abstract interface pointer.
  - Is it useful to break these up in order to faithfully obey the SRP?  I don't know but my guess in no ... or at least that we have "bigger fish to fry".
  - What we do need is that anytime our abstract interfaces get augmented thos new functions should be part of an existing responsibility, not introducing a new one.
2. OCP: Open Closed Princple
  - We are still in the R&D stages for this so abstract interfaces will be modified.  We just need a good reason to do so.  
  - For example if the extend and structure neutral interface in order get some Lattice_3D feature working we need ask the question: What does this change mean for atoms and molcules?
  - For lattice we created src/BasisSet/Band_FT_IBS.C which does almost the same thing as the generic src/BasisSet/Orbital_DFT_IBS.C that works for atoms and molecules.  If we are able to merge Orbital_DFT_IBS and Band_FT_IBS then that proves that there was no need to create a new separate interface just for lattice basis sets.
3. LSP: Liskov Substitution Principle
4. ISP: Interface Segregation Principle
  - There is always an urge to add getters and setters.  For setters, always ask: Is this something I can set and construction time?  For getters always ask: What are we going to do with the getter data?  Why not ask the owning class to do that task instead (maintain encapsulation).  
5. DIP: Dependency Inversion Princple
  - THis is possibly the most powerful concept, and often eneables adherance to the previous 4 princples.
  - This applies equally well for C++ classes, C++ modules (DAG enforced by compliler) and libraries (DAG enforced by linker).
  - example: Lib A depends on lib B.  B needs a way to send info to A.  Create an abstract interface AI (ha ha)in B, Classes in A derive from B::AI, when instances of the A classes are passed to B (as an AI*) B can then call back into A. 
  - src/ChargeDensity/ChargeDensity.C tStatic_CC, tDynamic_CC are a working example that the Hamiltonian library classes derive from and pass back into the ChargeDensity library.
  - As well as being a dependency inversion, this probably has a "Gang of Four" patter name.

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
