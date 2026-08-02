# Cleanup candidates (Claude's running list for the SOLID/OOD/CleanCode session)

Things noticed in passing while adding features — flagged here instead of fixed inline, so the
refactoring session can batch them.  (User keeps the master list; merge freely.)

## Fit-basis / mesh seams (from the 2026-08-01/02 symmetry+XC campaign)
- **`FIT_SF_ABS::SymmetrizeRaster`** — symmetry op living on a fit-basis interface.  Removable once a
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
