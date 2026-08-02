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
- **`GPWParams.reduceBZ` (bool)** — the §3 assert-switch as a bool on one basis family's params; should
  become the shared `SymmetryPolicy` once PW/APW/LAPW factories read `Lattice_3D::GetSpaceGroup()` +
  a policy, not a per-factory flag.

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
