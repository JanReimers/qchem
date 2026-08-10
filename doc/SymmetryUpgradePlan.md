# Symmetry Upgrade Plan — PW & GPW SCF

*Plan first, code second.* Scope: exploit lattice + point (and eventually magnetic)
symmetry across **all** the redundant work in a PW/GPW SCF, not just the k-mesh.
Companion to `doc/SpaceGroupPlan.md` (the detector + Tier A/B roadmap),
`doc/GPWPlan1.md` (the "SPACE-GROUP STREAM/COLLOCATION REDUCTION" plan section this
consolidates), `doc/SpinNativeDFTPlan.md` (spin-native XC), and
`doc/SymmetryRefactorPlan.md`.  SOLID/OOD debt encountered while executing this plan is
tracked in **`doc/CleanupCandidates.md`** (keep it growing; batch-fix in dedicated
refactor sessions).

**SESSION CLOSE (2026-08-02, second session) — WHERE THINGS STAND:** §7 steps 1–5 DONE —
step 5 (T2) CLOSED OUT this session: the MIXED special-orbit rule + NNLS landed
(`MakeInvariantAngularMesh(ops, L, avoid)`), the production-L growth re-measure came in at
the promised ~2× (1.97× mesh / 1.95× angular at L=29 vs GL-29), and the **`reduceBZ`
Auto→uniform carve-out is RETIRED** (Auto+reduceBZ now builds the mixed-rule site-adapted
invariant Becke mesh; Al re-pinned route-matched −2.1174805, reduced==full ~1e-5; Si
diamond IBZ passes on the old pin).  TWO DISCOVERIES en route (details in the W2c block
below): (i) plain LS weight-solving FAILS outright at L=29 (residual 1e-11 but a ~−5e-4
negative-weight floor persists to any pool size) — NNLS is a *requirement*, not an
optimization; (ii) the Becke builder's ε-borderline tail drops flip inconsistently across
orbit partners at scale (212 orphans at nR=40/L=29) — cured by the orbit-consistency
filter in `CreateSiteAdaptedBeckeMesh`, gated at production-L.  Step 6 (T3 streams) partly
done (T3.0–T3.2, T3.4 op-action; T3.4b open); **step 7's tier-4b prerequisite is now DONE
(2026-08-04, 64a17443 — §4 STATUS block; all three gates green)**, so step 7 (MnO/MnO₂
magnetic materials) is UNBLOCKED — remaining prerequisite for an AFM ordering is the
spin-polarized SAD seed with per-site moments (doc/SCFSeedingPlan.md §10); 8 (SALCs) not
started.  The XC route was redesigned along the way
(fit/grid separation: `VxcFit` ⊥ `MeshParams`; `Delta_XC` + `XC_GridEngine(mesh, fold)` +
the `CreateXCQuadrature` basis factory).
**Open next increments, in rough value order:** ~~(b) the free-run ROTATED-Lebedev
experiment~~ **DONE 2026-08-02** (§6a: alignment confirmed as the whole story, 5–10×
recovered by `angRot`; default stays GL-29 — the audited Lebedev tables stop at degree 11);
(c) T3 stream fold + the §5 commensurable raster menu (also delivers the raster fold that
retires `FIT_SF_ABS::SymmetrizeRaster`; CAMPAIGN SCOPED — see the §7 step 6 block: T3.0
design memo first, incl. the does-route-(b)-even-need-§5 question); (d) I3's one-functional
E/H design for the (PlaneWave, Becke) cell; ~~(e) tier 4b (Na-doublet/O₂-triplet boxes)~~
**DONE 2026-08-04** (§4 STATUS) → step 7 next, needing spin-SAD (SCFSeedingPlan §10);
~~(f) higher-degree Lebedev tables + the `AngularKind::Gauss`→`Lebedev` rename
(user, 2026-08-02)~~ **DONE same day** — rename + canonical Lebedev–Laikov import for
orders 86–434 (degrees 15–35; 74/230/266 excluded, negative weights), audit-gated;
MEASURED Leb-302 ≈ GL-29 at 67% of the dirs (§6a block).  The FREE-RUN DEFAULT FLIP to
Leb-302 is the remaining follow-on, blocked on the `nAngular` count-vs-degree semantics
(CleanupCandidates ergonomics pass).
**Cleanup debt found while feature-building goes to `doc/CleanupCandidates.md`** (standing
practice, user 2026-08-02) — keep it growing rather than fixing inline.  Detailed
per-increment records follow below.

**STATUS (2026-08-01):** §7 **steps 1+2 BUILT** — `qchem.Symmetry.Lattice_3D.Fold`
(`FoldPoints`/`FoldGrid`/`FoldGVectors`, per-member op edges, `SymOp {W|τ,σ}` with the
tier-4a `SpinAction` enum), the `SpaceGroup::FoldKMesh/FoldGVectors/FoldPoints` member
wrappers, and `ReduceToIBZ` re-expressed on `FoldGrid` (bit-identical, gated by
`Fold.KMeshFoldBitIdenticalToReduceToIBZ`).  Unit tests: `src/Symmetry/tests/L_Fold.C`.
§7 **step 3 BUILT** — `SymmetryPolicy` (impose-on-assert; `GPWParams.reduceBZ` [renamed
`imposeSymmetry` 2026-08-02] documented as the assert-switch, resolved ONCE in
`DetectPointOps` which now always detects and feeds policy + all op faces from one place), `SymmetryDefects` per-op order-parameter diagnostic
(GMap, G-space so commensurability-free), and the never-silent `[symmetry]` run line
(free run: max defect + which-ops-broke; imposed run: the release-audit assertion).
Gates: Al symmorphic + Si non-symmorphic IBZ still exact; unit negative control fires on a
broken ρ̃ (stabilizer ops stay clean = the which-op readout) AND demonstrates the
projector-blindness (post-`SymmetrizeGMap` defect ≡ 0 — §3's imposed-run caveat, measured).
MEASURED calibration: a converged (Δρ≈2.5e-6) free Si Γ run shows max defect ~1e-5 — the
through-SCF defect tier tracks the CONVERGENCE tolerance (per §8), so the run-line BROKEN
alarm sits at 1e-3.
§7 **step 4 (T1) BUILT** — `EvaluateSymmetricGMap` in GMap (fold the G-list via
`FoldGVectors`, evaluate the field at star reps, scatter members through the exact glide
identity f(Um)=e^{+2πi m·τ}f(m)); consumed by BOTH GPW structure-factor sweeps
(`MakeLocalPP` + the `MakeLocalPPLong` custom G-ball), §3-policy-gated (ops empty on free
runs; the imposed ops now live ONCE on `GPW_Evaluator::SetSymmetryOps`, shared with the
Vxc-raster star-average).  Exact for the static V_loc by construction (the ops are detected
from the very atoms the sum runs over).  Gates: unit reduced==full at 1e-13 on the
non-symmorphic FF×structure-factor field WITH rep-only call count; Al symmorphic + Si
diamond non-symmorphic IBZ energies unchanged through SCF.  Fold edge-op recording
optimized (BFS-discovery capture; the O(members×ops) search now only runs for non-closed
hand-made op subsets).  Full suite 638/638.
§7 **step 5 (T2) GROUNDWORK BUILT** — the geometric layer, no XC wiring yet:
`FoldPointsPeriodic` + `CountUnmappedPeriodic` (torus-metric point fold + per-op invariance
checker, bucket-hashed) and the pointwise star-average `SymmetrizeValues` (orbit mean == the
exact projector on an invariant set under a full group) in qcSymmetry; the Mesh-typed layer
in `qchem.BasisSet.Lattice_3D.Internal.SymmetrizeMesh` (**qcMesh deliberately stays a pure
quadrature leaf, uncoupled from qcSymmetry — user decision 2026-08-01**): `FoldMesh`,
`Symmetrize(Mesh)` (summed member weights — exact for a symmetric integrand even with
unequal star weights), `UnmatchedCounts`, and `MakeInvariant` (group-average any mesh into
an op-invariant one at preserved Σw and polynomial exactness — the generic route to the
site-adapted grid; Symmetrize folds it back to ~the original point count).  Unit-gated
(SymmetrizeMeshUT: invariance cured, 1e-13 fold exactness, orbit-mean projector idempotent,
torus boundary matching).  Full suite 643/643.
**T2 XC-wiring DESIGN RESOLVED (§6a, 2026-08-01):** the site-group ↔ angular-grid
connection IS the invariance precondition (user question answered: NOT independent), split
into **W1** (verification: `MakeInvariant` mesh + ρ star-average in `BeckeXC_Engine::Rho`;
the W1 adjoint is PROVED already-exact — P self-adjoint on orbit-symmetric weights + v
symmetric ⇒ the existing `Matrix()` quadrature is the exact derivative untouched; gate on a
coarse mesh, carve-out stays) and **W2** (production: per-atom-orbit site-adapted minimal
grids via `SpaceGroup::SiteStabilizer` [LANDED, diamond T_d=24 gated] — invariance at
today's point count, then the ≈|site-group| fold on ρ GEMM + functional evals; H assembly
full-mesh until the T3-shared rep transform; retires the carve-out; = the AngularMathPlan
revival).  **OWNERSHIP LANDED (§6a last bullet): qcStructure→qcSymmetry edge;
`Lattice_3D::GetSpaceGroup()` (detection + adapter, lattice-common); `qchem.SymmetrizeMesh`
re-homed to src/Structure/Lattice_3D; SpaceGroup gained ReciprocalMatrix/To{Fractional,
Cartesian}.**
**W1 BUILT + GATED (final shape — REVISED by user review, 2026-08-01):** the Becke route
has NO fit basis at all — it is a pure QUADRATURE, and the first cut's `BeckeFit_IBS`
(a zero-function pseudo-basis) + `cFIT_SF_Quadrature` face were a null-object smell,
deleted same-day.  Final design: the finished quadrature is assembled BEFORE the engine —
the Hamiltonian's Becke branch builds the mesh (`Structure::CreateIntegrationMesh`), and on
a §3-imposed run (the basis's imposed op set non-empty; direct face recovered via the new
`Symmetry::Lattice_3D::DirectOf` so the U=Wᵀ convention stays with the group)
group-averages it INVARIANT (`MakeInvariant`) and prepares the orbit fold (`FoldMesh`) —
the Structure→UnitCell cast inside this lattice-specific branch is a checked precondition
per the user's cast rule.  `BeckeXC_Engine(mesh, fold)` receives both at construction;
per iteration it applies the fold's orbit-mean to ρ (`SymmetrizeValues`) — everything
symmetry-shaped happens before the ctor, and NO fit-basis interface grew anything.  GATES:
`BeckeXC_IBZ_SiDiamond` — Becke+IBZ == Becke+full to **7e-5 Ha** (non-symmorphic, coarse
GL-9 recipe both arms; invariant mesh 1034→6120 pts = ~6× growth, 135 orbits); the free
arm's `[symmetry]` line measures max defect 6e-3 / 40 ops broken — the fixed-orientation
coarse-grid rotating-error channel, removed by the imposed arm's projector, exactly as
predicted.  Γ free-run equivalence unchanged (`BeckeXCMatchesUniformXC_SiGamma`).  The
Auto→uniform reduceBZ carve-out STAYS (W1 growth ~6×; explicit Becke opts in) until W2's
site-adapted no-growth grids.
**THE FIT/GRID SEPARATION (user, 2026-08-01):** the "Becke engine" is really the
DELTA-FUNCTION fit basis realization of FittedVxc — coefficients ARE the grid-point values,
H by direct quadrature — and the two user choices are ORTHOGONAL: `VxcFit {Auto, PlaneWave,
Delta}` (which basis represents v_xc) × `MeshParams` (which real-space grid).  Renames:
`BeckeXC_Engine`→`XC_GridEngine`, `PW_XC_Becke`→`Delta_XC` (nothing plane-wave about it);
`Auto` = the historical pairing (Delta on Becke, PlaneWave on the raster), so zero behavior
change.  The LOW-LEVEL MESH WORK moved OUT of `Ham_PW_DFT::BuildTerms` into the basis-side
factory `CreateXCQuadrature(st, mp)` (sibling of `CreateVxcFitBasisSet`, default = plain
mesh + no fold; the GPW override group-averages invariant + folds using its OWN cell and
ctor-injected ops — no casts, no `DirectOf` plumbing, the Ham branch is four lines).  I2 GATE
(`DeltaFitUniformGridMatchesPWFit_SiGamma`): the (Delta, uniform) cross cell reproduces the
PW fit to **1e-5 Ha** — the band-limit residual at this cutoff, isolated.  The (PlaneWave,
Becke) cell (I3) is asserted out until its ONE-FUNCTIONAL E/H pairing is designed: the
projection sum is trivial (the mesh drops out once c_j exist — user), but H must stay the
exact derivative of the quadratured E (the user's GDM-after-DIIS practice is the runtime
audit that would expose any mismatch).  Full suite 646/646.
**W2a BUILT — the invariant angular quadrature core** (`MakeInvariantAngularMesh(ops, L)`
in `qchem.SymmetrizeMesh`): generic Fibonacci-sphere seeds (stride-permuted so early orbits
COVER the sphere — clustered seeds left the moment system degenerate), full-size site-group
orbits only (special/bond directions skipped by construction — the Lebedev cube-corner
lesson), orbit weights from the ORTHONORMAL real-Ylm moment conditions (stable normalized
associated-Legendre recurrence, module-private; small ridged-Cholesky LS), orbits added
until residual < 1e-9 AND all weights positive.  Gated: O_h invariant + degree-9 exact on
every monomial vs closed forms, positive weights, Σw=4π, deterministic; trivial group
degenerates cleanly.  COST NOTE: generic orbits price ~(L+1)² directions vs GL's ~(L+1)²/2 —
the robustness premium; vs W1's measured ~6× group-average this is the ~2× route, and the
fold precondition is exact by construction.
**W2b BUILT — the site-adapted builder** (`UnitCell::CreateSiteAdaptedBeckeMesh(mp, ops)`):
per atom ORBIT of the imposed ops, the rep atom carries its site-group-invariant angular set
(stabilizer → Cartesian A·W·A⁻¹ → `MakeInvariantAngularMesh`), partners carry the op-rotated
copies (the atom-fold edge ops); the fuzzy-Voronoi partition runs on the final points
(geometric ⇒ weights inherit invariance).  `MakePeriodicBeckeMesh` gained an optional
per-atom angular-set parameter; GPW's `CreateXCQuadrature` routes imposed+Becke through it
(imposed+uniform keeps the W1 `MakeInvariant` path).  GATES: diamond mesh invariant under
all 48 ops BY CONSTRUCTION (`UnmatchedCounts≡0`; the plain shared-orientation mesh fails
30+ ops), Σw within 5% of plain; through-SCF Becke+IBZ == Becke+full at **5e-4** (different
rules now: degree-L-exact generic vs GL — the honest quadrature-difference class; W1's
same-rule comparison was 7e-5).  MEASURED at L=9/T_d: 240 dirs/atom = exactly the 10
invariant-harmonic orbits × 24; mesh 4718 pts / 100 orbits vs W1's 6120 (23% smaller, still
~4.6× the free mesh — the low-L generic-orbit premium; the ~2× asymptote needs production
L≈29, and a MIXED rule admitting special orbits that avoid atom-specific bad axes is the
noted refinement).  **Carve-out decision: Auto→uniform under reduceBZ STAYS** until the
growth is re-measured at the production recipe (L=29) and/or the mixed rule lands; explicit
Becke+IBZ is verified on both W1 and W2b routes.  Full suite 649/649.
**W2c BUILT (2026-08-02, second session) — the MIXED RULE + NNLS + the carve-out retirement:**
`MakeInvariantAngularMesh(ops, L, avoid)` now seeds its orbit pool with the site group's
SPECIAL orbits first (op axes via the ±1-eigenvector of each op, both ±n — opposite special
orbits are distinct without inversion, so diamond's ANTI-bond ⟨-1-1-1⟩ tetrahedron survives
the filter that kills the bonded one; mirror-plane |G|/2 semi-special orbits from projected
seeds), filtered against the ATOM's actual bond directions (`CreateSiteAdaptedBeckeMesh`
computes nearest-neighbour unit vectors, first shell ×1.05), then generic Fibonacci orbits.
Weights: **Lawson-Hanson NNLS** with per-DIRECTION cost-biased column selection (entering
column maximizes gradient/orbit-size — the Lebedev-efficiency bias) replaced the plain LS —
MEASURED NECESSARY: at L=29/T_d the LS residual converges to 1e-11 by K=46 but the min
weight plateaus at ~−5e-4 to Kmax however many orbits are added; NNLS guarantees w≥0 and
emits a sparse support.  MEASURED: T_d L=29 → 878 dirs vs GL-29's 450 = **1.95×** (the ~2×
asymptote CONFIRMED; W2b's coarse-L premium is gone — T_d L=9 → 76 dirs, L=5 → 18 = GL
parity; O_h L=9 → 4 orbits/50 dirs, Lebedev-9 class); full production mesh (nR=40, L=29,
diamond) 49172 vs 25007 pts = **1.97×**; Al O_h site 830 dirs (⟨110⟩ bonds filtered).
**THE ε-TAIL ORBIT-CONSISTENCY BUG (found by the production-L invariance gate):** the raw
site-adapted point set is op-covariant by construction, but `MakePeriodicBeckeMesh`'s
ε-borderline DROP decisions (the `<eps` screens + the `w>0` keep) are computed from
bit-different rotated distances — at ~70k raw points 212 orphaned points broke
`UnmatchedCounts≡0`.  Cure: post-filter dropping every orbit-INCOMPLETE point
(complete ⇔ |orbit|×|site stabilizer| == |ops|; only ~ε-weight points can be incomplete, so
this stays inside the ε-converged-series contract).  GATES: SymmetrizeMeshUT
`ProductionL29Growth` (growth < 2.5× pinned) + `ProductionAngularRecipeInvariantAfterTailDrop`
(UnmatchedCounts≡0 at nR=15/L=29, a recipe measured to exercise the borderline channel — 14
orphans cured) + the L=9 through-SCF gate re-measured (ΔE=−2.0e-3, the RULE-difference
class GL-9 vs minimal-mixed-76-dir, gate 3e-3, collapsing with L — passes 2e-3 at L=17,
production L on the comparison floor).  **CARVE-OUT RETIRED:** `ResolveXCMesh` Auto is
Becke unconditionally; `AlFCCMetalIBZExact` re-pinned ON THE BECKE ROUTE (reduced
−2.1174805 == full −2.11749 to ~1e-5 — folding exact on the Becke route; uniform-route
pair was −2.1169707, route difference ~5e-4); `SiDiamondIBZ_NonSymmorphic` passes
unchanged (−7.77846 pin, Becke-vs-uniform route difference < 2e-3 for Si).
The SCF-level broken-seed negative control lands with the §8 harness (needs a
symmetry-breaking SeedStrategy).  NOTE (§4/§9): non-collinear added —
`SpinAction{None,Flip}` is documented as the collinear collapse of the general spin
rotation (SO(3)/SU(2)); representation choice (2×2 spinor vs (ρ,m)) is an open §9 question.

This doc answers three things the user asked:
1. A **gap analysis** — where are we already using symmetry, where are we leaving it on
   the table.
2. **Critical OOD/SOLID feedback** on `Mesh SpaceGroup::Symmetrize(const Mesh&)` and the
   "can we encapsulate the ops?" question.
3. How the moving parts since the original request — **spontaneous symmetry breaking**,
   the **spin-native/Shubnikov** prerequisite, the **test-material ladder** — sequence
   into one plan.

---

## 0. What already works (the baseline this builds on)

- **k-star / IBZ reduction** — `SpaceGroup::Detect` → `ReduceToIBZ(N,shift,ops)` folds the
  Monkhorst-Pack mesh to the irreducible wedge; per-k star weights ride in `D`.
- **Density star-average** — the reduced-k density is symmetrized so the reduced run
  reproduces the full mesh *exactly*: Hartree in G-space (`SymmetrizeGMap`, glide phase
  `e^{+2πi(Um)·τ}`), XC on the non-negative real-space raster (`SymmetrizeRaster`, τ-voxel
  shift). **Non-symmorphic (diamond Fd-3m) is exact** (`SiDiamondIBZ_NonSymmorphic`
  −7.77847 vs full −7.77846).

So the *density-consistency* layer of symmetry is done. **Everything below is about the
redundancy still being computed at full multiplicity**: the {G} sum, the {r} quadrature,
and — the big one — the GPW collocation stream cache.

---

## 1. Gap analysis

Every mesh/loop carries a weight list `{w_i}`; a totally-symmetric integrand lets us keep
one representative per star and fold `w_i → w_i · n_star`. Status of each surface:

| Surface | Symmetry used today | Redundancy left | Lever | Precondition |
|---|---|---|---|---|
| **k-mesh {k}** | ✅ full IBZ fold + star-avg | — | done | — |
| **density ρ̃(G) / ρ(r)** | ✅ star-averaged (G phase + raster shift) | — | done | commensurate raster (met by FFT-shift) |
| **{G} sum** (structure factor, `⟨G\|V\|G'⟩`, Hartree G-loop) | ❌ full {G} summed | fold by \|G\|-shell × point group | ≈ \|point group\| on the G-loops (small absolute cost — T1 is the *prototype*, not a payoff) | integrand totally symmetric |
| **{r} quadrature** (uniform + Becke XC mesh) | ❌ full mesh integrated | **Becke mesh** folds by site-group order; the UNIFORM raster is FFT-bound (the full raster must exist for the transform + the stream gathers index it), so folding buys only the pointwise functional eval there — minor | ≈ site-group order per atom (Becke); small (uniform) | ρ symmetric on the mesh + raster commensurability |
| **GPW stream cache** (`(pair, offset R)` → raster pts) | ⚠️ only Hermitian `j≥i` ×2 | site-symmetry of `(pair, R)`; ~513M(64b)+ pts | **5–20×** (up to ×48 in O_h) — the headline.  NB the cache is built ONCE and REPLAYED per iteration, so route (a) is a MEMORY+BUILD lever only; only route (b) also cuts per-iteration work (§6) | commensurate raster on **every ladder level** |
| **H block-diagonalization (crystal SALCs)** | ❌ dense H per k | irrep blocks per k | matrix-size³ per block | Tier B (SALC), separate track |
| **spin / magnetic group** | ❌ unpol only (`PW_XC*`) | Shubnikov ops (spatial ⋉ spin-flip) | AFM/ordering cells | spin-native XC interface first |

**Reading of the gap:** the k-side is done; the *quadrature and collocation* side is
almost entirely unexploited. The stream cache is where the memory and the per-iteration
bandwidth actually are (the `(pair,R)` lists are the ~513M cached points — built ONCE per
ladder shape, then REPLAYED as the scatter/gather of every SCF iteration), so it is the
highest-value target — but it is also the one with the hardest precondition (pointwise
index permutation needs a τ-commensurate raster on every ladder level, §5).

The user's proposed algorithm — *group by |G| (|r|), collapse symmetry-equivalents,
`w → w·n_star`* — is exactly right, with two refinements:
- **|G|-bucketing is a prefilter, not the test.** Equal magnitude is necessary but not
  sufficient (accidental degeneracies put distinct stars in one shell). The *definitive*
  test is the group action `does op map i→j`. Bucket by |G| to make the search O(N·shell)
  instead of O(N²); then apply ops within a bucket. This is precisely what `ReduceToIBZ`'s
  orbit-closure already does — generalize it (§2), don't reinvent it.
- **The equality is "just faster" only at the symmetric fixed point.** Collapsing the mesh
  replaces `f(G_member)` by `f(G_rep)`; that is exact only if `f` truly has the group. For
  a converged ground-state ρ that respects the space group, yes — reduced == full to
  ~grid/SCF precision (the reordering-level ~1e-13 holds for the *weight bookkeeping*; the
  *physical* equality is only as good as the density's actual symmetry, which needs a
  commensurate grid). The moment the density does **not** have the group, mesh reduction
  stops being "just faster" and starts **imposing** symmetry — which is the whole
  spontaneous-symmetry-breaking story (§3). This is the single most important caveat in
  this plan.

---

## 2. Where the symmetrization logic belongs (the OOD/SOLID answer)

### 2a. Decision — `FoldPoints` primitive in qcSymmetry + a client-side `Symmetrize(Mesh)`

The folding is group theory and belongs in qcSymmetry, but the fold's *product* is a
**partition**, not a rebuilt container — so the primitive stays low-level and pure, and the
Mesh-typed convenience is a **client-side helper**. Two decisions, resolved with the user
(the container decision REVISED in the 2026-08-01 review):

- **Unify on the FOLD, not on a universal container** *(replaces the earlier "everything
  becomes a Mesh" direction).*  Each point set keeps its NATIVE container — {G} stays
  `vector<ivec3_t>` (the hot-loop currency of the evaluators; Mesh-ifying it would put
  double storage + `lround` into hot G-loops for a type-clarity nicety the integer fold
  path already delivers), the {r}/Becke quadrature stays `qcMesh::Mesh`, the k-mesh keeps
  the existing IBZ machinery.  What is SHARED across {k}/{G}/{r}/streams is the fold
  *algorithm* and its `Fold` partition.
- **Phases are EDGE data of the fold, not point data — so no ivec3_t Mesh and no phase
  channel is needed anywhere.**  The glide phase `e^{+2πi(Um)·τ}` depends on WHICH op maps
  rep→member: it lives on the `(member, opIndex)` edge `Fold::members` already carries.
  Consumers use it two ways: (i) *expansion* (the `SymmetrizeGMap` scatter) computes it on
  the fly from the op; (ii) *reduced sums* collapse it — for an integrand with a known
  transformation law, `Σ_members phase` folds into ONE complex star-weight `W_rep` per
  representative (symmorphic: the plain integer `n_star`), and the natural home for
  `(m, W_rep)` is the EXISTING G-space currency (`ΔG_Map` / an `(ivec3_t,dcmplx)` pair
  list) — not Mesh.  `qcMesh::Mesh` keeps its single invariant — real points, real
  quadrature weights, `Integral(f)=Σw·f` — and its fold is Mesh→Mesh (representatives,
  `w·n_star`), which is exactly the `Symmetrize(Mesh)` wrapper.  (A `Mesh<dcmplx>`/cMesh
  stays parked for the `LatticeSum` `(R, e^{ik·R})` use case, where a complex multiplier
  genuinely IS the weight semantics; the symmetry fold does not force it.)
- **The fold primitive speaks primitive vectors and returns a partition** — kept in
  qcSymmetry, no `qcMesh` edge. A **client-side helper** `Mesh Symmetrize(const Mesh&)`
  composes `FoldPoints` + Mesh-rebuild where Mesh already lives (qcMesh/qcBasisSet). This
  keeps qcSymmetry the pure group-theory leaf (mirrors `DetectPointOps` sitting *above* it)
  while giving clients the one-call ergonomics.

```cpp
// qcSymmetry (pure qcMath/qcCommon leaf).
namespace qchem::Symmetry::Lattice_3D {

struct Fold {                    // pure combinatorics; no typed container, no weights
    std::vector<int> owner;      // raw index  -> representative slot
    std::vector<int> repRaw;     // rep slot   -> raw index of the representative
    std::vector<int> starSize;   // rep slot   -> multiplicity
    // rep slot -> [(memberRawIndex, opIndex)]: WHICH op maps rep->member.
    // opIndex is the hook the consumer needs to apply the glide phase (G) or monomial sign (stream).
    std::vector<std::vector<std::pair<int,int>>> members;
};

// The OP TYPE is the full {W|τ} (review fix: bare Matrix3D cannot fold a NON-SYMMORPHIC raster --
// τ ACTS on {r} points, r -> W·r + τ; for {G}/{k} folds τ enters only as the edge PHASE).  The spin
// action σ (§4) is the third member when the Shubnikov axis lands.
struct SymOp { Matrix3D<double> W; rvec3_t tau; /* later: SpinAction sigma; */ };

// (1) The free CORE — takes an EXPLICIT op set. Kept because tests fold under hand-made ops,
//     and the imposed-SUBGROUP policy (§3) folds under an arbitrary op SUBSET, not "all ops".
Fold FoldPoints(const std::vector<rvec3_t>& pts,
                const std::vector<SymOp>& ops, double tol);              // tolerance path ({r})
Fold FoldGrid  (ivec3_t N, rvec3_t shift,
                const std::vector<SymOp>& ops);                          // EXACT integer path ({G},{k})

// (2) SpaceGroup MEMBER wrappers — the ergonomic, ENCAPSULATED surface: the caller passes NO ops,
//     the SpaceGroup forwards its own IMPOSED op set (the natural home for the §3 policy).
class SpaceGroup { /* ... */
    Fold FoldKMesh   (ivec3_t N, rvec3_t shift) const;               // reciprocal ops, TR-on  (the k-star)
    Fold FoldGVectors(const std::vector<ivec3_t>& m) const;          // reciprocal, EXACT integer path
    Fold FoldPoints  (const std::vector<rvec3_t>& pts, double tol) const;  // direct, tolerance path ({r})
};
}

// Client-side (where Mesh lives) — the ergonomic wrapper the user wants:
//   Mesh Symmetrize(const Mesh& in, const SpaceGroup& sg, ...);
//   = build Fold from in.Points(); rebuild a compact Mesh (rep points, folded weights).
```

The member wrappers just forward to the free core with the SpaceGroup's own (imposed) ops,
so the ops accessor stops being needed at fold call sites — the §2c encapsulation, realized.

- **qcSymmetry owns the group's action on GEOMETRIC objects** — orbit folding, star weights,
  subgroup selection, op generation. Speaks `rvec3_t`/`ivec3_t`. No Mesh, no basis, no density.
- **The client rebuilds the typed container** — `Symmetrize(Mesh)` for the pure-quadrature
  case; the G-space/stream paths consume the `Fold` directly (they need the per-member op —
  below). The `e^{+2πi(Um)·τ}` phase (`SymmetrizeGMap`) and raster gather (`SymmetrizeRaster`)
  already live here.

One primitive (`Fold`) reused for k, G, r, and stream offsets — consistent with the split
already used for k and density.

### 2b. Refinements that survive (design pins for the fold)

1. **Preserve integer-exactness for the lattice folds — but do NOT template `Mesh`.** {G}
   and MP-k folds are *exact lattice combinatorics* (`|Um| = |m|`, zero tolerance, integer
   permutation); {r}/Becke folds are *tolerance-based geometry*. The exactness belongs to
   the **fold entry point** (`FoldGrid`/`FoldGVectors` = the integer path; `FoldPoints` = the
   tolerance path), NOT to the coordinate storage type:
   - **Exact doubles suffice — no `Mesh<T>` template.** G-indices are small integers (|m| ≤
     a few hundred) and the ops are integer matrices, so `U·m` is bit-exact in `double`
     (nowhere near 2⁵³) and `lround` on an exact-integer double is exact.  (With the §2a
     containers-stay-native resolution, {G} never enters Mesh at all — the exact path lives
     in `FoldGrid`/`FoldGVectors` over `ivec3_t` directly, which also removes the old
     "store indices, not Cartesian" trap by construction.)
2. **The primitive returns the partition, `Symmetrize(Mesh)→Mesh` is the wrapper.**
   Consumers need *which op* maps rep→member — for the glide phase (G), the ±monomial sign
   (streams), the diagnostic (`D`). So `Fold` carries `members`-with-`opIndex`;
   `Symmetrize(Mesh)` is the thin wrapper that applies the fold and drops the op-index,
   correct only for a pure quadrature fold (the XC {r} mesh). G-space and streams take the `Fold`.
3. **Phases never enter `Weights()`** *(review resolution, replacing the earlier
   "weights double as phase factors" idea)*: a weight that is sometimes quadrature,
   sometimes phase, is a runtime-semantics trap (`Integral(f)` over phase-weights is
   silently wrong — against the compile-time-over-runtime pin).  Phases are fold EDGE data
   (§2a): expanded on the fly from the op, or collapsed into the consumer's complex
   star-weights `(m, W_rep)` in the existing G-space containers.  Mesh weights stay real
   quadrature weights, full stop.
4. **A stream is more than a Mesh — and the pair action is a TRIPLE.** The offset geometry
   folds as `R→W·R`, but the stream reduction also carries (i) the **atom permutation** —
   "the ops fix the atom assignment" holds only for 1-atom-per-species cells like NaF;
   diamond Si's two equivalent atoms swap under half the group, and LiMn₂O₄ has many-atom
   orbits, so basis functions permute across atoms — and (ii) the **monomial ±sign** layer.
   The general pair action is (atom-permutation × monomial-sign × offset-map); the
   pair/monomial rep theory stays with the basis (fed the op), designed for the triple from
   the start.

### 2c. Encapsulating the ops (deferred — but "give me all your ops" is the lazy interface)

For the **reductions**, the ops are already hidden: `FoldPoints`/`FoldGrid` and the
`Symmetrize*` calls never hand out matrices. That is the CLAUDE.md "answer high-level
questions, don't expose internal data" bias and lets the imposed-group policy (§3) be
enforced in one place.

Three clients *currently* reach for the raw matrices — crystal SALCs (Tier B), the stream
Cartesian-monomial ±sign, and the order-parameter diagnostic (ops applied to `D`). The
first instinct is "expose an ops accessor for them." **But (user) that is usually the lazy
interface**: in each case there is very likely a *higher-level* question to ask the space
group instead of "give me all your ops" — e.g. "give me the representation of my shell set
under the group", "give me the ±sign map for this monomial pair", "give me the
symmetry-defect of this density". Finding those higher-level operations (so the ops
accessor can shrink or disappear) is **deferred design work**, noted here so the accessor
is treated as a placeholder, not the intended surface. The boundary that guides it:
*SpaceGroup owns the group and its action on geometric objects; the higher-level query is
how a client asks the group to act on ITS objects without importing the matrices.*

---

## 3. The imposed-group policy & spontaneous symmetry breaking

This is the through-line that makes the whole upgrade *safe*, and the user's overarching
concern. Every reduction in §1 (mesh, stream, k, density-average) **imposes** the group it
folds under — exact iff the density has that group. So:

- **The reduction group is a RUN-LEVEL POLICY**, "the group to impose" — a subgroup of the
  geometric space group, **not** hard-wired to the lattice's full group.  **Default
  RESOLVED (review): impose-on-assert.**  Reduction under the full detected group is an
  *opt-in speed contract* the caller asserts (or it comes bundled with the mandatory
  release-audit below); the out-of-the-box default for an unknown system imposes nothing.
  One policy object feeds IBZ, {G}/{r} folds, stream reduction, and the density
  star-average alike.
- **Order-parameter diagnostic — for FREE runs.** On an un-imposed (or subgroup-imposed)
  run, monitor `D`'s symmetry under the **full** group `G`: the residual `‖D − P_G D‖`
  (or per-op defect) reports any lowering the SCF found — which ops broke, the order
  parameter.  A run always tells the user the symmetry it actually found.
  **CRITICAL LIMIT (review): this diagnostic CANNOT audit an imposed run** — the reduction
  projects every iterate, so `‖D − P_G D‖ ≈ 0` *by construction*, precisely on the runs
  where imposition might be hiding a broken ground state.  First-order signals cannot save
  it either: the symmetric solution is a genuine stationary point (`[F,D]=0` there even
  when unstable) — SSB instability is SECOND-order (a negative orbital-Hessian eigenvalue).
- **The release-check is the AUDIT for imposed runs** *(promoted from optional workflow to
  the backstop)*: converge under imposed `G` (fast, symmetric), then release to a
  subgroup/trivial group with a symmetry-broken **seed** (`SeedStrategy`) and re-converge;
  an energy drop ⇒ SSB is real ⇒ reported.  The density-continuation machinery exists.
  §8's negative control exercises exactly this loop.  (An orbital-Hessian / MO-stability
  lowest-eigenvalue check is the future cheap alternative to the full release.)
- **The two ordering workflows need exactly this** (`[[battery north-star]]`): (i) an
  SSB-search run imposes a subgroup (or trivial) in an ordering-commensurate supercell with
  a broken seed + `+U` (LDA/GGA delocalization error otherwise suppresses disproportionation);
  (ii) the cluster-expansion / MC route imposes each candidate ordering's stabilizer
  subgroup and compares energies — same code, enumerated groups.
- **The trap this design closes:** imposing `G` for speed on a system that *wants* to break
  it yields the wrong (too-high-symmetry, too-high-energy) answer *silently*.  Hence:
  impose-on-assert as the default, the release-audit bundled with any imposition, and the
  diagnostic on free runs — imposition is never silent by construction.

`SpaceGroup`'s SRP surface thus owns: **group/subgroup computation + selection**, geometric
orbit generation (points, `(pair,offset)`), raster-commensurability checks (§5), and
star-averaging — while the basis/density own how ops act on their objects and the SCF driver
owns the policy + the staged workflow.

---

## 4. Spin-native / Shubnikov interface shape (prerequisite, cheap)

**Do this before the reduction machinery, at the interface level only.** The magnetic
(Shubnikov) group composes each spatial op with a spin action (identity / spin-flip via
time reversal). If the field-symmetrization signatures are drawn against unpol ρ (a scalar
`rvec_t`), they get **re-cut** when spin arrives (ρ → {ρ↑,ρ↓}; time reversal swaps
channels). So:

- **Symmetrization signatures are spin-native from the start.** `SymmetrizeGMap` /
  `SymmetrizeRaster` / `FoldPoints`-consumers take a **spin-resolved** field (two channels,
  unpol = the ζ=0 collapse — per `feedback_spin_polarized_primary`), and each op carries an
  optional spin action `σ ∈ {none, flip}`. A pure-spatial op has `σ=none`; a Shubnikov op
  pairs a spatial op with `flip` (and the anti-unitary `ρ̃(−G)=ρ̃(G)*` conjugation on the
  G-side).
- **TWO TIERS (review fix — the gates below are NOT "no new physics"):**
  - **Tier 4a — the interface SHAPE (cheap, blocks nothing):** settle the two-channel
    signatures + the op spin-action enum on paper/in headers so the review doesn't
    encapsulate around the unpol special case.  `PW_XC`/`PW_XC_Becke` stay unpol; the
    molecular VWN5 is already spin-native (`SpinNativeDFTPlan` B1).  Gated by compilation +
    unit-level fold tests only — a doublet CANNOT run against today's unpol solid stack.
  - **Tier 4b — the polarized SOLID pipeline (real work, `SpinNativeDFTPlan` scope):**
    spin-resolved `D` through `Crystal_EC`, two-channel filling, spin-native XC evaluation
    through collocation/Becke.  THIS is what the runnable gates test:
    - **(a) Na pseudo-atom in a box, doublet** — the existing Na q1 PP, S=½, moment 1: the
      minimal end-to-end two-channel GPW run.
    - **(b) O₂ in a box, triplet** — O q6 GTH, cross-anchored against the molecular facade's
      spin-native triplet-O₂ gate (the spin sibling of `SiPseudoAtomInBoxMatchesFinite`).
  Only 4a is a prerequisite for the reduction machinery (§7 steps 2–6); 4b is a
  prerequisite only for the magnetic materials (§7 step 7) and can proceed in parallel.

  **Tier 4b STATUS (2026-08-04) — DONE, ALL THREE GATES GREEN (committed 64a17443).**  The seven
  seams are in: `Crystal_EC(nUp,nDown)` (single-count ctors = the ζ=0 collapse),
  `tPolarizedWF<dcmplx>` instantiated + the WF factory dispatching both lineages on
  `IsPolarized()`, `Polarized_CD → tPolarized_CD<T>` with the `FourierDensityBase<T>` face
  (↑+↓ sums; `tSpinDensity<T>` alongside), the spin-native Becke XC pair (`Delta_XC_Pol`
  channel-separable Dirac + `Delta_VcorrPol` coupled VWN5, sharing the pair's one
  `XC_GridEngine` via the new spin-resolved `RhoPol` — {↑,↓} cached as a PAIR under one
  density serial, spin-agnostic seed = the ρ/2 collapse), `Ham_PW_DFT(..., polarized)`
  (Delta route only; the polarized PLANE-WAVE fit throws until designed), facade
  `CalcOptions::ppValence` (Na's GTH default is semicore q9; q1 = the valence_lowq bases),
  and `GpwOptions::multiplicity` (driver convention: 1 = the EXPLICIT two-channel singlet).
  **Found + fixed along the way:** the GPW collocation-memo D-screen collided across spin
  channels — the shared integrate-back screened by the LAST collocated `D`, and a
  fully-polarized doublet's empty ↓ channel (D=0) blanked the whole `V_H` Fock block.  The
  screen is now the UNION (elementwise max magnitude) of every density the ladder has
  collocated (`CollocMemo::Dscr`; a magnitude screen may only widen — the no-cut pin).
  **COVERAGE GAP NAMED 2026-08-07 (the AFM-collapse post-mortem, §7 step 7 — read it):** all three gates
  below, and every polarized gate before MnO, pin their moment through the CHANNEL OCCUPATIONS (nUp≠nDown).
  A Fock that has lost its spin-dependence *entirely* still reproduces those moments, so all three stay
  green — which is exactly how a spin-blind density mixer survived to MnO.  MnO's AFM-II at nUp=nDn is the
  first order parameter this code has that the occupation numbers do NOT hold up.  Any future
  polarized-machinery gate wants an nUp=nDn arm, or an assert on something the occupations don't fix: the
  exchange splitting, or the ENERGY (the Mn atom does expose it — 68 mHa).
  Gates — ALL THREE GREEN: **ζ=0 collapse** (`PolarizedSingletMatchesUnpolarizedSiGamma`,
  two-channel singlet Si Γ = the unpolarized −7.11506 anchor to 0.12 mHa) ✅; **(b) O₂
  triplet in a box** vs the facade PP triplet: 26 mHa ✅ (needs the AUTO density cutoff —
  O q6 is hard; the borrowed Si Ecut=10 cost 0.7 Ha); **(a) Na doublet** ✅ — after a
  root-cause campaign whose verdict is a **SEED/BASIN pin, not a physics bug**: from the
  Uniform seed the LONE ↑ electron converges to a GENUINE self-consistent excited basin
  72 mHa above the minimum (diffuse 3s; DIIS honors it, GDM — a local descender — stays;
  box/grid/route-independent, so every health metric reads "converged").  The functional
  was proven correct everywhere: fixed-density probes (`DISABLED_NaFixedDensityTermProbe`)
  match analytic kinetic, the EXACT discrete G≠0 lattice-sum Hartree (0.383893/0.709411
  reproduced to 6 decimals at a resolved grid), and exact ζ=1 Dirac/VWN; and
  E_GPW[D*] = −0.1420 at the independent radial same-basis oracle's minimizer (oracle
  −0.1416; complete-basis −0.1922; facade −0.1371 — the residual is the facade's
  Dunlap-J/fitted-vxc tech, user's point).  **`SeedStrategy::IonicSAD` lands in-basin: 14
  iters, −0.141933, 4.8 mHa from the facade** — pinned in the gate with a did-E-move
  anchor.  PIN: one-electron (empty-minority-channel) GPW runs are uniquely basin-fragile
  — no partner electrons pull the density coreward (Na₂/O₂/F escape Uniform fine); seed
  them SAD-family.  (An EMPTY electron configuration segfaults in
  `cSCFAcceleratorDIIS::GetNProj` — the test driver now throws on a multiplicity/Nelec
  parity mismatch; minimal spin for odd Nelec = `twoS=Nelec%2`.)
- **NON-COLLINEAR (user addition, 2026-08-01): collinear is itself a collapse.**  The
  two-channel {ρ↑,ρ↓} formulation and the binary σ ∈ {none, flip} are the FIXED-AXIS
  (collinear) restriction of the general objects: a spin-space-group op carries a spin
  ROTATION (R_s ∈ SO(3) acting on the magnetization 3-vector m(r), equivalently SU(2) on
  the 2×2 spinor density ρ_αβ; time reversal = m→−m + conjugation), and spirals / canted
  AFM / general LiMn₂O₄ orderings need it.  Tier 4a's job therefore includes NOT baking
  "two scalar channels" into any signature deeper than the collinear tier — the spin
  action rides on the op (where `SpinAction` sits today, documented as the collinear
  collapse) and the field type is the swappable argument.  The REPRESENTATION choice
  (complex 2×2 spinor density vs (ρ, m) scalar+3-vector fields) is an OPEN question (§9)
  — it does not block the spatial reduction machinery.

---

## 5. Raster / stream commensurability (the hard precondition)

A **pointwise** index permutation (what streams and the {r} raster fold need) exists iff
every op maps the integer torus to itself: axes mixed by `W` need **equal N**, and each
translation component needs `N_i·τ_i ∈ ℤ` — diamond's τ=(¼,¼,¼) ⇒ **4 | N per axis, on
every ladder level.** The 5-smooth FFT padding does *not* guarantee this (45, 27, 15, 3 are
5-smooth but not ÷4; NaF's measured ladder had N=45 and N=3 levels). So stream-level
symmetry requires the raster menu to round up to **5-smooth multiples of
lcm(τ-denominators)** (45→48, etc. — modest).

- The **existing IBZ density symmetrization is exempt** — its FFT fractional shift `e^{iG·τ}`
  is exact for the band-limited interpolant at *any* τ. The constraint is specifically for
  **sparse compact streams**, which can't take an FFT shift without densifying.
- Diamond Si already forces this path in the suite (`SiDiamondIBZ_NonSymmorphic`: AutoGrid
  yields odd N, 4∤N is generic). A trigonal/hexagonal **screw** material (quartz P3₁21,
  τ=c/3; P6₁-family, τ=c/6) is the extension gate for non-power-of-2 denominators.
- **k≠Γ:** the offset lists also carry Bloch phases `e^{ik·R_n}`, and ops act on k too —
  stream reduction at general k must stay consistent with the IBZ k-star bookkeeping (one
  more thing `SpaceGroup` owns rather than the evaluator).

---

## 6. The three reduction targets, concretely

Each = (geometric fold in qcSymmetry) + (typed-data application in qcBasisSet) + (a
reduced==full test on a commensurate grid + a symmetry-broken negative control).

**(T1) {G} sum.** Fold `GetGVectors` under the point group (`FoldPoints` on the `ivec3_t`
list; |G|-bucketed). The Coulomb kernel `4π/|G|²` is constant within a star (`|Um|=|m|`), so
it folds cleanly. Apply in the structure-factor loop and the Hartree G-sum: evaluate at
representatives, weight by `n_star`. Cheapest to prototype (small, self-contained, no raster
commensurability issue since G-indices permute exactly).

**(T2) {r} quadrature mesh — in practice the BECKE mesh** (review fix: the uniform raster
is FFT-bound — the full raster must exist for the forward transform and the stream gathers
index it, so folding buys only the minor pointwise functional eval there).  The **Becke
angular grid** folds by site-group order and merges with the *site-group-adapted angular
grid* item (orbits must avoid special/bond directions — the Lebedev cube-corner lesson; the
site group fixes both how many and which directions, subsuming the GL-29 vs GL-17 default
question).  Requires ρ symmetric on the mesh + raster commensurability (§5).  **The
site-adapted (group-INVARIANT) mesh is also the PRECONDITION for pointwise Becke
star-averaging** — ops must map mesh points to mesh points, which today's GL grids don't —
i.e., T2 is what UNLOCKS the Becke+IBZ verification and retires the `reduceBZ`
uniform-raster carve-out.  Also the candidate cure for the **Becke ×
degenerate-open-shell oscillation** (a site-symmetric quadrature removes the
rotating-error channel — though note §3: a MINIMAL site-adapted grid presumes symmetric ρ;
deliberately symmetry-broken runs need the generic orientation-robust fallback grid).

**(T3) GPW stream cache — the headline 5–20×.** Fold `(pair, offset R)` under site symmetry:
`(pair, R) → (pair′, W·R)` + a raster-index permutation; for a 1-atom-per-species primitive
cell the ops fix the atom assignment, so most of the win is in the diffuse pairs'
hundreds-of-offsets lists collapsing by site-symmetry order. **Two routes** (from the
GPWPlan1 plan): (a) *orbit-expansion replay* — store irreducible streams, replay members
through per-op raster permutations (exact, no symmetry imposed, but cache-hostile gather);
(b) *collocate irreducible reps with orbit weights + star-average ρ once per iteration*
(reuses the existing `SymmetrizeGMap`/`SymmetrizeRaster`), mirror on integrate-back with the
representation transform of `h` — cuts the per-iteration scatter/gather by the same factor,
but is exact only for `D` symmetric under the group (**it imposes** — §3). Polarization
components: Cartesian monomials under cubic ops are signed axis permutations (pair→±pair,
no shell mixing) — the basis's job, fed ops.

### 6a. T2 wiring design (2026-08-01 — the site-group ↔ angular-grid question resolved)

**The user's question — is the XC wiring independent of the (uncoded) connection between an
atom site's point symmetry and its Becke angular grid? — NO: that connection IS the T2
invariance precondition.**  But it splits cleanly into a verification increment that
satisfies the precondition *generically* and a production increment that satisfies it
*minimally*:

- **The code today**: `MakePeriodicBeckeMesh` (UnitCell.C) builds ONE angular grid
  (`MakeAngular(mp)`), in ONE fixed orientation, shared by every atom — generically
  non-invariant under any nontrivial op (`UnmatchedCounts` now measures this).  The precise
  invariance condition: the collection {atom shell, grid orientation} must be op-COVARIANT —
  i.e. each representative atom's angular set invariant under its SITE group
  (`SpaceGroup::SiteStabilizer`, landed), symmetry-partner atoms carrying the op-rotated
  copies.  Ops fixing an atom then map its shell to itself; ops moving it map its shell onto
  its partner's.

- **W1 — verification wiring (generic invariance, no new grid design).**  Feed the engine
  `MakeInvariant(Becke mesh)` under the §3-imposed group: invariant with ORBIT-SYMMETRIC
  weights by construction (each member accumulates the same group-sum — and each rotated
  copy of a Becke quadrature is itself a valid Becke quadrature of the same cell, so
  exactness is preserved).  Insert the pointwise star-average (`FoldMesh` +
  `SymmetrizeValues`, geometry-fixed fold cached with the mesh) into `BeckeXC_Engine::Rho`
  — once per iteration, policy-gated.  NO quadrature fold in W1.
  **The adjoint is already exact — proof:** \f$E[D]=\sum_g w_g\,\epsilon((P\rho)_g)\f$ gives
  \f$dE/dD = \Phi^\dagger\,\mathrm{diag}(P^\dagger(w\circ v))\,\Phi\f$; on an invariant mesh
  with orbit-symmetric weights \f$P\f$ is self-adjoint in the \f$w\f$-inner product, so
  \f$P^\dagger(w\circ v)=w\circ(Pv)\f$, and \f$v=v(\rho_\mathrm{sym})\f$ is already
  symmetric ⇒ \f$Pv=v\f$: the EXISTING `Matrix()` quadrature is the exact derivative,
  untouched.  The only insertion point is the ρ side — the Becke sibling of the
  ∫ρV==Tr(Dh) machine-exactness pin.
  **Cost honesty:** `MakeInvariant` grows the mesh (≤|G|×, Φ-table memory with it), so W1
  runs the GATE on a COARSE mesh — legitimate, because the gate (Becke+IBZ == Becke+full,
  ~1e-7 through-SCF class) compares two runs on the SAME mesh, so grid quality cancels.
  The `reduceBZ` Auto→uniform carve-out (`ResolveXCMesh`) STAYS until W2.  [RETIRED with
  W2c, 2026-08-02.]

- **W2 — production wiring (the minimal site-adapted grids; the real T2 payoff).**  The
  bespoke connection: per atom ORBIT, the representative gets an angular set built as
  unions of site-group orbits of directions (generic directions, avoiding special/bond axes
  — the Lebedev cube-corner lesson) with weights from the angular moment conditions;
  partners get the op-images.  Invariance at ~today's point count (no blow-up), and only
  then does `Symmetrize` deliver the ≈|site-group| fold on the per-iteration ρ GEMM +
  functional evaluations.  The H assembly (\f$\Phi^\dagger\mathrm{diag}(wv)\Phi\f$) stays
  full-mesh until the representation transform of \f$h\f$ lands (the T3-shared machinery) —
  the honest statement of T2's lever.  W2 also subsumes the GL-29 vs GL-17 default question
  and is the candidate cure for the Becke × degenerate-open-shell oscillation (§6 T2).
  This is the AngularMathPlan revival and needs its own grid-design increment.
  **ROTATION INSIGHT (user question, 2026-08-02) — two separable facts:**
  (i) *Free runs:* rotating an efficient Lebedev grid OFF the bond axes is a legitimate,
  nearly-free accuracy fix — quadrature exactness is rotation-invariant; the measured 5–10×
  ρ-weighted loss was pure ALIGNMENT (the ⟨111⟩ orbit on diamond's bonds).  A one-knob
  experiment on the existing `GPW_BECKE_ANG=lebedev` path.
  **[EXPERIMENT RUN 2026-08-02 — CONFIRMED, default unchanged.]**  Knob: `MeshParams.angRot`
  (rigid Rodrigues rotation about the fixed generic axis (1,2,3)/√14, applied in `MakeAngular`;
  `GPW_BECKE_ROT` env; 0 = bit-identical).  Probe (`DISABLED_RotatedLebedevXCProbe_SiGamma`,
  converged Si Γ density, all vs the refined GL reference nR=80/L=29): unrotated Lebedev-50 has
  points EXACTLY on the bonds (0.000°); 0.4 rad moves the closest to 8.6° and recovers **7.6×
  on E_xc** (−7.0e-3 → +9.2e-4), **4.9× on ∫ρ** (+1.9e-2 → −3.8e-3), 2.4× on V_xc elements —
  the loss really was pure alignment.  GL-29 is rotation-insensitive (control: 1.050e-4 →
  1.075e-4).  BUT rotated Lebedev-50 (degree 11 — the audited tables stop there) still trails
  GL-29 ~9× on E_xc at this recipe, so the FREE-RUN DEFAULT STAYS GL-29; a Lebedev default
  would need verified higher-degree tables (then rotated by default).  Increment (b) CLOSED.
  **[TABLES LANDED same day — increment (f):]** canonical Lebedev–Laikov GEN_OH parameters
  imported for orders 86/110/146/170/194/302/350/434 (degrees 15–35; orders 74/230/266
  EXCLUDED by the generator audit — negative weights), each gated by the full monomial sweep
  to its claimed degree + positivity + Σw=4π (`CanonicalLebedevOrdersHaveClaimedDegree`,
  moment error < 1e-10).  A/B on the same converged Si density: **Leb-302 ≈ GL-29 on every
  metric at 302/450 = 67% of the directions** (dExc +7.3e-5 vs +1.05e-4, ρ-lost −5.4e-4 vs
  −8.9e-4, Vxc 1.77e-5 vs 1.75e-5), and at degree 29 even the UNROTATED bond-aligned grid is
  fine on Si — the alignment poison is a low-degree phenomenon (Leb-50rot stays 9× worse).
  **DEFAULT FLIP DEFERRED:** `nAngular` means COUNT for Lebedev but DEGREE for GL, and the
  imposed-run site-adapted builder consumes it as the degree — a naive `{Lebedev,302}`
  default would ask the NNLS builder for degree 302.  The flip lands with the
  `MeshParams` angular-ergonomics pass (CleanupCandidates: degree-typed interface), then
  Becke free runs get ~33% angular savings for free.
  (ii) *Imposed runs:* rotation destroys exactly what T2 needs — a rotated grid is invariant
  under the CONJUGATED group R·G_s·R⁻¹, not G_s, and the site group's special axes ARE the
  bond axes by definition (site symmetry is generated by the neighbours); for a non-axial
  G_s (e.g. T_d) the normalizer is finite, so no continuous rotation preserves invariance.
  Hence the generic-orbit premium.
  **THE MIXED RULE (the sharpened refinement):** a G_s-invariant grid may admit special
  orbits — and only SOME families lie on bonds (diamond T_d: ⟨111⟩ = the poison; the
  ⟨100⟩ 6-pt and ⟨110⟩ 12-pt orbits point BETWEEN bonds and are safe).  Extend
  `MakeInvariantAngularMesh`: seed the pool with the group's special-orbit representatives,
  filter against the ATOM-SPECIFIC bad-direction list (actual bond vectors from the
  structure, not the abstract group), same moment-condition solve — recovering much of
  Lebedev's small-orbit efficiency while keeping invariance-by-construction.  For AXIAL
  site groups (C_n/C_nv, the common low-symmetry sites) the normalizer is continuous:
  genuine rotation freedom about the axis to steer residual special directions away from
  neighbours.  [BUILT as W2c 2026-08-02 (see the status block): special-orbit seeding +
  bond filter + NNLS all landed; the axial-group CONTINUOUS steering (rotating the
  mirror-plane seeds about a C_n axis) remains open — today's axial sites just get the
  projected-seed defaults.]

- **OWNERSHIP RESOLVED (user + review, 2026-08-01) — the four-way Ham/Structure/Mesh/SpaceGroup
  standoff dissolves on two decisions:**
  1. **qcStructure is where Mesh and Symmetry meet without depending on each other.**  New
     library edge qcStructure→qcSymmetry; `Lattice_3D` OWNS its lazily-detected `SpaceGroup`
     (`GetSpaceGroup()` — the UnitCell→{A, sites} adapter moved out of the GPW factory,
     common to every lattice basis flavour PW/GPW/APW/LAPW), and the Mesh-typed fold helpers
     live in `qchem.SymmetrizeMesh` (src/Structure/Lattice_3D/).  qcMesh and qcSymmetry both
     stay pure leaves; `SpaceGroup` gained `ReciprocalMatrix`/`ToFractional`/`ToCartesian`
     so downstream calls need no separate cell plumbing.  WHAT a run imposes stays the §3
     run-level policy (today: `GPWParams.imposeSymmetry` — renamed from `reduceBZ` 2026-08-02,
     user request, since the switch imposes the whole group — resolved once in the factory).
  2. **The Becke-route tension was a design flaw, not a plumbing problem (user diagnosis):**
     `PW_XC_Becke`/`BeckeXC_Engine` conflate the term with a raw mesh, where every other
     fitted term holds a FIT BASIS that owns its integration mesh.  W1 therefore lands as a
     `BeckeFit_IBS` (cFIT_SF_ABS sibling of `PlaneWaveFit_IBS`): it owns the (invariant)
     Becke mesh + the fold, overrides the EXISTING `SymmetrizeRaster` virtual with the
     orbit-mean (`SymmetrizeValues`), absorbs the engine's Φ-table/GEMM machinery as
     internals, and gets its ops ctor-injected by the orbital basis exactly like
     `PlaneWaveFit_IBS` does.  ρ then flows through the SAME
     `FourierDensity::GetRhoOnGrid(fit)` + composite `SymmetrizeRaster` path as the uniform
     route — the Ham never sees a mesh, ops, a UnitCell, or a cast.

### 6b. T3.0 design memo (2026-08-02) — route (b) written out; the §5 verdict

**Scope:** the §7 step-6 T3.0 increment.  Route (b) — collocate irreducible `(pair, R)`
streams with orbit weights, star-average ρ once, mirror on integrate-back with the
representation transform of \f$h\f$ — written out to the adjoint level, and the open
question answered: **route (b) needs NO §5 commensurable raster menu.**  §5 remains a
precondition only for route (a) and for the (optional) voxel-fold retirement of
`FIT_SF_ABS::SymmetrizeRaster` on non-symmorphic crystals.

**1. Objects.**  A stream triple \f$t=(i\le j,\,n)\f$ carries \f$S_t\f$ = (wrapped raster
indices, analytic pair values) on the pair's ladder level; today's replay scatters
\f$c_t S_t\f$ with \f$c_t=\mathrm{fold}\cdot\mathrm{Re}[D_{ij}\overline{\varphi(n)}]\f$ and
gathers \f$b_t=S_t\cdot V\f$.  Triples are already Hermitian-canonicalized:
\f$(i,j,n)\sim(j,i,-n)\f$ (the same continuous product field, lattice-translated — on the
raster an EXACT integer voxel shift by \f$n\circ N\f$).

**2. The op action on triples.**  For an imposed \f$g=\{W|\tau\}\f$ (direct, fractional;
Cartesian \f$R_g=A\,W\,A^{-1}\f$), the basis-function map is \f$i\to(i',s_i,L_i)\f$: the
atom match \f$W f_i+\tau=f_{i'}+L_i\f$ (integer \f$L_i\f$ — the SpaceGroup detection
guarantees a partner) plus the monomial map.  When \f$R_g\f$ is a SIGNED AXIS PERMUTATION
(every op of a cubic-lattice crystal — diamond, rocksalt, FCC Al) a Cartesian monomial maps
to \f$\pm\f$ one monomial: \f$s_i=\pm1\f$, no shell mixing.  Non-cubic lattices get a per-op
runtime check; an op failing it is dropped from the stream fold (folding merely reduced,
never wrong — general Wigner shell-mixing is a later increment, NOT drawn into these
interfaces).  The triple action is then
\f[ g\cdot(i,j,n)=(i',\,j',\,W n + L_j - L_i),\qquad \sigma_t=s_i s_j, \f]
composed with Hermitian canonicalization when \f$i'>j'\f$ (edge datum = op index + flip
flag).  KEY FACTORIZATION: the pair map is \f$n\f$-INDEPENDENT, so the fold splits into
pair orbits \f$\times\f$ within-pair offset orbits under the pair's stabilizer — for a
1-atom-per-species cell the pair orbits are small (ops fix the atom assignment, permuting
only shell components) and the win is the diffuse pairs' hundreds-of-offsets lists
collapsing, exactly as §6 T3 scoped.

**3. Reduced collocation.**  With \f$D\f$ symmetric (\f$D_{i'j'}=\sigma_t D_{ij}\f$ — what
imposition asserts, §3) the signs cancel between weight and stream:
\f$c_{gt}S_{gt}=c_t\,O_g S_t\f$ as continuous fields, where
\f$O_g=\f$ (FFT \f$\tau\f$-shift) \f$\circ\f$ (voxel permutation \f$W\f$).  Hence
\f[ \rho_{\rm full} \;=\; P\Big[\sum_{\rm reps}\tfrac{|G|}{|\mathrm{stab}_t|}\,c_t\,S_t\Big],
\qquad P=\tfrac1{|G|}\sum_g O_g . \f]
An orbit with an odd stabilizer element (\f$g t=t,\ \sigma=-1\f$) has \f$D_{\rm rep}=0\f$
under the imposed symmetry and is annihilated by \f$P\f$ anyway — skip it at build (free
extra reduction, never a correctness burden).  On the raster the identity is exact up to
the BAND-LIMIT class (a member's samples vs the band-limited image of the rep's samples
differ by the beyond-\f$E_{cut}\f$ aliasing tail of the compact product — the same class as
the collocation error itself): the §8 through-SCF ~1e-7 tier.  On a commensurate grid
(τ=0, or \f$N\tau\in\mathbb Z\f$) \f$O_g\f$ is an exact index map and reduced==full lands
in the 1e-13 reordering tier — the unit-gate arrangement.

**4. Where \f$P\f$ is applied — NOT per level, and at NO new site.**  \f$W\f$ is an integer
matrix in fractional coordinates, so the voxel map \f$W\cdot g\f$ is exact at ANY \f$N\f$
(axis-mixing ops need equal per-axis \f$N\f$ — every `AutoGrid` level is cubic; keep the
`SymmetrizeRaster` guard per level); \f$\tau\f$ is exact via the FFT shift theorem at any
\f$N\f$ (§5-exempt — the IBZ exemption, per level).  \f$O_g\f$ commutes with the spectral
level transfers (G-shells map to same-\f$|G|\f$ shells), so the projector applies ONCE at
the sites that ALREADY run per iteration on imposed runs: `SymmetrizeGMap` on the combined
\f$\tilde\rho\f$ (`CompositeCD::GetRepulsion3C/Overlap3C`) and `SymmetrizeRaster` on the
raw XC raster (`CompositeCD::GetRhoOnGrid`).  The reduced scatter only restricts to reps
and re-weights by \f$|G|/|\mathrm{stab}|\f$.  Total charge is preserved even before
projection (\f$\int O_g f=\int f\f$), so the Tr(DS) bookkeeping is untouched.

**5. The adjoint — the §9 proof obligation, discharged.**  \f$E=E[P\rho_{\rm red}(D)]\f$
\f$\Rightarrow\f$ \f$dE/dD=(\text{reduced gather})\circ P^\dagger\f$ applied to \f$v\f$,
and \f$P^\dagger=P\f$ (each \f$O_g\f$ is orthogonal in the uniform raster inner product:
permutation \f$\times\f$ unit-modulus G-space phase; \f$O_g^\dagger=O_{g^{-1}}\f$ is in the
set).  So the EXACT derivative is: **symmetrize \f$v\f$ once, dense** (G-side:
`SymmetrizeGMap` on \f$\tilde V\f$ BEFORE the per-level band restriction; raster side:
`SymmetrizeRaster` on the raw \f$v_{xc}\f$) — the only NEW per-iteration dense step,
FFT-scale — **then gather rep streams only, then fill partners by the rep transform**
\f[ h_{i'j'} = \sigma_t\,h_{ij}\quad(\Gamma), \f]
with within-pair offset orbits carrying the signed multiplicity
\f$\sum_{\rm members}\sigma\f$ (\f$b(gn)=\sigma_g b(n)\f$ against a symmetric \f$V\f$).
With that recipe \f$\mathrm{Tr}(Dh)==\int\rho V\f$ to machine precision on the shared
REDUCED operator — the T3 sibling of the W1 adjoint proof, and the §9 "prove it in the
design" item.  No sparse stream is ever index-permuted, in either direction.

**6. D-aware kill weights (the §9 active-set item).**  The reduced replay decides each
kill ONCE per orbit, on \f$|c_{\rm rep}|\cdot\tfrac{|G|}{|\mathrm{stab}|}\cdot
\mathrm{maxv}_{\rm rep}\f$, and BOTH directions replay the same rep streams — the shared
active set is orbit-consistent BY CONSTRUCTION (the W2c ε-tail orphan mechanism cannot
arise: partner streams are never built, partner \f$h\f$ never gathered).  Freak-mode note
for the unit gate: reduced-vs-full comparisons can in principle differ by an
ε-BOUNDARY screen flip (partner arithmetic is ULP-close, not bit-equal, in `ForPairBox`'s
\f$|val|<\varepsilon\f$ cut); flips need \f$|val|\f$ within ULPs of \f$\varepsilon\f$ —
essentially never, but if the 1e-13 gate ever flakes, that is the mechanism (fix =
orbit-consistent screen decisions recorded in the fold, the W2c pattern).
**[MEASURED at T3.1 — the mechanism that ACTUALLY bites is the fp32 STREAM TIER, not screen
flips:** orbit-partner streams are rounded INDEPENDENTLY, so fp32-tier values (~6e-8
relative) break partner congruence at ~1e-8 — the a=10.26 diamond gate cell (474M-pt
demand, 172/528 pairs fp32) showed a 5.7e-9 symmetry defect in the FULL path itself,
ε-stable, before any fold code ran.  The reordering-tier gate therefore also requires
full fp64 tier coverage (gate cells sized to fit the 150M-pt budget).  Production
corollary, T3.2 bonus: the reduced build costs ÷~|orbit|, so MORE pairs fit the fp64
tier — imposed runs get an ACCURACY win on top of the memory/build win.]**

**7. k≠Γ boundary (T3.4).**  Per k-block, fold only under the LITTLE GROUP of \f$k\f$
(\f$gk\equiv k\f$ mod reciprocal lattice — at Γ the full group, where the entire win
lives); the cross-k projection already rides the IBZ star bookkeeping (the same
`SymmetrizeGMap` sites, summed over IBZ reps).  Member weights at general \f$k\f$ pick up
\f$e^{ik\cdot(L_j-L_i)}\f$-class phases — written out at T3.4, owned by SpaceGroup per §5.

**8. THE §5 VERDICT.**  Route (b) needs no commensurable raster menu: linear parts are
exact voxel maps at any cubic \f$N\f$, glides ride the FFT shift, streams are never
permuted.  §5 remains the precondition ONLY for (i) route (a)'s orbit-expansion replay
(sparse member streams BY index permutation — the no-impose fallback), and (ii) the TODO
cleanup replacing `FIT_SF_ABS::SymmetrizeRaster`'s FFT-shift by a precomputed voxel fold
(a τ-acting direct-grid `FoldGrid` needs \f$N\tau\in\mathbb Z\f$; symmorphic crystals can
take it today).  RECOMMENDATION: land T3.1+T3.2 without touching the raster menu; T3.3
(menu + route (a)) stays parked until a free-run memory-pressure case demands it; the
`SymmetrizeRaster` retirement follows the same schedule (the FFT-shift form is exact and
cheap — retiring it is ergonomics, not correctness).

**9. Ownership.**  The op action on basis functions (\f$i',s_i,L_i\f$ per op) is the
BASIS's job (the PG evaluator layer owns centers + monomials), fed `DirectOp`s + the cell;
the `(pair,R)` fold is stream bookkeeping and lives beside the stream cache (cached like
`itsStreamCaches`, keyed on the ops set; geometry-fixed across iterations and k-blocks).
Policy stays §3: ops arrive non-empty only under `imposeSymmetry`, via
`GPW_Evaluator::SetSymmetryOps` → threaded into `CollocateDensity`/`IntegratePotential`.
The BUILD also shrinks to reps (the memory/build win = budget ÷ ~orbit factor); the
two-pass parallel build's tier walk is unchanged, only its pair list shrinks; the
`[stream cache]`/`grids.stream` readout grows a reps + orbit-factor line.

**10. Unit gates (T3.1, the §8 tiers).**  (a) SYMMORPHIC 1e-13: Al FCC one-shot
collocate/integrate on a frozen symmetric \f$D\f$ — every τ=0, \f$O_g\f$ a pure voxel
permutation, reduced==full in the reordering tier.  (b) NON-SYMMORPHIC 1e-13: Si diamond
on a test-forced \f$4|N\f$ ladder (explicit Ecut), same one-shot; PLUS the production
(odd-\f$N\f$) grid at grid-class tolerance to MEASURE the band-limit gap item 3 predicts.
(c) NEGATIVE CONTROL: a symmetry-broken \f$D\f$ ⇒ reduced ≠ full and the §3 diagnostic
fires (never silent).

---

## 7. Sequencing

Ordered so each step is gated by a cheap, already-shaped test, and no interface is drawn
against a special case it will outgrow:

1. **✓ DONE — Spin-native interface SHAPE** (§4 tier 4a) — signatures + op spin-action enum only;
   gated by compilation + unit-level fold tests.  *Prevents a redo; blocks nothing.*
   (§4 tier 4b — the polarized solid pipeline with gates (a) Na doublet / (b) O₂ triplet —
   **✓ DONE 2026-08-04 (64a17443), all three gates green — see the §4 STATUS block.**
   Step 7's remaining seed-side prerequisite is the spin-polarized SAD with per-site
   moments, doc/SCFSeedingPlan.md §10.)
2. **✓ DONE — `Fold` primitive in qcSymmetry** (§2b) — generalize `ReduceToIBZ` → `FoldPoints`/
   `FoldGrid` with per-member op index; re-express the k-fold on it (bit-identical). Unit
   tests only.
3. **✓ DONE — Imposed-group policy object + order-parameter diagnostic** (§3) — run-level policy;
   diagnostic reports symmetry lowering. Gated by: symmorphic Al IBZ still exact; the unit-level
   broken density fires the diagnostic (the SCF-level broken-SEED control awaits a
   symmetry-breaking SeedStrategy, §8 harness).
4. **✓ DONE — (T1) {G} fold** — smallest real reduction; reduced==full on a symmetric density.
5. **✓ DONE — (T2) {r} + Becke angular fold**: the site-adapted invariant mesh + Becke+IBZ
   verification landed (§6a W1/W2a/W2b); the production-L growth re-measure (1.97× at L=29),
   the mixed special-orbit rule + NNLS (W2c), and the `reduceBZ` carve-out RETIREMENT
   (Auto+reduceBZ = Becke, Al re-pinned route-matched) closed it out 2026-08-02.
6. **(T3) stream cache reduction** — the 5–20× lever; needs §5's commensurable raster menu
   on every ladder level. Route (b) default (imposes, gated by the policy + diagnostic),
   route (a) as the no-impose fallback.
   **T3 CAMPAIGN SCOPING (2026-08-02, code recon done — anchors for the next session):**
   - *The cache* (`PG_Cart_MnD/Evaluator.C` ~L550-844): per `(pair i≤j, offset n)` a
     `PairOffsetStream{n, idx[], val[]/val32[], maxv}` — WRAPPED LINEAR RASTER INDICES + the
     analytic pair values, pure geometry, shared across iterations AND k-blocks (Bloch phase
     is a `cellphase_t` callback applied only at contraction: density `Re[D_ij·conj(φ(n))]`,
     integrate-back `+φ(n)`).  Budgets 150M fp64 / 850M fp32 pts; ≤4 caches keyed on ladder
     shape.  The torus wrap is `((g%N)+N)%N` per axis (`ForPairBox` L544) — the map any §5
     index permutation must commute with.  TODAY'S ONLY SYMMETRY: the Hermitian j≥i ×2 fold;
     `SetSymmetryOps` exists on the evaluator but feeds only T1 + the density raster fold.
   - *The raster menu hook* (`PW/Imp/Evaluator.C FFTGrid()` L64): per-axis `Next5Smooth` of
     the cubic `AutoGrid` — 5-smooth but NOT lcm(τ-denominator)-commensurate (§5's gap);
     ladder levels get their own `N_L` in `GPW BuildLevels`.
   - **T3.0 design memo DONE 2026-08-02 → §6b:** route (b) needs NO §5 menu (W integer =
     exact voxel map at any cubic N; τ rides the FFT shift; streams never permuted);
     adjoint discharged (symmetrize v once dense, gather reps, h-transform partners ⇒
     Tr(Dh)==∫ρV machine-exact on the reduced operator); D-aware kill orbit-consistent by
     construction; §5 demoted to route-(a)-only (+ the SymmetrizeRaster voxel-fold
     cleanup).  T3.3 parked pending a free-run memory-pressure case.
   - *The original open question, for the record:* route (b) may not
     need §5 commensurability at all.  Reduced scatter of rep streams gives a partial ρ
     whose group-average can ride the EXISTING dense-raster FFT τ-shift symmetrization
     (§5-exempt, per level); the integrate-back gather of rep pairs against a SYMMETRIC V
     needs only the h-side rep transform (matrix-element bookkeeping + Cartesian signed
     axis permutations), not a raster permutation.  Only route (a) (orbit-expansion replay
     of sparse streams) demonstrably needs the §5 menu.  If the memo confirms this, §5
     becomes a route-(a)-only increment and route (b) can land first without raster churn.
   - *Increments:* **T3.0** the design memo (route (b) adjoint written out: D-aware kill
     weights orbit-consistent, rep transform of h incl. polarization signs, per-level ρ
     star-average; decide the §5 question).  **T3.1 DONE 2026-08-02** the `(pair, R)` orbit
     fold primitive: `NR_Evaluator::SetStreamSymmetryOps(DirectOps, A)` (LatticeSum1E face,
     forwarded by PG_Cart::Orbital_IBS) builds per-op basis maps (signed-axis-perm check on
     R=AWA⁻¹, atom/monomial/radial matching → (i′, s_i, L_i)), the pair fold (rep/image
     edges σ+flip, odd-self-edge ⇒ dead, diagonal Hermitian-twin involution) and per-rep-pair
     offset-orbit multiplicities; CollocateDensity replays reps × (pairMult·within) with the
     member-rule D-kill, IntegratePotential gathers reps × within and fills images by
     h_{i′j′}=σ·h_ij.  GATES (GPW_UT `StreamFoldReducedMatchesFull_*`): FCC Si symmorphic
     ρ 2.3e-14 / h 8.9e-16; diamond Si non-symmorphic (4|N=16 grid) ρ 5.0e-14 / h 4.6e-15;
     adjoint at D-kill class; generic-detune negative control fires.  DISCOVERY: the fp32
     stream tier breaks partner congruence at ~1e-8 (see the §6b item-6 measured block) —
     gate cells sized for full fp64 coverage.  **T3.2 DONE 2026-08-02** route (b) reduced
     build+replay behind `imposeSymmetry`: `EnsureStreams` builds ONLY rep pairs/offsets under
     an armed fold (`SetStreamSymmetryOps` is idempotent per op set and CLEARS the stream
     caches on any change — a reduced cache must never serve an unfolded replay); the GPW
     factory arms the fold on IMPOSED **Γ-ONLY** runs, **OPT-IN via `GPW_STREAM_FOLD=1`**
     (read fresh so one process can A/B; multi-k IBZ runs keep full streams until T3.4's
     little-group + per-block arming).  **DEFAULT-ON RETRACTED 2026-08-03 — the open-shell
     finding:** the fold imposes STRICTLY MORE than the historical `imposeSymmetry` (the ρ
     star-average tolerates a symmetry-broken iterate D by projecting it pointwise; the
     reduced replay reads only orbit-rep D elements, asserting D itself symmetric).  A
     DEGENERATE OPEN SHELL breaks that permanently: the imposed Si pseudo-atom-in-a-box p²
     gate flips from the benign rotating-ρ mode into charge-transfer sloshing (~0.26 Ha off)
     with the fold armed — and only MARGINALLY (it passed one full-suite run before failing
     deterministically the next day: a bistable oscillator, not a tolerance issue).
     Default-on returns with an auto-arm criterion — gapped/closed-shell detection or Fermi
     smearing (the same cure as the Becke × open-shell channel) — the T3.4-adjacent item.
     MEASURED on the diamond unit-gate cell:
     528 pairs → 40 reps, 15000 offsets → 164, 72.7M → 1.02M cached pts (**71× build/memory**;
     the plan's 5–20× was conservative for high-symmetry cells).  Through-SCF gate
     (`StreamFoldImposedGamma_SiDiamond`): imposed Γ diamond, fold-off vs fold-on totals agree
     to **1.0e-6 Ha** (band-limit class on the non-commensurate production ladder N=15/8/24 —
     the §6b no-§5 verdict exercised in production), both arms converge, charge exact.
     DISCOVERY en route: the framework integral cache is content-keyed, so the fold state MUST
     enter `GPW_Evaluator::IDFragment` (`|sfold=n`) — without it the A/B's second run silently
     replayed the first run's unfolded cached tensor closures (bit-identical totals = the
     tell).  V-side note: v is NOT explicitly re-symmetrized before the reduced gather —
     V_H is exactly symmetric (built from the projected ρ̃) and Becke-route H_xc bypasses
     streams entirely (mesh GEMM), so the only unpaired piece is the uniform-route v_xc's
     glide-part band-limit defect, the same class the unfolded imposed route already carries;
     revisit only if a GDM/virial pin ever resolves it.  Anchor re-pin: NOT needed — no
     pre-existing suite test is imposed+Γ-only; the A/B gate is the pin.  **T3.3** §5 commensurable menu +
     route (a) expansion replay as the no-impose fallback (only if T3.0 says route (b)
     doesn't subsume the need).  **T3.4 (op-action machinery) DONE 2026-08-03**: k≠Γ op
     action consistent with the IBZ star bookkeeping — `SpaceGroup::LittleGroup(kFrac)` /
     `LittleGroupDirectOps` (SpaceGroup-owned, \f$(W^{-1})^\top k\equiv k\f$ mod 1);
     `SetStreamSymmetryOps(ops, A, kFrac)` folds under the little group only (guard drops
     movers), with the DERIVED simplification that replay multiplicities stay plain
     integers (\f$e^{2\pi ik\cdot Wn}\equiv e^{2\pi ik\cdot n}\f$ for little-group ops) —
     \f$k\f$ enters ONLY the edge factor \f$\zeta=\sigma e^{2\pi ik\cdot(L_j-L_i)}\f$: the
     dead rule (non-flip self-edge \f$\zeta\ne1\f$ ⇒ \f$D=0\f$; flip self-edges at k≠Γ pin
     only D's phase — never dead, conservatively not folded under) and the image fill
     \f$h_{i'j'}=\zeta h_{ij}\f$ (flip ⇒ conjugate of the whole product).  GATE
     (`StreamFoldReducedMatchesFull_SiDiamond_HalfK`): diamond at k=(½,½,½), little group
     order 12, COMPLEX frozen D=S^k, ρ 6.3e-14 / h 4.2e-15 (reordering tier; half-integer
     k ⇒ exact ±1 phases), 528→92 rep pairs (8.4×).  **T3.4b OPEN — multi-k production
     plumbing:** per-block arming needs either (i) per-fold stream caches with a UNION-of-
     reps build (win only on high-symmetry MP meshes; generic wedge points have trivial
     little groups → union → no reduction) or (ii) the STAR-SUMMED joint scatter (fold
     across the whole k-star against \f$\sum_k w_k\,\mathrm{Re}[D^k\ldots]\f$ — the real
     multi-k win, but a composite-level refactor: blocks currently collocate independently).
     Couple T3.4b to the AUTO-ARM criterion (the open-shell finding above) — both gate
     turning the fold on by default.
7. **MnO rocksalt AFM-II** (2 f.u., moments along [111]) — first *real* d-electron magnet:
   **CAMPAIGN IN PROGRESS (2026-08-04→06).**  Prerequisite DONE: the spin-polarized SAD seed
   (SCFSeedingPlan §10 increments A c08e81f8 + B 6a694c2a — spin tables incl. the Mn²⁺=3d⁵
   Atom_EC TM-cation fix, `Atom::itsSpinFlip`, `PolarizedSeedCD` + the `tSpinResolved_CD` face,
   RhoPol third branch; gates: AFM staggering algebra, ζ=0 collapse, Na re-pin; 665/665).
   The gate itself (`MnO_AFM2_RhombohedralGamma`: rhombohedral 2-f.u. cell a=8.40, q7/q6
   valence_lowq_sr Mn/O windows, IonicSAD AFM seed, free run, kT=5e-3) surfaced THREE findings:
   - **Lattice-sum truncation BUG (fixed en route):** `ForImageOffsets` enumerated to
     `rr+maxCellEdge`, assuming in-cell |dij| ≤ cell edge — FALSE for oblique cells (here |dij|
     up to 2.1× the edge): whole image shells silently dropped, Bloch overlap INDEFINITE
     (λmin −0.205 even with the trusted SIPP_SR Si control).  Fix = the exact triangle bound
     `rr+|dij|` (+ the GPW image-list sibling); kept sets provably identical where enumeration
     was already complete; suite 665/665 unmoved.
   - **Conditioning:** the dense 4-atom oblique cell saturates the s space — vet whack-a-mole
     at λ~1e-8 until Mn s thinned to 7× 0.10–24 (d 8× 0.18–36 sits at the atom floor).  Even
     the passing basis has cond(S)~7e8, and plain Cholesky EXPLODES (E~1e9 Ha at iteration 1);
     `CholeskyPivoted`+`orthoTol=1e-4` (the NaF-class recipe) restores physical energetics
     (E₁=−454.8 Ha ✓ vs 2Mn(−189.4)+2O(−13.9)+Madelung).
   - **COST (the fold motivation, measured):** the free-run magnetic cell wants 3.4e9 stream
     points vs the ~1e9 (8.6 GB) budget → 77% of pair-offset streams dropped → EVERY iteration
     pays analytic collocate+integrate sweeps: ~1 h/iteration + ~1 h build on a 14 GB box (the
     unbounded first attempt was killed at 30 h).  A converged free-run gate is compute-infeasible
     at this basis size — THE concrete motivation for extending the T3 stream fold to
     Shubnikov-imposed AFM runs (SpinAction::Flip consumption + the magnetic little group +
     T3.4b), which would make iterations GEMM replays again.  `GPW_MNO_NMAX`/`GPW_MNO_VERBOSE`
     bound/instrument hand-runs; CP2K oracle deck ready (IntegrationTests/CP2K/mno_afm2_gpw_sr.inp
     + Mn/O VALENCE-LOWQ-BASIS transcriptions).
   - **THE BLOCKER (2026-08-06): the OCCUPIED-d nonlocal PP is WRONG (~12.6x ≈ 4π-scale).**
     MnO is the FIRST system whose d-channel KB projectors are OCCUPIED (CsI's l=2 gate touches
     only virtuals — an l=2 defect was invisible to every prior anchor).  Evidence: our Mn q7
     pseudo-ATOM converges to −189.4 Ha vs the CP2K ATOM oracle **−14.243986** (l≤1 control
     species agree modulo the known atom-convention constant: F q7 ours −20.90 / CP2K −24.05
     with the CRYSTAL match at 0.19 mHa; O q6 −13.94 / −15.75); the MnO crystal (run 3, pivoted
     ortho) descends through −456 vs the CP2K crystal oracle **−61.470570** (72 Broyden steps —
     hard for CP2K too; Mulliken site moments **Mn ±4.654 μB**, O ~0, net charges ±1.69 — the
     AFM-II target physics).  RULED OUT: the parameter tabulation (byte-identical to CP2K's
     GTH_POTENTIALS), Qli/ProjR (Bessel-pair verified), LegendreP, RealYlm.  Both independent KB
     consumers (atomic mesh route + crystal analytic route) err TOGETHER ⇒ the defect sits in the
     shared `SeparablePotential` model layer — **RESOLVED FURTHER (2026-08-06, the probe
     `A_PP_Probe.DISABLED_MnOxygenAngularMeshBisect` — measured numbers + full diagnosis in its
     header): TWO SEPARATE consumer defects, model layer fully exonerated.**  (i) ATOMIC route:
     the atomic radial basis's 3-D face carries NO angular factor, so `PP_NonLocal`'s mesh overlap
     ⟨χ|βY_lm⟩ is structurally wrong against it — the default nAngular=1 rule manufactures the
     −189 catastrophe, exact angular resolution ZEROES every l≥1 projector (Mn → −8.72, its whole
     d nonlocal lost), and l=0 projectors carry a 4π-class overcount (F/O +3.15/+1.8 shallow; the
     l≤1 pinned atomic PP anchors sit on this route).  Fix: a per-l RADIAL KB assembly for atomic
     bases (analytic via `BetaGaussian`).  (ii) CRYSTAL route (−456 vs −61.4706): its own l=2
     defect in the GPW analytic KB Cartesian expansion — audit `Math.Angular`'s RAW SphericalShell
     coefficients + their unit-self-overlap normalization at the consumer.  Per-l s/p/d/f oracle
     gates (user request) land with the fix.  The gate + the Increment A Mn table entries (generated
     on the broken atomic route) both refresh after; gate DISABLED_ until then.
     Run logs: doc/logs/mno_afm2_run{1_partial,3_pivoted,4_postKBfix}.log.
   - **THE KB FIX LANDED 2026-08-06 (c2d86ec9) — and it was the ATOMIC route only.**  Root cause
     (measured by the new per-l oracle): an atomic block stores purely RADIAL χ with the irrep's
     Y_lm implicit ("fake radial" op(r), user), so `PP_NonLocal`'s 3-D mesh ⟨χ|βY_lm⟩ was
     meaningless — l=0 leaked into EVERY l block, every l≥1 projector integrated to ~1e-33.  Fix =
     `BasisSet::ImplicitAngular_IBS` + a per-l radial assembly.  ALL FOUR channels (s,p,d,**f**)
     now reproduce the analytic reference at 0.999996 with cross-l exactly 0
     (`A_PP.PerLKleinmanBylanderOracle`), and every atom matches its CP2K ATOM oracle: Mn −14.230
     (−14.2440), O −15.744 (−15.748), F −24.028 (−24.046), Si −3.7426 (−3.7470); with good bases
     Si Slater/Medium is 33 μHa and O Slater/High 255 μHa from CP2K, O matching TERM BY TERM.  Two
     "oracles" quoted in old test comments turned out to be broken-route artifacts; Na q1 (l=0+l=1)
     was silently wrong too and is now 100 μHa from CP2K.  Suite 666/666.  The removal of the
     underlying fake op(r) is CleanupCandidates **R1.0** (user direction: consumers move to
     `MakeOverlap`/`MakeOverlap3C` overloads, then op(r) becomes honest, then this capability retires).
   - **THE CRYSTAL SIDE IS **NOT** THE SAME BUG (2026-08-06).**  Both real-space routes are now
     oracle-matched on the very same Mn q7 PP (atomic radial −14.230 unpolarized; molecular
     Cartesian −14.6681 vs atomic polarized −14.6583 — they agree to 10 mHa, the 0.42 Ha vs CP2K
     being Hund polarization since the CP2K ATOM deck is restricted).  Two crystal-side findings
     replace the old "l=2 audit" guess:
     (a) **Basis conditioning is the leading explanation for MnO's ~356 Ha over-binding.**  CARTESIAN
     d carries the s contaminant, so the Mn window (7s + 8 d-shells × 6 components = 55 functions)
     is RANK-DEFICIENT before any physics — λmin 1.15e-07, cond 8.2e7, the GPW vet ABORTS on a
     single Mn atom in a box; the molecular facade survives only by dropping 5 near-null modes.  A
     near-null direction the SCF can occupy is textbook variational collapse, and MnO's 154-function
     cell (cond 7e8) sits in that regime.  SPHERICAL d (no contaminant) is the natural cure but is
     UNAVAILABLE on the GPW path (throws: no `Molecule::LatticeSum1E` — the parked S3b spherical
     work).  NEXT: regenerate a well-conditioned Cartesian Mn window (far fewer d shells, s/d
     exponents well separated — the SIPP→SIPP_SR "ill-conditioning is a BASIS problem" lesson),
     re-vet, then re-open the gate.  Probe: `GPW_SCF.DISABLED_MnAtomInBoxDChannelProbe`.
     (a2) **CURED 2026-08-06 (user's insight) + THE FIRST REAL MnO RUN.**  Cartesian d's s-contaminant
     was the redundancy: keep the d set, cut the s window to TWO (diffuse 4s tail 0.10 + one tight 24)
     — Mn box λmin 1.15e-07→3.0e-03, cond 8.2e7→2.1e3, at a cost of 2 mHa; trimming *d* instead costs
     0.55 Ha (the d set is the physics).  A SECOND, INTER-SITE redundancy then showed up in the cell
     (still 1.4e-07): the classic SR trim (Mn d≥0.38, O s≥0.39 p≥0.46) gives 122 functions at
     λmin 3.2e-03 / cond 2.4e3.  New GATE `GPW_SCF.MnAtomInBoxDChannel` = the first OCCUPIED-d species
     validated end to end through the crystal path (GPW −14.6380 vs facade −14.626).  NB CP2K solves
     our own shell list SPHERICALLY (55 Cartesian vs 47 spherical functions in its log) so it never
     sees the contaminant — the oracle comparison was apples-to-oranges in exactly the sore spot.
     **MnO run 6 (cured basis + fixed KB + regenerated seed tables) RUNS AND SETTLES** — no collapse,
     charge exact 26.000000, fingerprint DENSITY-DEGENERATE (E settled, ρ rotates), 30 iters:
     `Etot=−45.643  (Ekin 78.73  Een −84.45  Eee 29.14  Exc −14.88  Enn −57.34  E_alphaZ +3.17)`.
     TWO REMAINING DEFECTS, both now cleanly stated:
       * **THE AFM ORDER COLLAPSED**: site moments came out EXACTLY 0.0000000000 on BOTH Mn — the run
         relaxed to the non-magnetic ζ=0 solution.  This looked like the §10 trap (ζ=0 is a stationary
         point of the polarized functional), with the Kerker mixer, Fermi smearing and the missing
         constraint as co-suspects.  **SOLVED 2026-08-07 — it was none of the physics: the ρ̃ MIXERS ARE
         SPIN-BLIND.**  The order-parameter instrument (below) brackets it to a single step: the seed's
         diagonalized density carries `m_stag = +0.3656`, iteration 1 reads **EXACTLY +0.000000**, and it
         never returns.  Mechanism, structural and complete: `KerkerMixer`/`PulayMixer` hold ONE
         `FourierMixCD` — a bare ρ̃(G) map, the ↑+↓ TOTAL, with no spin channels — and `FockDensity()`
         hands *that* to every Fock build from iteration 1 on.  `XC_GridEngine::RhoPol` then finds neither
         a `cPolarized_CD` nor a `cSpinResolved_CD` and takes its third branch, the ρ↑=ρ↓=ρ/2 "spin-agnostic
         seed" collapse, so v_xc^↑ ≡ v_xc^↓, F↑ ≡ F↓, and the run is UNPOLARIZED from iteration 1.  The
         staggered seed never had a chance; ζ=0 was not a basin the SCF fell into, it was the only state
         the Fock could represent.  **A/B PROOF** (same run, linear D-mixing, which `MixIn`s the polarized
         density itself): `m_stag` = 0.203, 0.085, 0.136, 0.026 … and ε↑−ε↓ splits by up to 45 mHa —
         alive, sloshing with the energy (linear mixing is unstable here: E swings −51/−55/−44/−28, which
         is exactly why Kerker was adopted).  STOP-GAP LANDED: `MakeDensityMixer` now REFUSES the ρ̃ mixers
         on a polarized/spin-resolved density and falls back to linear D-mixing **loudly** (the file's own
         "[Mixer] DISABLED" idiom) — slower and sloshier, but never silently non-magnetic.  Diagnostic
         valve `QCHEM_SPINBLIND_KERKER=1` restores the broken arm to re-measure the collapse.  BLAST RADIUS
         is MnO alone: no other Kerker/Pulay caller sets multiplicity≥1 (NaF/Si are closed shell) and the
         molecular path never builds a ρ̃ mixer — but this would have bitten EVERY magnetic solid.  Gate
         `GPW_SCF.PolarizedRunKeepsItsSpin` (Mn q7 sextet box asking for Kerker, 38 s).  **WHICH OBSERVABLE
         HAS TEETH is itself a finding**: the on-site MOMENT does not.  At nUp=6/nDown=1 the moment is
         pinned by the CHANNEL OCCUPATIONS and survives a spin-blind Fock intact (0.64 vs 0.72) — an atom
         cannot expose this bug at all.  Only an order that must be SELF-CONSISTENTLY SUSTAINED can, which
         is precisely MnO's staggering at nUp=nDn, where nothing but v_xc^↑≠v_xc^↓ holds it up.  What the
         atom DOES expose is the ENERGY: without the exchange splitting it lands 68 mHa HIGH (−14.57 vs the
         physical −14.638), so the gate reports the moment and asserts the energy.
       * **THE REAL FIX LANDED 2026-08-07: `PolarizedDensityMixer`** (user's design) — ONE ordinary mixer
         per spin channel plus `PolarizedMixCD`, the view that gives their pair both faces the framework
         consumes (`FourierDensity` = the ↑+↓ total for Hartree; `cSpinResolved_CD` = the channels for
         `RhoPol`).  A COMPOSITION, not a second implementation: the leaves come from the same
         `MakeGSpaceMixer` the unpolarized path uses, so Kerker **and** Pulay are supported without the
         composite knowing which it holds.  **CHANNEL-BASIS QUESTION SETTLED against CP2K's source**
         (`qs_gspace_mixing.F`): CP2K transforms to (ρ_total, m) in its mixing DRIVER and then runs the
         same per-channel loop with the SAME Kerker factor and α on both — so it *does* Kerker-damp m, and
         runs Pulay/Broyden on it with independent per-channel histories.  For a LINEAR operator that is
         algebraically identical to per-channel (ρ↑,ρ↓) mixing (Kerker is linear in the residual ⇒
         m_mix = m_in + αK(m_out−m_in) either way), so this composite reproduces CP2K's Kerker EXACTLY,
         with no proof burden.  The basis only becomes a real choice for the NONLINEAR history mixers, and
         it is isolated in this one class.  (CP2K also floors its damping at kmin=0.1 and applies it only
         in a G band — ours is unbounded below; and its G=0 factor stays 1.0, matching our own
         deliberate "G=0 mixes at full α" convention, which a stale comment on `KerkerMixer` denied.)
       * **A LATENT NaN, exposed the moment the channels became real**: `VWN_Correlation::Zeta` was
         unclamped, and a band-limited ρ̃ rings slightly NEGATIVE in the tails — so a genuinely spin-resolved
         feed hands the functional ρ_σ<0 with the total still positive, i.e. |ζ|>1.  That is FATAL, not
         merely inaccurate: f(ζ) evaluates (1−ζ)^{4/3} via `std::pow`, NaN for a negative base, and ONE NaN
         mesh point makes the whole v_c matrix non-Hermitian (Blaze throws "Invalid assignment to diagonal
         matrix element" out of the Fock build).  Unreachable before the fix, because ζ was identically 0.
         Cured by clamping ζ to its physical [−1,1] at the ONE place both `GetEpsC` and `GetVc` read it —
         the same self-defence as the existing ρ>0 guards.  Gate `LDA_XC.NegativeChannelClampsToFullPolarization`.
       * **RESULT (MnO run 9, `doc/logs/mno_afm2_run9_polarizedmixer.log`): the structural collapse is GONE,
         the AFM ground state is still NOT reached.**  m_stag no longer reads EXACTLY 0 at iteration 1
         (+0.0046), it *regrows* to +0.107 by iteration 3 — the self-consistent exchange splitting doing its
         job — and thereafter it tracks the energy's sign, decaying to ~1e-4 as the run settles into the SAME
         −45/−46 limit cycle as run 6.  So the remaining barrier is no longer the mixer.  (Run 9 was
         OOM-killed at iteration 24 under a 12 GB cap — constant two-channel overhead on a run already at
         the ceiling, NOT a leak: an RSS trace of the two-mixer Mn box is flat at ~295 MB, identical to the
         single-mixer control.)
       * **WHAT THE REMAINING BARRIER ACTUALLY IS (run 10 sweep, 2026-08-07): OCCUPATION HOPPING, and the
         MAGNETIC BRANCH IS 4–6 Ha BETTER BOUND.**  Every iteration carries the `cfg *` flag and an `m`
         (partial-occupancy) frontier whose gap wanders over an order of magnitude (0.026–0.27 Ha): the
         occupied configuration re-shuffles EVERY step on the near-degenerate d manifold, and each reshuffle
         throws E by several Ha.  The run alternates between a MAGNETIC branch (−53.87 m=0.097; −51.95
         m=0.106) and a NON-magnetic one (−47.64 m=0.0025; −48.13 m=0.0028) — the magnetic branch is
         consistently the better bound.  So run 6's converged −45.64 was effectively the AVERAGE of the two,
         and "under-bound by 15.8 Ha" + "the AFM order collapsed" are plausibly ONE defect, not two: the
         magnetic branch alone reaches −53.9, closing roughly half the deficit to CP2K's −61.47.  Re-test
         once the run can hold a configuration.  No density-mixing knob can fix a run whose occupations hop
         every step ⇒ the next experiment is **kT / MOM-timing** (`RunGpwAnnealed` exists for exactly this),
         not mixing.
       * **α AND G₀ MEASURED (`MNO_ALPHA` / `MNO_KERKER_G0`, 16-iteration arms; logs `run10_alpha025_sweep_g*`).**
         α: 0.45 was the aggressive half of the borrowed NaF recipe (the NaF *Becke* gate); **α=0.25** (the
         NaF *Γ* slosh recipe, user) descends MONOTONELY to −53.87 instead of overshooting −51.6→−44.2, and
         settles ~4 Ha lower.  G₀ at α=0.25: **0.3 → E flat at −45.3 from iteration 5, moment extinguished
         (converges smoothly to the WRONG, non-magnetic answer — the worst outcome wearing the best-behaved
         run's clothes); 1.0 → oscillates, peak m=0.106, magnetic branch −53.87; 3.0 → slow monotone descent,
         STILL descending and the ONLY arm whose moment is still GROWING at the budget (0.028 and rising).**
         That ordering is the selectivity argument: at G₀=1 the Kerker factor is 0.61 on the AFM mode
         (|G|=1.24, the shortest odd-(h+k+l) mode of B=(4π/a)(I−J/4)) vs 0.30 on the lowest charge mode
         (|G|=0.65), so Kerker already favours the magnetism 2:1; G₀=3 sharpens that to 3.3:1 but crushes
         both (effective AFM step 0.146·0.25=0.037 — hence the slowness).  **The limit of the trend is not a
         bigger G₀, it is NOT KERKERING m AT ALL** — the user's original physics instinct: Kerker's
         G²/(G²+G₀²) models the Hartree restoring force against long-wavelength CHARGE fluctuations, which
         the magnetization does not have.
       * **(ρ,m) BASIS BUILT AND MEASURED — THE HYPOTHESIS IS REFUTED AS A CURE (run 11, 2026-08-08).**
         `PolarizedDensityMixer` gained the (ρ,m) channel basis (`QCHEM_MIX_RHO_M=1`): Kerker leaf on ρ,
         PLAIN LINEAR leaf on m (= Kerker at G₀=0, so no new mixer type was needed).  It is a 2×2 COUPLED
         mix in spin space — provably unreachable by any per-(ρ↑,ρ↓) leaf pair, which is exactly why the G₀
         sweep could not test it (a uniform filter on (ρ↑,ρ↓) IS that filter on (ρ,m)).  A/B at identical
         α=0.25 and G₀-on-ρ: peak m_stag 0.0618@it7 vs 0.1064@it7; final −0.0029 vs 7.3e-5; final E −49.08
         vs −49.12.  So undamping m leaves ~40× more residual moment (the effect is REAL) but LOWERS the
         peak and rescues nothing — both arms fall into the SAME two-branch oscillation at the same energy.
         **Kerker damping of m was not what killed the AFM order.**  Every mixing-side candidate is now
         fixed or eliminated: spin-blind ρ̃ mixer (a real bug, fixed), DIIS-per-spin (checked in source — it
         is JOINT: one B summed over both spins, one coefficient vector held by reference), G₀ selectivity
         (swept), Kerker-on-m (built + measured).  **THE REMAINING DEFECT IS THE OCCUPATION HOPPING.**  The
         (ρ,m) basis stays an opt-in capability (default unchanged): it is the physically-motivated
         construction and the one place we deliberately diverge from CP2K, but it is not the MnO cure.
         The `m_stag` assert stays the acceptance test.  Log `run11_rho_m_basis`.
       * **THE OCCUPATION IS THE DEFECT — CONFIRMED BY A CLEAN A/B, AND MOM IS WORTH ~9 Ha (runs 12/13,
         2026-08-08).**  Two instruments were meant to hold the configuration and NEITHER had ever fired:
         - **`UseMOM` was INERT in every MnO run to date.**  `tIrrepWF::FillOrbitals` gives `SmearingkT>0`
           precedence, and that branch consults the MOM reference ONLY when `MOMSmearPenalty>0` (default 0).
           With `kT=5e-3` set alongside `UseMOM=true`, MOM never changed a single occupation.  Reading a
           recipe off the options struct is not evidence that it ran.
         - **Even reached, the capture was one iteration too late.**  `SetMOM` is pushed in by `Iterate`,
           i.e. AFTER `Init`'s seed fill, so `itsUseMOM` is false when the seed is filled and the earliest
           self-capture is iteration 1's occupied set — after m_stag has already fallen 0.366 → 0.0046.  The
           delayed-IMOM design is built for a seed you want to settle AWAY from (the NaF diving ghost); MnO's
           seed IS the physics.  Cure: `GpwOptions::momFromSeed` self-adopts the seed's own occupied subspace
           through the existing grid-continuation face (`AdoptMOMReference`), before iteration 1 runs.
         - **The 0h MOM guard would have silently undone the experiment.**  It releases the reference after 3
           consecutive non-aufbau (hole) iterations, on the NaF-calibrated premise that a hole means a wrong
           reference.  That premise fails here: until the self-consistent exchange splitting opens, the
           magnetic branch's occupied states sit ABOVE empty ones in the raw ε ordering, so the hole is the
           SIGNATURE of the branch being held.  `MNO_MOM_HOLD` makes the persistence sweepable.
         **THE A/B (identical α=0.45, kT=0, everything else equal; only MOM differs):**
         | | m_stag trajectory | E (20 iters) | [F,D] | symmetry defect |
         |---|---|---|---|---|
         | run 13 aufbau (control) | dies at iteration 1, flat ±0.01 | −45.56, limit cycle | **stuck at 1.43** | 0.329 |
         | run 12 seed-pinned MOM | **0.19–0.30 throughout, final 0.29** | **−54.73** | 4.11 → 0.05–0.17 | 0.021 |
         The control's converged ε↑−ε↓ is ~1e-3 across the board: no exchange splitting, i.e. the textbook
         non-magnetic collapse.  **So holding the occupation buys BOTH the order AND ~9 Ha of binding, which
         is the standing "one defect, not two" hypothesis confirmed** — it closes 9 of the 15.8 Ha deficit to
         CP2K's −61.47, and it is the first MnO run whose [F,D] descends at all.  Logs `run12_momseed_cold`,
         `run13_cold_aufbau_control`.
       * **BUT run 12 IS NOT YET THE ANSWER — three things are wrong with it, and the instruments said so:**
         (i) **The moment is NOT staggered**: sites read (+0.0006, −0.579), i.e. ONE Mn carries all of it.
         `m_stag = ½(m1−m2)` is DEGENERATE between (+m,−m) and (0,−2m), so the trajectory line cheerfully
         reported "order SURVIVED" on a physically wrong state.  The readout now prints `m_net = m1+m2`
         beside it and flags the collapse; the `EXPECT_NEAR(m1,−m2)` gate was always honest and still fails.
         (ii) **A DEEP HOLE**: the run ends non-aufbau with an EMPTY level at ε=−1.29 Ha, ~2 Ha below the
         highest occupied one, and six levels below the printed frontier are empty in both channels.
         (iii) It still oscillates (±1 Ha, 5 sign flips in the last 8).
         So the next question is not "does MOM help" (settled) but **WHICH STATES the seed reference is
         refusing to occupy** — and whether refusing them is the cure or a second, opposite error.
       * **IT IS A SECOND, OPPOSITE ERROR: THE HARD PIN DEPLETES ONE Mn's d SHELL (run 14, the spectrum the
         display fix exposed).**  The doubly-empty levels are **exactly FIVE**, clustered inside 0.07 Ha at
         ε = −1.39…−1.32, sitting 0.8 Ha BELOW the highest occupied level — an l=2 manifold, and four of the
         five carry ε↑−ε↓ = 0.00000000 to eight digits.  **Exact spin degeneracy is the tell: an EMPTY shell
         feels no exchange splitting.**  Put together with the site readout (Mn1 −0.002, Mn2 −0.425), the
         state is unambiguous: **Mn1's 3d shell is EMPTY — a d⁵/d⁰ disproportionation, not antiferromagnetism.**
         Mn1 has no d electrons, hence no moment; its unscreened d levels sink to −1.4; and having sunk, they
         have no overlap with the seed reference, so the pin keeps them empty — a SELF-REINFORCING HOLE.
         That also disposes of run 12's tempting −54.73: a configuration with a 2.4 Ha hole is not
         variationally admissible, so the 9 Ha is not (yet) real binding.  **Both failures are now one
         sentence: the aufbau run cannot HOLD a configuration, and the hard MOM pin cannot LET GO of one.**
       * **WHAT THE MECHANISM PRESCRIBES**: the pin must be strong enough to break the ~0.03 Ha near-degenerate
         ties that drive the hopping, yet weak enough to be OVERRIDDEN by a state that dives below the
         frontier.  That is exactly `MOMSmearPenalty` Λ under Fermi smearing: the shift is Λ(1−s)², so any
         dive deeper than Λ re-occupies the state no matter how foreign its character.  Λ ≈ 0.3 Ha brackets
         the two scales (≫0.03 tie, ≪2.4 dive) — a two-sided argument, not a fitted knob.  Both extremes are
         already measured and both fail: Λ=0 (every run through 11) is the aufbau collapse, Λ=∞ (the cold pin,
         runs 12/14) is the d⁰ disproportionation.
       * **Λ SWEPT, AND THE Λ SCALE IS SET BY THE SCORE DISTRIBUTION, NOT BY THE TIE (runs 15/16).**
         **Λ=0.3 is indistinguishable from no MOM at all** (run 15: m_stag dead by iteration 1, E in the −46
         limit cycle, [F,D] pinned at 1.36).  **Λ=1.5 transforms the run**: it descends to **−60.08 Ha at
         iteration 11 — 1.4 Ha from CP2K's −61.4706**, against the −45.6 the aufbau run cannot leave, and the
         deep hole collapses from 2.4 Ha to 0.21 Ha, i.e. **the d⁰ disproportionation is CURED** (sites
         −0.215 / −0.419: both Mn magnetized now, though same-sign, so still not staggered).
         **WHY the 5× matters — the new `QCHEM_MOM_SCORES` instrument (one-shot, off by default) settles it.**
         At the first fill the up-spin scores run 0.952 → 0.687 with a cut gap of only **0.0147**, so the
         shift SEPARATING the last occupied from the first virtual is Λ·[(1−0.687)²−(1−0.701)²] ≈ 0.0086·Λ —
         **2.6 mHa at Λ=0.3**, against a raw-ε spread of 0.2–1 Ha.  The masked-Fermi shift is a hopelessly
         weak discriminator AT THE CUT.  What Λ=1.5 actually buys is the shift on the FOREIGN states
         (s≈0.49 ⇒ Λ(0.51)² ≈ 0.39 Ha), which is comparable to a dive.  **So the two MOM paths are not two
         strengths of one instrument, they are two different instruments**: the cold path RANKS (a 1.5%
         score difference decides the fill), the smeared path SHIFTS (only differences worth more than Λ in
         energy decide anything).  Λ must be scaled to the physical-vs-foreign score GAP, never to the tie.
       * **STILL NOT CONVERGED — the honest reading of run 16.**  Its fingerprint is **DIVERGING** (±5.9 Ha
         over the last 8 iterations), so **−60.08 is a TRANSIT, not a result**, and m_stag still flips sign
         iteration to iteration.  What run 16 does establish is that the ENERGY SCALE is now reachable: the
         "under-bound by 15.8 Ha" gap was, as hypothesised, mostly the occupation.  With the occupation
         instrument fixed, the instability is back to being a STEP-SIZE question — and the α evidence banked
         in run 10 was gathered under the collapsed aufbau dynamics, so it does not transfer and α must be
         re-swept on top of Λ.  Logs `run15_maskedfermi_L0.3`, `run16_maskedfermi_L1.5`.
       * **★ THE AFM-II STATE IS REACHED (run 17, 2026-08-08): Λ=1.5 AT α=0.25, 40 iterations.**  Re-swept as
         the previous bullet demanded — α was the only change from run 16 — and it lands:
         | | value | reference |
         |---|---|---|
         | E | **−60.1414** | CP2K −61.4706 → **1.33 Ha**, was 15.8 |
         | m(Mn1), m(Mn2) | **+0.625, −0.657** | **m_net = −0.032: STAGGERED** (the flag does not fire) |
         | m_stag | **0.641**, monotone from iteration 15, STILL RISING at 40 | seed 0.366 — it GREW past the seed |
         | fingerprint | **DENSITY-DEGENERATE (E settled, ρ rotates — benign)** | Eamp(last 8) 0.057 Ha, relAmp 9.5e-4 |
         | [F,D] | 4.11 → 0.015–0.075 | was STUCK at 1.43 |
         The spectrum now carries a real, large exchange splitting (ε↑−ε↓ = −0.086, −0.159, −0.247 …) instead
         of the control's uniform ~1e-3.  **This is the first MnO run that is simultaneously magnetically
         ordered in the AFM-II pattern, at the right energy scale, and settling** — and it is the SEED's own
         staggering, grown rather than merely survived.  Log `run17_L1.5_alpha025`.
       * **RUN 18 — the same recipe to 90 iterations: E = −60.9247, i.e. 0.55 Ha from CP2K (from 15.8).**
         The trajectory has three phases and all three are informative:
         - **iterations 1–58, Kerker+DIIS**: settles hard at E=−60.14, m_stag 0.66, Δρ 4.3e-4, [F,D] 7.8e-3 —
           and **the `cfg` flag goes BLANK for the first time in the campaign.  The occupation hopping that
           opened this investigation is GONE**: the run holds one configuration iteration after iteration.
         - **~iteration 70, the Ladder hands off to GDM** (direct minimisation, no mixing): E drops FURTHER
           to −60.92, so the DIIS fixed point was not the bottom.  Δρ turns noisy (GDM owns the density
           update, so the mixer's Δρ is not its convergence measure) and m_stag settles lower, 0.52.
         - **final**: m(Mn1)=+0.506, m(Mn2)=−0.529, **m_net=−0.022 — STAGGERED**, order SURVIVED, peak
           m_stag 0.668@64, fingerprint "E settled, ρ rotates — benign".
         Log `run18_L1.5_alpha025_long`.
       * **THE GDM LEG IS NOT DOWNHILL — AND IT IS THE OCCUPATION SEAM, NOT E[ρ] (2026-08-09, user's read).**
         User's observation on run 18: GDM goes haywire for ~5 iterations (64→68: −60.13, −45.60, −45.46,
         −42.15, **+0.44** — a 61 Ha excursion, [F,D] to 118) and is not strictly downhill in general; the
         usual reading is a non-variational E[ρ].  **That reading is already excluded for this stack** by work
         banked above `DISABLED_NaFImposedGDMSmearProbe`: ROUND 3 (user-measured, 2026-08-03) found STRICTLY
         NEGATIVE dE for imposed multi-k GDM on SR2, and the Si-diamond IBZ proxy is flawless imposed+GDM.
         E[ρ] is variational; the inconsistency is in the OCCUPATION SEAM, and the code says exactly where:
         - `tSCFIrrepAcceleratorGDM::OrbitalsAt` rotates a **FIXED occupied block** of size `itsNocc` (fixed at
           construction) and returns `e` = [occupied εs, then virtual εs] **in block order, NOT sorted**.
         - but `tSCFIterator::DirectMinStep` evaluates each trial point via `MoveOrbitals(t,false,mergeTol)`,
           and `tCompositeWF::MoveOrbitals` calls **`FillOrbitals(mergeTol)`** — so **every backtrack re-decides
           WHICH orbitals are occupied**, by aufbau, by Fermi μ, or (now) by MOM overlap.
         - Whenever that re-fill disagrees with the geodesic's own occupied block, **the energy being line-
           searched is not the energy along the geodesic.**  The direction was derived for one manifold; E is
           measured on another.  No monotonicity is owed, and none is observed.
         **The `h` flag PROVES the disagreement on every GDM row of run 18**: a hole IS "the fill did not take
         the lowest `no`".  So this run cannot be on-path even in principle.
         **It also explains the whole banked pattern at once**: Si (clean gap ⇒ aufbau fill always reproduces
         the geodesic block) is perfect; NaF SR (diffuse near-degenerate frontier ⇒ fill can disagree) walks
         uphill +23–56 mHa; NaF SR2 (diffuse directions removed) is strictly downhill; MnO (MOM at Λ=1.5 AND
         kT=5e-3 both re-deciding the fill at every trial point, with a standing hole) is the worst case.
         NB the earlier "MOM is ELIMINATED" round-2 verdict on NaF stands only for the COLD path it tested;
         no NaF arm ever exercised **masked-Fermi MOM through a GDM leg**, which is what run 18 did.
         **THE FIX IS THE LINE SEARCH, NOT THE FUNCTIONAL** — two independent defects, both local:
         (i) **Evaluate the trial at FIXED occupations.**  The line search must measure E along the geodesic,
             i.e. hold the occupied block the direction was built for and refill only on commit.  A refill
             inside the search makes E(t) discontinuous in t.
         (ii) **`DirectMinStep` COMMITS UNCONDITIONALLY.**  After 12 backtracks it runs
             `MoveOrbitals(t,true,mergeTol)` whether or not any t lowered E — `found` feeds ONLY the
             `GPW_GDMTRACE` print.  So a failed line search still takes a step, and the method is
             non-monotone BY CONSTRUCTION.  This is precisely the "actual-E-checked step acceptance (clamps
             the uphill)" that the ROUND 3 note listed as the required engine fix and which was never
             implemented.  A failed search should decline the step (or degrade to the mixed fallback the bail
             path already builds), not commit it.
         Neither is an E[ρ] problem, and neither is MnO-specific.
       * **MEASURED (run 20, `GPW_GDMTRACE=1` across the hand-off) — AND IT IS A DISCONTINUITY, NOT A GRADIENT
         ERROR.**  Six geodesic steps: 4 DESCENT, 2 FALLBACK.  The two FALLBACKs are QUALITATIVELY different
         and that is the whole answer:
         | step | verdict | best(E_t − E_cur) over 12 backtracks |
         |---|---|---|
         | iteration 65 | FALLBACK t=2.44e-04 k=12 | **+1.45e+01 Ha** |
         | iteration 76 | FALLBACK t=2.44e-04 k=12 | **+3.22e-05 Ha** |
         **THE ARITHMETIC CLOSES EXACTLY**: iteration 64 ends at E=−60.125186, the search reports its best
         trial is +14.53 Ha uphill, the code commits anyway, and iteration 65 reads −45.595130 — a jump of
         **+14.530 Ha**.  The reported blow-up IS the committed non-descent step, to five figures.
         **Why +14.5 Ha refutes the non-variational reading.**  As t→0 the geodesic tends to the IDENTITY, so
         E(t) → E(0) = E_cur *continuously, provided the occupation is held*.  A floor of +14.5 Ha still
         standing at t=2.44e-04 is therefore not a small gradient inconsistency — **it is a DISCONTINUITY in
         t**, and the only thing in the step that can jump discontinuously is the OCCUPATION, which
         `MoveOrbitals → FillOrbitals` re-decides at every trial point.  A non-variational E[ρ] cannot produce
         this signature; a re-fill that lands on a different branch produces exactly it.
         **And the run itself measures the honest noise floor**: the SECOND FALLBACK, at convergence, is
         +3.22e-05 Ha.  THAT is where a genuine E/H inconsistency would live (fit-vs-quadrature v_xc, grid
         charge leak) — **32 μHa, i.e. five orders of magnitude below the excursion.**  So E[ρ] is variational
         to ~1e-5 Ha here, exactly as the SR2 and Si results already said.
         The rest of the "5 haywire iterations" is a cascade from that ONE step: iterations 66–68 print NO
         `[gdm]` line at all, i.e. `BuildFockAndComputeSteps` FAILED and they silently took the bail path's
         diagonalize+mix from a wrecked density ([F,D] 44.7 → 118), recovering only by 69–71.
       * **THE OBVIOUS FIX IS BLOCKED, and that is worth knowing before someone tries it.**  "On `!found`,
         degrade to the mixed step" cannot reuse the existing bail path: that path is safe only because
         `ComputeStep` FAILED there, so `tSCFIrrepAcceleratorGDM::NextOrbitals` diagonalizes.  After a failed
         LINE SEARCH `ComputeStep` has SUCCEEDED, and `NextOrbitals` then returns `OrbitalsAt(1.0,true)` —
         **the full step, committed**, i.e. the very t=1 the search rejected first.  So declining a step needs
         an accelerator-level DECLINE/RESTART (drop the CG history, force the next call to diagonalize) or a
         trust-radius shrink-and-retry.  That is a `doc/SCFStrategyPlan.md` seam decision (the accelerator
         reports its mode; the iterator selects the driver), not a patch to `DirectMinStep` — left for the
         user rather than done unilaterally.  Log `run20_gdmtrace`.
       * **THE BAIL-OUT LANDED, AND IT IS NECESSARY BUT NOT SUFFICIENT (run 21, the A/B against run 20).**
         `tSCFAccelerator::RejectStep` + the retry loop in `DirectMinStep` behave exactly as specified:
         - **The rejected step is no longer taken.**  Iteration 65 now holds **E=−60.127** with m_stag 0.666,
           where run 20 committed the +14.5 Ha step and read **−45.595** with m_stag 0.116.
         - **The retries prove the diagnosis.**  Seven FALLBACKs report best(E_t−E_cur) = **+1.45e+01 EVERY
           TIME**, identical to three figures across six trust-radius shrinks.  A step that is merely too
           LONG improves as the trust radius falls; one that is scale-INVARIANT is a discontinuity.  The
           retries are, in effect, a controlled experiment that re-confirms the occupation reading.
         - Then EXHAUSTED → mixed fallback, as designed.
         **BUT THE RUN ENDS WORSE, NOT BETTER: −29.96, DIVERGING, m_stag decayed to 0.061** (run 20 ended
         −60.84).  The reason is visible in the trace and it is the SAME defect: iteration 66, the EXHAUSTED
         **mixed fallback itself**, drops to −45.37 — its `ρ_mix` column reads `Ker 0.25`, so this is the
         diagonalize-and-refill path, not the geodesic.  **Every path that REFILLS can jump branch**: the
         geodesic trial, the geodesic commit, the mixed fallback, the plain diagonalize.  Refusing the bad
         geodesic step just moves the branch jump one iteration downstream.
         Do not read run 20's −60.84 as the better answer either: it got there by committing a +14.5 Ha step
         and happening to recover.  Both runs are non-monotone and neither converged; what changed is that
         the failure is now HONEST (declined and reported) instead of silently committed.
       * **SO THE NEXT FIX IS THE OCCUPATION HOLD, and it is a SEMANTIC decision, not a patch.**  GDM rotates
         a FIXED occupied block (`itsNocc`); for the line search to mean anything, the fill must reproduce
         that block.  Plain aufbau ALREADY DOES — `TakeElectrons(ne)` fills in STORED order and
         `UpdateOrbitals` does not sort, so the stored order IS the geodesic's [occupied | virtual] order.
         **It is specifically `TakeElectronsFermi` (μ solved over the whole spectrum) and MOM
         (`TakeElectrons(ne,priority)`, ranked by overlap) that deviate** — which is exactly why Si and
         NaF-SR2 (cold, aufbau) are clean and MnO (kT=5e-3 AND Λ=1.5) is the worst case.  Two ways out:
         (a) **hold the block through the GDM leg** (fill in stored order for trial AND commit, so E(t) is
             the geodesic energy) — GDM then minimises E at integer occupation regardless of caller
             settings, and smearing/MOM are SUSPENDED for that leg, which must be announced;
         (b) **enforce the precondition** — the ladder does not hand off to GDM while smearing or MOM is
             active (the `RunGpwAnnealed` accSchedule pattern, `{"DIIS","GDM"} × {kT,0}`, already says this
             in its own comment: "GDM does not support Fermi smearing, so its stage must run kT=0").
         (a) makes GDM self-sufficient; (b) keeps one functional per leg but leaves MnO's GDM leg unusable
         while MOM is what holds its configuration.  Note either way the reported E jumps at the hand-off
         when fractional occupations become integer — that is honest, not a bug.
         Log `run21_bailout`.  New coverage: `UTSCFAccelerator` (5 tests on the contract; GDM had none).
       * **★ ROOT CAUSE, AND IT IS OLDER THAN ANY OF THIS: GDM HAD AN UNSTATED PRECONDITION (2026-08-09).**
         Chasing the +14.5 Ha through the held fill produced a trace line that named it outright:
         `[gdm] convention shift: E(t=0) held = −45.6 vs previous-iteration E = −60.1 (offset 14.5 Ha)`.
         The +14.5 Ha was never a bad direction, never a discontinuity *within* the search, and never a
         non-variational E[ρ] — **it is the energy gap between the MOM/smeared occupation and the LEADING-
         BLOCK occupation at the very same orbitals.**  Because every block GDM uses — `Cocc`, `Cvir`, the
         o–v gradient, the geodesic — is carved out of `itsCp` BY INDEX:
             `Cocc = submatrix(itsCp, 0,0, n,no)`
         i.e. **GDM assumes the occupied states are its leading `itsNocc` columns.**  Aufbau makes that true.
         MOM makes it false BY DESIGN (it occupies by overlap, deliberately not the lowest); Fermi smearing
         makes it false too.  `UseFD` receives the true D′ but only ever took its commutator norm, so nothing
         noticed: **the gradient was built for one manifold while the density occupied another** — silently,
         and since long before this session.  The minimizer was diligently optimising a configuration the run
         had already rejected.
         `UseFD` now checks it (`Tr(D′P_block)` over the leading `itsNocc` columns = `itsNocc` exactly under
         aufbau) and GDM DECLINES when it fails — gating `ComputeStep` as well as `Ready()`, since
         `NextOrbitals()` takes the step whenever `ComputeStep` succeeds.  **The measured violation on MnO is
         enormous and exactly the AFM configuration**: `Tr(D′P_block)` = **11 of 13** in ↑ and **7.00 of 13**
         in ↓ — six minority-channel electrons, essentially the whole d⁵ shell, outside the block GDM meant
         to rotate.  Preferred over gating on "is MOM on" because it tests the property GDM actually needs,
         so it covers smearing and any future occupation rule without knowing the caller's policy.
       * **★ RUN 24 — THE BEST-BEHAVED MnO RUN OF THE CAMPAIGN.**  With the GDM rung declining to
         diagonalizing steps:
         | | run 20 (none) | run 21 (bail-out) | **run 24 (precondition)** |
         |---|---|---|---|
         | E final | −60.84 *(after a +14.5 Ha excursion)* | −29.96 | **−60.1341** |
         | fingerprint | DIVERGING | DIVERGING | **FIT-FLOOR STALL** |
         | m_stag | 0.517 | 0.061 | **0.6606** |
         | sites / m_net | — | not staggered | **+0.652 / −0.669, m_net −0.017 (STAGGERED)** |
         Iterations 65→76 descend **monotonically** (Δρ 5.3e-3 → 1.4e-4, [F,D] 1.3e-2 → 2.4e-4) and the run
         stops because it hits the **fit/grid floor** — an honest, different problem — not because it went
         unstable.  The oscillation that has dogged MnO since run 3 is GONE, and the moment is the largest
         yet.  Do NOT compare against run 18's −60.92: that came from a wandering GDM leg optimising the
         wrong manifold.  **−60.134 is the trustworthy number** (CP2K −61.4706 ⇒ 1.34 Ha).  Log
         `run24_precondition`.
       * **★ THE OCCUPATIONS ARE THE WHOLE STORY — AND THE ANSWER IS +U, NOT A BETTER RECIPE (user's read
         + run 26, 2026-08-10).**  USER, reading run 24's spectrum: "for spin up I see 3 empty but deep
         orbitals, and for spin down I see 7".  Confirmed exactly — empty levels below each channel's OWN
         HOMO: ↑ n=6,7,12 (three); ↓ n=4,5,6,7,8,11,12 (**seven**), the ↓ channel occupying +0.378…+0.442
         while leaving −0.184 empty, ~0.6 Ha above what it refuses to fill.  My reading of that was that Λ=1.5
         had outlived its job and was pinning an EXCITED state.  **Run 26 refutes it.**  Two stages, kT fixed,
         Λ annealed 1.5 → 0 (the new per-stage penalty schedule, `MNO_ANNEAL_PENALTY`):
         | stage | Λ | E | m_stag | spectrum |
         |---|---|---|---|---|
         | 1 | 1.5 | **−60.14** | 0.657 | non-aufbau (the 3/7 deep empties) |
         | 2 | 0 (RELEASED) | **−49.53** | **0.0001, DIED at iteration 44** | clean aufbau, ε↑−ε↓ ≈ ±2e-4 |
         **Releasing the constraint raises the energy by 10.6 Ha and extinguishes the moment.**  So the deep
         empty levels are NOT the signature of a trapped excited state: the constrained state is far BETTER
         bound, and it is the one near CP2K (−61.4706).  What the release finds is a NON-MAGNETIC METAL — no
         exchange splitting anywhere (±2e-4), a fractionally-occupied frontier (0.98/0.98/0.94/0.18), and the
         benign DENSITY-DEGENERATE rotation.  Two distinct self-consistent branches, and plain LSDA aufbau
         cannot reach the good one.
       * **THE UNCOMFORTABLE PART, STATED PLAINLY.**  A converged NON-AUFBAU state is formally a saddle of the
         KS functional, not a minimum — so the magnetic branch as currently obtained is not yet a defensible
         ground state, even though it sits 10.6 Ha below the aufbau state we CAN reach and 1.34 Ha from CP2K.
         Those two facts together are the finding: **LSDA aufbau cannot reach the magnetic state, and the
         magnetic state is not aufbau-stable under LSDA.**  MOM is not corrupting the answer, it is SELECTING
         a branch the functional will not select for itself — which is exactly the LDA/GGA delocalisation
         error V1.29 names, and exactly what +U exists to fix: it pushes the occupied d down and the empty d
         up until the magnetic configuration IS the aufbau one.  So "+U" moves from follow-on to THE next
         step, and the acceptance test becomes "the magnetic solution converges aufbau-stable WITHOUT MOM",
         not "m_stag survived".  Log `run26_lambda_release`.
       * **⛔ THE +U CONCLUSION ABOVE IS RETRACTED, SAME DAY — CP2K IS NOT DOING +U (user's question).**
         Checked the banked deck `IntegrationTests/CP2K/mno_afm2_gpw_sr.inp`: **no `&DFT_PLUS_U`, no Hubbard
         keyword anywhere**, plain `LDA_X + LDA_C_VWN`.  What CP2K uses to reach −61.4706 is a `&BS`
         BROKEN-SYMMETRY guess on two distinct Mn KINDs (±5 d occupation, applied to the ATOMIC guess only),
         **Fermi smearing at 5.0E-3 au — exactly our kT**, Broyden mixing (α=0.2, β=1.5, NBUFFER 8),
         `ADDED_MOS 20`, plain diagonalization, and **no MOM and no constraint of any kind** — same
         functional, basis, PP, 144 Ry cutoff and kT as ours.
         **So the magnetic state IS reachable and aufbau-stable under LSDA, and our inability to reach it is
         OUR defect, not the functional's.**  The inference above was drawn from our own failure mode without
         checking whether the oracle needed the cure — backwards, with a same-recipe oracle sitting in the
         repo.  A failure to converge is evidence about the code until an oracle says otherwise.
       * **WHAT THE COMPARISON POINTS AT INSTEAD.**  Almost nothing differs except the SCF machinery, so the
         defect is tightly bracketed, and the event to instrument is the one measured since run 7 and never
         explained: **m_stag falls 0.366 → 0.0046 in ONE Fock build** from a provably staggered seed, before
         any mixing history exists.  At α=0.25 on a linear m channel that implies an OUTPUT density
         anti-staggered at m≈−1.08 — a sign-flipped response ~3× the input, not a decay.  CP2K's first step
         evidently does not do that.  Cheapest tests first:
         (a) **the polarized XC response at the seed density** — is v_xc^↑−v_xc^↓ the right sign and size?
             A single-point check on the seed, no SCF at all, and it is the one candidate that would explain
             a SIGN FLIP rather than a decay;
         (b) frontier coverage under smearing (CP2K's `ADDED_MOS 20`);
         (c) the mixer (Broyden vs Kerker+Pulay) — weakest candidate, since run 8's linear-mixing control
             already kept m_stag alive;
         (d) seed magnitude: ours is a per-site flip of the SAD d occupation, CP2K's `&BS` is the same idea —
             compare the two initial moments directly.
         **The acceptance test is unchanged (aufbau-stable WITHOUT MOM); only the mechanism changes — this is
         a bug hunt against a same-recipe oracle, not a functional upgrade.**
       * **WHAT IS STILL OPEN** (do not read the above as "done"):
         (a) **It does not pass the convergence gate**: the run ends on the GDM leg with lastΔρ=3.2e-2 against
             MinΔρ=1e-5.  E and the moment are stable to ~1e-3, but the gate measures ρ, and the DIIS→GDM
             hand-off resets what Δρ means.  A gate-passing MnO run needs either a longer DIIS-only leg or a
             convergence measure the hand-off does not invalidate.
         (b) **It is still NON-AUFBAU**: a 0.97 Ha hole persists (empty ε=−0.364 below occupied ε=+0.604).
             A held configuration that never becomes self-consistently aufbau is an excited state; the
             remaining 1.33 Ha may well be exactly this.  The question for the next session is whether the
             hole closes as the moment saturates, or whether Λ is still holding something it should release.
         (c) `cfg *` still flags most iterations, though that is expected under MOM (the occupied SET can be
             stable while its energy ORDER re-shuffles) — the `cfg` flag keys off orbital index, so under MOM
             it is no longer a hopping diagnostic.  A character-keyed configuration string would be.
         (d) The α×Λ surface has two measured points at Λ=1.5 (α=0.45 diverges, α=0.25 lands) and the kT and
             ANNEAL knobs are now wired but UNSWEPT — `MNO_ANNEAL` runs a descending schedule with density
             AND (under `momFromSeed`) occupied-character continuation across stages.
       * **A DISPLAY DEFECT THAT HID EXACTLY THIS (fixed 2026-08-08).**  `tPolarizedWF::DisplayEigen` dropped
         every level empty in BOTH channels, wherever it sat — so the −1.29 Ha empty level was invisible in a
         table that happily printed a +0.75 Ha virtual.  Worse, the rule was silently INCONSISTENT with kT:
         under Fermi smearing no occupation is exactly 0.0, so nothing was ever dropped and the smeared runs
         looked complete; at kT=0 the deep empties vanished.  Same table, same run, visibility depending on a
         convergence knob.  Now: a doubly-empty level is dropped only ABOVE the frontier, never below it (a
         hole is the whole point of looking).  The `tUnPolarizedWF` twin had the same idiom as a `break`,
         which truncated the table AT the anomaly; fixed the same way.  Also: the polarized table printed
         `setprecision(0)` occupations, rounding a smeared 0.996 to "1/1" — it advertised a clean integer
         configuration for a run whose own trace column was flagging partial occupancy every iteration.  It
         now shows fractional occupations, like the unpolarized twin already did.
       * **THE INSTRUMENT (general, landed with the diagnosis)**: `tSCFIterator::SetOrderParameter(name,
         probe)` — a caller-supplied named scalar measured on the WORKING density every iteration, rendered
         as an extra trace column, carried in `SCFProgress`, and printed once at "iteration 0 (seed)" as the
         reference every later value is judged against.  This is §9's "diagnostic metric" open question,
         answered by the concrete case: an explicit, caller-supplied ORDER PARAMETER, because which sites
         and which sign pattern is system knowledge the iterator cannot infer.  Zero cost when unset.  The
         GPW driver adds the post-mortem line (whole trajectory + the iteration at which the order died, at
         the 1%-of-seed threshold); MnO's probe is m_stag = ½[m(Mn1)−m(Mn2)] sampled 0.7 bohr off each Mn —
         and since the run fixes nUp=nDn, ferromagnetic drift is not an available escape, so m_stag is the
         WHOLE order.  Reusable verbatim for LiMn₂O₄ charge/spin ordering.
       * **UNDER-BOUND by 15.8 Ha** vs CP2K −61.4706 — far beyond basis incompleteness (our own free
         atoms sum to ≈−60.8, so a bound crystal must sit BELOW that).  CP2K's CRYSTAL decomposition is
         now banked for the comparison (same deck): `core-self −125.2197  coreH +42.3595  Hartree
         +35.7299  XC −14.3400  overlap ~0  => −61.4706`, against ours `Ekin 78.73  Een −84.45
         Eee 29.14  Exc −14.88  Enn −57.34  alphaZ +3.17 => −45.64`.  **CAUTION: the two PARTITIONINGS
         DIFFER** — CP2K folds the local PP's long-range part into its Hartree term and carries a
         Gaussian core SELF-energy (−125.22) that is cancelled inside that Hartree, whereas we split
         Enn (point-Zion Ewald) + E_alphaZ (G=0 alignment) + Een + Eee.  So the totals are the only
         apples-to-apples number until someone derives the mapping; the density-INDEPENDENT constants
         (ours Enn+alphaZ = −54.17) are where that derivation should start, since both KB routes are
         oracle-clean at the atom level and XC agrees to ~0.5 Ha (itself partly the ζ=0 collapse).
     Run logs: doc/logs/mno_afm2_run6_curedbasis.log (the settled run);
     `run7_orderparam_kerker` (m_stag dies at iteration 1 under the ρ̃ mixer) and
     `run8_linearmix_control` (m_stag alive under linear D-mixing) — the A/B that convicted the mixer.
     (b) **A real but unattributed l=2 discrepancy between the two CRYSTAL KB routes**: analytic vs
     mesh = 3.09e-2 relative on Mn (vs 2.4e-9 for the Si l=0,1 gate) — structural, not a scale factor
     (max|Va| == max|Vm| exactly).  Could still be the MESH arm being coarse on Mn's compact
     r_l=0.328 d projector; a densityEcut sweep attributes it.  3% cannot explain 356 Ha, so this is
     a separate (smaller) issue.  Gate `GPW.DISABLED_AnalyticSeparablePPMatchesMesh_DChannel`.
   the first genuine workout of the Shubnikov ops + imposed-subgroup policy + `+U` +
   order-parameter diagnostic together. The direct rehearsal for **LiMn₂O₄** charge/spin
   ordering (the north-star), which follows.  (Tier 4b DONE unblocks the two-channel
   machinery; the AFM-specific prerequisite is the SPIN-POLARIZED SAD seed with per-site
   moments — a spin-agnostic seed cannot express a staggered pattern, so the seed chooses
   the magnetic-ordering basin.  Design pinned in doc/SCFSeedingPlan.md §10: tables store
   the Hund pair in the UP-MAJORITY convention; assembly applies per-site flips (collinear
   configs) or SU(2) rotations (non-collinear) — the config is assembly-time data, never
   duplicated tables.  Also budget for `+U` and Fermi smearing on the d manifold.)
8. **Crystal SALCs (Tier B)** — H block-diagonalization per k; separate track
   (`SpaceGroupPlan` Increment 3), reuses `BuildSALCs`.

Steps 1–3 are prerequisites/scaffolding; 4–6 are the quadrature/compute wins; 7–8 are the
physics payoff.

---

## 8. Correctness harness (the invariant for every reduction)

For each of T1/T2/T3, two tests mirroring the existing IBZ pattern, with TWO tolerance
tiers (review fix — one number would make correct code "fail"):
- **reduced == full**, tiered:
  - **~1e-13 (reordering class)** for SINGLE-SHOT folds on a fixed symmetric input — T1's
    G-sum bookkeeping, a one-shot T3 collocate/integrate pair on a frozen symmetric D;
  - **grid/SCF class (~1e-7; the measured IBZ precedent is ~6e-8)** for anything taken
    THROUGH SCF — the converged density's symmetry defect is bounded by the convergence
    tolerance, and the through-SCF comparison inherits it (still needs the commensurate
    grid, §5).
- **negative control**: on a symmetry-broken density, reduced ≠ full, *and* the §3 audit
  fires — the FREE-run diagnostic reports the lowering, and the imposed-run RELEASE-CHECK
  recovers the broken (lower) solution — proving the reduction genuinely imposes symmetry
  and never does so silently.

---

## 9. Open questions / risks

- **Route (a) vs (b) for streams** — (b) is cheaper per iteration but imposes symmetry;
  (a) is exact-without-imposing but is a MEMORY+BUILD lever only (the cache is replayed,
  not rebuilt, per iteration — and (a)'s permuted gather likely makes the replay slower).
  Likely: (b) under the asserted-imposition policy, (a) as the no-impose fallback and the
  release-audit's ground truth.
- **The variational adjoint pin under route (b)** — today `∫ρV == Tr(Dh)` is machine-exact
  because collocate and integrate share the identical operator.  Route (b) inserts the
  star-average `P` on the density side and the representation transform on the h side; the
  pairing stays exact only if the h-side map is exactly `Pᵀ` composed with the SAME shared
  streams (`P` is self-adjoint in the grid inner product, so this is achievable — but it
  must be PROVED in the design, not discovered in a failing virial).  The D-aware
  active-set screen is part of the same proof: the representative's kill-weight must
  aggregate the whole orbit's `|Re[D_ij e^{-ik·R}]|` consistently on BOTH directions, or
  the shared-active-set invariant silently breaks.
- **Commensurability cost** — rounding every ladder level to 5-smooth multiples of
  lcm(τ-denominators) is modest for τ=¼ (÷4) but grows for screw axes (τ=c/6); measure
  before committing the stream fold on low-symmetry-denominator cells.
- **Diagnostic metric** — `‖D − P_G D‖` vs per-op invariance defect vs an explicit order
  parameter; pick one that is cheap per iteration and interpretable in the run report.
  **ANSWERED 2026-08-07 (by the MnO AFM collapse, §7 step 7): the EXPLICIT ORDER PARAMETER**, supplied by
  the caller — `tSCFIterator::SetOrderParameter(name, probe)`, one named scalar on the working density per
  iteration, as a trace column + `SCFProgress` field, with the seed's value printed as iteration 0.  Which
  sites and which sign pattern define the order is system knowledge the iterator cannot infer, so the probe
  is a functor, not a policy.  It paid for itself immediately: two point evaluations per iteration bracketed
  a "the SCF relaxes out of the magnetic basin" mystery to a single step and a spin-blind mixer.  The
  invariance-defect metrics remain the right instrument for the *symmetry* (§3) question — they answer
  "which ops broke", not "did the order survive", and the two coexist.
- **Spin op algebra** — the anti-unitary time-reversal part (`ρ̃(−G)=ρ̃(G)*` + spin flip)
  needs care in the G-space phase bookkeeping; settle it in the §4 interface, not later.
- **Non-collinear representation** — 2×2 spinor density (complex coefficients coupling
  up/down) vs (ρ, m) scalar + magnetization 3-vector.  They are equivalent
  (ρ_αβ = ½(ρ δ_αβ + m·σ_αβ)); the choice is about which is natural for XC (most
  non-collinear functionals are collinear functionals applied along the LOCAL m̂(r) axis →
  favors (ρ,|m|)) vs for the ops/densities (SU(2) rotations and k-space bookkeeping →
  favors the spinor).  Decide when 4b's successor (non-collinear pipeline) is scoped; the
  collinear two-channel tier is a strict subset of either.
