# Symmetry Upgrade Plan — PW & GPW SCF

*Plan first, code second.* Scope: exploit lattice + point (and eventually magnetic)
symmetry across **all** the redundant work in a PW/GPW SCF, not just the k-mesh.
Companion to `doc/SpaceGroupPlan.md` (the detector + Tier A/B roadmap),
`doc/GPWPlan1.md` (the "SPACE-GROUP STREAM/COLLOCATION REDUCTION" plan section this
consolidates), `doc/SpinNativeDFTPlan.md` (spin-native XC), and
`doc/SymmetryRefactorPlan.md`.

**STATUS (2026-08-01):** §7 **steps 1+2 BUILT** — `qchem.Symmetry.Lattice_3D.Fold`
(`FoldPoints`/`FoldGrid`/`FoldGVectors`, per-member op edges, `SymOp {W|τ,σ}` with the
tier-4a `SpinAction` enum), the `SpaceGroup::FoldKMesh/FoldGVectors/FoldPoints` member
wrappers, and `ReduceToIBZ` re-expressed on `FoldGrid` (bit-identical, gated by
`Fold.KMeshFoldBitIdenticalToReduceToIBZ`).  Unit tests: `src/Symmetry/tests/L_Fold.C`.
§7 **step 3 BUILT** — `SymmetryPolicy` (impose-on-assert; `GPWParams.reduceBZ` documented
as the assert-switch, resolved ONCE in `DetectPointOps` which now always detects and feeds
policy + all op faces from one place), `SymmetryDefects` per-op order-parameter diagnostic
(GMap, G-space so commensurability-free), and the never-silent `[symmetry]` run line
(free run: max defect + which-ops-broke; imposed run: the release-audit assertion).
Gates: Al symmorphic + Si non-symmorphic IBZ still exact; unit negative control fires on a
broken ρ̃ (stabilizer ops stay clean = the which-op readout) AND demonstrates the
projector-blindness (post-`SymmetrizeGMap` defect ≡ 0 — §3's imposed-run caveat, measured).
MEASURED calibration: a converged (Δρ≈2.5e-6) free Si Γ run shows max defect ~1e-5 — the
through-SCF defect tier tracks the CONVERGENCE tolerance (per §8), so the run-line BROKEN
alarm sits at 1e-3.  Full suite 637/637.  **Next = §7 step 4 (T1 {G} fold)**; the SCF-level
broken-seed negative control lands with the §8 harness (needs a symmetry-breaking
SeedStrategy).  NOTE (§4/§9): non-collinear added — `SpinAction{None,Flip}` is documented
as the collinear collapse of the general spin rotation (SO(3)/SU(2)); representation choice
(2×2 spinor vs (ρ,m)) is an open §9 question.

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

---

## 7. Sequencing

Ordered so each step is gated by a cheap, already-shaped test, and no interface is drawn
against a special case it will outgrow:

1. **Spin-native interface SHAPE** (§4 tier 4a) — signatures + op spin-action enum only;
   gated by compilation + unit-level fold tests.  *Prevents a redo; blocks nothing.*
   (§4 tier 4b — the polarized solid pipeline with gates (a) Na doublet / (b) O₂ triplet —
   is `SpinNativeDFTPlan` work that runs IN PARALLEL and is a prerequisite only for step 7.)
2. **`Fold` primitive in qcSymmetry** (§2b) — generalize `ReduceToIBZ` → `FoldPoints`/
   `FoldGrid` with per-member op index; re-express the k-fold on it (bit-identical). Unit
   tests only.
3. **Imposed-group policy object + order-parameter diagnostic** (§3) — run-level policy;
   diagnostic reports symmetry lowering. Gated by: symmorphic Al IBZ still exact; a
   deliberately symmetry-broken seed makes the diagnostic fire.
4. **(T1) {G} fold** — smallest real reduction; reduced==full on a symmetric density.
5. **(T2) {r} + Becke angular fold** with the commensurability menu (§5); folds in the
   site-group-adapted angular grid + the Becke+IBZ verification (currently the uniform-raster
   carve-out under `reduceBZ`).
6. **(T3) stream cache reduction** — the 5–20× lever; needs §5's commensurable raster menu
   on every ladder level. Route (b) default (imposes, gated by the policy + diagnostic),
   route (a) as the no-impose fallback.
7. **MnO rocksalt AFM-II** (2 f.u., moments along [111]) — first *real* d-electron magnet:
   the first genuine workout of the Shubnikov ops + imposed-subgroup policy + `+U` +
   order-parameter diagnostic together. The direct rehearsal for **LiMn₂O₄** charge/spin
   ordering (the north-star), which follows.
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
- **Spin op algebra** — the anti-unitary time-reversal part (`ρ̃(−G)=ρ̃(G)*` + spin flip)
  needs care in the G-space phase bookkeeping; settle it in the §4 interface, not later.
- **Non-collinear representation** — 2×2 spinor density (complex coefficients coupling
  up/down) vs (ρ, m) scalar + magnetization 3-vector.  They are equivalent
  (ρ_αβ = ½(ρ δ_αβ + m·σ_αβ)); the choice is about which is natural for XC (most
  non-collinear functionals are collinear functionals applied along the LOCAL m̂(r) axis →
  favors (ρ,|m|)) vs for the ops/densities (SU(2) rotations and k-space bookkeeping →
  favors the spinor).  Decide when 4b's successor (non-collinear pipeline) is scoped; the
  collinear two-channel tier is a strict subset of either.
