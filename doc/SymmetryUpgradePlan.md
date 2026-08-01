# Symmetry Upgrade Plan — PW & GPW SCF

*Plan first, code second.* Scope: exploit lattice + point (and eventually magnetic)
symmetry across **all** the redundant work in a PW/GPW SCF, not just the k-mesh.
Companion to `doc/SpaceGroupPlan.md` (the detector + Tier A/B roadmap),
`doc/GPWPlan1.md` (the "SPACE-GROUP STREAM/COLLOCATION REDUCTION" plan section this
consolidates), `doc/SpinNativeDFTPlan.md` (spin-native XC), and
`doc/SymmetryRefactorPlan.md`.

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
| **{G} sum** (structure factor, `⟨G\|V\|G'⟩`, Hartree G-loop) | ❌ full {G} summed | fold by \|G\|-shell × point group | ≈ \|point group\| on the G-loops | integrand totally symmetric |
| **{r} quadrature** (uniform + Becke XC mesh) | ❌ full mesh integrated | fold by star; **Becke angular** = the site-group-adapted grid item | ≈ site-group order per atom | ρ symmetric on the mesh + raster commensurability |
| **GPW stream cache** (`(pair, offset R)` → raster pts) | ⚠️ only Hermitian `j≥i` ×2 | site-symmetry of `(pair, R)`; ~513M(64b)+ pts | **5–20×** (up to ×48 in O_h) — the headline | commensurate raster on **every ladder level** |
| **H block-diagonalization (crystal SALCs)** | ❌ dense H per k | irrep blocks per k | matrix-size³ per block | Tier B (SALC), separate track |
| **spin / magnetic group** | ❌ unpol only (`PW_XC*`) | Shubnikov ops (spatial ⋉ spin-flip) | AFM/ordering cells | spin-native XC interface first |

**Reading of the gap:** the k-side is done; the *quadrature and collocation* side is
almost entirely unexploited. The stream cache is where the compute actually is (the
`(pair,R)` lists are the ~513M cached points, rebuilt identically every SCF iteration), so
it is the highest-value target — but it is also the one with the hardest precondition
(pointwise index permutation needs a τ-commensurate raster on every ladder level, §5).

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
Mesh-typed convenience is a **client-side helper**. Two decisions, resolved with the user:

- **`Mesh` becomes the shared representation for {G} balls and {r} grids** (today {G} is a
  bare `vector<ivec3_t>` and the raster is a bare `rvec_t`). They are all *weighted point
  sets* — exactly what `qcMesh::Mesh` is (SoA `Points()`+`Weights()`, a POD value type). So
  {k}/{G}/{r}/rasters unify onto one type, and the currently-separate `KMesh`/`IBZMesh`/
  raster shapes collapse.
- **The fold primitive still speaks primitive vectors and returns a partition** — kept in
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

// (1) The free CORE — takes an EXPLICIT op set. Kept because tests fold under hand-made ops,
//     and the imposed-SUBGROUP policy (§3) folds under an arbitrary op SUBSET, not "all ops".
Fold FoldPoints(const std::vector<rvec3_t>& pts,
                const std::vector<Matrix3D<double>>& ops, double tol);   // tolerance path ({r})
Fold FoldGrid  (ivec3_t N, rvec3_t shift,
                const std::vector<Matrix3D<double>>& ops);               // EXACT integer path ({G},{k})

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
     (nowhere near 2⁵³) and `lround` on an exact-integer double is exact. Templating `Mesh`
     on the point type is a large blast radius (widely-used value type) for a type-clarity
     nicety exact-integer doubles already deliver; if the index type is ever wanted visible,
     a thin `IntPointView` over the double storage beats templating the class.
   - **Store the integer INDICES `m` in the {G} Mesh**, not the Cartesian `G = B·m` (derive
     that only where the kernel `4π/|G|²` needs it). Storing Cartesian would be the trap —
     `R_cart·G` is genuinely real and forces the tolerance path onto an exact problem.
2. **The primitive returns the partition, `Symmetrize(Mesh)→Mesh` is the wrapper.**
   Consumers need *which op* maps rep→member — for the glide phase (G), the ±monomial sign
   (streams), the diagnostic (`D`). So `Fold` carries `members`-with-`opIndex`;
   `Symmetrize(Mesh)` is the thin wrapper that applies the fold and drops the op-index,
   correct only for a pure quadrature fold (the XC {r} mesh). G-space and streams take the `Fold`.
3. **Mesh weight semantics generalize** (user): the `Weights()` array — originally quadrature
   weights — doubles as **phase factors or phase × quadrature-weight**. That is the natural
   vehicle for the fold's reweighting (star multiplicity) *and* the glide phase in one place.
   Caveat: the glide phase is **complex**, but `qcMesh::Mesh` weights are real (`rvec_t`); a
   complex-weight Mesh (or a paired real-qweight + complex-phase) is needed for the G-space
   case — a clean fit with the project's real/complex-parity bias (everything available Pol
   *and* real/complex, not a bolted-on special case).
4. **A stream is more than a Mesh.** Mesh-ifying the offset geometry (`R→W·R`) is right, but
   the stream reduction still carries the `(pair, monomial-sign)` layer, which is not mesh
   data — so "streams become a Mesh" covers their geometric part only; the pair/monomial rep
   theory stays with the basis (fed the op).

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
  geometric space group, **not** hard-wired to the lattice's full group. Default = full
  detected space group; user/driver may set a subgroup or the trivial group. One policy
  object feeds IBZ, {G}/{r} folds, stream reduction, and the density star-average alike.
- **Order-parameter diagnostic — never silent.** While imposing `H`, monitor `D`'s symmetry
  under the **full** group `G`: measure the residual `‖D − P_G D‖` (or per-op invariance
  defect) after convergence. If the density is less symmetric than `G`, the run **reports**
  the lowering (which ops broke, the order parameter), never hides it behind a reduced mesh.
  A run always tells the user the symmetry it actually found.
- **Staged converge-then-release.** Converge under full `G` (fast, symmetric), then release
  to a subgroup / trivial group with a symmetry-broken **seed** (`SeedStrategy`) and
  re-converge; if the energy drops, SSB is real and reported. The density-continuation
  machinery for this exists.
- **The two ordering workflows need exactly this** (`[[battery north-star]]`): (i) an
  SSB-search run imposes a subgroup (or trivial) in an ordering-commensurate supercell with
  a broken seed + `+U` (LDA/GGA delocalization error otherwise suppresses disproportionation);
  (ii) the cluster-expansion / MC route imposes each candidate ordering's stabilizer
  subgroup and compares energies — same code, enumerated groups.
- **The trap to design against:** imposing `G` for speed on a system that *wants* to break
  it yields the wrong (too-high-symmetry, too-high-energy) answer *silently*. Therefore the
  safe default for an *unknown* system is to impose `G` only when asserted, or to run the
  release-check; the diagnostic is the backstop that makes any imposition honest.

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
- **No new physics required to fix the shape** — `PW_XC`/`PW_XC_Becke` stay unpol for now;
  the molecular VWN5 correlation is already spin-native (`SpinNativeDFTPlan` B1). This is
  purely: settle the two-channel signature + the op spin-action enum so the review doesn't
  encapsulate around the unpol special case.
- **Interface-shape tests (no new physics):**
  - **(a) Na pseudo-atom in a box, doublet** — the existing Na q1 PP, S=½, moment 1: the
    minimal end-to-end two-channel GPW run (two channels through collocation + XC).
  - **(b) O₂ in a box, triplet** — O q6 GTH, cross-anchored against the molecular facade's
    spin-native triplet-O₂ gate (the spin sibling of `SiPseudoAtomInBoxMatchesFinite`).

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

**(T2) {r} quadrature mesh.** Fold the uniform/Becke mesh points under the site/point group;
the **Becke angular grid** folds by site-group order and merges with the *site-group-adapted
angular grid* item (orbits must avoid special/bond directions — the Lebedev cube-corner
lesson; the site group fixes both how many and which directions, subsuming the GL-29 vs
GL-17 default question). Requires ρ symmetric on the mesh + raster commensurability (§5).
Also the candidate cure for the **Becke × degenerate-open-shell oscillation** (a
site-symmetric quadrature removes the rotating-error channel).

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

1. **Spin-native interface shape** (§4) — signatures + op spin-action enum; gates (a) Na
   doublet, (b) O₂ triplet. *No new physics; prevents a redo.*
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

For each of T1/T2/T3, two tests mirroring the existing IBZ pattern:
- **reduced == full** to ~1e-13 (reordering) for a density that *has* the imposed group, on
  a **commensurate** grid (the tolerance is only as tight as the density's actual symmetry,
  which is why the commensurate grid matters — §5).
- **negative control**: on a symmetry-broken density, reduced ≠ full *and* the
  order-parameter diagnostic fires — proving the reduction genuinely imposes symmetry and
  never does so silently.

---

## 9. Open questions / risks

- **Route (a) vs (b) for streams** — (b) is cheaper per iteration but imposes symmetry;
  (a) is exact-without-imposing but has a cache-hostile gather. Likely: (b) as the default
  under the policy, (a) as the audited fallback and the diagnostic's ground truth.
- **Commensurability cost** — rounding every ladder level to 5-smooth multiples of
  lcm(τ-denominators) is modest for τ=¼ (÷4) but grows for screw axes (τ=c/6); measure
  before committing the stream fold on low-symmetry-denominator cells.
- **Diagnostic metric** — `‖D − P_G D‖` vs per-op invariance defect vs an explicit order
  parameter; pick one that is cheap per iteration and interpretable in the run report.
- **Spin op algebra** — the anti-unitary time-reversal part (`ρ̃(−G)=ρ̃(G)*` + spin flip)
  needs care in the G-space phase bookkeeping; settle it in the §4 interface, not later.
