# Symmetry Upgrade Plan — PW & GPW SCF

*Plan first, code second.* Scope: exploit lattice + point (and eventually magnetic)
symmetry across **all** the redundant work in a PW/GPW SCF, not just the k-mesh.
Companion to `doc/SpaceGroupPlan.md` (the detector + Tier A/B roadmap),
`doc/GPWPlan1.md` (the "SPACE-GROUP STREAM/COLLOCATION REDUCTION" plan section this
consolidates), `doc/SpinNativeDFTPlan.md` (spin-native XC), and
`doc/SymmetryRefactorPlan.md`.  SOLID/OOD debt encountered while executing this plan is
tracked in **`doc/CleanupCandidates.md`** (keep it growing; batch-fix in dedicated
refactor sessions).

**WHERE WE LEFT OFF (2026-08-11, late).**  §7 steps 1–6 are DONE; step 7 (MnO AFM-II): **THE OPEN
DEFECT IS FOUND AND FIXED.**  Point-in-time status blocks from earlier sessions, and the full MnO
campaign narrative, are in **`doc/SymmetryUpgradeHistory.md`** — this section is the standing summary,
kept short on purpose.

**THE DEFECT (found by the planned v_xc probe, first try): the SEED's real-space ρ(r) had NO LATTICE
IMAGES.**  `SeedCD::operator()(r)` summed the per-atom radials over HOME-CELL atom positions only,
while the XC mesh samples the seed at points WRAPPED into the home cell (`kpt = r − A·n0`,
`MakePeriodicBeckeMesh`; the uniform raster's corner regions likewise).  Where the nearest atom is a
lattice IMAGE the raster read ~0 — and for the CORNER Mn, 7 of its density hump's 8 octants belong to
images.  Under the AFM staggering the erased hump is the MAJORITY channel of exactly one sublattice:
site-specific, rigid-translation-sensitive (`MNO_SHIFT`), uniform-raster-reproducing, spin-swapping —
every recorded symptom, including the "impossible" site-dependence of a pointwise LSDA v_xc (the
functional was fine; its ρ INPUT was site-wrong) and the positive iteration-1 Etot in every MnO log.
The earlier "seed exonerated" verdict and this defect were simultaneously true because the seed-gate
probe points all sat within ~2 bohr of a HOME atom, where an image-less sum is accurate.

**THE FIX + GATE:** `SeedCD` (and so `PolarizedSeedCD`) now evaluates the periodic image sum — FINITE
and EXACT because every radial table has compact support (`RadialDensity::Rmax`, ρ≡0 beyond it): the
only truncation is the table's own support, no distance cut (the no-cut pin).  Gate
`GPW_SCF.MnOSeedVxcMirrorOnBeckeMesh` tests ρ and v_xc mirrors POINT-BY-POINT on the actual Becke
raster via the wrapped-fractional partner map: before the fix 2420/3878 points broken, max v_xc defect
1.08 Ha at a point 0.58 bohr from an Mn1 IMAGE; after, 2e−15 / 6e−15.  **PIN: a periodic density's
real-space face must be correct at ANY point, because mesh consumers sample it WRAPPED — hand-picked
probe points near home atoms cannot exonerate it.**

**MEASURED POST-FIX (3-iteration bounded run, DEFAULT knobs — α=0.45, kT=5e-3, MOM still inert):**
iteration-1 Etot **−52.80** (was +7.373 in every log ever); E: −52.80 → −59.69 → **−61.394** = 0.077 Ha
from CP2K's −61.4706 after THREE iterations (the whole prior campaign's best: −61.370 at run 29 with
MOM+Λ scaffolding); m1/m2 = +0.43/−0.40 → +0.62/−0.58, m_stag GROWING (0.41→0.32→0.60), mirror held
(m_net=0.038), N↑=N↓ exact; still descending at the cap (Δρ=1.3e-3), no oscillation flips.  The free-run
§3 diagnostic reports exactly 6/12 ops broken — the grey ops that map +m→−m, i.e. the instrument seeing
the AFM order as designed.  Log: `doc/logs/mno_afm2_run30_postfix_bounded.log`.

**NEXT ACTION**: a full free run to convergence on the post-fix ground: defaults first (does it now
converge and pass the gate UNAIDED — no MOM, no anneal?), then re-derive which of the run-29 recipe
knobs (α, kT, MOM/Λ) are still load-bearing.  The occupation two-cycle and the 4–6 Ha branch gap were
measured AGAINST the broken v_xc — re-measure before believing either still exists.

**LANDED 2026-08-11, all green:** the JOINT Pulay history (one B summed over both spin channels, one
coefficient vector — the run-27 ejections are GONE: run 29 best −61.370, m_stag 0.50–0.66); ONE μ per
electron RESERVOIR (two μ were the Lagrange multipliers of a constraint nobody asked for — μ↑−μ↓ was
27 mHa); the seed fill no longer runs a μ-solve at kT=0 (it was losing charge on EVERY global-μ metal:
Al Σw·n = 2.25 against 3).  Gates: `JointPulay.*`, `SharedFermiLevelLetsTheMomentRelax`,
`MnOSeedSublatticesAreEqualAndOpposite`.

**INSTRUMENTS ADDED — and why they were missing:** `GPW_MNO_SITES=1` prints m1 and m2 SEPARATELY from
iteration 0, because `m_stag = ½(m1−m2)` reads +0.366 for both (+0.37,−0.37) and (0,−0.73) — the order
parameter was blind to the exact failure it existed to watch.  `svMin` is now a trace column (DIIS
conditioning had never been printed).  `GPW_BECKE_ATOMS` gives the per-atom mesh breakdown.  Also fixed:
`DisplayEigen` was FABRICATING a channel's energy when a level existed in only one channel (printing
ϵ↑−ϵ↓ = 0.00000000 exactly), which manufactured MnO's spin-up "hole".

**RUN 31 (2026-08-11, post-fix free run at defaults) — DIIS solves it; the GDM rung was the saboteur,
now structurally fixed.**  The DIIS leg descended monotonically to **−61.4138 by iteration 12** (past
run 29's scaffolded best, 0.057 Ha from CP2K) — then the ladder's tail hand-off advanced onto GDM,
which correctly DECLINED (Fermi-smeared D' is outside its integer-occupation manifold:
Tr(D'P_block)=12.80 of 13), and the loop degraded to bare damped-Kerker steps still labelled GDM — the
unstable iteration DIIS had been taming — un-converging the run to a −60.49 limit cycle ([F,D]=0.35).
USER CALL (live, watching the log): "for GDM we should not be seeing any Ker or relax."  FIX: the new
STANDING-PRECONDITION query `tSCFAccelerator::Engageable()` — GDM answers from D' idempotency
(smearing; knowable on every UseFD with no orbitals) AND the leading-block check (MOM); the LADDER now
VETOES any hand-off onto a non-engageable rung and RETREATS off one whose precondition fails while
active.  `Ker … GDM` rows are now structurally impossible: the ladder never sits on a rung it cannot
drive.  Gates: `SCFAcceleratorGDM.ASmearedDensityIsNotEngageable`,
`SCFAcceleratorLadder.HandOffIsVetoedAndRetreatsOnANonEngageableRung`.  Log:
`doc/logs/mno_afm2_run31_postfix_free.log`; run 32 (the fix, same defaults) =
`doc/logs/mno_afm2_run32_ladder_engageable.log`.

**RUN 32/33 (2026-08-11) — the SECOND saboteur: a DEPENDENT DIIS HISTORY the ABSOLUTE SVTol cannot
see.**  Run 32 (veto active, pure DIIS) hit the SAME −60.49 limit cycle at iteration ~27: svMin sat
PINNED at 1.0–1.1e-9 — near-singular relative to B's own |[F,D]|² ≈ 2.5e-5 scale, "just slightly above
the SVTol limit" (user, live) — while Nproj stacked to 8; a wild depth-8 extrapolation then amplified
the charge slosh.  This is the 2026-08-10 svMin-column prediction biting in production.  FIX:
`DIISParams::SVTolRel` — prune the oldest entry while svMin < SVTolRel·max_i B_ii (a test of
DEPENDENCE, not smallness); default OFF (molecular door unchanged), the SOLID door sets 1e-4.  Gate
`SCFAcceleratorDIIS.ARelativeConditioningPrunePurgesADependentHistory` (its control arm reproduces the
MnO shape).  RUN 33 (both fixes): NO collapse through the run-32 failure window — Nproj rides 3–7,
svMin healthy, E −61.4138±1 mHa, m_stag 0.648, gap 0.18 Ha open.  BUT Δρ wobbles at 0.8–2e-4 with cfg
flips: at kT=5e-3 the fractional frontier keeps re-shuffling, and that noise floor is outside a linear
Fock extrapolator's model — DIIS is at ITS floor, not the run's.  Logs:
`doc/logs/mno_afm2_run{32_ladder_engageable,33_svtolrel}.log`.

**RUN 34 (2026-08-11, `MNO_ANNEAL="5e-3,0"`) — kT=0 IS NOT REACHABLE BY A FILL: the frontier is a
genuinely TIED manifold.**  Stage 1 (kT=5e-3) reproduced the healthy plateau (A=−61.4140, m_stag
0.6482, 30 iters).  Stage 2's very FIRST kT=0 fill of the SAME orbitals jumped E to −60.63 and
[F,D] to 0.53 — the smeared state's occupation is FRACTIONAL within a tied frontier set (the 'm' flag
all through stage 1, beside an open 0.18 Ha gap), and integer aufbau must break the tie.  Every
resolution is bad: the run oscillates between an AUFBAU branch at −60.38 and a NON-AUFBAU ('h', hole)
branch at −61.16..−61.30 — i.e. at kT=0 a held non-aufbau configuration beats aufbau by ~0.9 Ha, and
neither is the ensemble.  DIIS never engaged ([F,D] 0.53 > FDMax) and GDM was never reached — the
machinery reported honestly at every step; the STATE is the problem.  m_stag survived throughout
(0.65).  VERDICT: at Γ-only sampling the kT=5e-3 fractional ENSEMBLE is the converged description of
this cell; the tie is (at least partly) a k-SAMPLING artifact — dispersion would split the tied
manifold.  Log: `doc/logs/mno_afm2_run34_anneal_kt0.log`.

**DECISION (user, 2026-08-11 late): the next semi-major step is SHUBNIKOV + imposeSymmetry** — impose
the MAGNETIC group (which preserves the staggering; the grey group would erase it), with the
imposed-SUBGROUP option "the ops common to two candidate orderings" so the SCF still chooses (§3).
Rationale: the per-iteration star-average under the magnetic group projects the frontier TIE-NOISE
out of ρ (the same mechanism that cures the Becke × degenerate-open-shell oscillation), so this is
the one route with a claim on making the Δρ gate REACHABLE — and it unlocks the imposed-run cost
levers (site-adapted mesh, T1, T3 folds) on the magnetic cell.  The free-run ensemble (runs 30–34)
is the §8 A/B reference.  Increments: S1 group construction (qcSymmetry, pure) → S2 channel-pair
star-average (σ=Flip = spatial op ∘ ↑↓ swap in SymmetrizeGMap/Raster; real densities, no conjugation
at collinear tier) → S3 mesh (spatial parts only — for AFM-II they span the grey group, the
site-adapted builder is unchanged) → S4 gates (imposed-M == free ensemble; grey-group imposition
must ERASE m_stag as the negative control; release-audit; the Δρ-floor measurement).

**S1 DONE (2026-08-11):** `AtomSite.spin` (collinear ±1/0 decoration; 0 = non-magnetic — which
species are magnetic is SEED knowledge, the caller decorates), `SpaceGroup::ShubnikovOps(decorated)`
→ `SymOp{W|τ,σ}` (the §4-tier-4a σ slot, now populated), and `CommonOps({M₁,M₂,…})` (same (W,τ,σ)
TRIPLE in every candidate — an op that is Flip for AFM-II and None for FM is NOT common).  KEY
MECHANISM: `Detect` keeps ONE τ-coset per W (primitive-cell assumption), but the magnetically
DOUBLED cell has TWO chemical-lattice cosets, and the extra one is the ANTI-TRANSLATION
{E|½½½}·Flip = the m1=−m2 sublattice mirror itself — so `ShubnikovOps` re-enumerates every valid τ
per detected W against the decorated basis (Detect and all grey consumers untouched).  Gates
(UTSymmetry `Shubnikov.*`, 5): MnO AFM-II 12 grey → 24 Shubnikov = 12 None + 12 Flip incl. the
anti-translation; group CLOSURE over all 576 products; undecorated/FM ⇒ all None; CommonOps(AFM,FM)
= the 12 sublattice-preserving ops; an incompatible (ferri) decoration drops ops.

**S2 DONE (2026-08-11): the channel-pair star-average primitives.**  THE DIAGONALIZATION: the pair
(ρ↑,ρ↓) splits into the TOTAL (EVEN under Flip — the plain average over ALL of M's spatial parts,
both chemical-lattice cosets included: the Hartree-side contract) and the MAGNETIZATION m = ρ↑−ρ↓
(ODD — character χ(op) = −1 on σ=Flip ops).  Landed: (i) `SymmetrizeValuesSigned(Fold, sigmas, vals)`
— the χ-weighted orbit mean for the Becke-mesh {r} path — plus the ODD-FIELD AUDIT
`FlipFixedPointsPeriodic`: a point some Flip op maps onto ITSELF carries m≡0 exactly (MnO: the O
sites, fixed by the sublattice-swapping inversion — m(O)=0 IS the physics) and the single-edge orbit
mean cannot see that, so callers zero flagged points (geometry-fixed, once per mesh); (ii)
`SymmetrizeGMap(rg, sops, oddUnderFlip)` — the G-space scatter with χ — EXACT projector with no
audit needed (the scatter accumulates ALL ops, forbidden components average to zero by themselves);
(iii) `MagneticSymmetryDefects(up, dn, sops)` — the S4 free-run instrument: σ=Flip ops compare
ACROSS channels (the anti-translation row IS the m1=−m2 mirror), None ops within; "which magnetic
ops broke" op by op; (iv) `ReciprocalOf(SymOp)` (U=Wᵀ, τ, σ — the direct→G-frame converter).
GATES: `Shubnikov.SignedFoldFixesTheStaggeredPatternAndGreyErasesIt` (the staggered pattern is a
FIXED POINT of the signed projector, a perturbed one projects onto the mirror ±0.675, the plain grey
average ERASES it to 0.025 — the unit-level form of why ShubnikovOps exists; the audit flags exactly
the O sites); `GPW.ShubnikovGMapOddFixesStaggeredMagnetizationAndGreyKillsIt` (analytic point-moment
structure factors: odd average = identity on the staggered m̃, even average kills it; the AFM total
is even and fixed); `GPW.MagneticSymmetryDefectsSeparateMirrorKeepersFromBreakers` (exact pair: all
24 ops clean; a shrunk sublattice moment — the run-31 disease — fires EXACTLY the 12 Flip ops).

**STILL OPEN, ranked:** (1) S3 — the wiring: WHO DECORATES (flip bits live on Atom, magnetic-or-not
is seed/library knowledge — assemble where the seed resolves species, thread to the factory's policy
resolution); the imposed-run plumbing (XCQuadrature/engine carry the σ list beside the fold + the
zero-flags; the composite's ρ̃ path takes the coset-complete spatial ops; RhoPol does the
(ρ,m)-split symmetrize-recombine); the uniform-raster pair overload on the FIT face rides along
(σ-carrying ops replace DirectOps in the fit-basis ctor).  Then S4 gates per the sketch.  (2) cheap
parallel arms for the 57 mHa: `MNO_SHARED_MU=1` + gentle anneal — user already probing shared-μ;
(3) the k-MESH free run (affordable after imposition); (4) free-energy-stationarity convergence gate
if imposition alone does not close the Δρ floor; (5) doc/CleanupCandidates.md D9–D12.
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
   **CAMPAIGN 2026-08-04 → 2026-08-11 — the full narrative, every run and every refuted hypothesis, is in
   `doc/SymmetryUpgradeHistory.md` §B.**  Current state, the open defect and the next action are in the
   WHERE WE LEFT OFF section at the top of this document.  The durable pins the campaign produced:
   - **the ρ̃ mixers must be spin-native** — a single-map mixer collapses v_xc to the ζ=0 branch from
     iteration 1 (fixed a228218d; negative control `GPW_SCF.PolarizedRunKeepsItsSpin`);
   - **a HISTORY may never be split per channel** — a filter is channel-diagonal and may be, an
     extrapolator is not: independent fits synthesise a moment that never occurred on the trajectory.
     One B summed over channels, one c (2026-08-11).  Same bug CLASS as the spin-blind filter, one level up;
   - **one μ per electron RESERVOIR** — two μ mean you have specified (nUp,nDn) and the moment is a
     CONSTRAINT; one μ means you specified N and the spins fall where they may.  Anything else is
     unphysical, and the constraint is never needed to FIND a ground state (the Mn S=5/2 atom reproduces
     the Hund fill exactly from N alone);
   - **an order parameter must not be blind to its own failure mode** — ½(m1−m2) cannot tell (+m,−m) from
     (0,−2m); report the site moments SEPARATELY;
   - **CP2K constrains nothing** on this deck (Fermi smearing, one μ, `&BS` seeds the ATOMIC guess only),
     so a constrained run is not oracle-comparable.

   **Why MnO is the right test material:** the first genuine workout of the Shubnikov ops + imposed-subgroup policy + `+U` +
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
