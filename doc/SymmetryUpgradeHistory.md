# Symmetry upgrade history — the closed record

Split out of `doc/SymmetryUpgradePlan.md` on 2026-08-11, when the working doc passed 1800 lines and the
MnO campaign narrative alone was ~800 of them.

NOTHING WAS TRIMMED, per the `doc/CleanupHistory.md` convention: material is moved here in FULL and the
plan keeps a summary plus a pointer, so every cross-reference ("run 27", "the T2 groundwork", "the
spin-blind mixer") still resolves and the reader can see what happened without the plan carrying it.

The same reason applies here as there, and this campaign proved it twice: the ejections that the joint
Pulay history cured in run 29 were diagnosed from the run-27 trace recorded below, and the "+U is the
next step" conclusion was RETRACTED on evidence that was only available because the original reasoning
had been written down in full.

---

## A. Session-close / status blocks (2026-08-01 → 2026-08-08)

Moved from the plan's preamble.  These are point-in-time snapshots; the plan now carries ONE
"where we left off" section instead.

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


---

## B. The MnO rocksalt AFM-II campaign, in full (2026-08-04 → 2026-08-11)

Moved from §7 step 7.  The plan keeps the step-7 header, a summary of the state, and the durable pins.
Read this for the reasoning, the refuted hypotheses (they are as valuable as the confirmed ones), and
the per-run numbers.

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
       * **★ THE MISSING INGREDIENT IS DENSITY HISTORY — MnO HOLDS ITS MOMENT WITH NO MOM AT ALL (user's
         correction + run 27, 2026-08-10).**  USER, on my listing the mixer as the weakest candidate because
         "run 8's linear-mixing control already kept m_stag alive": *"if linear mixing keeps m_stag alive, and
         Kerker kills it then I would conclude the opposite ... Pulay-Broyden might keep the moment."*  That
         is the correct inference and mine was backwards — a mixer swap that changes whether the order
         survives IMPLICATES the mixing path, it does not exonerate it.
         Acting on it exposed something never tested: **every MnO run to date used `PulayDepth=0`** — a damped
         Kerker step with NO density history — while the CP2K oracle runs `BROYDEN_MIXING` with `NBUFFER 8`.
         (And CP2K's `BETA 1.5` IS its Kerker damping, i.e. it damps HARDER than our G0=1.0 and still keeps
         the moment, which points away from "Kerker suppresses the moment" and squarely at the missing
         history.)  Run 27 — `PulayDepth=8`, `PulayStart=5`, **MOM OFF**, kT=5e-3, α=0.25:
         | | best E | m_stag | MOM |
         |---|---|---|---|
         | run 24 (Kerker, no history) | −60.134 | 0.66 | **REQUIRED** (Λ=1.5) |
         | run 26 stage 2 (MOM released) | −49.53 | 0.0001 (died) | none |
         | **run 27 (Pulay history)** | **−61.29** | **0.378, holds unaided** | **NONE** |
         **0.22 Ha from CP2K, and the moment survives with no constraint** — so MOM was never the cure, it was
         compensating for a mixer with no memory.
       * **WHAT RUN 27 IS NOT: converged.**  Eamp over the last 8 is 12.2 Ha.  The trace is not noise but an
         ATTRACTOR WITH INTERMITTENT EJECTIONS — it sits at E≈−61.25, m≈0.345, [F,D]≈3.2e-2, then is thrown to
         −47…−55 with m collapsing and [F,D]→1.0, then climbs back (iterations 57→59, 67→69).  The signature
         is Pulay history POLLUTION: once the history spans both branches the extrapolation lands on the
         non-magnetic one.  Note the shape is identical to the GDM defect fixed earlier — **a bad extrapolated
         step being COMMITTED rather than rejected** — which suggests the same cure (an actual-E or
         residual-spike check that discards the extrapolation and resets the history) rather than a knob.
         Sweeps running from α (CP2K uses 0.2) and depth.  Log `run27_pulay_nomom`.
       * **★ THE EJECTIONS ARE A SPIN-BLIND *HISTORY*: WE SOLVE PULAY SEPARATELY PER SPIN SECTOR (user,
         2026-08-10).**  USER: *"Should we be storing separate PB history for each spin sector?"* — and
         *"it is already concatenating by spatial irrep ... spin is just another irrep."*  **We are, and it
         is wrong.**  `PulayMixer` owns its own `std::deque<ΔG_Map> itsIns/itsOuts/itsResiduals`, and
         `PolarizedDensityMixer` builds ONE PER CHANNEL, so `PulayDepth=8` runs **two independent B-solves
         and gets c↑ ≠ c↓**.  The distinction the design misses:
         - a **FILTER** (Kerker's G²/(G²+G₀²)) is linear and channel-diagonal, so per-channel application is
           *identical* to joint application — the leaf design was correct for it;
         - an **EXTRAPOLATOR** is not.  Its coefficients come from a least-squares fit over the residual
           history, so independent fits give a mixed state (Σcᵢ↑ρ↑ᵢ, Σcᵢ↓ρ↓ᵢ) that **never occurred on the
           trajectory**.  Each channel still conserves charge (Σc=1 holds per channel) but **the MOMENT
           becomes an arbitrary synthesised combination of history moments** — which is exactly an ejection.
         **This is the same bug CLASS as the spin-blind ρ̃ mixer fixed at 041ddff3, one level up: that one was
         spin-blind in the FILTER, this one is spin-blind in the HISTORY.**  Two corroborations: our own Fock
         DIIS is already JOINT ("one B summed over both spins, one coefficient vector"), so the accelerator
         gets it right and the density mixer does not; and CP2K mixes the joint density with one Broyden
         history and keeps the moment.
         **The design-level defect** is `PolarizedDensityMixer`'s stated virtue — *"pure forwarding: the
         leaves are ordinary mixers … without this class knowing which it holds"*.  That abstraction is valid
         **only for MEMORYLESS leaves**; it treats a stateless filter and a stateful extrapolator as
         interchangeable.  The same objection applies to the `TotalAndMoment` basis (separate ρ and m
         histories would synthesise inconsistent (ρ,m) pairs).
       * **THE FIX, AND THE ONE OBSTACLE.**  Per the user: one extrapolation whose field is the DIRECT SUM
         over channels — spin as just another irrep, exactly as the Fock DIIS already sums B over its irreps.
         `qchem.Math.DIIS` needs no change at all; it already documents "the caller owns the history + builds
         B with ITS OWN inner product", and the joint inner product is simply the sum over channels:
         `B(i,j) = Σ_σ MapInnerRe(res_σ[i], res_σ[j])`, one `c`, applied to every channel.
         **Obstacle**: `PolarizedDensityMixer::Mix` calls `itsUp->Mix()` then `itsDn->Mix()` in sequence, so ↑
         would need `c` before ↓ has contributed.  A joint solve therefore needs residual STAGING separated
         from coefficient APPLICATION — i.e. `tFieldMixer` gains `CarriesHistory()` / `StageResidual() ->
         partial B` / `ApplyJoint(c)`, with memoryless leaves defaulting to the current single-phase path.
         That is the honest shape: the wrapper must know whether its leaves carry MEMORY, because that is
         precisely the property that decides whether splitting is legitimate.
       * **✅ THE JOINT HISTORY LANDED (2026-08-10), in the shape above.**  `src/ChargeDensity/DensityMixer.C`:
         - `tFieldMixer::History()` returns a `tFieldExtrapolator*` — **the capability query, not a bool**:
           the memory-carrying face IS the staging protocol, so a caller cannot ask the question without
           being handed the only legitimate way to act on the answer (`CarriesHistory()` returning a bool
           would have left the two-phase call site to be assembled by hand at every future call site).
           Memoryless filters return `nullptr` and keep the single-phase `MixField`.
         - `tFieldExtrapolator` = `StageResidual(GField) -> {resid, my block of B}` + `ApplyJoint(c)`; the
           step is DEFERRED between the two, which is what lets ↓ contribute before ↑ commits.
         - `MixJointly(channels, fields)` is the whole protocol in nine lines: stage all, sum the blocks,
           ONE `DIIS::Bordered`/`Coefficients`, apply the one `c` to every channel.  `qchem.Math.DIIS` was
           not touched, exactly as predicted.
         - **The unpolarized path is now literally the one-channel case** — `PulayMixer::MixField` is
           `MixJointly({this},{field})` — so the two paths cannot drift.
         - `PolarizedDensityMixer` forms BOTH channel fields itself and mixes through the field face in
           EITHER basis (it used to forward whole densities in `SpinChannels` mode and fields only in
           `TotalAndMoment`); one constructor now, `(a, b, fit, recip, basis)`.  So the joint history covers
           `(ρ,m)` and `(ρ↑,ρ↓)` by the same code — as required, since the basis choice never rescued it.
         - Gate: `UTChargeDensity`'s `JointPulay.*` (3 cases, new file `tests/JointPulay.C`) — drives a
           hand-built asymmetric two-channel trajectory through both paths against an INDEPENDENT reference
           (its own map algebra + its own DIIS calls) and pins (i) joint ⇒ one `c`, so the extrapolated
           moment is `Σcᵢmᵢ` with the same weights as the density; (ii) split ⇒ `c↑≠c↓` and a moment that is
           off the trajectory (the ejection in miniature); (iii) only extrapolators answer `History()`.
         **STILL OPEN — the STRUCTURAL half of the end state below**: `PulayMixer` still fuses the Kerker
         filter with the history (it takes its own Kerker step inside `ApplyJoint`), and
         `PolarizedDensityMixer` is still a `tDensityMixer`.  The DEFECT is fixed and the fusion can no
         longer split a history, but "preconditioner → extrapolator" is not yet two objects.  Logged in
         doc/CleanupCandidates.md.
       * **★ RUN 29 — THE EJECTIONS ARE GONE (2026-08-10, `doc/logs/mno_afm2_run29_joint_pulay.log`, 45 min).**
         Run 27's recipe EXACTLY (α=0.25, kT=5e-3, PulayDepth=8/Start=5, MOM OFF), joint history the only
         variable, banner confirmed live.
         | | run 27 (split history) | run 29 (joint) |
         |---|---|---|
         | best E | −61.29 | **−61.370** (iter 37; −61.35 at 29/63/72) |
         | m_stag | 0.378 | **0.50–0.66** (peak 0.6625 @50) |
         | ejections | −47…−55 with m→0.1, [F,D]→1.0 (57→59, 67→69) | **none**; those same iterations read −60.99/−61.04 and −61.14/−61.34 |
         | worst excursion after iter 40 | (repeated) | −60.24 |
         0.10 Ha from CP2K at the best iterations.  **The predicted failure mode is retired**; what is left
         is a different seam, below.
       * **WHAT RUN 29 EXPOSES INSTEAD — three things, none of them the mixer.**
         (a) **A TWO-CYCLE IN THE OCCUPATION, and the final eigen table names it.**  Iterations 50–55 read
             m_stag = 0.663/0.511/0.659/0.510/0.603/0.510 — the odd iterations return to 0.510 to three
             decimals.  The diagnostic that settles WHICH seam: **[F,D] falls 8e-3 → 2.2e-3 while Δρ sits at
             ~1e-1 and never shrinks** — orbitals settling with a flat density residual means the FILL is
             what moves.  The final table shows the tie directly: levels 6, 11 and 13 carry IDENTICAL ε↑=ε↓
             with occupations 0 vs 1, and 14/15 sit at 0.50–0.52.  At kT=5e-3 against ties this close μ picks
             a branch and re-picks it next iteration instead of averaging them.  **kT (the `MNO_ANNEAL`
             schedule, still unswept — item (d)) is the lever; α is not**: α damps the INPUT, it does not
             choose the OUTPUT branch, and run 28 already measured that α changes frequency, not possibility.
             (User's reading, 2026-08-10, and it is the right one.)
         (b) **THE LADDER HANDED OFF ON A COINCIDENCE.**  E77=−61.207178 and E78=−61.207208 — two iterations
             of the two-cycle happening to land on the same value — gave |dE/E|=4.91e-07, and the ladder
             advanced to GDM on that.  Δρ at that iteration was **1.10e-1**.  GDM then DECLINED to engage
             twice (`Tr(D'P_block)=12.554` — a smeared occupation is not its leading-13 block) and the
             endgame fell to −59.97 with Δρ=0.122, so the run is scored DIVERGING (Eamp last-8 = 1.41 Ha)
             even though its DIIS leg reached −61.37.  **The rung-advance gate reads dE ALONE**, which is
             meaningless on an oscillating trace while the mixer's own residual says plainly it is near
             nothing.  Gate it on Δρ (or [F,D]) as well, or the ladder will keep firing into oscillations.
             This is the mechanism behind open item (a) below, now caught in the act.
         (a2) **★ THE UNPHYSICAL DOUBLE μ — FOUND IN RUN 29's FINAL TABLE, FIXED SAME DAY (user, 2026-08-10).**
             USER, reading two adjacent rows: *"Spin down at ε=0.26435301 should be >0.51 … it looks like the
             two spin sectors have different chemical potentials."*  They do.  Inverting f=1/(1+e^((ε−μ)/kT))
             on each row at kT=5e-3 gives **μ↑ ≈ 0.29184, μ↓ ≈ 0.26455 — Δμ = 27 mHa** (cross-check: row 15's
             ε↑=0.29179 sits essentially AT μ↑ and reads 0.50).  A lower level cannot be less occupied when
             both are filled from one reservoir, so the numbers are self-consistent under two μ and
             impossible under one.
             **CAUSE — structural, not a bug in the solver**: `Crystal_EC` conserves nUp and nDn separately,
             and EVERY fill path solved μ inside a spin channel.  `FillOrbitalsAufbau` loops over
             `itsSpinWFs`; `FillOrbitalsGlobalFermi` did too, so its "global" only ever meant *global across
             k*.  Two reservoirs, two μ, always — and Δμ is an unrelieved driving force to move charge ↓→↑
             that the constraint holds open, i.e. **the converged state is not the free minimum**.  It also
             makes MnO's AFM constraint IMPOSED rather than found (m_net ≡ 0 by construction), while the
             CP2K oracle constrains nothing (Fermi smearing, one μ; `&BS` seeds the ATOMIC guess only).
             **FIX**: μ is a property of an electron RESERVOIR, so the fill now groups blocks BY RESERVOIR —
             key = (k unless k is shared) × (spin unless spin is shared) — and solves one μ per group
             (`FillOrbitalsGlobalFermi` → `FillOrbitalsSharedFermi`).  The BZ weight rides in the capacity
             `g` only when the reservoir actually spans k, which is exactly when `GetN` is a whole-mesh total
             rather than a per-block count; the metal path is therefore unchanged, and its six gates confirm
             it.  New `ElectronConfiguration::SpinsShareFermiLevel()` (default FALSE — a fixed-multiplicity
             run is constrained ON PURPOSE and must not drift), `Crystal_EC(..., spinsShareFermi)` making
             (nUp,nDown) a SEED PARTITION of a conserved total, and `GpwOptions::spinsShareFermi`
             (`MNO_SHARED_MU=1`).
             **GATE — a behavioural A/B, not an inspection**: `GPW_SCF.SharedFermiLevelLetsTheMomentRelax`
             seeds Si-Γ (closed-shell, gapped, non-magnetic) as a TRIPLET and runs it both ways:
             | | E | |
             |---|---|---|
             | shared μ | **−7.11512** | lands on the singlet anchor −7.11506 — the moment relaxed away |
             | per-channel μ | −7.08181 | stuck 33 mHa above: nUp−nDn=2 is conserved by construction |
             (TRAP, for the next person: the reservoir key CANNOT be an `Irrep` with a null `sym` for "no
             spatial identity" — `Irrep::operator<` dereferences `sym`.  It keys on `sym->SequenceIndex()`.)
             **NOT YET MEASURED ON MnO**: turning it on frees m_net, so the moment may relieve the 27 mHa by
             collapsing.  That is the honest test and m_stag is the instrument.
         (c) **THE MOMENT IS NO LONGER STAGGERED AT THE END**: Mn1=+0.3058, Mn2=−0.5849, m_net=−0.2791 —
             "the moment sits on ONE sublattice", and the symmetry audit reports 10 of 12 ops broken
             (defect 0.27).  NB this is the LOCAL probe 0.7 bohr off each Mn, not an integrated moment, so
             it does not contradict the fixed nUp=nDn; it says the two sublattices have developed different
             moment MAGNITUDES.  Whether that is the GDM leg's damage or predates it is unmeasured — the
             probe history (banked in the log) is per-iteration and would answer it.
       * **★★ THE ROOT DEFECT, LOCALISED (2026-08-11): THE FIRST FOCK BUILD BREAKS THE SUBLATTICE MIRROR.**
         User: *"with one magnetic species on one Wyckoff site"* ferrimagnetism is not physical, and *"something
         is pinning those numbers"*.  Both correct.  The chain of measurements:
         - **The seed is PERFECT.**  m1=+0.737419, m2=−0.737419 (sum 6.7e−16), N↑=N↓=13.0046.  Verified single-
           point, BATCHED (bit-identical), and for the mirror ρ↑(r)=ρ↓(r+t) at points INSIDE and OUTSIDE the
           home cell (≤2.7e−7, the seed's own G-truncation).  Gate: `GPW_SCF.MnOSeedSublatticesAreEqualAndOpposite`.
         - **ONE Fock build destroys it**: m1 → −0.000127 while m2 → −0.731313.  ε↑ sits ~30 mHa BELOW ε↓ at
           EVERY level and μ↑−μ↓ = −50 mHa at the first fill.  This IS the never-explained iteration-1 collapse
           (m_stag 0.366 → 0.0046) — same run, same number.
         - **It is SITE-specific, not channel-specific.**  `MNO_SWAP_SUBLATTICE` (seed the flip on the other Mn)
           gives the exact global-spin-flip mirror — so spin-flip equivariance is FINE — but the moment still
           dies on the atom at (0,0,0) and survives on the one at (½,½,½), whichever spin each carried.
         - **Rigid translation — an EXACT symmetry — changes the answer.**  `MNO_SHIFT=0.13` moves both Mn off
           the cell corner: m1 goes −0.0001 → **+0.4257** (seed value +0.737) and the iteration-1 energy goes
           **+7.373 → −48.53** (Ekin invariant at 87.5, the analytic control).  **The long-standing "the seed is
           so terrible that Etot starts >0" annoyance was never the seed — it is this defect, visible at
           iteration 1 in every MnO log ever recorded.**
         **EXONERATED, each by measurement — do not re-test these:**
         (i) the seed (above); (ii) the Becke PARTITION — per-atom `Sum(w)` identical to 1e−8 across all four
         atoms and 4×74.09339 = det(A) exactly (new `GPW_BECKE_ATOMS` instrument); (iii) the atom-centred mesh
         ENTIRELY — `MNO_XC_UNIFORM=1` runs the XC on the uniform raster and reproduces the defect unchanged
         (m1=−0.00046, m2=−0.7297, Etot=+7.306), so the corner/wrapping hypothesis is REFUTED; (iv) Kinetic /
         Vloc / Vnl — `DISABLED_TermTranslationInvariance` still passes (d = 0, 2.3e−8, 0), so the 2026-07-09
         wrapped-tail fix holds; (v) the Φ tables — Bloch-summed basis, and MnO is Γ-only so wrapping carries
         no phase.
         **⇒ WHERE IT MUST BE.**  The only spin-dependent term in a Fock build is v_xc (kinetic, external,
         local PP, KB and Hartree are all spin-blind and identical for both channels).  The defect survives
         changing the XC QUADRATURE, so it is in the polarized XC TERM'S SPIN HANDLING —
         `DeltaFittedVxcPol` / `DeltaFittedVcorrPol` / `XC_GridEngine::RhoPol` — and it is SITE-dependent,
         which for a pointwise LSDA functional should be impossible.  Note the shape it must have: a uniform
         ~30 mHa offset between the channels' whole spectra, NOT a level-by-level rearrangement.  Start by
         evaluating v_xc^↑ and v_xc^↓ on the SEED and testing v_xc^↑(r) = v_xc^↓(r+t) pointwise; that is
         §7's candidate (a), now reached by elimination rather than by guess.
         **ALSO REFUTED (user's hypothesis, worth recording because it was the best remaining mechanism):
         "cell imaging is simply not working inside Vxc".**  It would have explained all three facts at once
         -- both quadratures live in the home cell so both would break identically, and only a
         boundary-straddling atom would care.  But `GPW_Evaluator::Eval` DOES sum images
         (\f$\sum_R e^{ikR}\chi_i(r-R)\f$), its keep test `|r-Rc-ctr| <= cellRad+maxReach` is exact for r in
         the cell, and the image radius is already hardened FOR THIS CELL:
         `max(2*maxReach+2*maxCellEdge, maxReach+2*cellRad)` with the comment *"for an OBLIQUE cell (MnO
         rhombohedral) 2*cellRad exceeds 2*maxCellEdge, so the historical formula under-enumerated"*.
         `CellImages::Periodic` is the default and MnO uses it, so the `Rcut=0 => itsRc={0} => Eval == the raw
         orbital` trap named in `DISABLED_TermTranslationInvariance` does not apply.  Corroborating: the runs'
         own `ρ_lost/N ≈ 1.6e-5` says the mesh integrates the density essentially exactly, which a
         tail-losing Φ could not manage.
         **WHY NOTHING CAUGHT IT**: `m_stag = ½(m1−m2)` cannot distinguish (+0.37,−0.37) from (0,−0.73) — it
         reads +0.366 for both — so the campaign's own order parameter was blind to the failure it was built
         to watch.  `GPW_MNO_SITES=1` now prints m1 and m2 separately, from iteration 0.
       * **★ THE END STATE FOR `PolarizedDensityMixer` — DEMOTE IT, DO NOT DELETE IT (design Q&A, user asked
         "is it now a (provably?) dead concept or are there cases where it could be useful?", 2026-08-10).**
         Its valid domain collapses to exactly one thing, so the class survives but its SHAPE should not.
         - **Redundant where the leaves are the SAME.**  A memoryless filter is diagonal in channel space, so
           applying the same Kerker to each channel is identical to applying it jointly.  Once a joint
           extrapolator exists, `SpinChannels` with two identical leaves buys nothing — an implementation
           detail, not a concept.
         - **Irreplaceable where the leaves DIFFER.**  That is `TotalAndMoment`, and the code already carries
           the proof: Kerker on ρ paired with plain-linear on m is a 2×2 COUPLED mix in spin space,
           "provably NOT reachable by any per-(ρ↑,ρ↓) leaf pair".  Its physics is independent of anything
           found here — Kerker models the Hartree restoring force against CHARGE fluctuations and the
           magnetisation has no such force — so that argument survives the joint-history fix untouched.
         - **But the split must never extend to HISTORY, in ANY basis.**  Separate ρ and m histories
           synthesise an inconsistent (ρ,m) exactly as separate ↑/↓ histories synthesise an inconsistent
           moment.  `TotalAndMoment` + Pulay leaves is as wrong as `SpinChannels` + Pulay leaves; the basis
           choice does not rescue it.
         **⇒ THE END STATE: demote `PolarizedDensityMixer` from a MIXER to a channel-basis PRECONDITIONER.**
         It stops owning the mixing and becomes the filter stage feeding ONE joint extrapolator:
             `residual → channel-basis preconditioner (may differ per channel) → joint extrapolator (one B,
              one c, all channels)`
         That is precisely the VASP/QE/CP2K architecture — Kerker-preconditioned Broyden over a single
         history — and note **CP2K's `BETA 1.5` IS its Kerker preconditioner sitting in front of one Broyden
         history**.  So the two concepts COMPOSE rather than compete; we had them fused into one object, and
         that fusion is what let the history get split.
         Two consequences that make this more than a tidy-up:
         (i) **It retires the `tDensityMixer` inheritance on that class.**  Being a `tDensityMixer` is what
             made "pure forwarding, without knowing which leaf it holds" look reasonable; as a preconditioner
             it never sees history at all, and the abstraction becomes HONEST instead of merely convenient.
         (ii) **Run 11's verdict on the (ρ,m) basis may not survive the change and should be re-measured.**
             It was recorded as "REFUTED AS A CURE" — but measured with `PulayDepth=0`, i.e. with no
             extrapolator at all.  Preconditioning matters far more once something is learning an inverse
             Jacobian from the residual metric the preconditioner shapes.  Do not treat run 11 as settled.
         **Forward-looking:** non-collinear / Shubnikov work replaces two channels with a 2×2 spin density
         matrix.  "Which channel combinations get which filter" GENERALISES directly; "two independent
         leaves" does not.  A second reason the concept survives and the current shape does not.
       * **PREDICTION CONFIRMED (run 28, α=0.2 — CP2K's own value).**  If the ejections are synthesis rather
         than step size, α changes their FREQUENCY and not their possibility.  Measured: still ejecting —
         7 excursions above −56 Ha in the first 52 iterations, including −60.98 → −49.92 at iteration 52 with
         m_stag 0.396 → 0.126.  So the joint-history fix is the principled answer and α is not.
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

---

## C. Run 30 — the first-Fock mirror break: found, fixed, measured (2026-08-11, late)

**The elimination held and the "NEXT ACTION" probe landed on the defect first try.**  The plan's
contradiction — a pointwise LSDA v_xc failing SITE-dependently — resolved one level below the
functional: v_xc's ρ INPUT was site-wrong.

* **The defect.**  `SeedCD::operator()(r)` evaluated the real-space seed as
  Σ_atoms ρ_atom(|r−R|) over HOME-CELL atom positions only — **no lattice images** — while every mesh
  consumer samples the seed at points **wrapped into the home cell** (`MakePeriodicBeckeMesh` stores
  `kpt = r − A·n0`; the uniform raster's cell-corner midpoints likewise).  At a wrapped point whose
  nearest atom is a lattice IMAGE, the image-less sum reads ~0.  For the CORNER Mn, 7 of its density
  hump's 8 octants belong to images; the cell-centre Mn has no such loss.  Under the AFM staggering the
  erased hump is the MAJORITY channel of exactly one sublattice — which is why the moment died at the
  corner atom whichever spin it carried, why a rigid translation (an exact symmetry) changed the answer,
  why `MNO_XC_UNIFORM=1` reproduced it, and why every MnO log ever recorded started with a POSITIVE
  iteration-1 Etot.
* **Why the seed had been "exonerated".**  `MnOSeedSublatticesAreEqualAndOpposite` probes points within
  ~2 bohr of a HOME atom (including points outside the cell — near the home atom's own unwrapped
  position), where the image-less sum is accurate.  Its mirror comparison r ↔ r+t paired two such
  points.  Nothing ever compared ρ(r) with ρ(wrapped(r)) — the identity the mesh actually relies on.
  Lesson recorded as a pin in the plan: hand-picked probe points near home atoms cannot exonerate a
  periodic density's real-space face.
* **The probe** (the planned v_xc^↑(r) = v_xc^↓(r+t) pointwise test): gate
  `GPW_SCF.MnOSeedVxcMirrorOnBeckeMesh` builds the run's seed + a shrunk production Becke mesh
  (nR=20, GL-11), pulls the channel rasters through `XC_GridEngine::RhoPol` exactly as the Fock build
  does, maps every mesh point to its t-translate partner by hashed wrapped-fractional lookup (the two Mn
  blocks are exact t-translates — no interpolation), and compares ρ and v_xc (Slater + VWN via the same
  functional faces as `DeltaFittedV*Pol`).  BEFORE the fix: 2420/3878 points broken at >1e−6,
  max|ρ↑(r)−ρ↓(r+t)| = 0.761, max v_xc defect = **1.08 Ha**, worst point 0.58 bohr from an **Mn1
  IMAGE** (raster ρ↑ = 1.8e−4 where the true density is 0.76).  AFTER: 2.3e−15 / 6.0e−15.  Orphan
  points (partner tail-dropped by the free builder's eps-borderline keep flips): 6, max w·ρ = 4e−41 —
  the gate bounds their weighted density, not their count.
* **The fix.**  `SeedCD::operator()`/`Gradient` now sum lattice images: per atom, the image window is
  |f_j − n_j| ≤ Rmax·sqrt((AᵀA)⁻¹_jj) (the interplanar-spacing relation the Becke builder already
  uses), with a ball test before each spline call.  FINITE and EXACT — `RadialDensity` clamps to 0
  beyond its last node (`Rmax()` accessor added), so the only truncation is the table's own compact
  support; no distance cut (the no-cut pin).  `PolarizedSeedCD` fixes for free (it composes two channel
  SeedCDs); the molecular (double/AO) SAD path is open-boundary by design and untouched.
* **Run 30 (bounded, 3 iterations, DEFAULT knobs — α=0.45, kT=5e-3, MOM inert until iter 10, no
  anneal):** iteration-1 Etot **−52.804** (every pre-fix log: +7.373); E −52.80 → −59.69 → **−61.3938**
  = **0.077 Ha from CP2K's −61.4706 after THREE iterations** — the whole pre-fix campaign's best was
  −61.370 (run 29, with joint-Pulay + MOM + Λ=1.5 + α=0.25 + seed-pinning).  Site moments
  +0.43/−0.40 → +0.62/−0.58, m_stag GROWING 0.41 → 0.32 → 0.60, mirror held (m_net = 0.038), N↑=N↓
  exact, Δρ = 1.3e−3 and still descending at the cap, zero oscillation flips.  The free-run §3
  diagnostic read **6/12 ops broken** — precisely the grey ops mapping +m → −m: the order-parameter
  instrument seeing the AFM order, as designed (log: doc/logs/mno_afm2_run30_postfix_bounded.log).
  Energy decomposition at the cap:
  Ekin 75.53  Een −81.08  Eee 12.47  Exc −14.13  Enn −57.34  αZ +3.17.  (Note the run built stream
  caches and iterated at ~20 s — the pre-fix "9 min per analytic sweep" cost note did not apply here.)
* **What this retires vs. what it re-opens.**  The "UNDER-BOUND by 15.8 Ha" mystery, the positive-E
  start, the iteration-1 m_stag collapse, and the occupied-d-PP suspicion from the 2026-08-06 blocker
  note are all products of this one defect (the KB d-channel was already oracle-clean at the atom
  level).  The occupation TWO-CYCLE, the 4–6 Ha branch gap, and the MOM/anneal recipe were all measured
  AGAINST the broken v_xc — none of them should be trusted until re-measured on the fixed ground.  The
  next session's first move is a full free run at defaults.

---

## D. Runs 31–34 — the ladder/DIIS instrumentation session (2026-08-11, evening)

The post-v_xc-fix free run exposed two accelerator defects and one physics verdict, in sequence, with
the user watching the logs live (`doc/logs/mno_afm2_run{31,32,33,34}*.log`).

* **Run 31 (defaults):** DIIS descended to −61.4138 by iteration 12 (past run 29's scaffolded best);
  the ladder's tail hand-off then advanced onto GDM, which correctly DECLINED (smeared D′ outside its
  integer-occupation manifold — Tr(D'P_block)=12.80 of 13), and the loop degraded to bare damped-Kerker
  steps under the GDM tag, un-converging to a −60.49 limit cycle.  USER: "for GDM we should not be
  seeing any Ker or relax."  FIX: `tSCFAccelerator::Engageable()` (standing precondition ≠ per-iteration
  readiness; GDM = D′ idempotency [smearing, orbital-free, knowable pre-activation] + leading-block
  [MOM, known post-seed]); the LADDER now VETOES advances onto non-engageable rungs and RETREATS off
  ones that fail while active.  Gates: `ASmearedDensityIsNotEngageable`,
  `HandOffIsVetoedAndRetreatsOnANonEngageableRung`.
* **Run 32 (veto active, pure DIIS):** the SAME collapse at iteration ~27 — so GDM was an accelerant,
  not the cause.  The tell: svMin PINNED at 1.0–1.1e-9 for 10 iterations ("B is near singular but just
  slightly above the SVTol limit" — user, live) against B's own scale |[F,D]|²≈2.5e-5, Nproj stacked at
  7–8: the DEPENDENT-history failure the user predicted on 2026-08-10 when the svMin column was added.
  A wild depth-8 extrapolation amplified the charge slosh.  FIX: `DIISParams::SVTolRel` — prune while
  svMin < SVTolRel·max_i B_ii (dependence, not smallness); default OFF, solid door 1e-4.  Gate:
  `ARelativeConditioningPrunePurgesADependentHistory` (control arm reproduces the MnO shape; the rig
  needs a MULTI-direction perturbation family — [F,D] is linear in F, so a one-direction family is
  exactly rank-2 and the absolute gate already prunes it).
* **Run 33 (both fixes):** NO collapse through the run-32 failure window — Nproj rides 3–7, svMin
  healthy, E −61.4138±1 mHa, m_stag 0.648, gap 0.18 Ha.  But Δρ wobbles 0.8–2e-4 with cfg flips:
  DIIS is at ITS floor (frontier tie-noise), not the run's.  USER: "DIIS is not effective to finish
  the job."
* **Run 34 (`MNO_ANNEAL="5e-3,0"`):** the anneal-to-kT=0 finisher REFUTED for this cell — the very
  first kT=0 fill of the converged orbitals jumped E by +0.77 Ha and [F,D] to 0.53: stage 1's
  occupation is FRACTIONAL within a tied frontier set (the 'm' flag beside an open 0.18 Ha gap), and
  integer aufbau must break the tie.  The kT=0 stage oscillates between an aufbau branch (−60.38) and
  a BETTER non-aufbau hole branch (−61.16..−61.30) — neither is the ensemble; DIIS/GDM never engaged
  (honestly).  VERDICT: at Γ-only sampling the kT=5e-3 fractional ensemble IS the converged
  description; the tie is (at least partly) a k-sampling artifact.  Routes forward in the plan's
  STILL-OPEN: gentle-anneal extrapolation + shared-μ arm (cheap), the k-mesh run (right), a
  free-energy-stationarity convergence gate for smeared ensembles (strategy).

Standing comparison: our smeared plateau A=−61.4140 (E−TS; −TS=−0.0159) vs the CP2K smeared oracle
−61.4706 — 57 mHa, down from 15.8 Ha at the session's start.

---

## E. The Shubnikov S-track lands: S1-S4 in one day (2026-08-11, evening/night)

User decision after runs 30-34: the next semi-major step is imposeSymmetry under the MAGNETIC group,
with an imposed-subgroup option permitting either of two candidate orderings.  Four increments, each
committed green the same day (S1 5ad45c06, S2 d9fe59fe, S3 f8fd4fa0, S4 with this entry):

* **S1** — `SpaceGroup::ShubnikovOps(decorated)` + `CommonOps` (see §D of the plan-side record).  KEY:
  Detect keeps one τ-coset per W; the magnetically doubled cell's EXTRA coset is the anti-translation
  {E|½½½}·Flip = the m1=−m2 mirror — ShubnikovOps re-enumerates every coset itself.
* **S2** — the (ρ,m) diagonalization: total EVEN under Flip, m ODD (χ=±1).  SymmetrizeValuesSigned +
  the FlipFixedPointsPeriodic audit (a Flip-fixed point carries m≡0 exactly — MnO's O sites);
  SymmetrizeGMap(…, oddUnderFlip) exact by scatter; MagneticSymmetryDefects (Flip rows compare ACROSS
  channels).
* **S3** — the wiring: MagneticDecoration by the seed's own species rule (IonicSADTargets factored out
  = one resolution), GPWParams::siteSpins, factory resolution to the Shubnikov group,
  XCQuadrature/engine σ+flags, RhoPol's (ρ,m) split.
* **S4** — the through-SCF measurements, with two live-caught findings:
  - **Run 36**: the DETECTED grey group of the MnO cell is sublattice-PRESERVING (FindTau's first-coset
    rule: every W admits τ=0 fixing both Mn) — so "grey erases m" is structurally unreachable in
    production; the erasure mechanics are proven at the S2/S3 unit gates.
  - **Run 37** (user, live: "m_stag = 0.00000"): the S3 choice to fill the LEGACY op faces
    coset-complete was WRONG — the composite star-average acts PER CHANNEL and the ρ̃ MIXER consumes
    the channels separately through those faces; a sublattice swap is not a symmetry of one channel,
    so m was annihilated machine-exactly at iteration 1.  FIX: legacy faces = the σ=None subgroup
    only; the flip content is enforced at PAIR level (the engine).  PIN: a per-channel face may only
    be projected under ops that are symmetries OF THAT CHANNEL — the spin-blind-mixer bug class, one
    level down.  The Mn₂ toy gate had passed because it runs no ρ̃ mixer.
  - **★ Run 38**: THE FIRST CONVERGED MnO AFM-II of the campaign — 41 iterations, DEFAULT knobs,
    lastΔρ 9.9e-6 PASSES the 1e-5 gate (the magnetic star-average closed the tie-noise floor),
    E=−61.41455 (0.5 mHa below the free plateau = the release-audit), m_stag 0.6513.  FM arm also
    converged and sits 38 mHa BELOW — the ordering is reversed at Γ-only+LSDA(no U), with no CP2K FM
    oracle banked yet; the k-mesh + (+U) routes own that question.
  - Gates added: `MnOImposedShubnikovKeepsTheSeedStaggering` (seed-level end-to-end incl. the grey
    erasure control), `ImposedShubnikovHoldsAFMThroughSCF_Mn2Box` (through-SCF, ~35 s), the
    MNO_IMPOSE knob (0 free / 1 Shubnikov / 2 detected-grey), and the free-run MAGNETIC
    `[symmetry]` line (MagneticSymmetryDefects on the channel pair).

Score at close: MnO AFM-II converged+staggered at defaults, 56.1 mHa from the CP2K AFM oracle —
the campaign started the day 15.8 Ha under-bound with the moment dying at iteration 1.
