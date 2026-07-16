# GPW (Gaussian And Plane Waves) — Plan & State

**Self-contained orientation for a fresh session.** GPW = Gaussian orbitals on a periodic lattice, with
the electron density collocated on a real-space grid, FFT→G-space Poisson for Hartree, XC on the grid,
integrated back against the Gaussians to form the KS matrix (CP2K / Lippert–Hutter). It is the north-star
that makes ab-initio solids → battery voltage curves possible.

This doc supersedes the GPW sections of `doc/MolecularPP_HarmonizationRound2.md`.

**The doc is split into two major sections: [DONE](#done) (compact timeline + the still-load-bearing
records) and [TODO](#todo--next) (what's left, in priority order), then the durable invariants + pointers.
Full archived narratives live in `doc/GPWHistory.md` — read THIS file to orient; open the history only for
archaeology.**

---

# DONE

Everything here is committed on `main`; the GPW suites (`GPW_UT`, `GPW_SCF_UT`) are green.  GPW is a **new
evaluator, not a new IBS** — it satisfies the plane-wave concepts and the whole `Ham_PW_DFT` KS stack drives
it verbatim.  **Full per-increment narratives: `doc/GPWHistory.md`** — below is the compact timeline, then the
still-load-bearing records in full (naming, the CP2K recipe, the C+D analytic-rewrite state, and the §0a
runtime close-out incl. the CP2K NaF oracle + convergence findings).

## Compact timeline (details in doc/GPWHistory.md)
- **1E at Γ** (`ab2c6a76`): Bloch lattice sums delegated to the molecular basis via the engine-neutral
  `Molecule::LatticeSum1E` seam (new edge qcLattice_BS→qcMolecule_BS); home cell == finite matrices <1e-12.
- **DFT tier by collocation** (`cc123b3b`,`63fbf70c`): GPW fills the PW `Repulsion3C`/`Overlap3C` tensors →
  the entire PW_Hartree/PW_XC/IrrepCD stack reused; Coulomb factorised through G-space (weight × 4π/G²).
- **First-light periodic SCF** (`dcef8528`,`db314e6a`): `Integrals_Pseudo<dcmplx>` realised (G-space local PP
  — box-independent, PW G=0 convention; KB via qcMesh) → the real `cSCFIterator`; atom-in-box == finite DFT.
- **General-k + multi-k plumbing** (`b2a29249`): `e^{ik·R}` through the stack; one `GPW_IBS` per BZ k with
  weights.  SIPP→SIPP_SR conditioning lesson (ill-conditioning is a BASIS problem); N3/N5 removed (`10ad6e29`).
- **Bulk over-binding root-caused + FIXED == CP2K** (`a4c94ec5`,`95e8f4a8`): the 16 Ha translation-variance
  was the KB bra using the RAW home orbital (fix: the Bloch orbital); the FFT-raster suspicion was a red
  herring.  Γ SR/2a −7.11505 == CP2K −7.11506.
- **Multi-k validated vs CP2K; complex-k FIXED** (`5fe61aeb`,`1980d6ef`,`745d03ff`): Γ-centred 2×2×2 matches
  CP2K grid-for-grid; the first genuinely-complex k exposed two GPW-evaluator bugs (density ket-conj slot;
  KB image phase must be conj) — fixed; shifted 2×2×2 −7.86673 == CP2K default −7.86744.  Also
  `GetTotalCharge` Tr(D Sᵀ)→Tr(D S).
- **NaF convergence campaign — correctness closed** (2026-07-12): auto `densityEcut` (basis-derived floor,
  <0/0/>0 convention), `∫ρ_grid−N` readout (`ReportGridCharge`), trajectory fingerprint, ionic-seed library
  (PW iters 35→17), Kerker ρ-mixing (`FourierMixCD`; Si-exact −8.24758; +DIIS tames the NaF charge-transfer
  limit cycle).  Grid under-resolution was the dominant cause; then real charge-transfer dynamics.
- **Runtime round 1** (`7708d2dc`,`05e44fab`): OverlapMatrix→zgemm 4×; OpenBLAS pinned 1 thread; magnitude
  screening on the 1E lattice sums (~4×, PSD at any enumerated reach).
- **Runtime round 2 = sampling multigrid DEAD END** (`c94269c8`..`38b63d7b`): sampling collocation rings,
  aliases at bulk (2.66 Ha), needs a hard Rcut → pivot to the analytic method.
- **Analytic kernels A/B/cross-cell** (`0d09a6d5`,`068b4e96`,`729b6355`): per-pair exp-tail boxes +
  modulo-wrap, exact adjoint, screened cross-cell offsets, `G_ERI3::apply` matrix-free seam.
- **§0a Si runtime leg** (`9ff982ba`): stream-cache lockout fix + coverage readout + same-D/phase-independent
  memos; Γ 157→31 s, multi-k 475→89 s, bit-consistent; shifted-MP complex-k gate ENABLED == CP2K to 0.2 mHa.
- **§0a NaF leg** (`b0f497c6`): ANALYTIC KB via the `⟨χ|g⟩` Gaussian seam (== mesh to 4.6e-11; the >33-min
  mesh setup wall dead) + fp32 stream tier; NaF end-to-end 2h15m.
- **§0a D-aware radii + CP2K NaF oracle** (`4c71450c`): eps/|coef| kill+shrink with shared-active-set
  integrate (adjoint stays exact); NaF 40m41s; CP2K same-basis oracle **−27.93128** (own q-tag-free basis).
- **§0a NaF convergence findings** (`35789164`): CP2K recipe machinery (no DIIS, E-gate, tuning knobs);
  α=0.025/G0=1 converges Ecut=40 (pinned anchor −27.73); the fine grid's unphysical attractor (E≈−39)
  captures ALL linear mixing → quasi-Newton mixing + XC consistency are the TODO leads.

## Naming (`5f609d2f`) — remember these
- `Overlap(f)` = ANY 1-electron `⟨i|f|j⟩` (f may be a potential); `Repulsion` = the 2-electron `1/r12`.
  So the reciprocal-space field→KS-matrix bridge is `Band_FT_IBS::MakeOverlap(f)` / evaluator
  `OverlapMatrix(f)` — **not** `MakePotential`/`PotentialMatrix`. `Make` = uncached.

## THE CP2K METHOD (Quickstep / Lippert–Hutter) — the authoritative GPW recipe (deep-dived from `~/Code/cp2k`)
Read the CP2K source (`src/grid/ref/grid_ref_{collocate,integrate}.c`, `qs_collocate_density.F`,
`qs_integrate_potential_product.F`, `pw_env/gaussian_gridlevels.F`, `task_list_methods.F`, `aobasis/ao_util.F`).
The recipe — every piece fixes a wall we hit:
1. **ANALYTIC collocation, NOT sampling.** Each primitive PRODUCT is ONE Gaussian: `p=z_a+z_b`, centre
   `R_p=(z_a R_a+z_b R_b)/p`, prefactor `exp(−z_a z_b/p·|R_ab|²)`, times a Cartesian polynomial (binomial
   re-expansion about `R_p` — CP2K's `cab_to_cxyz`; **we already have all this in `Ω`/`H2` in `GaussianRF.C`**).
   Evaluated analytically on grid points inside an exp-tail radius — never a sampled pre-summed orbital.
2. **No Gibbs ringing by construction.** The box ends where the poly×Gaussian `< eps_rho_rspace` (a smooth
   tail), so there is no truncation discontinuity. (This is the fix to the hard-`Rcut` ringing.)
3. **Integrate-back = exact adjoint** (same kernel, gather flag flipped): gathers **Hermite moments of V** over
   the same box. Only **V** is sampled (weighted by the analytic Gaussians), never the sharp orbital product —
   which is WHY it stays accurate on a coarse grid where naive sampling aliases.
4. **REL_CUTOFF multigrid, done right.** Each pair → the coarsest level with `cutoff ≥ p·rel_cutoff`
   (`gaussian_gridlevel`); V is transferred to ALL levels up front via FFT (spectral → no ringing). Analytic +
   matched grid → coarsening is accurate (unlike our sampling multigrid). This is the ~10–100× speed.
5. **Periodicity + screening, no hard cutoff.** Density is collocated from the DENSITY MATRIX `P` over
   NEIGHBOUR-LIST pairs `(i, j@cell R)` — a screened image sum (include only where `|⟨χ_i|χ_j^R⟩| > EPS_PGF_ORB`,
   default 1e-5) — with each compact box MODULO-WRAPPED onto the grid. So: a **screened** image sum (no hard
   Rcut → no ringing) PLUS the wrap (an atom at the cell edge tiles automatically). k-points: the grid density
   is always real/cell-periodic; ALL k-dependence lives in `P(R)=Σ_k w_k e^{ikR}` — collocation is k-agnostic.

## GPW ANALYTIC REWRITE COMPLETE -- Increments C + D LANDED, sampling machinery DELETED (2026-07-14)
The analytic collocate/integrate path IS the SCF path; the whole sampling stack is gone.  Si SR/Rcut=2a Gamma
== CP2K to 0.18 mHa; the Gamma anchor runs in ~40 s (warm) / ~2 min (cold).  14/14 GPW+GPW_SCF green.
- **Kernels multigrid + general-k** (`LatticeSum1E::CollocateDensity/IntegratePotential`): level-vector form
  (finest-first `N_L`/`ecut_L`; K=1 = single grid; the old `MakePotentialMatrix{,MG}` sampling faces DELETED,
  interface stays minimal per the Band_FT_IBS lesson).  Bloch phases enter via a `cellphase_t` CALLBACK
  (`e^{ik.R_n}` of the integer cell offset -- the k-convention stays lattice-side; the kernels enumerate offsets
  internally so a pre-built (Rs,phases) pair cannot work).  Phase conventions: density weight `Re[D_ij e^{-ik.R}]`
  (ket-conj), integrate-back `e^{+ik.R}` -- the same +/- pairing as the KB-projector fix.  2x Hermitian fold
  (loop j>=i, double off-diagonal weights -- the (j,i,-R) twin is exact).
- **REL_CUTOFF assignment + level guards.**  Pair -> coarsest level with `ecut_l >= kRelSafety *
  ecut_fine*(a_i+a_j)/(2 a_max)`; **kRelSafety=2 is LOAD-BEARING**: at 1x the charge-calibrated ratio leaves
  ~e^{-2.5} spectral tails per pair at its own cutoff -- fine for charge (1e-5) but Si Gamma sat 5 mHa below
  CP2K until doubled (CP2K's REL_CUTOFF default is ~3x stiffer than the auto floor).  `EnsureLevels` keeps a
  level only if its SPACING resolves the sharpest assignable pair (`h^2 p_max <= 1`, Poisson e^{-2pi^2}) -- this
  replaced a naive min-N floor after a degenerate Ecut=0.375 level came out as a 1x1x1 GRID whose sigma~2.4
  pairs lost percent-level charge (the crystal-gate 28.37-vs-26.56 failure).
- **Wiring (Band_FT_IBS + IrrepCD UNTOUCHED).**  `Repulsion3C/Overlap3CTensor` -> matrix-free
  `G_ERI3::apply=MakeCollocator(coulomb)` (per-level collocate -> FFT -> rho-tilde combined NESTED in G-space,
  Coulomb kernel folded); `OverlapMatrix(Vtilde)` -> per-level spectral restriction of V + analytic
  `IntegratePotential`.  `G_ERI3::weights` (dense tensor) DELETED.
- **Local PP stays FINE-level (every pair): the analytic sibling of the sampling-MG lesson.**  Routing V_loc
  through the ladder gave Si Gamma -10.72 vs -7.115 WITH charge exact -- a coarse-level pair sees V only
  through its level's {G}, and V_loc's mid-G content is a Ha-scale term.  V_H/V_xc (smooth) stay on the ladder.
  Static, so the one fine sweep is not per-iteration.
- **Two latent bugs fixed:** the seam's `Recip().GetCell().MakeReciprocalCell()` double-reciprocal
  RE-ORIENTS a non-cubic (FCC) primitive cell (evaluator now STORES `itsCell`); the box's fractional bounding
  used `reach/minEdge`, UNDER-COVERING a skewed cell by sqrt(3)/sqrt(2) (now exact per-axis `reach*||row(A^-1)||`
  -- this WAS the committed crystal gate's 2.4e-7 residual).
- **Perf: pair-box STREAM CACHE (the analytic PhiOnGrid successor).**  The per-(pair,offset) box streams
  (wrapped idx + value) are pure geometry -- identical across SCF iterations AND k-blocks (phases apply at
  contraction) -- so they are built once per ladder shape and replayed as gathers: Si Gamma SCF 51.5 min ->
  2m15s cold (23x), bit-identical.  Plus: prefactor-SHRUNK per-offset reach (a far cross-cell image costs a
  near-empty box), ellipsoid pre-screen before the exp/poly evals, incremental grid walk (hoists the
  per-point ToCartesian, was 14% of the profile).  perf now works locally (paranoid=1).
- **Anchors re-pinned on the CONSISTENT scheme (all in `GPW_SCF_UT.C`).**  The analytic collocation is ALWAYS
  screened-complete Bloch, so `Rcut=0` 1E matrices would MIX SCHEMES (observed: `Tr(D S_home)=8` while
  `int rho_grid = 16.55` -- the durable-pin violation).  The old fast Rcut=0 anchors are GONE:
  `SiliconGammaConverges` = SR/Rcut=2a/dE=20 == **CP2K -7.11506 to 0.18 mHa (-7.114883, 18 iters)** [the old
  DISABLED CP2K gate, now THE enabled anchor]; `SiliconMultiKPlumbing` = SR/2a 2x1x1 **-7.45137** (real
  dispersion; the "k-blocks==Gamma at Rcut=0" trick is gone); `SiPseudoAtomInBoxMatchesFinite` box grown
  a=11->16 (cross-cell prefactor 2.7e-2 -> 4.6e-4) == finite sipp to 0.023.  Kernel gates re-referenced to
  **Tr(D S^Bloch)** (charge now conserves to 8e-8; the old "3%" was the WRONG home-only reference, latent since
  the cross-cell commit) + a machine-precision multigrid seam-adjoint gate (Tr(D H) == <apply(D),V>).
- **[grid charge] readout through the live SCF: lost = -1.4e-6 e** -- collocation charge-exact in production.

## §0a RUNTIME CLOSE-OUT — COMPLETE (2026-07-15/16).  The full records:

**(0a) Si LEG DONE (2026-07-15) — Γ 157→31 s (5×), multi-k 475→89 s (5.3×), all anchors BIT-consistent;
complex-k REVALIDATED (shifted-MP gate ENABLED).**  The profile OVERTURNED the commit-message attribution:
the 6→14 min suite regression was NOT the AUTO enumeration radius (the analytic kernels enumerate offsets
per-pair internally, Rcut-independent; every O(|Rs|) consumer loop is cheap norms).  ~85% of the multi-k
anchor was the pair-box kernels re-evaluating analytically per iteration — the STREAM-CACHE BUDGET (added in
the same commit) was the regression:
- **EnsureStreams lockout bug**: after the FIRST over-budget pair, `budget=0` un-cached every later pair.
  Fixed to skip-and-continue packing; + a one-line `[stream cache]` coverage readout per build (pairs
  cached/total, pts cached/dropped) — the tuning instrument for NaF.
- **Budget 100M→150M pts** (~1.8 GB): Si SR demand is 104.9M — at 100M its 7 most-DIFFUSE pairs (the biggest
  boxes) re-evaluated every iteration × k-block ≈ the whole regression.  Si now caches 300/300.
- **Same-D collocation memo** (`GPW_Evaluator::CollocMemo`, shared by the Coulomb + overlap tensor closures):
  each iteration collocated the SAME D twice (RefreshRhoGrid + GetRepulsion3C, ~10% each in the profile);
  the second call now replays the level densities.  EXACT-equality keyed on D → bit-identical.
- **Phase-independent integrate-back memo** (`NR_Evaluator::IntegrateMemo`): h_ij(k)=w Σ_n e^{+ik·Rn} B_ij(n)
  with B k-INDEPENDENT — memoized on the EXACT (ladder shape, scale, V_L), so the static local-PP sweep
  (~10% PER k-block) is paid once per geometry and the per-iteration KS fields once per V instead of per k.
  Contraction order == direct evaluation order → bit-identical on hit (field equality is exact per-element;
  NEVER blaze relaxed equal).
- **Complex-k through the analytic kernels: VALIDATED.**  `SR_2x2x2ShiftedMP_vs_CP2K` (8 k-blocks, genuinely
  complex phases) ENABLED as a regression gate at AUTO Rcut: **−7.86724 vs CP2K −7.86744 (0.20 mHa)**, charge 8,
  CONVERGED Δρ=4.5e-8, ~2.5 min (the memos make the 8 k-blocks share the static sweeps).  Γ-centred 2×2×2
  stays disabled (redundant coverage).  196/196 UTMain green.
**(0a) NaF leg (2026-07-15, same day): ANALYTIC KB + fp32 stream tier LANDED — the setup wall is dead;
per-iteration collocation volume is now the whole NaF story.**
- **ANALYTIC KB ASSEMBLY (the big one).**  The measured NaF setup wall was `MakeSeparablePP`'s mesh
  quadrature — `Eval` (the truncated-Bloch orbital sum) over a 358k-point eCut=160 mesh ≈ billions of exp
  calls: the mesh-path run burned **>33 min without finishing setup**.  CP2K never touches a grid here: GTH
  projectors are polynomial×Gaussian, so ⟨χ|β Y_lm⟩ is analytic.  Now ours is too: qcPseudopotential grew the
  OPTIONAL capability face `SeparablePotential_Gaussian::BetaGaussian` (the radial's CLOSED Gaussian form
  Σ_t c_t r^{l+2n_t} e^{−α_t r²}; HGH/GaussianProjector/MultiSpecies implement it), and the molecular seam
  grew `LatticeSum1E::MakeOverlap(Rs, phases, GaussianFunction)` — b_i = Σ_R phases[R]⟨χ_i|g(·−R)⟩ with
  g = {centre, α, Cartesian-monomial terms}: PURE Gaussian language (user pin: the basis interface talks
  integrals-over-functions; no Fourier/potential vocabulary).  GPW expands β·Y_lm → monomial Gaussians
  (`YlmCartesian` pins `Math::SphericalShell` to the mesh path's own `RealYlm` convention numerically;
  `MultiplyR2` folds the r^{2n} powers) and calls the seam per radial term.  Models without the face keep the
  mesh path (contract intact).  **Gate `GPW.AnalyticSeparablePPMatchesMesh`: analytic == mesh to 4.6e-11**
  (SR/AUTO complete enumeration; at an UNDER-enumerated Rcut the two paths truncate differently — the mesh's
  Bloch orbital reaches χ-image×β-image separations up to 2·Rcut, the analytic single sum stops at Rcut — the
  "two schemes" pin again; measured 9.3e-2 for diffuse SIPP at 1.5a, so the gate pins the COMPLETE setting).
  All four SCF anchors byte-stable (the Si mesh KB was already converged; the win is runtime).
- **fp32 STREAM TIER (the coverage lever).**  Stream budgets are now TWO-TIER: fp64 150M pts (bit-identical
  replay; all Si shapes live here → every anchor/kernel gate unchanged) + fp32 700M pts (~5.6 GB; overflow
  pairs store float values instead of falling to on-the-fly; ~6e-8 relative replay noise, invisible at NaF's
  anchor scales; the collocate/integrate ADJOINT stays machine-exact — both directions replay the SAME
  stream).  NaF coverage 16% → 89%.
- **NaF end-to-end (charge 8.0000000000, 60-iter cap): 2h15m, peak RSS 8.2 GB** on the analytic-KB +
  fp32-tier build.  The remaining cost was PER-ITERATION collocate/integrate volume: 850M cached pts
  replayed ~5 sweeps/iter + 314 small pairs (102M pts, first-fit packing victims) on-the-fly each sweep +
  the one-time scale-6 static-PP fine sweep.  CP2K calibration on this box: Si Γ 3.6 s (ours 31 s),
  Si 2×2×2 shifted 32 s (ours 149 s).

**(0a) D-AWARE RADII + FULL PACKING + THE CP2K NaF ORACLE (2026-07-15, later the same day; 198/198 green).**
- **D-aware density-magnitude screening (CP2K's eps/|coef| radii), `kDensityEps=1e-10`.**  What lands on the
  grid is c·χχ (c = fold·Re[D e^{−ik·R}]), so the tolerance a box must honour is eps/|c|:  (a) each cached
  stream stores its max|value| and replay SKIPS a (pair, offset) whole when |c|·maxv < eps (one compare);
  (b) on-the-fly boxes get the CONTINUOUS shrink — eps/|c| threaded into `ForPairBox` (clamped so |c|>1
  never grows past the geometry screen);  (c) `IntegratePotential` gains an OPTIONAL `screenD` (the seam
  already speaks `chmat_t` densities): the SAME |c|·maxv criterion keeps the IDENTICAL active set in both
  directions, so the collocate/integrate ADJOINT stays machine-exact on the shared truncated operator (the
  variationality ledger's property).  GPW passes its `CollocMemo` D (the iteration's own density); screened
  calls bypass the V-keyed B-memo (cheap by construction); the static PP keeps memo + full sweep.  A pure
  magnitude screen (smooth tails, no Gibbs).  Machine gates UNCHANGED (charge 8.5e-8/2.1e-7, adjoint exact,
  analytic-KB 4.6e-11); Si anchors within pins (Γ/shifted identical to print; multi-k −7.45133 vs −7.45134,
  trajectory 14→17 iters — kills drop 1e-10-level terms, not bit-identical by design).
- **fp32 budget 700M→850M**: NaF now caches 528/528 pairs (0 dropped; 76 fp64 + 452 fp32), peak RSS 9.0 GB.
- **NaF re-time: 2h15m → 40m41s (3.3×).**  Setup (stream build + static-PP sweep) is now a large fixed
  share; the D-aware kills are WEAK while the density sloshes (large |D| everywhere) and strengthen as it
  settles — so the next multiplier is convergence itself.
- **CP2K NaF ORACLE (doc/CP2Kresults.md): Etot = −27.93128 Ha** on OUR transcribed low-q SR basis
  (`naf_gpw_sr_diag.inp`: q-tag-free own basis fixes the q1-vs-q9 abort; damped Broyden α=0.2 +
  diagonalization).  CP2K's ENERGY settles to 1e-6 by ~130 iterations while its DENSITY limit-cycles forever
  (RMS 0.03–0.12) — the SAME charge-transfer cycle we see (its OT run never settled E at all, −25.7↔+253):
  the disease is the system+basis (overlap cond ≈ 8e3), not either implementation.  CP2K's grid also leaks
  2.0e-4 e at 320 Ry (our readout's class).  **OUR Kerker(G0=1)+DIIS at relax 0.3 does NOT settle E in 60
  iterations — iteration 60 lands essentially randomly (−24.03, +887.55 across two runs; charge exactly 8
  throughout).**
**(0a) NaF CONVERGENCE increment (2026-07-16): the linear-mixing axis is EXHAUSTED — the production grid
needs QUASI-NEWTON DENSITY MIXING (the one CP2K ingredient we lack).**
- **Recipe machinery landed** (`DISABLED_NaFRocksaltGamma`): `tSCFAcceleratorNull<dcmplx>` (NO DIIS — the
  mid-cycle Fock extrapolations ARE the +900 Ha spikes: they land exactly on the Nproj=8 iterations),
  fixed-α Kerker, exit on the relative-E gate `MinΔE` with `MinΔρ=1e30` (CP2K's density never converges
  either — E-flat is the physical criterion), env tuning knobs `NAF_{ECUT,ALPHA,KERKER_G0,NMAX}`.
- **α scan at Ecut=40** (cheap grid): α=0.2/0.1 → ±75 Ha period-~48 limit cycles that pass THROUGH the fixed
  point; α=0.05 → contained ±1 Ha, not decaying; **α=0.025 → converges** (~−27.75, ±0.04 residual wobble;
  −27.7304 at the pinned 200-iteration endpoint).  G0=1.0 is the sweet spot — BOTH 2.5 and 0.5 destabilize
  (the Kerker screen must match the charge-transfer mode, not smother or under-damp).  The Ecut=40 answer
  sits 0.2 Ha above the 320-Ry oracle — the leaky-grid gap (Ecut=40 loses >5 e⁻ of F's collocated density).
- **The FINE (auto=160) grid grows a second, UNPHYSICAL attractor** at E≈−39 (Exc≈−143, ∫ρ_grid swinging
  5.1↔7.7 vs Tr(DS)=8): the mid-slosh D loads the sharpest F pairs beyond the grid calibration; the XC of
  that spiky/locally-negative ρ feeds back; the state is self-consistent garbage.  It captures plain damped
  Kerker at EVERY α (0.2 → 0.01 all dive in, sliding past −26 on the way).  Damping sets the rate, not the
  destination — a wrong basin needs a different METHOD.  CP2K converges the SAME map with BROYDEN
  (quasi-Newton, 8-step history, α=0.2).
- **Test now pins the CONVERGING regime** (Ecut=40/α=0.025/200 iters → −27.73 ± 5e-2; both bad attractors
  land ~+65 / ~−39, far outside): a true mixing-regression anchor until the production grid converges.

---

# TODO / NEXT

**Orientation (2026-07-16).**  §0 (A–D + the runtime close-out) is DONE — the analytic multigrid path IS
the SCF path; Si == CP2K at Γ / 2×1×1 / shifted 2×2×2; NaF runs end-to-end in ~40 min against a same-basis
CP2K oracle (**−27.93128**, `doc/CP2Kresults.md`).  The open problem: NaF's PRODUCTION-grid SCF is captured
by an unphysical attractor (E≈−39) that plain damped mixing cannot escape at any α.  The two lead
increments attack that — **0b first** (it deletes the bad basin itself), then **0c** (robust navigation) —
followed by the runtime follow-ups (0d) and the standing queue (1)–(5).

## 0b. XC CONSISTENCY + ρ-FLOOR — one discrete functional (deletes the bad basin; unblocks GDM/OT; the GGA prereq)
**The ringing/variationality ledger (user pin, 2026-07-14), where each error source stands:**
- Gibbs PROPER (hard SPATIAL truncation) — GONE: pair boxes end in smooth exp tails; cross-cell offsets +
  1E sums magnitude-screened; auto-Rcut derives the enumeration from the basis, so S is PSD to eps << λ_min.
- SPECTRAL truncation (density-grid Ecut, per-level G-spheres) — PRESENT but EXPONENTIALLY controlled
  (e^{−ecut/2p} via kRelSafety/relCutoffScale), NOT Gibbs-like, and rendered variationally SAFE by
  ADJOINT-EXACTNESS (the machine-precision seam gate Tr(D H)==⟨collocate(D),V⟩: the KS matrix is the exact
  gradient of the DISCRETIZED energy).  Truncation alone does not break GDM; energy/gradient INCONSISTENCY does.
- **The ONE remaining inconsistency = XC**, and it is now KNOWN to be load-bearing: the NaF fine-grid
  garbage attractor (E≈−39, Exc≈−143) lives off the E_xc/H_xc representation FORK — the MATRIX carries
  v_xc FITTED onto the finite {G} (band-limited), while the ENERGY takes the ¾-virial `E_x=¾⟨ρ|v_x⟩` +
  `FittedEpsXc` correlation on the raw grid.  "E rewards density spikes; the band-limited H never resists"
  → a self-consistent garbage cycle BELOW the physical state, which also means variational minimizers
  (GDM/OT) would seek it out BY DESIGN.  (CP2K has no such state: one pointwise representation + floors;
  its E settles from an atomic guess with a plain quasi-Newton mixer.)

**The fix = the Hartree precedent** (the machinery already exists):
1. Evaluate v_xc POINTWISE on the grid density — no {G} fit — and feed it through the per-level spectral
   restriction + analytic `IntegratePotential` (adjoint-exact, exactly the V_H route).
2. Take E_xc = Σ w·ε_xc·ρ on the SAME grid — retiring the ¾-virial and `FittedEpsXc` on the periodic path
   (also the GGA prerequisite: the old "route E_xc through ∫ε_xc·ρ" TODO).
3. **ρ-FLOOR (the CP2K guard)**: ρ<ε → ε_xc,v_xc treated as 0 before the XC evaluation.  The raw density
   of a PSD D is pointwise ≥0, but the KERKER-MIXED auxiliary density is a FILTERED field and CAN go
   pointwise negative — that mixed field feeds V_xc every iteration.  Audit the VWN path on pathological
   input while there (only SlaterExchange has the ρ≤0→0 guard today).
Result: ONE discrete functional; every SCF state a genuine stationary point (H_xc = ∂E_xc/∂D exactly);
plausibly widens the Kerker α-window enough that existing mixing converges the production grid.
Instruments: the E(λ) line-search + FD potential-consistency probes; develop on the 31 s Si anchor,
validate on NaF vs the −27.93128 oracle.

## 0c. PULAY/BROYDEN ρ̃-MIXING behind the DIP mixer face (`tDensityMixer`) — user design, 2026-07-16
Mixing is today hardwired inside `tSCFIterator::Iterate` (the `KerkerG0>0 ? KerkerUpdate(relax) :
MixIn(1−relax)` branch + the inlined adaptive-α heuristics).  Extract the face and inject the concrete
from the top (SOLID DIP — the existing `tSCFAccelerator<T>*` ctor-injection precedent):
- **Face** `tDensityMixer<T>` (qcChargeDensity — it speaks ChargeDensity and needs FourierMixCD; no new
  lib edges): `double Mix(cd_t& cdInOut, const cd_t& cdFresh, const MixSignals&)` + `Reset()`;
  `MixSignals={E,[F,D]}` so adaptive policies live INSIDE concretes.
- **Concretes**: `NullMixer` (pass-through — what a GDM/OT-driven SCF wants: a minimizer must not fight a
  mixer); `LinearMixer(α₀)` (today's D-mixing + the adaptive-α policy moved in VERBATIM — molecular SCF
  bit-preserved); `KerkerMixer(α,G0)` (today's KerkerUpdate + its periodic-basis validation moved into
  construction); `PulayMixer(α,G0,m)` (NEW: last-m (ρ̃_in, residual) history, small residual-norm LS,
  Kerker-preconditioned update — the VASP/QE/CP2K scheme; Broyden = a sibling behind the same face).
  Null/Linear T-generic; Kerker/Pulay dcmplx/periodic-only.
- **Plumbing**: `cSCFIterator` ctor gains the mixer pointer beside the accelerator; the Calculation facade
  constructs the concrete (options beside AcceleratorOptions); `SCFParams.KerkerG0/StartingRelaxRo` remain
  as facade DEFAULTS (no call-site break) and the iterator stops reading them.
- **Accelerant on top**: grid-continuation seeding (converge Ecut=40 → seed the fine grid — start in the
  right basin).  Gate: the NaF test on the production grid vs the −27.93128 oracle.
Convergence pays twice: fewer iterations AND stronger D-aware kills on a settled density.

## 0d. Runtime follow-ups (after 0b/0c)
- **OpenMP over pairs** (the parallel-execution TODO): collocate/integrate are embarrassingly parallel
  (per-thread ρ accumulators); CP2K's ssmp is threaded on top of everything above.  ~4× on this box.
- **Setup share** (stream-cache build + the scale-6 static-PP sweep) = the next profile target once
  per-iteration cost is mixing-limited.

Then the standing queue: **(1) DROP SR** (rank-reduction + auto-tol, below); **(2) low-q multi-species
bases → Si/NaF/CsI**; **(3) CP2K reference library**; **(4) IBZ**; **(5) cleanups**.

## 1. DROP SR — rank-reduction through the periodic stack + auto-tol
The `_SR` basis is a hand-tuned crutch (drop the most-diffuse primitive so the Bloch overlap is cleanly PD).
We PROVED (2026-07-13; record: doc/GPWHistory.md) that the FULL basis + screening + canonical Eigen/SVD ortho with tol in the
~1000× spectral gap gives a clean overlap transform (‖VᴴSV−I‖=6.6e-11) — BUT the SCF is **BLOCKED**: truncation
reduces the working dim (NaF 37→33) and the periodic stack (`Crystal_EC`/`cDM_CD`/collocation) assumes the full
`n` → "Matrix sizes do not match" (`DISABLED_NaFFullBasisEigenTol`). The MOLECULAR path handles rectangular V;
the PERIODIC path does not. So dropping SR = two pieces:
- **(a) Rank-reduction through the periodic stack** — let a truncated ortho (`V` is `n×(n−k)`) flow through
  `Crystal_EC` (band count `n−k`), `cDM_CD` (density still full `n×n` via `C=V·U'`), and the collocation;
  mirror the molecular path's rectangular-V handling. This is the real work and gates (b).
- **(b) The user-friendly automation** (agreed design; resolved-investigation record: doc/GPWHistory.md): **auto-Rcut**
  [**DONE `9714f58d`** — Rcut<0 = AUTOMATIC, radius from the basis, 3-mode convention] via a basis reach scalar (wall B — the lattice enumerates `CellsInSphere(MaxReach+span)`; exponents
  stay behind the molecular-basis wall, k-convention stays lattice-side), removing the `Rcut` param for one ε
  (CP2K `EPS_PGF_ORB`; CP2K sets no user Rcut). **Auto-tol** via `LASolver` GAP DETECTION (pure LA): force-drop
  `d[i]≤0`, scan the low region for the largest consecutive ratio; if `> R_threshold` (default **30**, exposed
  at the Calculation facade) it's a CLEAN gap → cut there, else fall back to the ε-tol + WARN. `orthoTol<0`=auto
  / `=0`=none / `>0`=explicit (mirrors `densityEcut`). **Auto-cut allowed but NEVER silent** — always `cerr` WARN
  (count + gap ratio + clean/ambiguous). Vision: collapse to ~one CP2K-like ε.

Until (a) lands, **SR stays** (dimension-preserving, cleanly PD, no truncation).

---

## 2. Low-q multi-species bases → Si/NaF/CsI cross-validation (PW + GPW + CP2K)

**Valence-basis GENERATOR — DONE** (`qchem.ValenceBasisGen`; full record: doc/GPWHistory.md): pseudo-atom
SCF → even-tempered valence blocks → `BasisSetData/valence_lowq.bsd` (F 8s+6p E=−21.10, Na 5s+2p E=−0.144;
enum `VALENCE_LOWQ`), tests `UnitTests/ValenceBasisGen_UT.C`.  Pinned lessons: validate against the physically
relevant CHARGE STATE (F⁻ for NaF); oracle GS-energy matching is the WRONG objective (N≈8 windows, refine later
from a NaF orbital-coefficient heat-map); keep per-l exponents DISJOINT (the shared-exponent Gaussian94 reader
bug is flagged in `PG_Cart/Imp/IrrepBasisSet.C` — flipping it re-pins every density-fit anchor).  NEXT: Cs/I.

Hand-roll SIPP-style **low-q valence Gaussian bases** for Na/F/Cs/I so GPW (and CP2K) can run NaF + CsI, then
triangulate our two codes against CP2K on Si/NaF/CsI. Unblocks **multi-species GPW** (the battery-oxide path,
[[project_battery_voltage_goal]]) and yields the CP2K runtimes. The CP2K reference library (§2) is the oracle.

**Why blocked today.** Our GTH PPs are low-q — verified in `gth_potentials.json` LDA: **Na q1, F q7, Cs q1,
I q7** (Na/Cs also ship q9 semicore; F/I only q7). CP2K ships only q9 semicore Gaussian bases for Na/Cs and
**no GTH basis for iodine**, so it aborts on the valence mismatch. The fix is a matched low-q Gaussian valence
basis — which **GPW needs anyway** (GPW = Gaussian orbitals), so the work is shared.

**Include PW? YES — it is the basis-INDEPENDENT anchor, nearly free.** Our plane-wave code needs NO Gaussian
basis (orbitals ARE plane waves; only PP + Ecut) and already has NaF −20.3293 (Ecut=6) / CsI −11.3868 (Ecut=4)
[`606a54ff`]. Converging its Ecut gives the complete-basis limit. Three-way triangulation:
- **GPW vs CP2K** (SAME Gaussian basis + PP + functional) → IMPLEMENTATION correctness (the tight gate).
- **GPW vs PW** (Gaussian basis vs complete) → BASIS quality (the gap = Gaussian incompleteness; GPW ≥ PW in
  energy, i.e. less bound, as an incomplete basis under-binds).
- **PW vs CP2K** (both → complete-basis as CP2K's basis grows + cutoffs converge) → cross-code sanity.
PW is the leg that separates "is our GPW code correct" from "is the Gaussian basis good enough."

**Basis recipe (mirror `sipp.bsd`/`sipp_sr.bsd`).** Uncontracted even-tempered valence (one primitive per .bsd
shell, `nprim=1 coeff=1`), + a `_SR` variant dropping the most-diffuse primitive(s) for Bloch conditioning
(the SIPP→SIPP_SR lesson: ill-conditioning is a BASIS problem, [[feedback_scf_accuracy_levels]]). Valence
shells (from the PP q):

| el | q (Zion) | valence | shells | notes |
|----|----|----|----|----|
| Na | 1 | 3s¹ | s (+p polar) | 1 val e⁻ (alkali) |
| F  | 7 | 2s²2p⁵ | s+p | tight 2p → hard atom, higher cutoff |
| Cs | 1 | 6s¹ | s (+p) | heavy, diffuse 6s |
| I  | 7 | 5s²5p⁵ | s+p | **no GTH Gaussian basis anywhere** — first one; soft, big r_loc |

Seed α_max from the GTH `r_loc`, α_min from the valence ⟨r⟩, ratio ~2.5–3 (SIPP s = 2.0/0.7/0.25). New files:
`BasisSetData/{na,f,cs,i}_lowq{,_sr}.bsd` + `BasisSetData` enum entries + the loader map (mirror sipp/sipp_sr).

**Validation loop (per element → per compound).**
1. Build the `.bsd` (+ SR variant).
2. Finite pseudo-ATOM cross-check (the `SiPseudoAtomInBoxMatchesFinite` pattern): `Calculation(atom,
   {.basis=…, .pseudopotential=true})` converges, and GPW-in-box == that finite molecular DFT. Converge the
   basis by adding/tightening functions — NOT against Slater/High (different basis, a loose oracle: SIPP Si
   −3.759 vs Slater/High −3.337).
3. Transcribe the `.bsd` → CP2K `BASIS_SET` format (`El NAME`, nset, per-set `n lmin lmax nexp nshell` +
   exponent/coeff — the `UnitTests/CP2K/SIPP-SR-BASIS` pattern) + a CP2K deck (mirror `si_fcc_gpw*.inp`,
   `POTENTIAL GTH-PADE-q{1,7}`).
4. **Compounds:** NaF (rocksalt FCC), CsI (CsCl simple-cubic). Run **PW, GPW, CP2K**. Record Etot + runtime in
   `doc/CP2Kresults.md`; add did-E-move anchors: GPW → `GPW_SCF`, PW → `PlaneWaveDFTUT`.

**Multi-species GPW plumbing (small — the bases are the real work).** `Ham_PW_DFT` already has the multi-
species ctor (`{{"Na",1},{"F",7}}`, PW path `606a54ff`) and it drives GPW verbatim, so GPW multi-species =
thread the species→q map through `RunGPW`/`GPWFactory` in place of the single `element`/`q=4`. Ewald + the G=0
alignment are already per-atom (Zion per species); `MultiSpecies_Local/SeparablePotential` routers exist.
**DONE — multi-species GPW FIRST LIGHT (2026-07-11): NaF rocksalt Γ converges** (multi-species `Ham_PW_DFT`
ctor `{{"Na",1},{"F",7}}` on the generated `valence_lowq` basis, Na 5s2p + F 8s6p): 22 iters, **charge=8
conserved**, Etot=−25.086 (Enn=−14.00 = ionic Madelung, matches PW). Grid-underconverged (`densityEcut=40`,
Rcut=0) so not yet comparable to PW −20.3293. Gate `GPW_SCF.DISABLED_NaFRocksaltGamma` (~140 s: F's tight
40-a.u. exponent forces a fine density grid). Rcut=2a + SR basis (PSD overlap) → Etot=−23.556 (removes ~1.5 Ha
of the Rcut=0 over-binding).

**NaF cross-validation PINS (2026-07-11; full record: doc/GPWHistory.md):** both codes agree the answer FOR
THIS GAUSSIAN BASIS is ≈ −23.6 (CP2K transiently passes −23.64 vs our −23.556); the ~3.3 Ha gap to PW's
complete-basis −20.3293 is Gaussian-basis INCOMPLETENESS (the "GPW vs PW = basis quality" leg).  The SCF
instability root is the near-singular Bloch overlap METRIC (min eig 7.5e-4, cond≈8000 at SR), NOT occupation;
magnitude screening fixes the TRUNCATION artifacts but not intrinsic over-completeness → SR stays until §1.

**Gates / deliverables.** `doc/CP2Kresults.md` rows Si/NaF/CsI × {PW, GPW, CP2K} (Etot + runtime); `GPW_SCF`
NaF/CsI converge (charge, Etot) == CP2K same-basis; the GPW−PW gap documented (basis quality). **Pitfalls:**
iodine is the first GTH Gaussian basis for the element (validate its pseudo-atom carefully); F's tight 2p is
the hardest (needs the highest cutoff, per the PW NaF vs CsI experience — F set the cutoff, not the heavy I).

## 3. CP2K reference library (the oracle for §2) — BUILT; growing it
CP2K's Quickstep **is** the reference GPW implementation (Lippert–Hutter); its per-term breakdown points
straight at a bug (as this session's hand-rolled breakdown did: Een ×15.7 → local PP → the raster).
I can run CP2K directly: `~/Code/cp2k/build/bin/cp2k.ssmp`, decks in `~/Code/cp2k-runs/`.
- **DONE — CP2K 2026.1 built** (serial ssmp, gcc 15.2) at `~/Code/cp2k` (sibling to qchem6, outside the git
  tree). Toolchain: OpenBLAS+FFTW+libxc+libxsmm+DBCSR, no MPI/libint. Build: `tools/toolchain/build_cp2k.sh`
  (CMake, NOT the old arch-file `make`). Run needs `source install/setup` +
  `LD_LIBRARY_PATH=install/lib`.
- **DONE — FCC-Si Γ reference (SIPP_SR, GTH-PADE-q4, LDA_X+VWN5):** **Etot = −7.11506 Ha, charge 8**,
  converged by `CUTOFF` 80 Ry (≈40 Ha). Breakdown: Core-H (kin+PP) +5.565, Hartree +10.380, XC −2.544;
  PP total −7.548 (local −8.489, nonlocal +0.941); core self-energy −20.516. (CP2K's GPW electrostatic split
  differs from ours — compare the TOTAL + the cleaner sub-terms kin/XC/nonlocal-PP.) **Γ gate — MET** (−7.11506).
  Also Si **2×2×2 = −7.86744 Ha** (`si_fcc_gpw_222.inp`). Results table: **`doc/CP2Kresults.md`**; decks:
  **`UnitTests/CP2K/`**.
- **PP already aligned:** our `src/Pseudopotential/Data/gth_potentials.json` IS the CP2K GTH-PADE database
  (Si GTH-PADE-q4 params match ours exactly — verified). **Basis: same exponents, transcribed to CP2K
  `BASIS_SET` format** (uncontracted → one set per primitive; see `UnitTests/CP2K/SIPP-SR-BASIS`).
- **NaF/CsI:** the hand-rolled low-q bases + decks are now **§1's plan** (was "blocked"; the plan resolves it).
- **Si 2×2×2 cross-checks DONE + validated:** `si_fcc_gpw_222.inp` (shifted MP, **−7.86744** == our GPW after
  the complex-k fix) + `si_fcc_gpw_222_gamma.inp` (Γ-centred, **−7.77846**, matches our GPW −7.7778).

### Parameters to line up (qchem ↔ CP2K) — keep this table current
| quantity | qchem (ours) | CP2K keyword | note / pitfall |
|---|---|---|---|
| method | GPW | `&DFT &QS METHOD GPW` | (CP2K default is GPW) |
| cell | FCC primitive, a=10.26 a.u. | `&CELL` (A/B/C vectors, `BOHR`) | match lattice vectors exactly; `PERIODIC XYZ` |
| atoms | Si (0,0,0),(¼,¼,¼) frac | `&COORD SCALED` | match fractional coords (the corner atom at 0 is the bug trigger — compare it deliberately) |
| pseudopotential | GTH-LDA q4 (Zion=4) | `POTENTIAL GTH-PADE-q4` | same params (ours from CP2K) |
| orbital basis | SIPP_SR (3s3p, uncontracted) | `BASIS_SET` (our exponents, CP2K format) | convert file; keep it uncontracted |
| exchange | Slater/Dirac Xα=2/3 | LIBXC `LDA_X` | equivalent |
| correlation | **VWN5** | LIBXC `LDA_C_VWN` (=VWN5) | **NOT `PADE`** (that's PZ correlation) — must force VWN5 |
| density cutoff | `densityEcut` (Ha) | `&MGRID CUTOFF` (**Ry**) | **1 Ha = 2 Ry**; ours 8–12 Ha = 16–24 Ry is ~10× too low (CP2K default 300–600 Ry) — see TODO 1 |
| multigrid | single grid | `&MGRID NGRIDS`, `REL_CUTOFF` (Ry) | start `NGRIDS 1` to match; align `REL_CUTOFF` later |
| k-points | `MakeKMesh(shift)` (MP; shift=0 Γ-centred, shift=½ classic MP) | `&KPOINTS SCHEME MONKHORST-PACK` | CP2K's MP is SHIFTED (k=±¼ for even N) — use `kShift=½` to match; its Γ-centred list needs `SCHEME GENERAL` (see `si_fcc_gpw_222_gamma.inp`). CP2K prints its k-list (`grep BRILLOUIN`). Complex-k FIXED (`745d03ff`); the shifted gate awaits revalidation through the analytic kernels (§0a). |
| Poisson | G-space periodic | `&POISSON PERIODIC` | |
| spin | closed shell (Si₂, 8 e⁻) | `LSD` off | |
| occupation | integer, no smearing | `&SCF SMEAR OFF` | insulator |
| lattice-sum reach | AUTO (`Rcut<0`: radius from the basis + magnitude screening; `9714f58d`) | `EPS_PGF_ORB` / neighbour lists (auto) | both sides parameter-free now |
| **energy breakdown** | Ekin, Een, Eee, Exc, Enn, Ealign | CP2K: Overlap, Core (kinetic+PP), Hartree, XC, Core-self/Ewald | **term SPLITS differ** — match the TOTAL first, then map sub-terms (kinetic, XC are the cleanest to compare) |

## 4. Symmorphic space groups + auto-IBZ (deferred qcSymmetry track)
Only the POINT-GROUP part matters for k-reduction. Cubic Si (O_h, centrosymmetric) → IBZ = 1/48 of the BZ,
no time-reversal factor. Validate bit-level: reduced-mesh-with-weights == full-mesh (against Step 2). The IBZ
is an *efficiency* layer, not a correctness requirement — hence it comes AFTER a working full-BZ reference.

## 5. Deferred cleanups (do once bulk works — "the working code is the definitive declaration")
- **Rigorous periodic external PP:** `MakeLocalPP`/`MakeSeparablePP` quadrature the HOME-CELL orbitals against
  the cell's OWN atoms (no periodic-image PP) — exact at Γ / large box, an approximation for a dense crystal.
  Sum the PP over lattice images (analogous to Ewald / the PW G-space assembly).
- **DRY the PP field adapters into `qcPseudopotential`:** `RealYlm`/`BetaYlmField` are byte-identical in
  `PP_{Local,NonLocal}.C` (molecular terms) and replicated in the GPW evaluator. Hoist into a public module in
  `qcPseudopotential` (below both libs). Pure refactor; verify `L_PP` + `A_PP` + `GPW_SCF` unchanged.
- **`cMesh` = `Mesh<dcmplx>` (user-directed):** the `(Rs, phases)` pair (a `{R}` + `{e^{ik·R}}` weighted point
  set) and the density/quadrature grids should collapse to a `template<class W=double> class Mesh` — the
  integration algorithm is identical for real/complex weights, only the weight TYPE differs (confirmed vs
  `src/Mesh/Quadrature.C`). Then a `FourierMesh_R` ({R}) and `FourierMesh_k` ({k} + real BZ weights, unifies
  with today's `KMesh`). A cross-cutting refactor (Quadrature.C + bit-identity across ~29 consumers);
  currently marked with `// future: one cMesh` comments.
- **GGA Vxc fit grid (`relCutoff`) — CORRECTNESS for GGA, guarded now (`44bebe88`):** GPW uses ONE absolute
  `densityEcut` grid for both ρ (Hartree) and v_xc, and `GPW_IBS::CreateCD/VxcFitBasisSet` IGNORE `mp.relCutoff`
  (the CP2K REL_CUTOFF the Hamiltonian derives from the functional's `GridCutoffFactor()`; `PlaneWave_IBS` DOES
  honor it, building its Vxc grid at `Ecut*relCutoff`). LDA relCutoff==1 so it's exact — but a GGA's ∇ρ wants a
  DENSER v_xc grid. Fix = build a separate Vxc grid at `densityEcut*relCutoff`, mirroring the PW Vxc line. A
  guard `assert(relCutoff<=1)` now fires loudly on a GGA-on-GPW attempt instead of silently using the LDA grid.
- **Multi-grids + whole-density collocation — DONE** (the C+D analytic rewrite, see the DONE entry).
- **Common `PW_Evaluator`/`GPW_Evaluator` base:** the shared FFT engine (`PeriodicGridEvaluator`) is already
  factored; a base earns its keep once general-k k-space logic is shared.

---

# Durable pins / invariants (carry into all GPW work)
- **PP-smoothness is GPW's enabler; GAPW is out of scope (first pass).** All-electron cores are too sharp;
  validate with a well-conditioned GTH valence basis, never all-electron.
- **Use well-conditioned bases for SCF.** Ill-conditioning is a BASIS problem, not a solver/code bug (SIPP
  diffuse → SIPP_SR; N3/N5 removed). "LASolver" symptoms are basis conditioning. `N3/N5` no longer exist.
- **GPW is a Coulomb/Hartree STRATEGY orthogonal to the orbital basis** — a third one beside exact-4-centre
  (`Vee`) and density-fitting (`FittedVee`). Same `⟨χ|V_H|χ⟩` out, different internals.
- **Never assume `orbital == fit`.** Any fit/aux basis comes from the orbital basis via `Create{CD,Vxc}
  FitBasisSet(...)` — the factory is the seam even when trivial.
- **Fit quality is measured by grid-convergence of ρ, NEVER by ΔE_total** (the fit is non-variational).
- **Spin-polarized is the native formulation**; unpolarized is the ζ=0 collapse. New periodic terms
  spin-native (`FittedVxcPol`/`FittedVcorrPol`).
- **Regression style:** periodic/GPW energies are "did-E-move" anchors (pin the converged value, no
  `Converged()` guard). Where a real-space-on-lattice quantity must equal its finite counterpart, assert
  bit-consistency (`L_PP`-style) rather than an absolute oracle.
- **Two self-consistent schemes — do NOT mix:** (A) complete-Bloch analytic single-sum matrices (what GPW
  has, correct as Rcut→∞); (B) truncated-Bloch collocation Gram matrices (always PSD). Scheme-B overlap +
  scheme-A analytic kinetic gave `Ekin=−300`. Stay in scheme A at a converged Rcut (overlap PSD there).

### Symmetry comes AFTER a working GPW (independent optimisation layer, does not gate GPW)
Symmorphic space groups → BZ reduction (irreducible wedge) → SALC with plane waves. None of these gate GPW.

---

# Pointers
- **doc/GPWHistory.md** — the full archived DONE narratives, resolved investigations (indefinite-S,
  conditioning, NaF diagnostics), dead-end records, and complete commit archaeology.
- Tests: `UnitTests/GPW_UT.C` (1E + Bloch invariants; analytic collocation/adjoint gates;
  `AnalyticSeparablePPMatchesMesh` == mesh KB to 4.6e-11),
  `UnitTests/GPW_SCF_UT.C` (enabled anchors: `SiliconGammaConverges` == CP2K −7.11506 ± 2 mHa,
  `SiliconMultiKPlumbing` −7.45134, `SR_2x2x2ShiftedMP_vs_CP2K` == CP2K −7.86744 ± 3 mHa (the complex-k gate),
  `SiPseudoAtomInBoxMatchesFinite`; DISABLED: NaF, the Γ-centred 2×2×2 gate (redundant), conditioning sweeps),
  `UnitTests/L_PP.C` (finite==lattice PP), `UnitTests/PlaneWaveDFTUT.C` (PW anchors).
- CP2K decks + results: `UnitTests/CP2K/`, `doc/CP2Kresults.md`; CP2K itself: `~/Code/cp2k/build/bin/cp2k.ssmp`.
- Recent commits: **`8dba0625`** (C+D analytic rewrite, sampling deleted), **`9714f58d`** (auto-Rcut,
  budgeted stream cache, sharp-field PP ladder), **`9ff982ba`** (§0a Si leg: lockout fix + memos, complex-k
  gate enabled), **`b0f497c6`** (analytic KB + fp32 tier), **`4c71450c`** (D-aware radii + CP2K NaF oracle),
  **`35789164`** (NaF convergence: recipe machinery + fine-grid attractor findings).  Older: doc/GPWHistory.md.
- Build/test: `cd build/Release && ninja UTMain && ./UnitTests/UTMain`.
