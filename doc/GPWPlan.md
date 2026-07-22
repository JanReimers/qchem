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
- **§0b XC-consistency: FALSIFIED by the FD probe** (`f82db70e`): new gate `GPW.XCPotentialConsistencyFD`
  proves H_xc == ∂E_xc/∂D to FD accuracy (h² scaling to 2e-10) in both the smooth and ρ<0-guard regimes —
  the LDA discrete functional was already exactly consistent.  Full record below.
- **§0b′ gated ladder-completion rung + NaF ROOT CAUSE** (`a218c69c`): the top rung (energy-calibration-
  gated; `RelCutoffSafety` seam accessor; order-free `PairLevel`) + the D=S⁻¹ probe pins the NaF 4.9-e
  grid-charge loss on the Rcut=2a ENUMERATION-SCHEME MISMATCH (grid-independent −2.25 e; fp32 + screens
  vindicated at 7e-9/3.4e-7 per unit |D|).
- **BANISH-Rcut** (`bf3d70ad`): "there is no cut in R" — `(Rs,phases)` deleted from the seams, series
  ε-converged per shell pair inside `LatticeSum1E` (`ForImageOffsets`), KB convention simplified to the
  plain phase oracle, `Rcut`/`collRcut` gone from the GPW surface (finite mode = `CellImages` enum).
  Anchors identical, multi-k 123→84 s; NaF scheme mismatch DEAD (iter-1 charge −4.9 e → −2.4e-6 e);
  true conditioning exposed (λ_min=1.03e-6).
- **SR2 basis + instability CLASSIFIED** (`3f77c96e`): `valence_lowq_sr2.bsd` (the spectrum fingered the
  Na p 0.05 triplet; λ_min→1.57e-3, NaF 6× faster) — but the departure spikes SURVIVE: α-independent,
  DIIS-resistant, smooth growing mode from a clean fixed point ≈−27.73 → hypothesis = near-degenerate
  HOMO/LUMO at Γ (giant response).  The OPEN problem; full records below.
- **NaF Γ-instability MECHANISM MEASURED — band-gap instrument** (2026-07-17): new `ReportBandGap` flag
  appends ε_HOMO/ε_LUMO/gap to the verbose SCF line.  Verdict: the fixed-point gap is HEALTHY (~0.35 Ha,
  wide-gap insulator) so the static-degeneracy hypothesis is FALSE; the real mechanism is a giant-response
  DIFFUSE VIRTUAL whose ε_LUMO dives 0.2–0.5 Ha during the charge-transfer slosh, transiently crossing the
  occupied manifold (gap → 1e-4) → aufbau occupies it → +5–7e3 Ha spike (period ~27).  Records in §0b″.
- **NaF Γ-instability CURED (occupation-swap disease) — MOM wired up** (2026-07-17): the crystal's within-irrep
  fill (`TakeElectrons` = energy order) never touched the parked cross-irrep MOM, so MOM was wired into the
  irrep fill: `TOrbitals::TakeElectrons(ne, priority)` + `SCFParams::UseMOM`/`MOMStartIter` (threaded via a new
  `tSCFWaveFunction::SetMOM`) +
  **delayed IMOM** (aufbau for ~10 fills, then capture {F 2s, F 2p} ONCE and hold — running MOM drifts,
  iter-0 IMOM anchors the raw seed → both catastrophic).  NaF Ecut=40 now CONVERGES −27.76 (Δρ 6e-4, 196
  iters, partial-occ 0, diving virtual banished to −45 Ha unoccupied); vs CP2K oracle −27.93 the 0.17 Ha is
  the grid.  One residual iter-19 MIXING spike remains → 0c Pulay.  198/198 green (`SCFParams::UseMOM` off by
  default).
- **§0c SCF-STRATEGY REFACTOR + PULAY DONE** (2026-07-18, full design `doc/SCFStrategyPlan.md`): the SCF
  convergence machinery is now a role-seam framework — density-mixer seam (`tDensityMixer`: Linear/Kerker;
  bit-identical extraction `f4f48431`), loop-driver virtual dispatch replacing the `WantsLineSearch` mode `if`
  (`388b33d3`), and ONE shared paper-faithful `qchem.Math.DIIS` engine serving BOTH Fock-DIIS and density-Pulay
  (`c41f06f9`+`a60a04de`).  **`PulayMixer`** (Kerker-preconditioned density-DIIS, priming via `PulayStart`)
  accelerates NaF Ecut=40 **196→63 iters** to the SAME −27.756.  Flexed via `scfrun` (which grew a molecular
  `--mol` mode + a fixed SCFParams misalignment): Boron 16 / Sc 21 / O2 triplet 13–15 iters; O2-HF-triplet
  display SEGV fixed (`78b8f66a`).
- **DIRECT FINE-GRID NaF FALLS INTO THE −39 BASIN** (2026-07-19, `30d0eb87`; MOM+Pulay, auto Ecut=160, 15m45s):
  "converges" (Δρ 2.9e-5) to E=+54.3 garbage — the Kerker descent goes STRAIGHT into −39, Pulay thrashes on it.
  So the production-grid failure is a DENSITY/GRID-basin problem, NOT occupation/mixing: MOM+Pulay necessary
  but NOT sufficient.  Next-session plan (grid-continuation seeding + basin removal + OpenMP, basin kept as a
  test fixture) recorded in the TODO §0e below.
- **§0e-PP CP2K local-PP split + Q1 grid speedup** (branch `gpw-0e-pp-local-split`: `94544683` split, `83d827b9`
  Q1; 202/202 green).  Local PP split at the `LocalPotential` form-factor level — LONG (softened-Coulomb → folded
  into `PW_Hartree`'s G-space Poisson) + SHORT (poly×Gaussian → external `PW_Pseudo`); `FormFactorLong` primary,
  base provides `FormFactor=Long+Short` (Design A); a matrix-identical ENERGY-RELOCATION refactor (Si Γ −7.11506
  + NaF −27.756 held).  **Q1 — the ~295 s NaF fine-grid `MakeLocalPP` setup wall is the `relCutoffScale`, over-set
  to 6 by the DENSITY SCREEN** (the increment-1 `−280`/`−259` was `OverlapMatrix`'s `screenD` zeroing off-diagonals
  of the FIXED `V_long`, NOT aliasing; unscreened, smooth==stiff to 4e-3 for soft Si).  Default 6→3 = ~2× (Ecut=160
  578 s→128 s @scale 2), all gates green (Si Γ now 31 s); env knobs `GPW_LOCALPP_SCALE`/`GPW_LOCALPP_FULL` for the
  later 2/4 verify.  The ANALYTIC V_local (short BUILT+finite-validated but dormant; long = the Ewald crux) is a
  SEPARATE accuracy upgrade → **TODO §0e-PP** (re-gates to converged CP2K −27.93).
- **§0e step 1 (grid-continuation) + step 2 (XC-collapse ROOT-CAUSED & FIXED)** (2026-07-20, branch
  `gpw-0e-pp-local-split`; commits `75e1d4c8`,`758b92a8`,`c816cb39`,`10a91a1e`,`1e13df74`,`4e84284c`,`b65e4185`,
  `a7769f81`; Si Γ −7.11485 bit-identical, adjoint machine-exact, 28 PW anchors green).
  **Step 1 — grid-continuation seeding**: explicit-density seed ctor (`tSCFIterator`, enum ctor delegates) +
  `AdoptMOMReference` (transfer the converged coarse WF's occupied subspace as the fine MOM reference — the
  density seed ALONE gave a wrong −23.3).  Avoids (not removes) the −39 basin.
  **Step 2 — the fine-grid Exc collapse (NaF pinned −24.4 vs oracle −27.93) ROOT-CAUSED**: F's tight density
  (product α≈80) UNDER-RESOLVED on the fit grid → the collocated ρ aliases into huge negative lobes (`negCharge
  −9.3 e`, `neg-frac 0.50` even at the converged fixed point; Si clean at 0.08 %) → the XC `ρ>0` guard's
  grid-sensitive interaction collapses Exc.  **FIXED** by (a) the **fit-grid thread-through** (`MakeRepulsion3C`/
  `MakeOverlap3C` build the tensor over the REQUESTED fit basis's grid, not the block's own), (b) the
  **`Overlap3C` ADJOINT** (`G_ERI3::applyAdjoint` + `ContractAdjointG_ERI3`: the KS matrix `⟨i|v|j⟩=Σ_k
  v-tilde(G_k)⟨i|e^{iG_k}|j⟩` is the BACKWARD contraction of the same tensor, carrying the fit grid — killing the
  grid-less `MakeOverlap(field)` that silently used the coarse grid), and (c) **density-fit densification** =
  making `cutoffFactor` big enough to resolve the product (**4→8**; measured F: `negCharge −9.3→−0.03 e`, clean
  SCF −26.198).  **ONE-GRID cleanup**: `itsGrid`→`itsFFT_R_G_Grids`, `GPW_CDFIT_SCALE`+second grid deleted,
  `CreateCD/VxcFitBasisSet` self-documenting (`{G}_ρ`=`DensityGrid()`, `{G}_vxc`=relCutoff·`{G}_ρ`).  Bugs fixed
  en route: `CollocMemo` grid-collision segfault (`c816cb39`), `RhoOnGrid` out-of-band aliasing guard (`10a91a1e`).
  **HONEST PICTURE**: the old −27.75 was an ALIASING COINCIDENCE; the resolved answer is −26.198, the gap to the
  oracle now the still-coarse **local-PP** base grid, NOT the density.  REMAINING → **TODO §0e step 2** (the
  grid-matched CP2K validation + the grid/exponent diagnostics; the local-PP resolution).

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
> ⚠ SUPERSEDED by §0b′ (same day): this whole sub-block was measured on the CORRUPTED map (Rcut=2a
> scheme mismatch).  The −39 attractor, the −27.73 pin, and the "quasi-Newton is the missing ingredient"
> conclusion are ALL corrupted-map artifacts — see §0b′ for the honest map (mismatch deleted; the real
> blocker is the Γ giant-response instability, TODO §0b″).  Kept for the archaeology only.
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

## §0b XC CONSISTENCY — RESOLVED BY FALSIFICATION (2026-07-16).  The full record:
**The fork does NOT exist; the LDA discrete functional is ALREADY exactly consistent.**  The probe is the
new gate `GPW.XCPotentialConsistencyFD`.
**The instrument came first (as this section prescribed) and overturned the premise.**  The probe replicates
the PW_XC chain verbatim at the evaluator seam (collocate → nested {G_L} combine → `RhoOnGrid`; pointwise
v_xc → raster `ForwardFFT` → per-level restriction → analytic `IntegratePotential`) and compares the central
FD `[E_xc(D+h dD)−E_xc(D−h dD)]/2h` against `Re Tr(H_xc dD)` on the FCC-Si crystal (cross-cell pairs + a
real ladder), with a bilinear Hartree control:
- **Positive-density regime: rel err 8.0e-8 (h=1e-3) → 2.0e-10 (h=1e-4) — exact h² scaling, i.e. pure FD
  truncation converging onto the analytic answer.  H_xc IS ∂E_xc/∂D.**  Hartree control 3e-10.
- **Indefinite-D regime (ρ_q<0 over part of the grid — the Kerker-mixed-field case): same h² scaling
  (5e-6 → 3e-8).  The ρ≤0→0 guards are CONSISTENT between E and H** (both SlaterExchange AND
  VWN_Correlation already guard `rho>0.0 ?` — the "only SlaterExchange has the guard" worry was stale).
Why the old fork description was wrong: `PW_XC::GetEnergy` already takes ∫ε_xc·ρ on the fit grid (the
¾-virial survives only as `ExFunctional::GetEpsXc`'s default, EXACT for Dirac; VWN overrides), and the
"band-limited fit" of v_xc is the fine-grid projection onto the fit ball — which is EXACTLY the gradient of
the grid-sum energy w.r.t. the ball-limited ρ̃ the energy itself uses.  One discrete functional, end to end;
`FittedEpsXc` is molecular-path-only and was never on the periodic route.

**Consequences (re-scope):**
- The NaF fine-grid attractor (E≈−39, Exc≈−143) is a **GENUINE basin of the (under-resolved) discretized
  functional**, not a consistency artifact: mid-slosh D loads the sharpest F-F pairs beyond the grid
  calibration → collocated ρ aliases (∫ρ_grid swings 5.1↔7.7 vs Tr(DS)=8, spiky/locally-negative) → E_xc
  is legitimately huge-negative WITHIN the discretization, and since H_xc is its exact gradient, the SCF
  map is self-consistent there.  A variational minimizer (GDM/OT) would find it too — the escape is not
  consistency but (a) never wandering into the basin (quasi-Newton mixing with small steps = what CP2K's
  Broyden does on the same map; grid-continuation seeding = start in the physical basin) and/or (b)
  removing the basin by resolving the sharp pairs (stiffer fine-grid calibration; CP2K leaks only 2e-4 e
  at the same 160 Ha — understand its EPS_RHO/REL_CUTOFF stiffness if (a) is not enough).
- **ρ-FLOOR: already effectively present for LDA** (both functionals zero at ρ≤0, verified consistent by
  probe 2).  An explicit ε-floor remains only as the **GGA prerequisite** (∇ρ/ρ powers diverge at tiny ρ)
  — fold it into the GGA increment together with the `relCutoff` Vxc-grid item (§5).

## §0b′ TOP RUNG + NaF ROOT CAUSE + BANISH-Rcut + SR2 (2026-07-16, one session).  The full records:
**Two separate things came out of this increment: the ladder-completion rung (LANDED, small-but-real energy
fix, decision pending on scope) and the ACTUAL root cause of the NaF grid-charge catastrophe (an
ENUMERATION-SCHEME MISMATCH — not grids, not precision).  The instruments: `GPW.SharpestPairChargeConservation`
+ `GPW.DISABLED_IllConditionedChargeProbe`.**

**(1) The top rung — LANDED (code in tree), measured, scope decision pending.**
`PairLevel`'s requirement `req = kRelSafety·ecut_fine·(αᵢ+αⱼ)/(2α_max)` is unsatisfiable for pairs with
αᵢ+αⱼ > α_max; one rung at `RelCutoffSafety()·ecut_fine` (appended LAST — `ecut_L[0]` STAYS the resolution
reference, selection made order-free; new seam accessor `LatticeSum1E::RelCutoffSafety`) completes the
ladder by construction.  The local-PP path (relCutoffScale=6) keeps the BASE sub-ladder (`itsNBaseLevels`)
— its stiffened rule would flood the doubled grid with mid pairs.  Machine gates (adjoint, FD-consistency,
charge) all carry over.  MEASURED: the rung is an ENERGY-tail fix ONLY —
- CHARGE is rung-INVARIANT (~1e-9 with or without): the G=0 coefficient survives ball truncation by
  construction, and pow2-padded rasters keep box sampling at ~e^{−50}.  (The gate documents this.)
- Si anchors (explicit Ecut=20 = 2.5× their auto floor): moves SUB-mHa (Γ −7.11485→−7.11482, shifted
  −7.86724→−7.86713 — all within existing gates, no re-pin forced), cost 1.6–4× (the global N³ work:
  Γ 29→48 s, shifted 167→430 s, atom-in-box 25→107 s).
- DECIDED (user, 2026-07-16): **GATED on the energy calibration** — the rung is added only when the
  reference grid sits below `RelCutoffSafety()·cutoffFactor·α_max` (every AUTO run gets it; the Si anchors'
  explicit Ecut=20 ≥ 16 skip it and return to baseline speed).  ALSO NOTED: the auto-floor
  `cutoffFactor=4` calibration ("Ecut=40 loses >5 e⁻ of F") is SAMPLING-ERA data (2026-07-12, pre-analytic-
  rewrite) — the analytic path conserves charge at ANY Ecut, so the production Ecut may be recalibratable
  DOWN by an ENERGY criterion (a large runtime lever that also shrinks the rung's cost).

**(2) NaF iteration-1 grid-charge loss ROOT-CAUSED = ENUMERATION-SCHEME MISMATCH (the "two schemes" pin,
violated by the NaF config itself).**  The probe (D=S⁻¹: PSD, Tr(DS)=n EXACT, entries ~1/λ_min — the
loading a mid-slosh SCF produces):
| error source | measured | per-unit-\|D\| |
|---|---|---|
| **Rcut=2a-truncated S vs screened-complete collocation** | **−2.247 e at \|D\|=450, GRID-INDEPENDENT** (identical Ecut=40 vs auto=160, across fp32 tiering) | 5e-3 |
| kScreenEps screening tails | −0.36 e at \|D\|=1.05e6 | 3.4e-7 |
| fp32 stream tier | ~7e-3 e at \|D\|=1.05e6 | 7e-9 |
- The collocation enumerates its cross-cell offsets INTERNALLY to the complete magnitude screen
  (VALENCE_LOWQ_SR α_min=0.0857 → pair reach ≈33 au), while the NaF SCF builds S over `Rcut=2a`=17.5 au —
  S/charge/diagonalization live in the TRUNCATED scheme, ρ̃/Hartree/XC in the COMPLETE one.  Mid-slosh D
  loads the near-null (diffuse) directions where truncated-S is most wrong → the e-scale ∫ρ−Tr(DS) swings
  (iter-1: 4.9 e), a corrupted SCF map, and (plausibly) the −39 basin.  NOT fixable by mixing (0c) or by
  grids (rung) — the map itself is inconsistent.
- At AUTO Rcut the mismatch vanishes (err/|D| ÷15000) BUT the complete-enumeration S is genuinely
  near-singular: **λ_min ~ 1e-6** (|S⁻¹|~1e6).  The 2a truncation was double-dutying as a conditioning
  crutch (the SR .bsd header even says "PD at a MODEST Rcut").  The precision machinery is VINDICATED
  (fp32 + screens hold their calibrations even at million-scale loading).
**→ BANISH-Rcut — THE REFACTOR LANDED SAME DAY (2026-07-16, in-tree; user directive, attempt #4 — this
time with the crutch measured to corrupt the map).  STATUS + measurements:**
- **The `(Rs, phases)` arguments are GONE from `LatticeSum1E`**: the 1E/KB builders take
  `(cellphase_t, UnitCell)` and sum their series to ε internally per shell pair via `ForImageOffsets`
  (the collocation kernels' own exact-threshold screen — 1E and collocation are now ONE scheme by
  construction).  New finite `MakeOverlap(g)` overload for the home-mode KB.  The KB phase convention
  SIMPLIFIED: with internal symmetric enumeration the historical `(−Rs, conj-phase)` artifact reduces to
  the PLAIN phase oracle (m=−n substitution) — validated by `AnalyticSeparablePPMatchesMesh` AND the
  complex-k shifted-MP anchor (−7.86724, bit-identical).
- **`Rcut`/`collRcut` DELETED from `GPW_Evaluator`/`GPW_IBS`/`GPWFactory`** (`itsR/itsPhase/rcutEff` gone;
  the only remaining image list is the INTERNAL ε-derived Eval/mesh-KB set).  The finite mode is now
  `CellImages::HomeCellOnly` (an `enum class` so a stray numeric can never silently select a mode); its 1E
  matrices are the finite molecule's own cached faces, widened — the home-cell gates run in MILLISECONDS.
- **Anchors: Si Γ −7.11485 / 2×1×1 −7.45133 / shifted −7.86724 — IDENTICAL to pre-refactor**; multi-k
  RUNTIME improved 123→84 s (the exact-threshold enumeration is leaner than the old conservative ball).
  14/14 GPW kernel gates green (adjoint, FD-consistency, charge, KB==mesh).
- **NaF AT COMPLETE ENUMERATION (the measurement this was for): the scheme mismatch is DEAD.**
  Iteration-1 diagonalized-density grid charge: **−4.9 e → −2.4e-6 e**.  True conditioning measured:
  λ_min(S)=1.03e-6, cond=6.0e6 — and **Cholesky survives** (the SCF runs).  The IONIC SEED still loses
  1.09 e (seed construction on the near-singular basis hits the |D|-amplified precision floors — one
  iteration only; the diagonalized densities are clean).  Setup share grew (~28 min to iteration 1 at
  auto Ecut: bigger streams + the 5-level ladder) — the 0d OpenMP/setup item.
**The Ecut=40 recipe measurement ON THE HONEST MAP (α=0.025/G0=1, 200 iters) — the verdict:**
- **The map is healthy and has a genuine fixed point ≈ −28.00**: after the seed transient the SCF descends
  SMOOTHLY and monotonically (−27.64→−27.9999 over ~29 iterations, repeatedly), grid charge −3e-3 clean
  throughout (no slosh, no mismatch).  Note −28.00 vs the corrupted-map "anchor" −27.73 and the CP2K
  320-Ry oracle −27.93128 — the old "0.2 Ha leaky-grid gap" attribution was itself a corrupted-map
  artifact; the honest Ecut=40↔oracle comparison awaits actual convergence + the production grid.
- **The ONE remaining disease is the NEAR-NULL OCCUPATION EVENT** — and (user challenge, answered) it is
  NOT the linear algebra: cond=6e6 costs Cholesky ~7 of 16 digits, V=S^{−1/2} amplifies 10³ — all exact
  enough.  The instability is the RAYLEIGH QUOTIENT of the near-null state: ε_null = vᴴFv/vᴴSv is a ratio
  of two near-zeros, sensitive as δε ≤ ‖δF‖/λ_min.  The LEGITIMATE per-iteration Fock update during the
  Kerker descent is ~1e-2 (Δρ≈2e-2) and projects strongly onto v (the null combination is built of the
  same diffuse functions the mixed V_H/v_xc fields move), so the spurious band sweeps up to 1e-2/1e-6 =
  1e4 Ha per iteration; when its trajectory carries it below the Fermi edge, AUFBAU OCCUPIES IT (a
  1/√λ≈10³-amplitude vector enters D) → E=+1e4, [F,D] 0.12→150.  Signature: DETERMINISTIC period ~29
  (six spikes, iters 45/74/103/131/160/185 — trajectory-driven, not noise), discontinuous [F,D] (an
  occupation swap), charge CLEAN throughout.  This is the classic QC near-linear-dependence collapse —
  molecular codes drop S-eigenvalues below 1e-6..1e-8 for exactly this reason; the criterion that matters
  is ‖δF‖/λ vs the gap, not cond(S).  Previously MASKED by the 2a truncation (λ_min 7.5e-4).  The old
  ±75 Ha limit cycles / the −39 attractor / the "+900 DIIS spikes" all belong to the corrupted map; the
  NaF energy pin is SUSPENDED in the test until the near-null fix lands.  (Fix menu: the basis trim /
  rank-reduction below; occupation control (level shift / MOM) would stabilize around the garbage band
  but leaves ε_null polluting the band structure — not the clean fix.  Verification instrument for the
  SR2 session: print the lowest band energies per iteration — the spurious level should dive across the
  Fermi edge one iteration before each spike.)
**SR2 TRIM — DONE same session (`valence_lowq_sr2.bsd`, enum VALENCE_LOWQ_SR2): drop Na p 0.05 + s 0.0857
(the SPECTRUM identified them: SR's three degenerate 1.03e-6 near-null modes = exactly the Na p 0.05
triplet; F kept intact for the anion).  λ_min 1.03e-6 → 1.57e-3 (cond 2715, Cholesky residual 4e-14);
NaF 200 iters 17 min → 3 min (the deleted diffuse shells owned the biggest boxes).  BUT THE SPIKES
SURVIVED — the conditioning/near-null diagnosis is DISPROVED as the mechanism (measurement-driven, round 3):**
- **α-INDEPENDENT**: 10/10/13 spikes at α=0.025/0.0125/0.00625 (period ~27; smooth descent to the SAME
  ≈−27.73 fixed point each cycle, then a SMOOTH climb-away over ~5 iters before the +5e3-scale blowup —
  a growing departure, not a discontinuous occupation swap).  Rules out the plain linear-mixing gain
  story UNLESS the response multiplier is ~1e4 (α·|λ|≫1 even at α=0.006).
- **DIIS (quasi-Newton in Fock space, `NAF_DIIS=1` knob) does NOT fix it** on the honest map — 51
  excursions, no smooth descents, endpoint −27.1±2 with `En>EMax` flapping.  (Its ban was for the
  corrupted map; on the honest map it fails DIFFERENTLY — fighting the same mode.)
- **The surviving hypothesis: a GIANT RESPONSE MODE from a near-degenerate HOMO/LUMO at Γ** — a tiny gap
  makes χ ~ 1/(ε_v−ε_c) huge: explains the α-independence at practical α, the smooth departure, CP2K's
  OWN eternal density limit-cycle on this same system (RMS 0.03–0.12 forever), and DIIS's failure.
  Γ-only NaF in this minimal ionic basis SHOULD be wide-gap — if the measured gap is tiny, that itself
  is the finding (basis? PP? Γ-only folding?).
Side effect of the refactor: the "(Rs,phases)→one cMesh" future note is MOOT for these seams (no weighted
point set crosses the interface — the stronger form of that cleanup); KMesh + quadrature meshes keep it.

---

## §0b″ NaF Γ-INSTABILITY — mechanism MEASURED + occupation-swap disease CURED by MOM (2026-07-17).  The full record:
**The classified facts (records in DONE §0b′): the honest, conditioned map descends smoothly to its fixed
point and departs via a GROWING mode — α-INDEPENDENT (10/10/13 spikes at α=0.025/0.0125/0.00625, period
~27, smooth climb-away over ~5 iters), NOT conditioning (SR2 λ_min=1.6e-3 shows the same spikes as SR
1.03e-6), NOT DIIS-fixable (`NAF_DIIS=1`: 51 excursions, `En>EMax` flapping).**

**1. BAND-GAP INSTRUMENT — DONE (2026-07-17, `ReportBandGap` flag on the verbose SCF line; extracts
ε_HOMO/ε_LUMO/gap from `wf->GetEnergyLevels()`, flags a partially-occupied frontier).  The hypothesis is
REFINED, not simply confirmed — the mechanism is now directly visualized (Ecut=40/α=0.025, `GPW_SCF_UT`):**
- **The FIXED-POINT gap is HEALTHY: ε_LUMO−ε_HOMO ≈ 0.33–0.37 Ha (~9–10 eV).**  NaF/Γ in this ionic basis
  IS a wide-gap insulator at convergence (the plateaus iters 30–37, 55–66 sit at gap ≈ 0.35).  So the
  *static* near-degenerate-HOMO/LUMO version of the hypothesis is **FALSE**.
- **Each spike is preceded ONE iteration earlier by ε_LUMO DIVING 0.2–0.5 Ha** — a diffuse virtual with a
  GIANT RESPONSE to the low-G charge-transfer slosh.  The gap collapses (iter 12 → 2.8e-2, ε_LUMO crashing
  +0.167 → −0.077; iter 68 → 1.2e-4 with ε_H/ε_L DEGENERATE) as the virtual crosses the occupied manifold;
  then AUFBAU fractionally occupies it (`[partial-occ HOMO]` fires exactly on the spike iters 14, 41) → a
  ~1/√λ diffuse vector enters D → E=+5e3…+7e3 Ha, [F,D] 0.09 → 130.  Deterministic period ~27 (spikes
  14/41/68 in one run) — matches the classified fingerprint exactly.
- **→ mechanism = a giant-response DIFFUSE VIRTUAL causing a periodic aufbau LEVEL-CROSSING, NOT a small
  static gap.**  (The "χ ~ 1/(ε_v−ε_c)" framing was close but the small denominator is TRANSIENT — created
  by the slosh, not intrinsic; the transition density onto the diffuse virtual is what makes the response
  giant.)  Also explains CP2K's eternal density limit-cycle (RMS 0.03–0.12) on this same system.
- **FRONTIER-WINDOW refinement (2-occ/4-virt window per iteration) — two sharper facts:**
  - **it is ONE ISOLATED hyper-responsive virtual, not a wide-band cluster.**  At the dive (iter 11→12)
    the LUMO crashes +0.167 → −0.077 (0.24 Ha in ONE step) while its virtual NEIGHBOURS (+0.42, +0.79)
    barely move, and at the plateau the LUMO sits ~0.25 Ha clear of the next virtual.  So the giant
    response is a *single* diffuse (Na-3s-like) conduction state that overlaps the charge-transfer region
    — NOT an over-complete diffuse-band cluster.  This argues **3b (physical-but-responsive), not 3a
    (basis ghost)** — and de-prioritises the §1 rank-reduction angle for this instability.
  - **the spike IS an OCCUPATION SWAP — this CORRECTS the §0b′ "growing mode, not a swap" note.**  At each
    spike the F 2p level drops from (6.0) to (4.0) electrons: the diving virtual captures 2 e out of the
    F 2p manifold (iters 14, 41).  The smooth dive (the "growing mode") TERMINATES in the aufbau swap —
    they are two phases of ONE event, not alternatives.  → **MOM (pin the {F 2s, F 2p} occupied subspace)
    is the direct fix**, and should be clean because it is an isolated single-state swap.

**2. MOM FIX — WIRED UP + VALIDATED (2026-07-17; `SCFParams::UseMOM`/`MOMStartIter`).**
- **NOT Fermi smearing / not the "gap≈0" branch** — the fixed-point gap is large, so there is no static
  degeneracy to smear; smearing would leave a residual fractional-occupation error at a wide-gap insulator.
- **The measured mechanism is a clean, isolated, single-state OCCUPATION SWAP (F 2p 6 e → 4 e)** — exactly
  what MOM prevents.  The parked MOM machinery (`tIrrepWF::MOMScores`/`CaptureMOMReference`) lived ONLY in
  the molecular cross-irrep aufbau (`tCompositeWF::FillOrbitalsAufbau`), which the crystal never runs (a
  crystal k-block is a fixed-EC single irrep filled by `TakeElectrons` = pure energy order).  So MOM was
  wired into the **within-irrep** fill: new `TOrbitals::TakeElectrons(ne, priority)` (occupy highest-overlap
  first), driven from `tIrrepWF::FillOrbitals`; the per-run knobs `SCFParams::UseMOM`/`MOMStartIter` threaded
  through `tSCFWaveFunction::SetMOM` (+ activation on a captured
  reference, NOT on the accelerator engaging — NaF's Null accelerator never engages).
- **The reference-capture POLICY is the whole game (two wrong variants measured + rejected):**
  RUNNING MOM (re-capture every iteration) DRIFTS — a spike corrupts the reference, MOM then locks a
  +0.74 Ha level occupied while a −50 Ha level stays empty → wrong −24.4.  IMOM-from-iter-0 anchors the RAW
  SEED (mid-transient, shapes still shifting) → catastrophe (+5 Ha occupied, −112 Ha empty).  **DELAYED
  IMOM WINS** (`MOMStartIter`, default 10): plain aufbau for ~10 fills to descend to the physical fixed
  point, THEN capture {F 2s, F 2p} ONCE and hold it fixed.
- **RESULT: NaF Ecut=40 now CONVERGES.**  Occupation swaps VANISH (partial-occ count 0), the diving virtual
  is banished (−45 Ha, UNOCCUPIED), and the SCF descends SMOOTHLY+MONOTONICALLY to **−27.76** (Δρ 6e-4 at
  150 iters; gap 0.50 Ha) — the physical fixed point the spiking run only ever visited transiently.  vs the
  CP2K oracle −27.93128 at 320 Ry, the ~0.17 Ha is the Ecut=40 grid.
- **ONE residual excursion survives (iter ~19) — but partial-occ 0, so it is NOT an occupation swap: a
  density-MIXING transient (the charge-transfer slosh).  → 0c Pulay/Broyden is the next lever** (damp the
  slosh; also accelerate the slow linear-Kerker tail; matches CP2K's Broyden on this map).  MOM and 0c are
  complementary: MOM stops the swap, Broyden stops the slosh.
- Open sub-question (de-prioritised): WHICH diffuse virtual dives?  A single Na-3s-like state, not an
  over-complete cluster — ties loosely to §1 but MOM makes it a spectator, so §1 stays a curiosity here.
3. **0c (Pulay/Broyden mixer face)** on the conditioned map, now the LEAD remaining item (kills the residual
   iter-19 mixing spike + accelerates the tail); its `MixSignals` trust-region signal (∫ρ_grid − Tr(DS))
   stays — now purely a precision/conditioning health meter.  Also probe: the ionic SEED's 1.09-e
   precision-floor loss (may already be gone with SR2's conditioning).

## §0c PULAY/BROYDEN ρ̃-MIXING — the mixer face + shared DIIS engine landed (2026-07-18; design in doc/SCFStrategyPlan.md).  The full record:
> **SUPERSEDED/EXPANDED by `doc/SCFStrategyPlan.md` (2026-07-18)** — the mixer is one seam of a four-role
> ISP model (orbital / occupation / density / loop) with a single shared extrapolator (DIIS≡Pulay, one
> paper-faithful engine on either the F or ρ residual stream) and an occupation seam that extends to Fermi
> smearing.  Read that doc for the design + increment plan; the sketch below is retained for context.

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

---

# TODO / NEXT

**Orientation (2026-07-19, end of session).**  Everything through §0c is **DONE** — §0 through SR2, §0b″
(band-gap instrument + MOM cure) and §0c (the SCF-strategy refactor: mixer seam, loop-driver, ONE shared DIIS
engine, and Kerker-preconditioned Pulay) now sit as full records in the [DONE](#done) section above; §0c
design in `doc/SCFStrategyPlan.md`.  NaF Ecut=40 converges (MOM+Pulay, 63 iters, −27.756).  **The ONE
remaining NaF problem is the PRODUCTION GRID (§0e below): the direct auto-Ecut=160 run falls into the −39
density/grid basin — MOM+Pulay are necessary but not sufficient, so grid-continuation seeding + basin removal
(+ OpenMP to make iteration bearable) is the next-session critical path.**  Then the runtime follow-ups (0d)
and the standing queue (1)–(5).

## 0e. NaF PRODUCTION GRID — the one remaining NaF problem (NEXT, critical path)

**DIRECT FINE-GRID RUN MEASURED — 2026-07-19 (MOM + Pulay depth6/start35, auto Ecut=160, 45527 G, 15m45s,
NMAX=100): FAILS to the unphysical basin; grid-continuation seeding is now the CRITICAL PATH, not just an
accelerant.**  The run "converges" (Δρ=2.9e-5, 90 iters) but to E=+54.3 (εH=92/εL=139, Eee=+152/Exc=−137 =
the aliased/negative-ρ garbage breakdown).  The trajectory is the smoking gun: the **Kerker priming descent
goes STRAIGHT into the −39 basin** (iters 20→34: −24→−39.85, smooth), then **Pulay engaging on that garbage
state thrashes** (+45/+102/…) to +54.  Verdict: the fine-grid failure is a DENSITY/GRID-basin problem, NOT
occupation (MOM keeps occ sane) and NOT mixing (Pulay only accelerates — it can't escape a basin, and on the
pathological −39 map it destabilises).  MOM+Pulay are necessary but NOT SUFFICIENT for the production grid.
→ NEXT: implement grid-continuation seeding (converge Ecut=40 physical −27.76 → seed the fine grid with THAT
ρ → start in the physical basin, never wander into −39); and/or (b) stiffen the fine-grid calibration to
REMOVE the basin (CP2K leaks only 2e-4 e at 160 Ha — understand its EPS_RHO/REL_CUTOFF stiffness).

**AGREED PLAN for the next session (user, 2026-07-19) — keep the −39 basin as a TEST FIXTURE; both fixes are
TOGGLEABLE options so the default (ionic seed + current grid) still exposes it, and each fix is verified with
the OTHER turned OFF:**
- **Step 0 — OpenMP over the collocate/integrate pairs — DONE (2026-07-19), but the fine-grid win is smaller
  than hoped; the real lever is now the SETUP (0d).**  `PG_Cart_MnD::NR_Evaluator::CollocateDensity`
  (per-thread private ρ accumulators + a `critical` reduce) and `IntegratePotential` (write-independent per
  pair — no reduction) are OpenMP-parallel over the flattened `(i,j)` pair list.  **Opt-in at runtime via the
  env knob `GPW_OMP_THREADS` (>1; default 1 = serial), NOT `OMP_NUM_THREADS`** — because the Si anchors and a
  threaded NaF run share one UTMain binary and cannot be separated by a global harness pin (the same reason as
  the NAF_*/GPW_ILLCOND_ECUT knobs), so no harness pin was needed (serial by default keeps the anchors
  byte-identical; 201/201 UTMain green).  Toolchain: this LLVM install ships no libomp → `-fopenmp=libgomp`
  (which honours the pragmas but does NOT define `_OPENMP`, so the code gates on our own **`QCHEM_OPENMP`**
  macro); no `<omp.h>` (private-buffer + critical pattern).  See [[project_openmp_runtime]].
  - **MEASURED (NaF fine grid, auto Ecut=160, 4 threads):** per-iteration collocate/integrate **~10.4 → ~6.1
    s/iter ≈ 1.7×** — the per-iteration cost is memory-BANDWIDTH-bound (scatter/gather replay over the cached
    streams), so 4 cores only buy ~1.7×.  Threading confirmed engaged (4 workers running simultaneously at
    Ecut=40 via /proc); charge/energy correct under threads.
  - **THE FINE-GRID WALL IS THE SETUP — AND IT IS `MakeLocalPP`, NOT `EnsureStreams`, AND IT DOES NOT
    PARALLELISE (profiled 2026-07-19, SR2 basis; a per-phase `std::chrono` breakdown + a threaded A/B).**  The
    `cSCFIterator` ctor is **~320 s** of the fine-grid run; inside it: overlap-S is instant, the `EnsureStreams`
    stream build is only **~25 s** (129.5M pts for SR2 — the old "~950M / 5-min" figure was a pre-SR2 basis),
    and the remaining **~290 s is `MakeLocalPP`** (the `relCutoffScale=6` static local-PP sweep in the
    iteration-0 Fock, `GPW/Imp/Evaluator.C:484`).  A 0d attempt to OpenMP `EnsureStreams` (parallel per-pair box
    eval + `critical` budget tiering) was implemented, verified byte-identical serial — and then **REVERTED
    because it gave ZERO speedup** (25 s → 25 s at 4 threads).  Same for `MakeLocalPP` through step 0's
    `IntegratePotential` path (ctor 329 s → 318 s = noise).  ROOT CAUSE: both are dominated by a **few
    ultra-diffuse pairs with enormous boxes** (`MakeLocalPP`'s own comment: "an ultra-diffuse pair's box on N=64
    × ~180 offsets stalls the setup for hours") — a **load imbalance** the biggest pair runs alone on one thread,
    so *per-pair* OpenMP cannot help (99 % CPU throughout).
  - **→ the real fine-grid lever is an ALGORITHMIC `MakeLocalPP` fix, not threading.**  Step 0's per-iteration
    ~1.7× stands (committed); the setup is a separate, algorithmic problem.

### 0e-PP. `MakeLocalPP` SETUP WALL — the CP2K local-PP split (analysis 2026-07-19; NEXT implementation)
**Root cause (measured per-pair):** the `relCutoffScale=6` sweep is **1.6e9 grid points / 290 s**, spread over
406 pairs (NOT a few giant ones — no load imbalance), because scale=6 drags the DIFFUSE pairs (e.g. F s
α=0.275, reach ~9 au × ~180 cell images) onto field-resolution grids.  The energy `∫χ²V_loc` is dominated by
the DEEP WELL near the nucleus (erf/r ~ `Zion/r_loc`, width `r_loc`), so resolving it needs `ecut~1/r_loc²`
for **every** contributing pair — and the diffuse pairs DO contribute (measured below).
- **DEAD END 1 — threading:** memory-bandwidth-bound, 290 s → 290 s on 4 threads (per-pair OR intra-pair
  can't help a bandwidth wall).
- **DEAD END 2 — reduced-exponent level rule** (`p_eff = p/(1+2p·r_loc²)`, parameter-free from r_loc+basis):
  FALSIFIED.  Si Γ **over-binds to −7.216** (vs −7.11485) — the diffuse pairs genuinely couple to the well, so
  coarsening them aliases.  And it doesn't even help cost: `p_eff≈p` for diffuse pairs, so the giant-box pairs
  aren't coarsened at all.  (Also measured: a FIXED long-range `scale` can't serve both elements — Si soft
  r_loc=0.44 ok at scale 3, F hard r_loc=0.2 gives NaF −23.6 at scale 3.)  **Conclusion: per-product grid
  integration of V_loc cannot be both correct and cheap by ANY level rule — the well must be sampled per pair.**
- **→ THE FIX = the CP2K split** (the analogue of Ewald's erf/erfc; `r_loc` is our α, fixed by the PP):
  - **LONG-RANGE `−Zion·erf(r/√2 r_loc)/r`** = a Gaussian core charge → fold into the **G-space Poisson**
    (`PW_Hartree`, one electrostatics term).  The deep well is sampled **once per atom** (the core-charge
    collocation), not per orbital pair → no giant boxes; energy `E_een_long = Σ_G ρ̃_elec(G)·V_long(G)`, exact
    and adjoint-consistent with the matrix (stays variational, `E=Tr(D·H)`).  CP2K: `rho_core` passed to
    `pw_poisson_solve` (`qs_ks_methods.F`); split in `qs_core_hamiltonian.F:54`.
  - **SHORT-RANGE `poly×Gaussian`** = compact, CONVERGENT lattice sum (no Ewald) → **analytic** via the
    `LatticeSum1E::MakeOverlap(GaussianFunction)` seam the analytic KB already uses.  No grid, no `scale`.
  - **Deletes the user knob** (the grad-student-first-day goal): the grid is derived from `r_loc` (PP) +
    basis exponents; the fine cutoff is already auto (`4·α_max`).  CP2K still exposes `CUTOFF`+`REL_CUTOFF` as
    user convergence knobs (the rite-of-passage) — this design is MORE automated.
  - **Implementation increments (each gated on Si Γ == −7.11506, then NaF == −27.756, then re-time):**
    (1) `PW_Hartree` owns the `LocalPotential` + structure; total field `V_H[ρ_elec]+V_long(G)`; energies
    `E_hartree=½Tr(D·V_H)` + `E_een_long=Tr(D·V_long)`; the G=0 alignment (`FormFactorG0`) moves with it.
    (2) short-range → analytic Gaussian seam; drop the `MakeLocalPP` grid sweep entirely.  (3) re-time NaF.
    Interface: expose the GTH split on `LocalPotential` (`FormFactorLong/Short`, a core-charge/`r_loc`
    accessor).  Decided (user 2026-07-19): fold into `PW_Hartree` — "cleaner physics, one Poisson solve."
  - **STATUS 2026-07-19 — increment 1 (the split) + Q1 (the grid speedup) are DONE (branch
    `gpw-0e-pp-local-split`; compact record in the [DONE](#done) timeline).**  The split is a `LocalPotential`
    form-factor property (`FormFactorLong` primary + base-provided `FormFactor=Long+Short`); `PW_Pseudo` does the
    SHORT local, `PW_Hartree(fb,st,loc)` folds the LONG `V_long` into its Fock matrix + owns its energy/alignment
    — a matrix-identical ENERGY-RELOCATION refactor (Si Γ −7.11506 + NaF −27.756 held, 202/202).  **Q1 corrected
    the plan:** the ~295 s wall is the `relCutoffScale`, NOT the per-pair sweep — it was over-set to 6 by the
    DENSITY-SCREEN bug (the `−280`/`−259` was `OverlapMatrix`'s `screenD` zeroing off-diagonals of the FIXED
    `V_long`, NOT aliasing; unscreened, smooth==stiff to 4e-3 for soft Si).  Default `relCutoffScale` 6→3 → ~2×
    (Ecut=160: 578 s → 128 s @scale 2), all gates green.  So implementation increment (2) above is RE-SCOPED:
    **the grid knob (Q1) is the perf fix; the analytic seam is a separate ACCURACY upgrade** — it re-gates to
    converged CP2K (band-limiting cancellation: grid short/long each ~0.5 Ha off, cancelling in the smooth full
    V_loc, so analytic-short + grid-long misses the gate by 0.55 — BOTH must go analytic together).
  - **STEPS (a)+(b) DONE (2026-07-22): the SHORT local PP is ANALYTIC in production; the grid sweeps are
    STANDALONE-exact under the ABSOLUTE κ rule.**  User-approved sequence (before mixed-radix FFT, so no
    sharp object remains on the raster to corrupt at smaller N):
    - **(a) The absolute pair→level rule** replaces `relCutoffScale` at the seam
      (`LatticeSum1E::IntegratePotential(..., absRelCutoff)`; `PairLevel`): req = κ·(αᵢ+αⱼ), coarsest
      satisfying, finest fallback — CP2K `gaussian_gridlevel` semantics.  THE INSIGHT: the absolute rule
      bounds EVERY pair's spectral tail by e^{−κ/2} UNIFORMLY, independent of the field's sharpness —
      which is what "REL_CUTOFF" really is and why CP2K's numeric-but-smooth V_long is sub-mHa.  κ=30 Ha
      (e^{−15}) default for the local-PP sweeps (`LocalPPRelCutoff`, env `GPW_LOCALPP_RELCUTOFF` for
      verification); the density-side collocate/integrate keep the RELATIVE rule (adjoint-paired).
      `MakeLocalPP` now runs the FULL ladder (top rung included); `relCutoffScale`/`GPW_LOCALPP_SCALE`/
      `GPW_LOCALPP_FULL` deleted.  Gate `GPW.LocalPPKappaSelfConverged`: κ=30 vs κ=60 → Full 7.6e-9 /
      Long 1.4e-9 / Short 1.6e-8 (the e^{−15} class on the nose).  KNOWN HOLE in that self-check: pairs
      SATURATED at the ladder top are κ-invariant by construction, so self-convergence is blind to them —
      harmless for the r_loc-soft LONG (sharp pairs meet a tiny field tail at the rung ball), visible for
      the SHORT at cheap ladders (the gate's 3.6e-3 cross-val residual at Ecut=10 = exactly this class).
    - **(b) `GPW_IBS::MakeLocalPotentialShort` → the ANALYTIC `MakeLocalPPShort`** (exact 3-centre
      Gaussian lattice sums; grid fallback for non-Gaussian models).  Safe ONLY after (a): the old grid
      short's ~0.5 Ha band-limit error CANCELLED the grid long's — exact-short + lenient-long missed the
      gate by 0.55 (the recorded trap).  **WIRING BUG CAUGHT by the new cross-val: the periodic G=0
      convention** — the grid sweep drops ΔG=0 (cell mean → Ealign via `FormFactorG0Short`) while the
      analytic sum integrates it: 5.7% disagreement → subtract V̄·S for a periodic Structure (the same
      `isFinite()` physics decision as `PW_Pseudo`) → 0.36% (the remaining = the grid REFERENCE's
      saturation at the cheap test ladder, see (a)).
    - **VERIFICATION:** Si Γ anchor −7.11526 IDENTICAL to 5 decimals between κ-ruled-grid-short and
      analytic-short (μHa-level agreement at production settings); vs CP2K −7.11506 the anchor moved
      −7.11482→−7.11526 (0.24 mHa below → 0.20 above — same distance class, the ±2 mHa gate holds, and
      V_loc discretization no longer leans on long/short cancellation).  Atom-in-box (finite branch)
      green; (a)-only sweep 199/199.  NEXT: re-run the NaF SR2 oracle config (expect ≈−24.4317 within
      ~mHa) + re-time the setup (the short sweep is deleted; the long sweep remains, κ-ruled).
  - **REMAINING TODO — analytic V_local LONG (the Ewald/core-charge crux). Branch `gpw-0e-pp-local-split`.**  Both pieces are
    EXISTING `GaussianRF` kernels (no new Boys function): short = `Overlap3C(χ_i,χ_j,g_short)`, long =
    `−Z_ion·Repulsion3C(χ_i,χ_j,g_core)` (the erf-Coulomb IS a normalized Gaussian core charge, exp `1/2r_loc²`).
    SHORT is BUILT + finite-validated but DORMANT (`LocalPotential_Gaussian::ShortRangeGaussian`,
    `LatticeSum1E::MakeLocalGaussian` = the 3-centre `Overlap3C` MATRIX sibling of the 2-centre `MakeOverlap(g)`
    VECTOR, `GPW_Evaluator::MakeLocalPPShort`).  **LONG is the crux:** the `Repulsion3C` lattice sum is
    conditionally convergent (erf→1/r Madelung tail) ⇒ needs a G-space/Ewald neutralizing background, NOT a
    real-space sum.  Both go analytic TOGETHER; the exact total re-gates NaF to converged CP2K −27.93 (a WIN
    over the Ecut=40 grid −27.756).  **DO AFTER the −39 basin fix** — the fine-grid SCF diverges regardless of
    V_local, so the energy can't be verified until it converges; then also use `GPW_LOCALPP_SCALE=2/4` to verify
    grid scale-convergence.
- **Steps 1 (grid-continuation seeding) + 2 (XC-collapse ROOT-CAUSED & FIXED): DONE 2026-07-20 — moved to the
  [DONE](#done) timeline** (full record there: the fit-grid thread-through, the `Overlap3C` adjoint, the
  density-fit densification → one-grid `cutoffFactor` 4→8, the `itsFFT_R_G_Grids` rename, the two bug fixes, and
  the HONEST PICTURE — the old −27.75 was an aliasing coincidence, the resolved answer −26.198, the residual gap
  now the coarse LOCAL-PP base grid).  What remains is the validation + local-PP work below.

  **★ NEXT SESSION — QUEUED (user, 2026-07-20). The definitive numerics check + the instruments for it:**
  - User story inserted here:  We use real and reciprical space grids for FFT, integrals (some using Parsevals theorem), and to define CD (rho) and Vxc PW fit basis sets.  Sometimes the code uses the same grid (for example the FFT G grid is automatically used for the CD fit basis ... is this fully justified ... I don't know!).  In the last session we discovered that the whole fit basis set was completely ignored inside GPW_IBS::Repulsion3C(const CDFitBasis* c) !! THis should be fixed now.  In general the evlautor classes should not be making high level policy decisions like "what is the CD fit basis set?". There is also a grid for integrating over V_local part of the PP.  Why we use a grid at this ttype of integral is unkown to me.  Is this the same as the FFT grid?  SHould it be? (I think no becuse we need to take rloc into account)  Also every time I see the term Ecut ... I immediatly think "Ecut for *what*??"  Anyway as you I can see I am very confused and losing confidence in the GPW code, so I decided it is time to cross check with CP2K all grids and the energy break down.  After that I would like see a table of all grid usages, how grid range spacing is decided, and if it is a user knob or decided by a sensible algo.  If we can get all of this right then there should be no fake, unphysical E=-39 basin of attraction to avoid. In general NaF is an excellent test case becuse it has very diffuse Na basis functions and very sharp F basis functions ... we need handle both of those with distinct optimized grids.  It forced us to the grid manabngemt near perfect.  I am totally open to any suggestions on high level strategy (CP2K cross check or other).
  1. **GRID-MATCHED CP2K VALIDATION** — does GPW's NaF GS energy == CP2K −27.93128 (we should also check the whole energy breakdown) when EVERY density-side grid is
     made IDENTICAL to CP2K's?  Match all four: (a) the **FFT grid** `itsFFT_R_G_Grids` = CP2K `CUTOFF` (N, Ecut);
     (b) the **ρ fit grid** `{G}_ρ` (CD fit) = CP2K density grid; (c) the **v_xc fit grid** `{G}_vxc` (Vxc fit) =
     CP2K XC grid; (d) the **V_local integration grid** = CP2K's local-PP grid (its `REL_CUTOFF` multigrid
     assignment).  PURPOSE: remove grid resolution as a variable — if the energies then AGREE, GPW's numerics ==
     CP2K (validated end-to-end); if they DON'T, it pins the remaining difference as the collocation METHOD (our
     Fourier round-trip aliases where CP2K's REAL-SPACE collocation stays graceful — the measured 2× : our
     8·α_max ≈ CP2K's 4·α_max).  Read CP2K's grids from its log (`&MGRID`: `NGRIDS`/`CUTOFF`/`REL_CUTOFF` + the
     per-level N it prints); force GPW to those exact N/Ecut (explicit `densityEcut`, a matched `relCutoff`, and a
     local-PP-grid override).
  2. **GRID DIAGNOSTIC PRINT (cout, at the START of every run)** — one essential line PER stored grid:
     `N=45 |Gmin|=0.01 |Gmax|=160` (FFT divisions N, min/max |G|).  Print ALL of them (FFT/ρ-fit/vxc-fit/local-PP)
     so we can SEE what GPW uses and line it up against CP2K.  (Essentials only, not the full {G} list.)
  3. **ORBITAL-BASIS EXPONENTS (cout, run start)** — print α_min and α_max of the orbital basis, so the
     α_max→grid policy (`cutoffFactor·α_max`) is visible and checkable.
  (Items 2–3 are the instruments for item 1.  Also still open from before: the FINE auto grid at `8·α_max` puts the
  LOCAL PP on the fine grid too — should close the Een gap toward −27.93; and the CP2K-vs-us real-space-collocation
  2× as a future efficiency lever.)

  **★ RUN 2026-07-21 — items 1–3 DONE; the energies DO NOT AGREE → the gap is the collocation METHOD,
  not grid settings.**  (The user-story questions above are answered in the new **`doc/GPWGrids.md`** —
  the requested table of every grid, its sizing rule, and knob-vs-algorithm status.)
  - **Instruments LANDED (items 2–3):** `GPW_Evaluator::ReportGrids` — run-start cout of the basis
    exponents (α_min/α_max/cutoffFactor → the auto floor) + one line per STORED grid (FFT reference,
    every ladder level incl. the top rung, `{G}_ρ`/`{G}_vxc` at their factories, the local-PP sub-ladder
    + relScale), each with N/Ecut/n_G/|G|min/|G|max; called once per run from the `GPW_BasisSet` ctor.
    Plus two grid-MATCHING knobs (env-gated verification instruments): `GPW_MGRID_ECUTS` (explicit
    sub-level cutoff list — CP2K's progression-3 ladder is unreachable by our factor-4 default;
    skips-with-warn entries ≥ the block's reference so a coarser same-process block keeps a valid ladder)
    and `GPW_RELCUTOFF` (Ha: switches `PairLevel` to CP2K's ABSOLUTE `gaussian_gridlevel` rule
    `req=(αᵢ+αⱼ)·REL_CUTOFF`, finest-as-fallback, ignoring the relative rule + relCutoffScale).
  - **CP2K oracle RESTORED on the new machine** (the migration lost `~/Code/cp2k`): conda-forge
    **CP2K 2026.1** in the `cp2k` env of `~/miniforge3` (recipe in `UnitTests/CP2K/README.md`); source
    shallow-clone at `~/Code/cp2k` (v2026.1 — for reading algorithms + `data/GTH_POTENTIALS`); run dir
    `~/Code/cp2k-runs/`.  Si Γ **−7.11505788** and NaF **−27.9312751** reproduce the recorded oracles
    exactly (NaF grid leak 1.95e-4 e, same class).  **CP2K's ACTUAL NaF grids (from the `PW_GRID|` log):
    160 Ha/36³, 53.3/24³, 17.8/12³, 5.926/8³; REL_CUTOFF 30 Ha; pair spread 3973/3149/3532/2342** —
    table in `doc/CP2Kresults.md`.  CP2K's N=36 is MIXED-RADIX (2²·3², FFTW-class); our radix-2-only
    FFT pads the SAME 160-Ha ball to **128³ = 45× the points** — a standing efficiency lever.
  - **THE MEASUREMENT** (grid-continuation test, `GC_FINE_ECUT=160 GPW_MGRID_ECUTS=53.33,17.78,5.93
    GPW_RELCUTOFF=30`, 4 threads; full log `~/Code/naf_gridmatched.log`): coarse seed (Ecut=40, stiff
    rule) converges 49 iters to −24.099 — the CP2K-stiff assignment REMOVES the aliasing that flattered
    the old −27.76 coarse number (honest-picture class).  Fine stage: same Ecut ladder + assignment as
    CP2K, DENSER raster (128³ vs 36³), MOM+Pulay: converges CLEANLY, 22 iters, charge 8.0000000000 —
    **Etot = −23.6739 vs CP2K −27.9313: Δ = 4.26 Ha.**  Breakdown (the splits differ — ours
    Ekin/Een/Eee/Exc/Enn/Ealign vs CP2K's compensating-core scheme — so compare the clean common terms):
    **Ekin ours 18.2570 vs CP2K 19.1408** (an ANALYTIC, grid-free term: 0.88 Ha means the converged
    DENSITIES differ — the discrete functionals have different fixed points, not just different
    bookkeeping); **Exc −4.4837 vs −3.7398** (ours over-negative: the Gibbs-lobe signature); the
    electrostatic remainder carries the rest (~5.9 Ha, density-shift and V_loc-discretization mixed).
  - **VERDICT + refined mechanism:** with Ecut/ladder/assignment matched, N finer on our side, basis+PP+
    functional identical, and both SCFs stable/charge-exact, the 4.26 Ha is the DISCRETIZATION METHOD.
    Sharpest identified difference: **we project the collocated ρ onto the {G}_ρ BALL before XC**
    (`RhoOnGrid` over the n_G=16145 ball coefficients — a hard spherical truncation of the sharp F
    products → Gibbs negative lobes → shallow XC), while **CP2K evaluates XC on the RAW collocated
    raster values** — it never ball-limits ρ (its 36³ raster even RETAINS corner G-content beyond the
    160-Ha sphere; its FFT/ball only transports the smooth Poisson/level-transfer fields).  Both rasters
    alias; ours is additionally Gibbs-truncated — the measured 2× (our clean 8·α_max vs CP2K's 4·α_max)
    is the price of that ball round-trip: a FITTING-SEAM design question, not raster resolution.  Second
    contributor: grid-integrated V_loc vs CP2K's analytic-short + core-charge-long (§0e-PP, unchanged).
    **NEXT LEADS (in order): (1) evaluate XC on the raw collocated ρ — skip the ball projection between
    collocation and XC (the single most CP2K-aligning change; re-examine `CreateVxcFitBasisSet`'s role);
    (2) analytic V_local (§0e-PP); (3) mixed-radix FFT (36³-class rasters) for the 45× raster cost.**
    USER (2026-07-21): agrees with 1/2/3; lead 1 started — plan in §0f below.

## 0f. LEAD 1 — XC on the RAW collocated ρ (STARTED 2026-07-21; user-approved)
**Hypothesis to kill or confirm (from the grid-matched verdict):** our XC eats `RhoOnGrid(ball ρ̃)` — the
BALL-projected density — and the hard spherical truncation of the sharp F products is the Gibbs/negative-lobe
source that CP2K (XC on the raw collocated raster, never ball-limited) does not have.
- **Increment 0 — the BALL-RESOLUTION PROBE (no refactor; existing knobs).**  Re-run the grid-matched NaF
  with ONLY the reference ball enlarged: `GC_FINE_ECUT=480` (3× — the ball then covers the raster-box corner
  content a 160-box keeps), same `GPW_MGRID_ECUTS=53.33,17.78,5.93` + `GPW_RELCUTOFF=30`, raster expected to
  stay 128³ (pow2 padding absorbs the √3).  PREDICTION if the hypothesis is right: E moves decisively from
  −23.67 toward −27.93 with the Exc term normalizing toward CP2K's −3.74 (the V_loc/Hartree band-limits ride
  along — the term breakdown separates the XC share from the Een share).  If E barely moves, the ball story
  is dead and the method difference is elsewhere (V_loc first suspect).
- **Increment 0 RESULT (2026-07-21, same day): the BALL HYPOTHESIS IS FALSIFIED — and the REAL story
  surfaced: the grid-matched runs converged to a MOM-PINNED EXCITED STATE.**
  - The 480-Ha ball (nG 16145→83659, raster unchanged 128³, coarse seed bit-identical −24.0988) gives
    **−23.67372 vs the 160-Ha ball's −23.67387 — 0.15 mHa, every term sub-mHa** (log
    `scratchpad/naf_ballprobe480.log`).  Tripling the ball moves NOTHING: the collocated ρ, v_xc, V_H
    and the fine-level V_loc are all ball-CONVERGED at Ecut=160 on this raster.  The Gibbs/ball-
    truncation mechanism is DEAD; do NOT build the raw-raster path on this motivation (increment 1
    below is SUSPENDED pending the aufbau probe).
  - **The `frontier ε(occ)` instrument then exposed the actual disease: BOTH grid-matched endpoints are
    NON-AUFBAU.**  Converged frontier: `−0.2767(4.0)  +0.2105(2.0) | −0.3631(0.0) ...` — an
    UNOCCUPIED level at −0.363 Ha sits BELOW both occupied frontier levels; MOM (reference transferred
    from the differently-discretized coarse stage via `AdoptMOMReference`) holds 2 e⁻ in a +0.21 Ha
    level with a −0.36 Ha hole underneath = a pinned EXCITED state, plausibly the bulk of the 4.26 Ha.
    The §0e★ "collocation METHOD" verdict is SUSPENDED until the aufbau ground state converges on the
    matched grid.  (Corollary: the step-1 fine fixed points −24.393 / −24.099 / −23.674 — all
    AdoptMOMReference runs — are suspect for the same reason.)
  - **Probe 3 RESULT — pure-aufbau fine stage converges CLEANLY to −24.4317** (22 iters, Δρ 2e-6, gap
    0.233 Ha, proper occupations F2s(2)+F2p(6), no spikes — the §0b″ swap disease did NOT return on the
    seeded start).  The pinned excited state was 0.76 Ha; **3.50 Ha to the −27.93 oracle remained** —
    which prompted the BASIS AUDIT: our runs use `VALENCE_LOWQ_SR2` (n=28) but the CP2K oracle deck ran
    `VALENCE-LOWQ-SR` (Na keeps s 0.0857 + p 0.05) — the oracle predates the 0b′ SR2 trim.  **THE
    APPLES-TO-APPLES RUN SETTLES IT (CP2K on a transcribed SR2, `naf_gpw_sr2_diag.inp`):**

    | | qchem GPW (aufbau, matched grids) | CP2K 2026.1 (SR2) |
    |---|---|---|
    | **Etot** | **−24.4316608** | **−24.4312134** |
    | Exc | −4.95597 | −4.95531 |
    | SCF | 22 iters, Δρ 2e-6, clean | **converged in 16 steps** (no limit cycle!) |

    **GPW == CP2K to 0.45 mHa (Exc to 0.7 mHa) — the implementation is VALIDATED end-to-end.**  The
    §0e★ 4.26 Ha decomposes EXACTLY: **0.76 Ha = the MOM-pinned excited state** (an SCF-strategy bug,
    below) **+ 3.50 Ha = the SR-vs-SR2 BASIS difference** (a comparison error: the oracle basis was
    never re-matched after the SR2 trim).  The "collocation METHOD" verdict is RETRACTED — the ball
    projection, the grid V_loc, and the raster convention are all inside 0.45 mHa at matched settings.
  - **The 3.50 Ha is REAL VARIATIONAL PHYSICS, not a numerical artifact: BOTH codes agree on BOTH
    bases** (SR ≈ −28.0 [our §0b′ honest-map fixed point ≈−28.00; CP2K −27.93 E-flat] vs SR2 ≈ −24.43
    [both codes, sub-mHa]).  So the near-null (λ~1e-6) SR modes the SR2 trim dropped were NOT
    physically null — Na's diffuse s 0.0857 + p 0.05 carry ~3.5 Ha of interstitial/F⁻ variational
    freedom.  **SR2 is a well-conditioned but physically POOR truncation**; the conditioning-vs-
    completeness tension lands squarely on the §1 rank-reduction track (run the FULL/SR basis, let the
    ORTHO transform drop the null space — not the basis) or on re-optimized diffuse exponents at
    λ_min~1e-3.  Also: CP2K's clean 16-step SR2 convergence == our clean 22-iter run confirms the
    ill-conditioning story of the SCF misery on BOTH sides.
  - **Two SCF-strategy lessons banked:** (a) `AdoptMOMReference` ACROSS a discretization change can pin
    an EXCITED state (the transferred occupied subspace need not span the new grid's aufbau ground
    space) — needs a guard: detect a persistent hole (unoccupied ε below an occupied ε at convergence)
    and either release MOM or re-capture; (b) the `ReportBandGap` εH/εL summary line MASKED the hole
    (it printed gap=0.67 from the lowest virtual ABOVE the HOMO index while a −0.36 Ha virtual sat
    BELOW) — the `frontier ε(occ)` window is the honest instrument; teach the gap line to flag
    non-aufbau (εL taken over ALL unoccupied, not just index-above).
  - **LEADS RE-SCOPED:** raw-raster XC (increment 1) = DEAD (falsified twice over).  Analytic V_local
    (§0e-PP) = demoted from accuracy-blocker to robustness/perf (its grid errors are inside 0.45 mHa
    at these settings).  Mixed-radix FFT = stands, pure efficiency (128³ vs CP2K's 36³ = 45× raster
    points for the same answer).  The NaF critical path is now **basis completeness under conditioning
    (§1 rank-reduction)** + the MOM guard.
- **Increment 1 — the raw-raster path (SUSPENDED by increment 0's falsification; kept for the record).**  Combine the per-level
  collocated densities onto the fine raster by FULL-BOX zero-pad upsampling (per-level FFT → embed level-L
  G-box into the fine G-box → accumulate → one inverse FFT): `ρ_raw` on the fine raster, no ball anywhere.
  XC (E_xc grid sum AND pointwise v_xc) evaluates on `ρ_raw`.  Integrate-back: FFT(v_xc) full box →
  per-level BOX truncation (the exact adjoint of zero-pad upsampling) → iFFT per level → the EXISTING
  analytic `IntegratePotential` seam.  Adjoint-exact by construction ⇒ H_xc = ∂E_xc/∂D; re-gate with
  `GPW.XCPotentialConsistencyFD`.  Hartree/Poisson STAYS on the ball (diagonal kernel, variational — the
  legitimate projection; `doc/GPWGrids.md` row 2).  Design decision en route: where the real-space-density
  seam lives (the ΔG_Map-speaking `G_ERI3`/`Band_FT_IBS` faces are ball-shaped; the raw path wants
  rvec_t rasters — likely a new capability on the fit-basis/G_FieldEvaluator side, since "what ρ does XC
  see" is the fit basis's policy question).
- **Increment 2 — re-calibrate `cutoffFactor` DOWN.**  If the raw path holds XC clean at 4·α_max (the CP2K
  operating point), retire the 8 → halve the auto grids (the 2× runtime lever measured in §0e step 2);
  re-anchor Si/NaF.

## 0d. Runtime follow-ups (after 0b/0c)
- **OpenMP over the per-iteration collocate/integrate pairs — DONE (step 0 above).**  Memory-bound → ~1.7×.
- **`MakeLocalPP` is the fine-grid SETUP wall (~290 s of a ~320 s ctor) and needs an ALGORITHMIC fix, not
  threading** (profiled 2026-07-19; full record in §0e step 0).  The `relCutoffScale=6` static local-PP sweep
  forces a few ultra-diffuse pairs onto huge fine-grid boxes → a load imbalance that per-pair OpenMP cannot
  touch (measured: no speedup).  Fix leads: (a) smarter sharp-field level assignment so an ultra-diffuse pair
  (whose own spectrum kills the field tail) stays on a deep coarse level — the sweep's own comment argues this;
  (b) intra-pair (over-offset) parallelism for the few giant pairs.  A per-pair-OpenMP `EnsureStreams` build
  (only ~25 s, and also load-imbalanced) was tried and reverted — no benefit.  CP2K's ssmp is threaded on top.

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

**NLCC vs semicore — decision point when the TM-oxide (Mn/Ni/Co, battery-track) bases are built.**  Our XC
is valence-only (E_xc[ρ_val], v_xc[ρ_val]) — CORRECT for the GTH-PADE set we ship (`gth_potentials.json`
has NO NLCC/core-charge entries; the core-valence XC linearization is absorbed at PP generation, and CP2K
runs the same PPs the same way, so all oracles are apples-to-apples).  The GTH remedy where linearization
fails (spin-polarized TM cores) is historically SEMICORE promotion (the Na q1→q9 pattern; sharp semicore
density → much higher grid cutoff), the alternative is NLCC-GTH (Willand 2013 style; CP2K supports an NLCC
section).  If NLCC is chosen: the core density is an analytic per-atom Gaussian → ONE more static
collocation onto the same grid (like the local-PP sweep), then ε_xc/v_xc evaluated at ρ_val+ρ_core in BOTH
the energy and the integrate-back field; ∂ρ_core/∂D=0 so H_xc stays the exact gradient and the
`XCPotentialConsistencyFD` gate covers it unchanged.  Forces add the core-motion term (forces increment).

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
  **PREREQUISITE NOW IN PLACE (2026-07-20): the fit-grid seam is honest.**  Previously `GPW`'s
  `MakeRepulsion3C(c)`/`MakeOverlap3C(c)` (the shared `EPW_Orbital_DFT_IBS` mixin) DROPPED the fit basis `c`
  and rebuilt the tensor from the block's own `itsGrid` — so a denser `CreateVxcFitBasisSet` grid would have
  been SILENTLY IGNORED (the policy factory and the tensor builder were two disconnected sources of truth for
  the density-fit `{G}`, reconciled only by both hard-coding `DensityGrid()`).  Now `GPW_IBS` overrides those
  two seams to build the tensor over the REQUESTED fit basis's grid (`c` IS-A `PW_Grid_Evaluator`;
  `GPW_Evaluator::Repulsion3CTensor(grid)`/`Overlap3CTensor(grid)` + a grid-parameterized `BuildLevels` ladder).
  Bit-identical while the factory wraps `DensityGrid()` (Si Γ −7.11485 / multi-k −7.45133 / adjoint
  machine-exact / all GPW gates green), and the block's own `OverlapMatrix`/`MakeLocalPP` (KS-assembly, not a
  requested table) keep `itsGrid`.  So densifying `CreateVxcFitBasisSet` will now ACTUALLY take effect for the
  collocated ρ̃ — the GGA increment can diverge the CD/Vxc grids without the tensor silently overriding it.
  (PW's own `relCutoff` Vxc path — `PlaneWaveDFT.ItemK_RelCutoffDensifiesAndConvergesVxc` — was left untouched,
  deliberately not lumped into the GPW-scoped fix; audit it separately if the shared mixin is ever unified.)
- **Multi-grids + whole-density collocation — DONE** (the C+D analytic rewrite, see the DONE entry).
- **Common `PW_Evaluator`/`GPW_Evaluator` base:** the shared FFT engine (`PeriodicGridEvaluator`) is already
  factored; a base earns its keep once general-k k-space logic is shared.

---

# Durable pins / invariants (carry into all GPW work)
- **THERE IS NO CUT — IN THE R DIRECTION (user pin, 2026-07-16).**  Real-space lattice sums are
  ε-CONVERGED SERIES for a FIXED operator: magnitude screening is the ONLY truncation mechanism; a radius
  must never appear as a parameter, member, or concept in any interface — not user-facing, not internal.
  A truncation radius yields a DIFFERENT operator, not "the operator to ε" (measured: the Rcut=2a NaF
  metric lost 2.25 e per mid-slosh loading), AND must never be a conditioning crutch (that job belongs to
  the basis or to rank-reduction).  The G DIRECTION is different in kind: the Ecut ball is a PROJECTION
  onto a finite auxiliary subspace — variational (adjoint-exact), exponentially controlled, systematically
  improvable — i.e. a legitimate resolution dial, not a cut.  End state: ONE knob per direction —
  ε in R (convergence tolerance), Ecut in G (projection resolution).
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
  `AnalyticSeparablePPMatchesMesh` == mesh KB to 4.6e-11; `XCPotentialConsistencyFD` — H_xc == ∂E_xc/∂D to
  FD accuracy in both the smooth and the ρ<0-guard regimes, the 0b falsification gate),
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
