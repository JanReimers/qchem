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

- **Grid instruments + grid-matched CP2K validation → GPW VALIDATED** (`aecbb410`,`a6560e3e`): run-start
  grid diagnostic (`ReportGrids`) + match knobs (`GPW_MGRID_ECUTS`/`GPW_RELCUTOFF`); CP2K restored via
  conda-forge; the 4.26 Ha "gap" decomposed = 0.76 MOM-pinned excited state + 3.50 REAL SR↔SR2 basis
  physics (ball/Gibbs hypothesis falsified by the 480-Ha probe); **NaF SR2 == CP2K to 0.45→0.19 mHa**;
  `doc/GPWGrids.md` = the grid inventory.  Full records: §0e★/§0f below.
- **§0e-PP (a)+(b): absolute κ rule + ANALYTIC short V_loc in production** (`5d963b04`): req=κ·(αᵢ+αⱼ)
  (CP2K `gaussian_gridlevel`; e^{−κ/2} pair tails, κ=30) replaces `relCutoffScale`; analytic 3-centre
  short (the periodic G=0 double-count caught by the new gate); Si Γ grid-vs-analytic identical to 5
  decimals.  Gate `GPW.LocalPPKappaSelfConverged`.
- **MIXED-RADIX RASTERS** (`1837b21e`,`aaaf1ea0`): PocketFFT submodule behind `qchem.FFT` (pow2 →
  radix-2 verbatim = bit-identical) + `FFTGrid()` pads to 5-SMOOTH N: 199/199 with ZERO re-pins,
  **NaF verification 94 min → 190 s (30×) at 0.4 µHa** (fine raster 128³→72³); CP2K gap 900×→33×.
  OpenMP pair-loops (`GPW_OMP_THREADS`) remain opt-in (~1.7×, bandwidth-bound).

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

## TRAPS — the distilled do-not-revisit list (full records: doc/GPWHistory.md)
1. **Aliasing-flattered energies.**  Every pre-2026-07-20 NaF number in the −27.7..−28.0 band (Ecut=40
   era) was an under-resolved-grid COINCIDENCE, not physics.  Resolved-grid truth: SR2 −24.4314 /
   SR −24.4324.  Never trust an energy whose grid fails the negCharge/XC probes.
2. **Oracles are not always oracles.**  The CP2K SR −27.93128 was a 3.50 Ha EPS_PGF_ORB screening
   artifact (truncated overlap metric × 1/λ amplification on a λ~1e-6 basis) — complete with a FAKE
   eternal density limit cycle.  Rule: CP2K oracle runs on ill-conditioned bases use tight EPS
   (`naf_gpw_sr_tight.inp`) and require REAL density convergence; also DIFF THE BASIS before comparing.
3. **MOM across a discretization change can pin an EXCITED state** (0.76 Ha on NaF; hole at −0.36 Ha
   under occupied levels).  The εH/εL gap line MASKS it — read the `frontier ε(occ)` window (0h guards).
4. **The long/short local-PP split pieces carry ~0.5 Ha grid errors that CANCEL in the sum** — never mix
   an exact piece with a band-limited partner (measured trap).  The absolute κ rule (e^{−κ/2} uniform
   pair tails) is what makes a piece standalone-exact — up to LADDER-TOP SATURATION (κ·p above the top:
   harmless only where the field is r_loc-soft — fine for LONG at production tops, why SHORT is analytic).
5. **Truncated metrics corrupt maps** (Rcut=2a lost 2.25 e; CP2K's screening = the same class).
   THERE IS NO CUT.  Also: overlap-null ≠ physics-null — but for THESE bases the near-null diffuse
   modes carry only ~1 mHa (SR↔SR2), so auto-dropping them is safe.
6. **Geometry-keyed caches + lattice image clones = unbounded growth** (the 11-GB OOMs: content-keyed
   MnD Ω/RNLM/H3 on ~2k images/pair; stream-cache pair builds materialised before tiering).  A lattice
   SERIES is consumed once — stream it (§5 LRU design retires the interim clear-based band-aid).
7. **Ops on this 14-GB box:** ONE heavy run at a time, always inside `systemd-run --user --scope
   -p MemoryMax=10G` (systemd-oomd kills by CGROUP — an unscoped run makes the DESKTOP APP take the
   blame); full-SR-class stream demand is 5.3B pts (`GPW_STREAM_BUDGET_PTS[_F32]` caps).
8. **Sampling collocation (pre-analytic era) aliases at bulk; hard Rcut rings (Gibbs)** — the analytic
   CP2K method (compact exp-tail boxes, screened image sums, modulo wrap) is the only scheme in the code.

# TODO / NEXT

**⇒ The FORWARD QUEUE now lives in `doc/GPWPlan1.md` (2026-07-23): param-struct graduation →
per-system display → Cache2/3 LRU → diffuse robustness + Fermi smearing → B_ij(R) with IBZ.  This
file remains the RECORD of the 2026-07 campaign (0.5(b)/(c)/(f), 0h, C=2, raster A/B) + the durable
pins; deep archive in `doc/GPWHistory.md`.**

**Orientation (2026-07-23).**  GPW is VALIDATED against CP2K at sub-mHa on every honest comparison (Si,
NaF-SR2 0.19 mHa, NaF-full-SR 0.10 mHa — after the CP2K SR "oracle" −27.93 was RETRACTED as its screening
artifact; TRAPS #2).  NaF production (SR2, matched grids) runs in 190 s.  The queue: §0.5 runtime —
b, c, f1, f2 DONE 2026-07-23 (raw-XC feed: collapse basin removed, C=3 NaF 0.2 mHa vs CP2K) — then
**0h SCF guards (pulled ahead: they gate the C 8→3 default flip)**, 0.5(a) raster policy as one combined
{policy × C} calibration, (e), (d), 0i analytic long, §1 diffuse-basis robustness (demoted to
automation), then the standing items (2)–(5).

## 0.5 RUNTIME IMPROVEMENTS (the consolidated performance queue)
Current standing (all converged, same machine): Si Γ 48 s vs CP2K 3.5 s; NaF SR2 190 s vs 5.8 s (~33×);
the full-SR diagnostic ran 8.45 h vs 88.6 s — its pathology is items (b)+(c) below.  In EXECUTION order
(dependency-aware leverage; letters are stable labels, kept for cross-references).  Prelude: the
test-side `ctest -j16` upgrade (separate session) lands before (b) so every step below gets fast
confirmation runs.
- **(b) DONE 2026-07-23 — free the coarse stage's stream caches after the seed handoff.**  The class fix
  is TWO halves, both needed because the fine shape's cache is built DURING the handoff (the coarse block
  collocates the seed on the fine fit grid) while the coarse caches still hold the budget:
  (1) RELEASE — new seam capability `Molecule::LatticeSum1E::ReleaseStreams(N_L, ecut_L)` (PG_Cart
  forwards to the MnD evaluator, which erases the shape's caches, refunding the GLOBAL budget), called by
  `~GPW_Evaluator` with its own ladder shape — so `bsC.reset()` genuinely frees;  (2) SELF-HEAL — a
  StreamCache records `droppedPts` + the budget headroom it was offered; an INCOMPLETE cache rebuilds
  when headroom has GROWN (a resident shape was released), a complete cache NEVER rebuilds (bit-stable
  replay preserved on every anchor path), an unchanged starved cache never churns.  Gate:
  `GPW.StreamCacheReleaseUnstarvesLaterGrid` (two grids on one shared molecular basis under a tight env
  budget: starve → release → self-heal → no-churn; note the static-field IntegrateMemo replays identical
  fields WITHOUT consulting stream caches — the gate varies the field per call).  The grid-continuation
  SCF test now tears the whole coarse stage down (iterator → EC → basis, all unique_ptr) right after the
  seed/MOM handoff.  Suite 564/564.  VALIDATION NUANCE (full GC run 2026-07-23): on SR2 the fine demand
  (65M pts) fits tier 2 even with coarse resident — the starvation needs the full-SR-scale demand (952M),
  so the unit gate (tight budget) is the regression anchor, not the SR2 run.  That run also showed
  `DISABLED_NaFGridContinuation`'s energy pins were STALE — RE-PINNED same day on the clean landscape:
  coarse Ecut=40 fixed point −23.6929 (E-flat converged, 515 iters; the old −27.76 anchor was the
  retracted aliasing-era value — the leaky coarse grid honestly UNDERBINDS by 0.74 Ha now), fine aufbau
  ground state −24.4337 ± 0.01 (22 iters, 2.5 mHa from CP2K's Ecut=160-class −24.4312); fine defaults
  flipped to PURE AUFBAU (the MOM-transfer path = the 0h repro, env-gated).
- **(c) DONE 2026-07-23 — stream-budget follow-ups.**  (1) Byte-aware build transient via STREAMING
  fp32 demotion: the moment a pair's count exceeds the open fp64 budget it can only land in tier 2, so
  the already-built offsets demote immediately and the rest build demoted — the transient is bounded by
  the fp32 STORAGE (+ one offset's fp64 box) instead of 12 B/pt over the whole pair (~10 GB at the 850M
  default; the full-SR OOM class).  Same doubles narrowed → replay bit-identical, tiering unchanged.
  (2) `[stream cache]` readout prints the EFFECTIVE (env-overridden) budgets (rode with (b)).
- **(f) BALL-PER-ROLE re-calibration (user question 2026-07-23): one ball currently serves three roles
  with three different natural calibrations** — Hartree ρ ball ≈ 2–3·α_max (charge-converged, measured);
  the ρ FED TO the XC nonlinearity ≈ 4–8·α_max (a POINTWISE non-negativity requirement, not spectral:
  Gibbs lobes + the ρ>0 guard = the Exc collapse; the C=8 calibration's real content); the v_xc OUTPUT
  ball ≈ ρ's/3 (LDA v_xc ~ ρ^{1/3}: the cube root of a Gaussian peak has exponent p/3 — SMOOTHER than ρ,
  opposite to the GGA ∇ρ lore).  The governing ball sets the raster (∝ Ecut^{3/2}), so C 8→4 ≈ 2.8×
  fewer raster points machine-wide.  EVIDENCE C=8 is over-conservative for production: the matched NaF
  runs at Ecut=160 (C=4, CP2K's own operating point) agree with CP2K to 0.1–0.2 mHa; the 480-ball probe
  bought 0.15 mHa.  **(f1) MEASURED 2026-07-23 — VERDICT: lowering the default C is REJECTED on the
  current XC path.**  Seeded pure-aufbau NaF SR2 (grid-continuation harness, GC_SEED_MOM=0), fine stage
  C-sweep:  C=8 (Ecut=320) → −24.4337, 22 iters, negCharge −9.5e−3 (clean; 2.5 mHa from the CP2K SR2
  −24.4312, itself an Ecut=160-class number);  C=4 (Ecut=160) → **XC-COLLAPSE from a SEEDED PHYSICAL
  START** (−41.8 at cap, negCharge −91, Exc −109 — the basin is reachable even without the ionic-seed
  descent);  C=3 (Ecut=120) → converges but DIRTY (−24.477 = 43 mHa overbound, negCharge −0.196).  So
  the "C=8 over-conservative, CP2K runs at C≈4" evidence does NOT transfer to our Fourier-round-trip
  XC feed — the earlier matched-160 agreement rode the forced FULL CP2K ladder + recipe, not the auto
  path.  The C 8→2-3 prize is REAL but gated ENTIRELY on (f2) below (the DM-ρ raw XC feed); re-run
  this exact sweep as (f2)'s acceptance gate.  CAUTION: the retired GPW_CDFIT_SCALE two-grid fork —
  any re-split must buy real money over one-grid simplicity.
  **(f2) BUILT + ACCEPTED 2026-07-23 (three increments, commits 84592979 / f2-inc2 / f2-inc3).**
  Mechanism as designed: (inc 1) `G_ERI3::applyRaw/applyRawAdjoint` — ρ_DM(r) on the integration raster
  (finest level RAW, others spectrally transferred, `TransferBand` zero-pad/truncation transposes; FD
  adjoint rel ~1e-8, gate `GPW.RawXCConsistencyFD`); (inc 2) `FourierDensity::GetRhoOnGrid`
  (empty-means-fallback) + `PW_XC` raw-first with the E/H pair always from ONE discrete functional (PW
  bit-identical); (inc 3) the SCF DYNAMICS — the ρ̃-mixing working density is raw-shadowed
  (`ApplySpectralFilter` full-box smooth Kerker on the raster; `FourierMixCD` carries it; Kerker+Pulay
  evolve it with the SAME algebra/coefficients; late-activation for SAD seeds; PW drops out empty) —
  inc 2 alone left the dynamics ball-fed (measured: C=4 still collapsed with a clean FINAL eval).
  **ACCEPTANCE (seeded pure-aufbau NaF SR2): the XC-COLLAPSE BASIN IS REMOVED — negCharge == 0.000 at
  C=8, 4, AND 3 (ρ_min ≥ −1e-3 e/bohr³ worst case).**  C=8 → −24.4325 (22 iters; 1.3 mHa from CP2K);
  **C=3 (Ecut=120) → −24.4310, 22 iters, 0.2 mHa from CP2K's −24.4312** — the C 8→3 grid unlock is REAL
  on accuracy.  BONUS: the raw feed removed the ball-Gibbs noise from the XC residual and the COARSE
  Ecut=40 stage now E-flat CONVERGES in 45 iters (ball era: 515).  CAVEAT → 0h: C=4 (Ecut=160)
  converges NO-collapse but lands a NEAR-GAPLESS wrong state under pure aufbau (gap ~2e-3 Ha, frontier
  occupation swaps, 95 mHa off) — an OCCUPATION-SEAM instability at that discretization, not XC; so the
  DEFAULT cutoffFactor stays 8 until the 0h occupation guard (or Fermi smearing) lands — the C 8→3 flip
  is then a one-line default + suite re-pin.  Historical design note (the original insight, kept for
  the record): feed XC the DM-ρ, pointwise NON-NEGATIVE by construction (PSD D ⇒ φᵀDφ ≥ 0) — the C=8
  cleanliness constraint dissolves.  **DEFAULT FLIPPED C 8→2 (2026-07-23, user-approved: C=2 resolves
  the density at its own exponent 2α_max — the natural constant now that the guard-rail compensation is
  gone).  En route the ladder's TOP-RUNG GATE was found borrowing cutoffFactor as its calibration
  (flipping C silently dropped the rung at explicit Ecuts → GPW.LocalPPKappaSelfConverged failed at
  0.087); the gate now has its own fixed constant kRungGateC=8 — DECOUPLED, all explicit-Ecut ladders
  bit-identical.  NaFRocksaltGamma re-anchored −24.4357±0.01 and re-tuned to α=0.25 (user measurement:
  converges 20 iters / 63 s — the ball-era α=0.025 was calibrated against limit cycles that no longer
  exist).  Suite 565/565.**  NOT via op(r) sampling (that is the deleted PhiOnGrid era, ~1e9 exp/iter) — the
  collocation STREAMS already materialise exactly those χᵢχⱼ(r) samples: the RAW D-weighted level
  densities before the FFT/ball combine ARE ρ_DM(r) to screening-ε (worst negatives ~1e-10, not the
  ball's Gibbs −0.77 e).  So: keep the fine level RAW for the XC feed (skip the ball truncation there
  only), spectrally upsample the coarse levels (their content is genuinely band-limited — benign),
  keep the BALL for Hartree/Poisson (variational, exact).  H_xc=∂E_xc/∂D via the existing box-gather
  adjoint with box-truncation replacing ball-restriction per level (the suspended §0f increment-1
  design, RESURRECTED with the right motivation: not an accuracy fix — the falsified role — but the
  C=8→2-3 unlock); re-gate `GPW.XCPotentialConsistencyFD` + negCharge probes.  This is CP2K's own
  arrangement (XC on raw collocated values).  BEFORE (a) because the raw-collocation XC feed is the
  same softening CP2K relies on to run clean on ball-only rasters — validating BallOnly against the
  current Fourier-round-trip XC path would measure a configuration (f2) deletes.  **DO TOGETHER with
  §5's fit-basis ctor ISP split (user 2026-07-23): the Vxc-fit `({G}_vxc, integration grid)`
  two-argument ctor is the natural vehicle for the per-role balls — if making that ctor explicit is
  easy while in here, pull it forward from §5 rather than riding the CreateCD/VxcFitBasisSet plumbing
  and re-touching the same seam later.**
- **(a) The raster POLICY — the ~8× lever.**  Full design in §0.5(a) just below.  AFTER (f): run ONE
  combined calibration matrix {RasterPolicy × cutoffFactor C} through the negCharge/XC probes — the
  two levers multiply (~8× × ~2.8× raster points) but share the single Gibbs-into-XC failure mode that
  (f2) dissolves, so the validation campaign is shared instead of run twice.
- **(e) Cache2/3 byte-budget LRU** (§5, user-approved): also the robustness fix; runtime-relevant because
  it retires the per-pair `ClearGeometryCaches()` rebuild cost on healthy bases (currently unmeasurable,
  but the LRU makes the policy principled).
- **(d) B_ij(R) k-independent 1E memo** (user design: cache B(R), never M(k) — "keep k out of the key"):
  `LatticeSum` currently folds `phase(n)` into the accumulation, burying the k-independence; storing the
  per-pair, per-offset reductions once makes every additional k-block's 1E build a ~ms phase contraction
  (the `IntegrateMemo` pattern, already shared across k-blocks via the one molecular basis).  The
  multi-k enabler (Si 2×2×2: 8 builds → 1 build + 8 contractions).  LAST here: orthogonal to the
  raster items (1E build, not raster-scaled) and only pays on multi-k runs — schedule with/into §4
  (IBZ); on the current Γ-dominated benchmarks its near-term leverage is small.
- (OpenMP pair loops remain opt-in `GPW_OMP_THREADS`, ~1.7× — bandwidth-bound; CP2K threads on top of
  everything, so parity ultimately needs the structural items above first.)


### 0.5(a) RASTER POLICY — the remaining ~8× raster factor vs CP2K (a designed choice, likely a knob)
At the SAME Ecut ball our raster is 72³ where CP2K's is 36³ (~8× the points).  This is not waste by
accident but a POLICY difference — and now that everything else is matched, it is the whole remaining
grid-cost gap:
- **Ours — ALIAS-FREE (difference-set) rasters:** `AutoGrid` = 4m+1 per axis (m = the ball's max index),
  so the raster resolves the full DIFFERENCE set \f$\{G-G'\}\f$: the product of ANY two ball waves
  (bandwidth 2m) is sampled exactly, the FFT of the collocated ρ gives EXACT ball coefficients, and the
  raster is a true quadrature for every \f$\langle G|f|G'\rangle\f$ with f in the ball.  Discretization
  error lives ONLY in the ball radius (Ecut) — N is never a physics dial.  This is why the 5-smooth flip
  re-pinned nothing.
- **CP2K — BALL-ONLY rasters:** N ≈ 2m+1-class (its 36 at the 160-Ha ball).  Products of two ball waves
  ALIAS on that raster — the fold-back into the ball is ACCEPTED as discretization error, controlled by
  converging CUTOFF (and softened by evaluating XC on the raw raster values).  The bet: the aliased
  product tails are the same \f$e^{-E_{cut}/2p}\f$ tails the CUTOFF calibration already budgets for, so
  paying 8× in points to capture them exactly at FIXED Ecut is wasteful — better to raise Ecut a little
  on a cheap raster if needed.
- **PROPOSED (user 2026-07-22: "sounds like a user knob"): a raster POLICY enum, not a numeric dial**
  (the no-grad-student-knobs rule): `RasterPolicy { AliasFree /*default*/, BallOnly }` at the factory
  surface, printed by `ReportGrids`.  AliasFree stays the default (correct-first; N provably not a
  physics variable).  BallOnly is the measured-efficiency option — ANOTHER ~8× on every raster-scaled
  cost (the NaF fine raster would drop to CP2K's own 36³-class, est. run ≪ 60 s).
- **A/B MEASURED 2026-07-23 — CP2K'S BET CONFIRMED on the raw-XC landscape.**  Mechanism landed:
  `RasterPolicy { AliasFree /*default*/, BallOnly }` on `PW_Evaluator`/`PW_Grid_Evaluator` (AutoGrid
  4m+1 vs 2m+1, both 5-smooth), threaded through every GPW grid (density + ladder + top rung), printed
  by `ReportGrids`; A/B instrument `GPW_RASTER_POLICY=ball` (factory-surface promotion rides the
  default decision).  RESULTS (raw XC + C=2 + guards; negCharge == 0.000 EVERYWHERE — the raw feed's
  pointwise samples are raster-independent, exactly the prediction; convergence behaviour unchanged):
  Si Γ −7.11514 (0.13 mHa from AliasFree, CLOSER to CP2K −7.11506);  NaF auto/Ecut=80 −24.4304
  (0.9 mHa);  NaF Ecut=320 −24.4311 (**0.1 mHa from CP2K −24.4312**);  the only real price is the
  C=1-class regime (Ecut=40: −43 mHa aliasing) — BELOW the C=2 default floor.  So at/above the default
  floor BallOnly costs ~1 mHa and buys another ~8× raster points (NaF total vs this morning:
  100³ → 25³-class ≈ 64× with the C flip).  DEFAULT-FLIP DECISION pending (user): flipping moves every
  raster-derived anchor (one more re-pin wave, the 5-smooth-flip precedent); calibration recorded
  either way.


## 0h. SCF-strategy guards — PULLED AHEAD of 0.5(a) (user-approved 2026-07-23): the C 8→3 default flip
## (~2.8× raster, measured clean by the (f2) sweep) is gated on THESE guards, not on more grid work; and
## (a)'s BallOnly A/B should run on a stable occupation footing.  Also in scope now: the C=4 near-gapless
## aufbau instability from the (f2) sweep (occupation swaps at gap~2e-3 — the smearing/occupation seam's
## first concrete customer).
- **MOM cross-grid guard — DONE 2026-07-23.**  `HomoLumo` fixed (εL over ALL unoccupied + `hole` flag +
  `[HOLE: non-aufbau]` marker — the old above-the-HOMO-index scan masked the diagnostic);
  `tSCFWaveFunction::ReleaseMOMReference()` (drop the reference + re-arm the delayed-IMOM capture);
  iterator guard: 3 consecutive hole iterations → WARN + release + convergence veto (≤2 releases/run);
  NEVER SILENT: any run ENDING non-aufbau warns on cerr regardless of recipe.  VERIFIED on the 0h repro
  (GC_SEED_MOM=1 GC_FINE_MOM_START=1): fine lands −24.43252 in 16 iters — identical to the aufbau pin to
  8 decimals (was −23.680, +0.75 pinned).  **DISCOVERY: the guard caught the COARSE stage's own recipe
  red-handed — its capture-at-fill-10 reference pinned a +0.75-class excited state in EVERY earlier
  measurement (the −23.68/−23.69 "coarse fixed points"); the honest raw-XC Ecut=40 fixed point is
  −24.4357, only 3.2 mHa from the Ecut=320 answer.  Under raw XC + the guard, the Ecut=40 (C=1!) grid is
  NEARLY CONVERGED on NaF — which reframes the C-default question entirely: the calibration sweep for
  the default C should now scan DOWN TO C~1-2, and grid continuation may be unnecessary for production
  tolerances.**
- **`ReportBandGap` hole-masking fix — DONE 2026-07-23** (folded into the guard work above).

## 0i. Analytic V_local LONG — DEMOTED to robustness/perf (was the §0e-PP crux)
With (a)+(b) landed, the SHORT is analytic in production and the grid LONG is standalone-exact to the
ladder top (κ rule) — NaF sits 0.19 mHa from CP2K, so the analytic LONG is no longer an accuracy
blocker.  Kept on the list for ROBUSTNESS (very hard PPs / deliberately cheap ladders, where ladder-top
saturation shows — the §0e-PP saturation corollary) and to delete the LAST V_loc grid sweep entirely.
The recorded crux, unchanged: the LONG's `−Z_ion·Repulsion3C(χᵢ,χⱼ,g_core)` lattice sum (the erf-Coulomb
IS a normalized Gaussian core charge, exp `1/2r_loc²`) is conditionally convergent (erf→1/r Madelung
tail) ⇒ needs the G-space/Ewald neutralizing background, NOT a real-space sum — i.e. CP2K's ρ_core-into-
the-Poisson-solve arrangement, which for us means folding the core charge into `PW_Hartree`'s existing
G-space solve rather than assembling a separate matrix.  Kernels all exist (`GaussianRF`, no new Boys
function); verification gates: Si Γ −7.11506 + NaF SR2 −24.4312 (both codes' clean oracle pair).

Then the standing queue: **(1) DROP SR** (rank-reduction + auto-tol, below); **(2) low-q multi-species
bases → Si/NaF/CsI**; **(3) CP2K reference library**; **(4) IBZ**; **(5) cleanups**.

## 1. ROBUST HANDLING OF DIFFUSE BASIS FUNCTIONS (was "DROP SR"; demoted from critical path 2026-07-23)
**Goal: a grad student can add diffuse functions AT WILL; the code detects near-null overlap modes,
drops them in the ORTHO TRANSFORM (never the basis), WARNS loudly, and everything downstream just
works.**  Now known SAFE (the dropped modes carry ~1 mHa — the SR↔SR2 measurement) and known NON-URGENT
(λ~1e-6 runs stably via Cholesky + the seeded-aufbau recipe — the §1 probe).  Essential steps:
1. **Rectangular V through the periodic stack** (the real work; unblocks `DISABLED_NaFFullBasisEigenTol`):
   a truncated ortho (V is n×(n−k)) must flow through `Crystal_EC` (band count n−k), `cDM_CD` (density
   stays full n×n via C=V·U′), and the collocation — mirror the molecular path's rectangular handling.
2. **Auto-tol via LASolver GAP DETECTION** (pure LA): force-drop d[i]≤0; scan the low spectrum for the
   largest consecutive ratio; ratio > R_threshold (default 30, exposed at the Calculation facade) = a
   clean gap → cut there; else ε-tol fallback + WARN.  `orthoTol<0`=auto / `=0`=none / `>0`=explicit
   (mirrors `densityEcut`).
3. **NEVER SILENT**: every auto-drop reports count + gap ratio + clean/ambiguous on `cerr`; a λ_min /
   condition-number line joins the run-start diagnostics (`ReportGrids` sibling) so conditioning is
   visible BEFORE it bites.
4. **Default-path SCF robustness** (with 0h): the seeded-aufbau/MOM-guard recipe that tames these bases
   lives in test env knobs today — promote the working policy into the facade defaults.
5. **Gates**: full-basis (VALENCE_LOWQ) NaF == the SR/SR2 answers ± the dropped-mode mHa; Si anchors
   untouched; the never-silent WARN asserted in a test.
Vision: collapse to ~one CP2K-like ε.  (Auto-Rcut half is DONE — enumeration is ε-complete in-seam.)


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
- **ENV-KNOB GRADUATION (user 2026-07-23): the env parameters were measurement-campaign quarantine —
  now that SCF/MOM/DIIS/ρ-mixing are understood, move the settled ones into typed params.**  Sort:
  (1) SCF POLICY (α, KerkerG0, Pulay depth/start, MOM start — already `SCFParams` fields; delete the
  `NAF_*`/`GC_*` env OVERRIDES and hard-code the tests' recipes, possibly as named `SCFParams` presets
  or nested structs); (2) VERIFICATION INSTRUMENTS (`GPW_MGRID_ECUTS`, `GPW_RELCUTOFF`,
  `GPW_ILLCOND_*` — CP2K-oracle matching only; keep env but mark as instruments, or fold into a test
  helper); (3) OPS/MEMORY VALVES (`GPW_STREAM_BUDGET_PTS[_F32]`, `GPW_OMP_THREADS` — emergency valves;
  env defensible, else a small runtime-config struct).  Do after the C=2 re-pin settles so the recipes
  being frozen are the final ones.
- **ISP-split the fit-basis ctors (user design, 2026-07-23 — final form).**  `PlaneWaveFit_IBS` is
  built from the whole `PW_Grid_Evaluator` (ball + raster + the orbital tier — the last is pure
  baggage).  The honest signatures differ PER ROLE:
  - **CD/ρ fit ← `{G}_ρ` (a `PW_Ball`: recip, k, Ecut, members) ALONE.**  The alias-free raster
    computes the ball coefficients EXACTLY and is DERIVABLE from the ball (`FFTGrid()` already does) —
    it is an implementation detail of the COLLOCATION ENGINE, not a property of the ρ fit basis.
    (The raster round-trip IS the fast evaluation of the analytic projections — per-pair boxes + one
    FFT vs ~1e10 closed-form Gaussian-FT evaluations at NaF scale — so it stays, but INSIDE.)
  - **Vxc fit ← (`{G}_vxc`, integration grid) as TWO arguments.**  v_xc(ρ(r)) is NOT band-limited
    (the nonlinearity), so its projection is a genuine QUADRATURE with a real accuracy choice — the
    integration grid is an independent degree of freedom (and the natural GGA densification hook),
    deservedly explicit at the seam.
  Makes ball-vs-raster structural instead of documentary (this week's confusion is the evidence), and
  dovetails with §0.5(f)'s per-role ball calibrations.  **SCHEDULING (user 2026-07-23): if the explicit
  Vxc-fit ctor is easy while implementing §0.5(f), do them TOGETHER.  ASSESSED during (f2) 2026-07-23:
  DEFERRED — the mechanical split is easy to WRITE but impossible to GATE today (no path distinguishes
  the two arguments: LDA has ball == grid, GGA is guarded out), so it would ship as an unexercised
  degenerate seam.  Land it with its FIRST REAL CONSUMER: the GGA relCutoff increment — which, note,
  post-(f2) densifies the raw feed's INTEGRATION GRID (∇ρ resolution), not the ball; the two-argument
  ctor is then load-bearing and testable.**
- **Cache2/Cache3 BYTE-BUDGET LRU + per-cache RAM report (user-approved 2026-07-23; the intended
  REPLACEMENT for the clear-based band-aid).**  The MnD geometry caches (Ω/RNLM/H3) currently stay
  correct on lattice paths only because the drivers call `ClearGeometryCaches()` per pair and the 3C
  kernel ships a duplicated `Overlap3CStream` — cache-then-clear is the WRONG shape (user).  The clean
  design: give Cache2/Cache3 a BYTE budget with LRU eviction (`RAMsize()` plumbing exists), and let the
  ALGORITHM (not the user) select the policy per scope — a lattice driver pushes a scoped tiny budget
  (size-1 preserves the `const&`-returning contract AND the within-triple component reuse at O(1)
  memory), molecular paths keep the generous default.  One mechanism with per-cache specificity
  (Ω/RNLM/H3); the `ClearGeometryCaches()` calls and `Overlap3CStream` duplication then RETIRE.  Also:
  extend the end-of-run `IntegralsCache RAM usage report` with per-cache byte sizes (growth visible in
  every log instead of discovered by OOM).
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
