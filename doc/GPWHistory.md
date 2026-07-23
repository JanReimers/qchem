# GPW History — archived DONE narratives & resolved investigations

Moved out of `doc/GPWPlan.md` (2026-07-14) to keep the plan lean: the PLAN carries the compact timeline,
the current state (the CP2K recipe + the analytic-rewrite record), the TODO, and the durable pins; THIS file
preserves the full per-increment narratives, diagnostic transcripts, dead-end records, and resolved
investigations for archaeology.  Nothing here should be needed to orient a new session.

---

## Increment 1 — periodic Gaussian 1E integrals at Γ (`ab2c6a76`)
- `GPW_Evaluator` (`src/BasisSet/Lattice_3D/Evaluators/GPW/`) satisfies `isPW_1E_Evaluator`; `GPW_IBS`
  (`src/BasisSet/Lattice_3D/GPW_IBS.C`) = `EPW_Orbital1E_IBS<GPW_Evaluator>` + identity. Scalar = **`dcmplx`**.
- Overlap/kinetic(`⟨p²⟩`)/nuclear are real-space Bloch lattice sums, **delegated** to the molecular Gaussian
  basis via the engine-neutral capability **`Molecule::LatticeSum1E`** (`src/BasisSet/Molecule/LatticeSum1E.C`),
  realised by `PG_Cart::Orbital_IBS` → `PG_Cart_MnD::NR_Evaluator` (analytic McMurchie–Davidson kernels +
  `GaussianRF::AtCenter`). GPW reaches it by an abstract→abstract cross-cast (no Gaussian internals cross into
  qcLattice_BS). New library edge `qcLattice_BS → qcMolecule_BS` (no cycle). libCint is the faster follow-up.
- Validated (`UnitTests/GPW_UT.C`): home cell `R={0}` reproduces the finite matrices `<1e-12`; images give
  textbook large-cell convergence.

## Increment 2 — DFT-tier collocation (Hartree/XC machinery) (`cc123b3b`, `63fbf70c`)
- GPW satisfies `isPW_DFT_Evaluator` and reuses the **entire** PW-DFT stack (`PW_Hartree`/`PW_XC`/`IrrepCD`)
  by filling the `Repulsion3C`/`Overlap3C` tensors with dense collocation weights `W_c(i,j)=(1/Ω)∫χ_iχ_j
  e^{-iG_c·r}`. `G_ERI3` gained `weights`; `ContractG_ERI3` branches. Tensor caching delegated to the framework.
- **Coulomb factorisation:** `W_c(i,j)` is a SINGLE-`r` integral (density side); the second electron + `1/r12`
  are the diagonal Poisson kernel `4π/|G_c|²`. Full Coulomb = weight × kernel, factorised through G-space.

## Increment 3 — first-light periodic SCF (`dcef8528`, `db314e6a`)
- Closed the last tier (external PP) by making **`GPW_IBS` realise `Integrals_Pseudo<dcmplx>`**, so `PW_Pseudo`
  and the **entire `Ham_PW_DFT` drive a GPW basis verbatim** through the real `cSCFIterator`. Zero new
  Hamiltonian code. `MakeLocalPP` = G-space form factor (ΔG=0 dropped, box-independent, PW alignment);
  `MakeSeparablePP` = KB projector via `qcMesh::Overlap` on `CreateIntegrationMesh`.
- **Validated (`UnitTests/GPW_SCF_UT.C`):** crystalline Si (Γ, FCC primitive, 8 val e⁻) converges, charge 8,
  **Etot = −8.24758**; Si pseudo-atom-in-box reproduces the finite SIPP molecular DFT to grid tolerance.

## Implementation 4 — general-k GPW (Step 1) + multi-k BZ plumbing (`b2a29249`)
- **General-k:** the `e^{ik·R}` Bloch phase runs through the stack. `Molecule::LatticeSum1E` now takes an
  adjacent `(Rs, phases)` pair (`cvec_t`) and returns `chmat_t` (Hermitian). `GPW_Evaluator` does complex
  Bloch `Eval`/collocation (`BuildWeights` **conjugates the i-slot** per `ρ=ΣD_ij χ_i*χ_j`, full n²), a
  complex Hermitian KS bridge, complex Bloch KB projector. New complex-input `PeriodicGridEvaluator::
  ForwardFFT(cvec_t)`. Phase = `exp(2πi k_frac·n)` (integer cell index — convention-safe). **{R} and
  {e^{ik·R}} are kept bundled/adjacent** (a future `qcMesh cMesh = Mesh<dcmplx>`).
- **Γ bit-identity held** (phase=1, conj no-op): the gapped Si-Γ crystal is unchanged.
- **Validated:** 4 matrix-level Bloch invariants in `GPW_UT` (k-invariance at Rcut=0; phase-is-live +
  Hermiticity at k≠0; Bloch translation law `χ^k(r+R0)=e^{ik·R0}χ^k(r)`; `S(−k)=conj(S(k))`).
- **Multi-k plumbing:** `GPW_BasisSet` iterates `lat.MakeKMesh()`, one `GPW_IBS` per k **with the BZ weight**
  (`BlochFactory(N,ik,kp.weight)` — a missing weight had given charge = Nk×Nelec). Gate
  `GPW_SCF.SiliconMultiKPlumbing` (2×1×1 Rcut=0 == Γ, charge 8, Etot −8.24758). A `collRcut` decouples the
  collocation reach from the overlap Rcut (feasibility for the diffuse basis).

## Basis conditioning: SIPP_SR + N3/N5 removal (`b2a29249`, `10ad6e29`)
- The diffuse SIPP test basis (Si s=0.09/p=0.06, RMS ~5 a.u.) goes near-linearly-dependent when Bloch-summed
  in a solid: **min eig(S(k)) = 4.3e-6 (SIPP) vs 0.0164 (SIPP_SR)** (drop the 2 most diffuse). `sipp_sr.bsd`'s
  overlap is PSD + converged at **Rcut=1.5a** (vs SIPP's 3a, still near-singular); with it the dispersive SCF
  is numerically STABLE (no divergence). **Lesson (durable): an ill-conditioned overlap is a BASIS problem,
  not a solver/code bug.** Consequently **N3/N5 were removed** from `BasisSetAccuracy` (now {Low,Medium,High});
  the UTAtom_BS tests that used them migrated behaviour-preservingly to inline-JSON `N3Basis/N5Basis` helpers.

## Bulk over-binding ROOT-CAUSED (`a4c94ec5`) — the atom-on-FFT-raster-node bug
- With a well-conditioned basis the dispersive-bulk SCF converges but to Etot ≈ −15 (≈ 2× PW −7.76). Ruled
  out in turn: **charge = 8** (not a double-count); the **density collocation is consistent** with the
  analytic overlap once images restore the corner atom's leaked density; **kinetic + separable-PP matrices
  unchanged** with images.
- **Cause:** `GPW_Evaluator::OverlapMatrix(V)` (local-PP + Hartree + XC integrate-back) quadratures on the
  **FFT raster `A·(i/N)`**, where a lattice-point atom (the FCC corner atom at 0) sits EXACTLY on a grid node
  → its sharp density peak is over-weighted against the deep PP well (Vloc trace −29→−52 with images; Een
  −1.06→−16.6). **Decisive:** shift all atoms by ⅛ cell → Etot −15.2→−8.4 (must be invariant). The
  **separable PP is immune** (it already uses the offset qcMesh MIDPOINT mesh `A·((i+½)/n)`); **raising
  densityEcut does NOT help** (r=0 is a node at every N).
- Landed: DISABLED diagnostics in `GPW_SCF_UT` (`SR_TranslationInvariance`, `PPMatrixTraceProbe`,
  `CollocationVsAnalyticOverlapWithImages`, `SR_CornerAtomVsDensityEcut`) + the cause documented in
  `OverlapMatrix`. (A partial fix — midpoint mesh for `OverlapMatrix` only — was explored + reverted:
  incomplete, and it moved the committed anchor.)

## Bulk over-binding FIXED — GPW bulk matches CP2K to 1e-5 (was TODO 1) (`95e8f4a8`)
The root cause was **one thing wearing two costumes: an incompletely-wrapped Bloch orbital.**
- **The real bug (KB nonlocal, the 16 Ha term):** `MakeSeparablePP` used the **raw home orbital `*itsOrb`** as
  the projector bra on the single-cell mesh. A boundary-straddling corner atom lost its wrapped tail → `b_i`
  ≈ half (corner trace 21 vs interior 37) → the nonlocal PP was translation-variant by ~16 Ha. **Fix: use the
  Bloch-summed orbital (`Eval`, precomputed on the mesh) as the bra.**
- **The FFT-raster `Vloc`/Hartree/XC term was a RED HERRING:** once the orbital is fully wrapped (`Rcut ≥ 2a`)
  its translation-variance also vanishes (the on-node over-weighting self-corrects when the full periodic
  density is present). So **the voxel-grid-shift (old TODO 1b, Option A) was reverted entirely** — simpler.
  Both terms go to Δ = 0.0000 at `Rcut ≥ 2a` (`GPW_SCF.DISABLED_TermTranslationInvariance`).
- **Validation vs CP2K (Γ, SIPP_SR, Rcut=2a):** Etot **−7.11505** (CP2K −7.11506), charge 8, Exc −2.544
  (CP2K −2.544). Nonlocal-PP term hits CP2K's +0.9406. **Committed anchors safe:** at `Rcut=0`, `Eval` = the
  raw orbital, so `SiliconGammaConverges` (−8.24758) and the atom-in-box are unchanged.
- **Perf:** cached `PhiOnGrid` (geometry-fixed; was recomputed every SCF iteration) → the CP2K gate dropped
  ~25× (1100 s → ~45 s at N=32/`densityEcut=20`). A GEMM quadrature was tried and reverted (not faster at
  n=13). Gate `GPW_SCF.DISABLED_SR_GammaRcut2a_CP2KReference` (N=32, −7.11467, ~0.4 mHa grid gap, tol 2e-3).
- **Test cleanup:** removed 11 obsolete over-binding-investigation diagnostics (`GPW_SCF_UT` 541→268 lines).

## Multi-k GPW dispersion VALIDATED vs CP2K (`5fe61aeb`)
Dispersive multi-k GPW runs (unblocked by the KB fix): Γ-centred 2×2×2 MP, SIPP_SR, Rcut=2a → charge 8, real
dispersion (Γ −7.11467 → 2×1×1 −7.451 → 2×2×2 −7.7778). **Grid-for-grid at the SAME Γ-centred mesh: our
−7.7778 vs CP2K −7.77846 (~0.7 mHa, the N=32 grid gap).** So the general-k GPW *physics* is validated. The
90 mHa vs CP2K's *default* −7.86744 is purely the **k-mesh CONVENTION** (Γ-centred vs the classic shifted MP,
k at ±¼ — confirmed from CP2K's own k-point list). Decks: `si_fcc_gpw_222.inp` (shifted) + `si_fcc_gpw_222_
gamma.inp` (Γ-centred). Test `GPW_SCF.DISABLED_SR_2x2x2GammaCentred_vs_CP2K`.

## Shifted Monkhorst-Pack support (`1980d6ef`) — and it EXPOSED the next bug
Threaded an optional fractional MP `shift` so `k=(ik+shift)/N` through `BlochQN`/`BlochFactory` →
`Lattice_3D::MakeKMesh` → `GPW_BasisSet`/`GPWFactory` → `RunGPW` (shift=0 = Γ-centred, backward-compatible;
shift=½ = CP2K's `k=±¼`). `GPW_BasisSet` recovers the integer index as `lround(k·N − shift)` (plain
`lround(k·N)` is wrong for shift=½). 186/186 green; Γ-centred anchors unchanged. **But running the shifted
mesh exposed two complex-Bloch-phase bugs — now FIXED, see next.**

## Complex-k GPW FIXED — CP2K default shifted 2×2×2 matches −7.86744 (was TODO 1) (`745d03ff`)
The shifted mesh (k at ±¼) is the **first genuinely-COMPLEX Bloch phase** (`e^{ik·R} ≠ ±1`), so D and every
k-block matrix are genuinely complex. It over-bound (single k=¼ block: Een → −18.9, Etot → −15.2, no
convergence). The plan's own localization was **WRONG** — it blamed the shared framework complex-D path
(`cSCFAcceleratorDIIS`/`Crystal_EC`/`cDM_CD`) and cleared "the density collocation" and "the GPW evaluator".
In fact **BOTH bugs were in the GPW evaluator** (`src/BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C`);
the framework complex-D path was correct all along (it had just NEVER been run at complex k — PW-DFT's multi-k
tests are all Γ-centred too, so this was its first genuine exercise).
- **Bug 1 — collocation density convention (`BuildWeights`).** The weight conjugated the **bra (i)** slot
  (`conj(Φ_i)·Φ_j`), making ρ̃ the TRANSPOSE-density `Σ D_ij χ_i* χ_j` — a *different real field* at complex k.
  The physical density is `Σ D_ij χ_i χ_j*` (= `IrrepCD::operator()` `trans(φ)·D·conj(φ)`, = the PW delta
  path, = `Σ_occ|ψ|²`). Fix: conjugate the **ket (j)** slot. The plan's `‖W₀·Ω − S(k)‖ = 4e-6` "rules out
  collocation" diagnostic was a red herring — it checks the overlap *integral*, not the *D-contraction slot*.
- **Bug 2 — KB projector image phase (`MakeSeparablePP`), the dominant over-binder.** The projector-image sum
  used `e^{+ik·R}`; the correct Bloch projection `b_i = ⟨χ_i^k|β_home⟩` tiles all-space (`∫_all f = Σ_R
  ∫_cell f(·+R)`) and the Bloch law `χ^k(r+R)=e^{ik·R}χ^k(r)` puts a **conjugated** `e^{−ik·R}` on the
  R-shifted projector. At complex k this **halved the nonlocal-PP trace** (`TrVnl` 42→22 at k=¼) → a spurious
  deep core level (−3.79) → over-bind. Fix: `ph = conj(itsPhaseC[r])`.
- **Bonus — `IrrepCD::GetTotalCharge`** used `sum(D % S)` (= `Tr(D Sᵀ)`), the exact anti-pattern its sibling
  `DM_Contract` documents; corrected to `sum(D % trans(S))` = `Tr(D S)`. No-op at real k.
- **All three are inert at Γ / half-integer k** (phase ±1 self-conjugate, real orbitals) → every committed
  anchor byte-identical (Si Γ −8.24758, atom-in-box −3.73567, Γ-centred 2×2×2 unchanged). **Validation:** the
  single k=¼ block now converges (17 iters, Etot −7.565, physical); **the full shifted 2×2×2 converges (21
  iters, charge 8) to Etot −7.86673 vs CP2K −7.86744** (0.71 mHa = the N=32 grid gap). Gate
  `GPW_SCF.DISABLED_SR_2x2x2ShiftedMP_vs_CP2K` now asserts −7.86744 (disabled: ~5 min SCF). **Multi-k GPW over
  the full BZ (any k) is now DONE and CP2K-validated at both Γ-centred and shifted meshes.**

## NaF convergence campaign — DONE 2026-07-12 (correctness resolved; runtime optimisation is now TODO §0)
**OUTCOME:** every correctness axis closed — auto-floor `densityEcut` (`04e31a8e`), `∫ρ_grid`/fingerprint diagnostics (`3a87dba6`), diffuse ionic seed (`e1f986da`, PW iters 35→17), and Kerker ρ-mixing (`d66b7c8e`, Si-exact, NaF+DIIS converges). The code is correct in principle; the only remaining blocker is RUNTIME → TODO §0 (profile-first). The diagnosis/execution record follows.
A long diagnostic session got GPW NaF to first light (charge 8) but neither our GPW nor CP2K CONVERGES
cleanly on the low-q valence_lowq basis. The picture below is the reprioritised, corrected understanding to
start from (several of this session's early claims were wrong and are struck through — see the corrections).

### DIAGNOSTICS LANDED 2026-07-12 (disentangling infrastructure — uncommitted)
The strategy: stop reading the SCF endpoint as the instrument — the SCF does double duty (*find* ρ AND *score*
ρ), which is what tangles conditioning / fit-noise / charge-transfer / grid. Split every question into
"is the functional right?" (fixed-ρ, **zero-SCF**) vs "does the iteration find the min?" (dynamics). Then each
confound gets an orthogonal probe (everything else held fixed). Progress so far:
- **PROBE 1 — conditioning axis CLOSED (red herring, PROVEN).** `DISABLED_NaFOverlapConditioningSweep` now also
  reports the orthogonaliser residual `‖VᴴSV − I‖` (V=S^-½ from the SAME LASolver the SCF uses). At the SR/Rcut=2a
  operating point (min eig 7.5e-4, **cond 8252**) the residual is **2.5e-14** — machine ε. cond(V)=√cond(S)≈90, so
  even Fock inversion amplifies to ~1e-12. **The metric is NOT what ails the NaF SCF.** And the probe cleanly
  separates the TWO things that both read as "conditioning": (A) near-singular-but-PSD = red herring (above);
  (B) INDEFINITE S at a sharp `Rcut` (full basis min eig −0.42→−0.60→**−0.11 even at 2a**, never PSD in-window) —
  Gibbs ringing of the truncated Bloch autocorrelation, fixed by **magnitude screening, NOT a bigger Rcut** (Rcut=2a
  does not rescue the full basis; only SR barely escapes). Confirms the OPEN INVESTIGATION section's conclusion.
- **`∫ρ_grid − N` READOUT — DONE** (the §0 ask below, "ADD the same readout"). New opt-in toggle
  `qchem::Hamiltonian::ReportGridCharge()` (mirrors `ReportOverlapConditioning`), printed by `PW_XC::RefreshRhoGrid`
  each iteration: `[grid charge] integral rho_grid=… Tr(DS)=… lost=…` — the electrons lost to grid truncation
  (== CP2K's "Electronic density on regular grids: <int> <error>"). Wired into `RunGPW` (gated on `verbose`) + the
  NaF test. Si/Rcut=2a/dE=20 loses only −1.2e-5 e⁻ (soft PP, grid converged); **NaF is where it lights up** (F's
  tight 40-a.u. exponent). All GPW+PW anchors byte-identical (toggle off by default).
- **E[ρ] SEAM investigated → PROBE 2 (fitted-vs-collocation ΔE) is DEGENERATE in GPW; DON'T build it.** In GPW the
  energy is ALREADY collocated: `PW_XC::GetEnergy` = `∫ε_xc(ρ)ρ` by direct grid quadrature; `PW_Hartree` = exact
  W-tensor × 4π/G². So there is no molecular-style "fit energy" to diff. The non-variationality is an energy/gradient
  **INCONSISTENCY**: the SCF minimises the *projected v_xc MATRIX* (`PW_XC::CalcMatrix` fits the *nonlinear*
  v_xc(ρ) onto finite {G}), but the ENERGY uses direct full-grid quadrature — two discretisations of a non-band-limited
  field, consistent only as the grid resolves (SAME ROOT as the deferred `relCutoff`/denser-v_xc-grid cleanup, §4).
  ⇒ the sharp variationality instruments are the **E(λ) line-search** or an **FD potential-consistency check**
  (`dE_xc/dρ` vs the assembled v_xc matrix), NOT fitted-vs-collocation. The "Fock-first" gotcha = the serial-keyed
  freshness cache (`newCD`/`itsFitVersion`) handing back a stale ρ-grid to a colliding `Version()`; `MixIn`/`ReScale`
  bump the serial so synthesized densities are safe.
- **PROBE — NaF grid-charge diagnostic RAN, and it is DECISIVE. The NaF problem is GROSS GRID UNDER-RESOLUTION,
  not conditioning or the fit.** At `densityEcut=40` the readout shows `Tr(DS)=8.000000` EVERY iteration (the
  density MATRIX conserves charge perfectly) while `∫ρ_grid` **oscillates 2.4 → 7.2 → 4.8 → 3.2 → 2.7**, settling
  ~**2.8 — the grid holds under 3 of 8 electrons (>5 e⁻ lost off-grid).** F's tight 40-a.u. exponent produces a
  density (product exponent ~80) the `densityEcut=40 Ha` grid cannot represent, so Hartree+XC are built on a
  density missing ~65% of its charge → garbage potential → the SCF oscillates. This CONFIRMS the plan's long-held
  suspicion ("F is the hard atom, wants a fine grid; our 16–24 Ry is ~10× too low vs CP2K's 300–600 Ry") with a
  HARD number. **The dominant NaF fix is densityEcut ≈ 200 Ha (CP2K's CUTOFF 400 Ry), or multi-grids (§4)** — NOT
  a better seed/mixing/basis (those are second-order until the grid holds the charge). Seed-independence is now
  LOW priority (both seeds hit the same grid wall).
- **Observer trajectory FINGERPRINT — DONE + validated** (`Fingerprint()` in `GPW_SCF_UT.C`, fed by `SetObserver`).
  Classifies a run in one line by its time-series signature: CONVERGED / DENSITY-DEGENERATE (E settled, ρ rotates —
  benign, the Γ open-shell atom) / OSCILLATING (charge-transfer sloshing) / FIT-FLOOR STALL / DIVERGING. Self-checks:
  Si Γ crystal → CONVERGED (relAmp 4e-8); Si atom-in-box → DENSITY-DEGENERATE (relAmp 3e-3) — correctly NOT flagged
  as sloshing. `RunGPW` gained a `seed` param (Uniform default; CoreGuess is the other real GPW seed — SAD/IonicSAD
  fall back to Uniform on dcmplx) so seed-independence is a one-liner when wanted.
- **NaF densityEcut LADDER RAN (SR/Rcut=2a, 40 → 120) — GRID CONFIRMED as the dominant cause, and it disentangles
  a SECOND axis.** Grid loss `|∫ρ_grid − 8|` fell **5.2 e⁻ (Ecut=40) → 0.24 e⁻ avg, with the well-resolved
  iterations hitting 7.997 (loss 2.7e-3) (Ecut=120)** — a ~20× improvement, exactly the grid-resolution signature.
  Two consequences: (1) Etot moved −23.556 → **−23.982** and the DIIS↔GDM gap shrank toward GDM's −23.936, so
  **part of the earlier "non-variational" DIIS/GDM discrepancy was GRID-TRUNCATION NOISE, not functional
  non-variationality** — the fit-noise floor was overstated. (2) At Ecut=120 the run STILL oscillates (∫ρ_grid
  swings 5.2 ↔ 7.99 across iters, iter-capped at 60) — but now the grid is mostly resolved (good iters → 7.997),
  so **the residual oscillation is exposed as GENUINE charge-transfer dynamics**, cleanly separated from the grid
  noise that masked it at Ecut=40. So the axes have peeled apart: grid (dominant, ~fixed by cutoff) → then
  charge-transfer mixing instability (the real remaining SCF problem, NOW the right target for seed/Kerker/smearing).
- **AUTO CUTOFF IMPLEMENTED (2026-07-12): `densityEcut` is now BASIS-derived, not a user burden.** GPW has NO
  orbital/wavefunction `Ecut` (Gaussians are analytic) — `densityEcut` is its ONLY grid cutoff, and it is a
  DENSITY-scale quantity (the density is the product of two orbitals, exponent `2·α_max`, so its constant already
  folds in the ×2 over a single-orbital cutoff — do NOT confuse it with an orbital `Ecut`). New
  `Molecule::LatticeSum1E::MaxExponent()` (scalar α_max summary; realised by `PG_Cart` → `NR_Evaluator`, walks the
  radials — no primitive escapes) drives a **three-mode `densityEcut`** in `GPW_Evaluator` (threaded through
  `GPW_IBS`/`GPW_BasisSet`/`GPWFactory` with a new `cutoffFactor` param):
    - **`< 0` = AUTOMATIC (recommended):** grid = `cutoffFactor·α_max` — the caller need not know the Hartree value.
    - **`= 0` = DFT tier OFF** (1E-only; unchanged).
    - **`> 0` = EXPLICIT:** honoured as given, but **`cerr` WARNING** if below `cutoffFactor·α_max` (respects the
      expert's insisted-on value, does NOT silently clamp — consistent with "don't hide the problem").
  `cutoffFactor` default 4, choose `C ≥ 4` for a finer grid (calibrated: F α_max=40 → floor 160, in the ladder's good
  regime between 120 "good" and CP2K's 200 "converged"). **Inert on every committed anchor** — SIPP α_max=2 → floor 8,
  below all Si `densityEcut` (≥10), so Si Γ −8.24758 / atom-in-box / GPW+PW 30 tests byte-identical (verified: no
  warning, all green). NaF test switched to AUTOMATIC (`densityEcut=-1` → 160). Burden removed: the user never guesses
  the Ha value. (A GGA `relCutoff>1` Vxc densification sits ABOVE this floor, still deferred.)
- **CHARGE-TRANSFER OSCILLATION LOCALISED (2026-07-12, user-run GPW NaF at auto-160 + IonicSAD):** the diffuse
  seed + auto-floor grid give a PERFECT iter-1 (`∫ρ_grid=7.9998`, loss 2e-4), then the SCF DYNAMICS sharpen the
  density → it ALIASES off the 160-grid (`∫ρ_grid → 6.32`, 1.68 e⁻ lost) → slosh. **Root = DIIS activates TOO
  EARLY:** `DIISParams.EMax=8.0` starts DIIS as soon as `[F,D]<8` (≈ immediately), extrapolating from a thin
  history → over-concentrates the density. **GPW-SPECIFIC** (PW NaF is fine at EMax=8, 17 iters): PW is natively
  band-limited and CANNOT alias; GPW's grid-collocated density can. So grid (auto-floor) + seed (diffuse) both
  WORK — the residual is pure extrapolation dynamics. **EXPERIMENT (a) DONE — DIIS-early REFUTED.** `EMax=0.5`
  (delay DIIS) did NOT help: with ρ sloshing, `[F,D]` never drops below 0.5 so DIIS never even activates, yet PURE
  relax mixing STILL limit-cycles — a STABLE period-~6 cycle (`∫ρ_grid` 7.999→4.35→7.27→2.40→6.05→3.13→repeat,
  amplitude 5.6 e⁻; `Tr(DS)=8` always, so the swing is entirely the COLLOCATED grid density sharpening+aliasing,
  not the DM). **⇒ the instability is the density-mixing fixed point itself (charge-transfer sloshing); linear/
  relax mixing cannot damp it — NOT a DIIS/seed/grid problem.** (Runtime 34.5 min — the speed item is real.)
  **(b) KERKER / high-G-damped mixing is now the CONFIRMED fix** (damps exactly the short-wavelength sharpening
  that aliases — the GPW-shaped cure), expectation: break the limit cycle. Strategy (user): defer known fixes
  (magnitude-screening, Kerker) until the baseline is UNDERSTOOD, so each next fix has a clear expectation — met.
  **KERKER FOUNDATION DONE (2026-07-12), iterator wiring NEXT.** New `ChargeDensity::FourierMixCD` — a G-space
  density holding a raw `ρ̃(G)` map + reciprocal lattice, presenting the `FourierDensity` face (so
  `DoSCFIteration(ham, mixedρ)` builds the next Fock from it — the SAME seam the SAD seed uses; verified
  `DoSCFIteration` takes `const tChargeDensity*`, not a `tDM_CD`, so ρ-mixing needs NO Fock-build rework). Static
  `KerkerMix(in, ρ̃_out, α, G0)` applies `ρ̃_in + α·G²/(G²+G0²)·(ρ̃_out−ρ̃_in)`. **Unit-validated
  (`KerkerMix.*`, 3 tests, no SCF):** charge conserved EXACTLY (f_K(0)=0 → G=0 never mixed, so ∫ρ stays N even
  when the output aliased to 2 e⁻ — directly kills the ∫ρ_grid-collapse symptom); low-G damped / high-G passed;
  G0→0 = plain linear mixing. **ITERATOR WIRING DONE (2026-07-12):** `SCFParams.KerkerG0` (0=off, default) gates an
  optional ρ-mixing branch in `tSCFIterator` (dcmplx-only via `if constexpr`): `KerkerSetup` builds the G-space
  fit basis (`Band_FT_IBS::CreateVxcFitBasisSet`) + the initial `FourierMixCD` from the seed; the loop drives the
  Fock from `FockDensity()` (the mixed ρ̃ when active, else the working D) and `KerkerUpdate` re-collocates ρ̃_out
  and folds it in. **DEFAULT PATH BYTE-IDENTICAL** (KerkerG0=0 → `itsMixedRho` null → linear D-mixing everywhere;
  verified: Si Γ −8.24758, molecular M_Calculation, PW Si all unchanged). The `DoSCFIteration(ham, const
  tChargeDensity*)` seam (DIP) is what let a ρ̃-only density substitute for D with NO Fock-build rework.
  **DEBUGGING (2026-07-12, on FAST Si — not the 34-min NaF):** first run was bit-identical to no-Kerker →
  KerkerSetup was silently bailing (Release has no asserts): the iterator stored the raw `st` from
  `Lattice_3D::GetStructure()` which returns a TEMPORARY (`make_shared`), so it DANGLED by `Iterate` time and the
  UnitCell cast failed. FIXED: deep-copy the cell in the ctor (`itsKerkerCell`), and the setup now reports LOUDLY
  (cerr) instead of a NDEBUG-silenced assert. Kerker now ACTIVATES. TWO issues remain, both found on Si in ~10 s:
  (1) **convergence-metric bug (the real one):** the loop gates on `‖D_out−D_out_prev‖`, but with ρ-mixing `ρ̃_mix`
  changes SLOWLY so `D_out` does too → it reports "converged" while `ρ̃_out≠ρ̃_mix` (NOT self-consistent) → Si stops
  early at −8.24662 vs −8.24758. **FIX: gate ρ-mixing on the SCF RESIDUAL `‖ρ̃_out−ρ̃_in‖`, not the D_out change.**
  (2) the real fixed-point bug: **Kerker FROZE G=0** (`f_K(0)=0`, its plane-wave charge-conservation feature). In
  PW that's right (ρ̃(0)=N/Ω, fixed); in GPW ρ̃ is a fit-basis PROJECTION whose (0,0,0) is SHAPE-dependent (that's
  why `Ω·ρ̃(0)=1.38≠8`), so freezing it stranded the XC's mean density at the seed value → wrong fixed point
  −8.24662, residual floored. **FIX: mix G=0 fully (`f_K(0)=1`) — the SCF diagonalization conserves charge, so no
  freeze needed.** Also added the ρ-residual convergence gate (`‖ρ̃_out−ρ̃_in‖`) and made `FourierMixCD` carry N
  explicitly (the shape-dependent ρ̃(0) can't give it). **RESULT: Si-Kerker now converges to −8.24758 EXACTLY
  (26 iters) = the D-mixing fixed point** — the fast correctness gate PASSES, machinery sound. KerkerMix unit
  tests updated (G=0-mixes), default path byte-identical (Si −8.24758). **NaF RESULT (pure Kerker G0=1.0): PARTIAL
  — machinery works, tuning not there yet.** It transformed the chaotic period-6 cycle into a clean PERIOD-2 flip
  (`∫ρ_grid` 7.999992 ↔ 2.414158; the 8.0 state is beautifully resolved, loss 8e-6 — Kerker IS producing clean
  densities) and moved the energy toward physical (−20.10 vs PW −20.33, up from the garbage −24.5). But it's a
  STABLE period-2 cycle, not converged — UNDERDAMPED. Two tuning causes (NOT correctness — Si proved that):
  (a) **G0=1.0 too weak** — NaF's charge-transfer mode is inter-atomic (Na–F ~4.4 a.u. → G~1.4), ABOVE G0, so
  barely damped; needs G0~2–3. (b) **α ran away to 1.0** — the `if(FD<FDold) relax*=1.5` growth wasn't guarded
  for Kerker (now FIXED: Kerker holds α=StartingRelaxRo). **NEXT (tuning, but 34-min/run): G0~2–3 + fixed
  α~0.1–0.2; and/or re-enable DIIS (excels at breaking period-2 flips — Kerker damps amplitude, DIIS kills the
  flip).** The period-2 flip + near-physical energy say we're close.
- **KERKER + DIIS: CONVERGING — charge-transfer axis essentially CRACKED (2026-07-12).** The combo broke the
  period-2 flip: the ρ-residual Δρ decreases MONOTONICALLY in the tail (iters 54–60: 1.5e-2 → 3.6e-3, clean
  ~×0.8/iter), energy settled at −24.01 (near the basis-limited GPW value ~−23.6/−23.9 for VALENCE_LOWQ_SR — NOT
  the complete-basis PW −20.33). Hit nmax=60 at Δρ=3.6e-3 (just above the 1e-3 gate — ~4 more iters would
  converge); two mid-run DIIS crashes (iters 50–53) it recovered from. **CONCLUSION (user): the code is correct
  IN PRINCIPLE — no bug making charge jump around; the oscillation was real SCF dynamics, now tamed.** Kerker is
  built, unit-tested (`KerkerMix.*`), Si-validated to the EXACT fixed point (−8.24758), default-safe (191/191),
  and it converges NaF. **STOP heuristic mixing trials here** — the 34-min loop makes tuning nmax/G0/α/DIIS too
  expensive. **The next lever is RUNTIME (round-1 DONE below).**

## RUNTIME round 1 (zgemm) + conditioning (magnitude screening) — DONE (2026-07-13)
**Profiled NaF with `perf`** (9.86M samples): the #1 hotspot is the per-iteration integrate-back
`GPW_Evaluator::OverlapMatrix(Vtilde)` at **43.8%** (dense `M_ij=w Σ_p conj(Φ_pi)V_p Φ_pj`), NOT the grid/
screening the plan had assumed. FFT ~18%, orbital-eval ~14% (mostly one-time `PhiOnGrid`), W-tensor 0.9%.
- **`OverlapMatrix` → OpenBLAS `zgemm`** (`7708d2dc`): `M=w·Φᴴ·(V.∗Φ)` via `cblas_zgemm(ConjTrans)`. Isolated
  bench (`GPW.DISABLED_BenchOverlapMatrix`, NaF scale): scalar 12.09 → **3.00 s/call @1-thread = 4.0×**. blaze's
  own product is scalar (`BLAZE_BLAS_MODE=0` never reaches our TUs) → the direct cblas call. Full NaF **28:14**
  (the untouched ~40% FFT/setup dilutes it → that's the §0 target). `VPhi` scratch reused (`506f6a11`, hygiene).
- **OpenBLAS adopted + `openblas_set_num_threads(1)`** pinned in `gtestmain.C`/`scfrun.C` (`7708d2dc`): OpenBLAS
  auto-sizes threads by load → non-reproducible last-ULP drift (an SCF E moved >2e-5) → pin for determinism
  (keeps 4.0×, mostly SIMD). 4 over-tight anchors loosened for OpenBLAS-vs-netlib roundoff. **193/193 green.**
- **Magnitude screening on the 1E lattice sums** (`05e44fab`): per-component reach `√(−ln ε/α_min)` (ε=1e-10) in
  `NR_Evaluator::LatticeSum`, shared identically across S/T/V_nuc (consistency = correctness for `HΨ=εSΨ`).
  **Sparse → ~4×** on the 1E sums (a generous Rcut is now free — the 2a-tuning pain is gone). Roadmap `ac272432`.
- **Conditioning FINDINGS** (`c015b038`, `c96db327`): the full valence_lowq basis is intrinsically OVER-COMPLETE
  when Bloch-summed (min eig→0⁻) — screening EXPOSES this but cannot fix it; its ~1e-6 null cluster sits in a
  clean ~1000× spectral gap. Canonical Eigen/SVD ortho with tol in the gap is clean at the OVERLAP level
  (‖VᴴSV−I‖=6.6e-11) but the SCF is **BLOCKED** by a periodic-stack rank mismatch (37→33 "Matrix sizes do not
  match"). ⇒ **SR stays the GPW conditioning answer**; dropping it → TODO §1. Agreed auto-tol/auto-Rcut design
  in the (resolved) OPEN INVESTIGATION section.

## RUNTIME round 2 — sampling patch/multi-grid (a DEAD END that redirected the whole approach) — 2026-07-13
Chased the `O(n²·Npts)` scaling with SAMPLING-based collocation (keep `PhiOnGrid` = the Bloch orbital sampled
on the grid; restrict the dense contraction to per-orbital/per-pair PATCHES; then per-exponent grid LEVELS).
It works but hit a ceiling, and the failures are the reason we pivoted to CP2K's ANALYTIC method (next section).
- **Measure-first (`GPW.DISABLED_PatchSparsityProbe`)**: on the 2-atom NaF per-point scatter is a dead end
  (0.37× MACs but every grid point occupied, scalar → loses the zgemm's SIMD); per-pair patches are ~break-even
  (0.25× MACs at exact ε, scalar ≈ the zgemm). Patch-PRUNING is asymptotic in atom count — the 2-atom cell has
  nothing to prune, so the lever is the GRID axis (multi-grid), not the pair axis.
- **Increment 1 (`c94269c8`)** — molecular-side patched integrate-back (`LatticeSum1E::MakePotentialMatrix`,
  per-orbital compressed χ columns, contract each pair on the support overlap). Bit-consistent vs dense.
- **Increment 2 (`d6079a68`,`8b69e4c8`,`cdb695c0`,`38b63d7b`)** — REL_CUTOFF-style multi-grid. KEY FINDING: an
  integrate-back-only multigrid is a DEAD END (Si Γ −21.4 vs −8.25) — coarsening a diffuse pair against the
  SHARP local PP is catastrophic. Split it: keep the STATIC local PP dense, route only the smooth DYNAMIC
  Hartree+XC through the ladder → Si Γ −8.24851 (grid tol), **3.4× on the hotspot** (`maxLevels=2` cap).
- **WHY IT'S A DEAD END**: (a) the sampling multigrid is only grid-TOLERANCE accurate (~10 mHa cost, argues
  against making it the default vs the CP2K-validation reference); (b) it BREAKS at genuine bulk — NaF Rcut=2a
  MG vs dense = **2.66 Ha** (the diffuse-pair coarse-grid sampling ALIASES with the image sum); (c) it still
  needs a hand-tuned `Rcut`/`collRcut` sphere on the collocation → **Gibbs ringing** (a hard truncation of the
  sampled Bloch sum, the root of the indefinite-overlap symptom). Sampling-then-quadrature cannot coarsen and
  cannot escape the hard cutoff. ⇒ the whole sampling collocation is the wrong foundation.

## GPW REWRITE to the CP2K method — Increments A, B + cross-cell (DONE, kernels validated) — 2026-07-14
The analytic collocate/integrate core, molecular-side (reusing `GaussianRF`/`Ω`; primitives stay encapsulated).
- **A (`0d09a6d5`)** `Molecule::LatticeSum1E::CollocateDensity(D,cell,N)` — `ρ=Σ_ij D_ij χ_iχ_j` collocated
  analytically per pair on compact exp-tail boxes, modulo-wrapped. New `UnitCell::ToFractional`. Gate
  `GPW.AnalyticCollocationConservesCharge`: `∫ρ=Tr(D S)` to 2e-8; corner atom wraps IDENTICALLY (no ringing).
- **B (`068b4e96`)** `IntegratePotential(V,cell,N)` — the exact adjoint (`ForPairBox` shared by scatter/gather).
  `⟨collocate(D),V⟩=⟨D,integrate(V)⟩` to machine precision (variational).
- **Cross-cell fix + `G_ERI3.apply` (`729b6355`)** — the periodic Γ density is a product of BLOCH orbitals
  `χ_i^G χ_j^G = Σ_R'' χ_i^0 χ_j^R''`, so collocation loops the SCREENED cross-cell offsets (`ForImageOffsets`,
  magnitude-screened on the product prefactor — NOT a hard Rcut), not just the home pair. (A single-grid SCF
  gave −16.77 vs −8.25 precisely because these were missing.) Gate
  `GPW.DISABLED_AnalyticCollocationCrystalChargeConservation`: Si crystal `∫ρ=Tr(D S^G)` to 2.4e-7.
  **Density seam per the design review (Band_FT_IBS + IrrepCD UNCHANGED)**: `G_ERI3` gains an optional
  matrix-free `apply: D→ρ̃` closure (a 3rd realization beside PW-delta and GPW-dense-weights); `ContractG_ERI3`
  dispatches to it. `GPW_Evaluator::MakeCollocator` builds it (collocate→FFT, Coulomb kernel folded). Built but
  NOT wired: single-grid analytic is impractically slow (SIPP's diffuse α=0.25 → full-grid boxes × the image
  sum; Si Γ ~60 min), so C (wire it) and D (multigrid) are COUPLED and land together — see TODO §0.

## (archived from TODO §0) the original C+D design record
**[DONE 2026-07-14 — see the "GPW ANALYTIC REWRITE COMPLETE" DONE entry.  What remains from this section is
only E (auto-Rcut + k in P(R) — the cellphase callback already carries any k; genuinely-complex k awaits the
shifted-MP gate re-validation) and the NaF re-timing.  Kept below as the design record.]**
The sampling collocation is a proven dead end (rings, can't coarsen, needs a hard Rcut — see DONE). The analytic
collocate/integrate KERNELS are DONE + validated to machine precision (Increments A/B/cross-cell, DONE). What
remains is to make them the SCF path — and they are only PRACTICAL with the multigrid (single-grid analytic is
~60 min for Si Γ because SIPP's diffuse products have full-grid boxes), so **C and D land together.**

**C — wire the analytic path in (density + integrate-back), delete the sampling machinery.** The seam is ready:
- **Density**: `Repulsion3CTensor()`/`Overlap3CTensor()` → `g.apply = MakeCollocator(coulomb)` (already built,
  reverted). `ContractG_ERI3` dispatches to `apply(D)` = collocate→FFT. **Band_FT_IBS + IrrepCD UNCHANGED.**
- **Integrate-back**: `OverlapMatrix(Vtilde)` → `IntegratePotential(V,cell,N)` (widen real→Hermitian at Γ);
  `MakeLocalPP` uses the same (the analytic integrate-back is accurate for the sharp PP — no coarsening error).
- **DELETE**: `PhiOnGrid`, `BuildWeights`, `DenseOverlapMatrix`, `PatchedOverlapMatrix`, `MultiGridOverlapMatrix`,
  the sampling MG level machinery, `Rcut`/`collRcut` for the collocation, `G_ERI3::weights`, and the superseded
  sampling tests (`PatchedIntegrateBackMatchesDense`, `MultiGridIntegrateBackMechanics`, `SiliconGammaMultiGrid`,
  `NaFMultiGrid*`, `BenchOverlapMatrix`, `PatchSparsityProbe`). Re-pin the anchors to the analytic values (they
  MOVE ~grid tolerance from the sampling values — Si Γ was −8.24758 sampling).

**D — REL_CUTOFF multigrid (what makes it fast AND keeps it accurate).** Assign each pair to the coarsest level
with `cutoff ≥ (α_i+α_j)·rel_cutoff`; collocate/integrate each pair on ITS level; transfer V to all levels via
FFT (spectral, no ringing); combine ρ̃ nested in G-space. Because the collocation is ANALYTIC (not sampled), a
diffuse product on a coarse grid matched to its width is accurate — the sampling multigrid's fatal flaw is gone.
This is CP2K's per-exponent multigrid; each pair's box becomes a ~constant small point count → the ~10–100× that
closes the gap (NaF ~1 min).

**E — auto-Rcut + k in P(R).** Replace the last `Rcut` param with EPS_PGF_ORB screening (the neighbour-list reach
is `√(−ln ε/α_min)`, already the collocation offset screen); general k via `P(R)=Σ_k w_k e^{ikR}` (collocation
is k-agnostic — the grid density is real/cell-periodic). Then a single ε drives everything, CP2K-like.

**Validation ladder**: re-pin Si Γ (analytic == sampling to grid tol) → the crystal charge/adjoint gates
un-DISABLED (fast once multigrid) → NaF vs CP2K same-basis (the real target). Kernels + `MakeCollocator` +
`G_ERI3::apply` + the cross-cell `ForImageOffsets` are all in place; C+D is the wiring + the level machinery.

### (ARCHIVED, pre-2026-07-13 — superseded by §0/§1 + DONE) original profiling / experiment notes
**Discipline (user-directed): measure before optimizing.** Start by PROFILING, then pick the fix.

**Profile, don't guess (ready-to-run entry points):**
1. **Profile our GPW NaF.**
   - **perf (PREFERRED — low overhead, profile the REAL 34-min run).** Blocked today by `kernel.perf_event_paranoid=4`
     (unusually high; ≥3 blocks unprivileged perf). Fix with sudo (persists via `/etc/sysctl.d/99-perf.conf`):
     `sudo sysctl kernel.perf_event_paranoid=1`. The build already keeps frame pointers + `-g`
     (`-O2 -g -fno-omit-frame-pointer`), so `--call-graph=fp` is enough (no heavy dwarf):
     ```
     perf record --call-graph=fp ./UnitTests/UTMain --gtest_also_run_disabled_tests \
       --gtest_filter='GPW_SCF.DISABLED_NaFRocksaltGamma'
     perf report -g 'graph,0.5,caller'
     ```
   - **callgrind (NO ROOT — the fallback that works today; EXACT counts, view in KCachegrind).** ~20–50× slower,
     so use a SMALLER-but-same-code case; the function-level hotspot RANKING transfers. Si GPW (fast, but Ecut~12 /
     Rcut=0 — under-weights NaF's fine-grid+image cost) for structure, or NaF at reduced `nmax` for the real terms:
     ```
     valgrind --tool=callgrind --callgrind-out-file=cg.out ./UnitTests/UTMain \
       --gtest_filter='GPW_SCF.SiliconGammaConverges'     # then: kcachegrind cg.out
     ```
   - gprof/`-pg` REJECTED (user: unreliable under modern inlining/opt).
2. **Profile/time CP2K** on the same NaF (`~/Code/cp2k/build/bin/cp2k.ssmp`, deck `UnitTests/CP2K/naf_gpw.inp`)
   for a wall-clock TARGET + its term breakdown (points at where the time *should* go).
3. (Optional, cheap insurance, SEPARATE from runtime) one `valgrind --tool=memcheck` pass to definitively close
   "is there a memory bug?" — expected clean (Si-validated, charge-conserved). **ASan/memcheck ≠ profiler:** they
   find memory bugs; `perf`/`callgrind` find runtime.

**Candidate hotspots + their DIFFERENT fixes (the profile picks among these — they are NOT the same optimization):**
- **[MY BET] The fine grid `densityEcut=160`** (auto-floored for F α_max=40): Npts ∝ Ecut^1.5, so ~11× the
  Si-scale grid. EVERYTHING per iteration scales with Npts — density collocation `ρ=ΣD_ij χ_iχ_j`, FFTs,
  integrate-back. **Fix = MULTI-GRIDS** (map F's tight primitive to its own fine grid, diffuse functions to coarse
  ones — CP2K's per-exponent multigrid; §4). Likely the dominant lever, and NOT magnitude-screening.
- **Image collocation (Rcut=2a)** — the `|R|≤Rcut` sphere re-summed at every grid point. **Fix = magnitude-
  screening** (per-pair `|⟨χ_i|χ_j^R⟩|>ε`, CP2K's `EPS_PGF_ORB`) — sparse + drops the arbitrary Rcut + fixes the
  indefinite-S (§OPEN INVESTIGATION). Real, but probably a one-time `PhiOnGrid`/overlap cost, not the per-iter one.
- **`PhiOnGrid` one-time build** — O(Npts × nOrb × nImages); already cached across iterations, but the build at
  Ecut=160 × Rcut=2a images could be large. Fix = magnitude-screen the images (above) + the multigrid (above).
- **The W-tensor / integrate-back** — O(nGfit × nAO²) storage+FFTs. Fix = whole-density collocation (§4).
**Deliverable:** a profile-backed ranking, then implement the top 1–2 fixes, re-time NaF vs CP2K. THEN return to
the (now-cheap) mixing tuning (G0/α/DIIS/nmax) to get the converged NaF number for `doc/CP2Kresults.md`.
- **SPEED: the GPW NaF run is VERY slow (≫ CP2K) — magnitude-screening is now a SPEED item, not just correctness.**
  The `|R|≤Rcut` sphere drags EVERY function (incl. tight ones that overlap nothing at 2a) out to Rcut=2a, and the
  collocation re-sums that whole image set at every grid point every SCF iteration. Per-pair `|⟨χ_i|χ_j^R⟩|>ε`
  screening (CP2K's `EPS_PGF_ORB`) makes cost scale with REAL overlaps (sparse) AND drops the SR/Rcut crutch AND
  fixes the indefinite-S correctness issue (§OPEN INVESTIGATION). **Back on the active list** (was deferred).
- **NEXT (reprioritised, in order):** (a) [RUNNING] confirm DIIS-delay tames the slosh; (b) Kerker/high-G-damped
  mixing (needed regardless); (c) magnitude-screen the overlap (SPEED + correctness + drops Rcut); (d) push the
  ladder to Ecut≈160/200 (or multi-grids) for a converged NaF Etot vs CP2K/PW; (e) multi-grids
  (§4) to make the high cutoff affordable (F's tight primitive → own fine grid, diffuse Na → coarse). Considered:
  a CP2K-style grid-charge RESCALE `ρ_grid *= N/∫ρ_grid` (cheap monopole guard + the same number CP2K prints) —
  but it is a COUNT fix, not a SHAPE fix (aliasing corrupts ρ̃(G≠0), which Hartree/XC see), so it does not replace
  cutoff convergence. **DEFERRED (user, 2026-07-12): do NOT add it yet — it would MASK the very `∫ρ_grid` swing
  (5.2 ↔ 7.99 at Ecut=120) we are using to TRACK the residual charge-transfer oscillation; pinning ∫ρ_grid=8 blinds
  the instrument before the problem is solved.** Add only AFTER grid is converged AND the oscillation is fixed, and
  even then CP2K-style with the rescale MAGNITUDE always printed (= the loss = the diagnostic; never silent). Note
  the periodic G=0-dropped Hartree is charge-BLIND (unlike molecular 1/r RI, which is why molecular DFT constrains
  ∫ρ̃=N and PW does not need a constrained fit).
- **IONIC SEED — generator DONE (2026-07-12), wiring NEXT.** Root cause of the useless IonicSAD: it scaled the
  NEUTRAL F valence density's AMPLITUDE ×8/7, keeping the COMPACT neutral shape (69 vs 58 iters vs Uniform,
  `PlaneWaveDFTUT.C:1473`) — a real F⁻ is spatially DIFFUSE. **Design (user): an OFFLINE-generated LIBRARY of seed
  densities (neutral + chemically-plausible ions), NEVER an atom SCF at lattice-run time** (a production run is a
  lookup — robust for a newbie; the SCF "surprises" are confined to the offline generator). Unified with the
  same atom-SCF machinery that makes the valence BASES: `qchem::ValenceBasisGen::GenerateSeedDensity(recipe)`
  runs the charge-state pseudo-atom SCF (`recipe.electrons`: neutral F 7, F⁻ 8, Na⁺ 0) and samples ρ(r) on a log
  mesh into a library entry (schema of `atomic_valence_densities.json`). **Validated (`ValenceBasisGen.
  FluorineSeedDensityAnionIsDiffuse`): F⁻ ⟨r⟩=1.62 vs neutral F 1.18 (37% more diffuse), charge 8.00 vs 7.00** —
  the anion diffuseness the seed needs, for free. **INCREMENT 2 DONE (2026-07-12) — and it WORKS:** (a) F⁻ entry
  (Nelec=8) captured into `atomic_valence_densities.json` (neutral F/Na/Si preserved); (b)
  `GetAtomicDensity(Z,functional,dbfile,Nval)` + `HasAtomicDensity` select a charge state by valence count
  (Nval<0 = neutral, backward-compatible); (c) `SeedCD`/`IonicSAD` now pull the library's CHARGE-STATE density
  (target `Nval-q`: F→F⁻ 8 e⁻ diffuse, Na→0 e⁻) with scale 1 — no more amplitude hack; falls back to
  neutral×scale only if the library lacks the ion. **Validated in PW (`PlaneWaveDFT.FrameworkNaFThroughSCFIterator`):
  IonicSAD HALVES the iterations — 17 vs Uniform 35** (was 69 vs 58 = WORSE with the compact seed), charge 8,
  same converged −20.3293 (seed-independence). Guarded `EXPECT_LT(I.iters,U.iters)`. GPW NaF test switched to
  IonicSAD (same machinery → same win; the slow auto-160 SCF not re-run this session). Full `-A_*` suite green
  (shared molecular SAD path unaffected — neutral lookups still first-match). Ionic-seed axis: DONE.

**Reprioritised diagnosis (what's actually going on):**
- **Overlap conditioning ≈ RED HERRING.** min eig(S)=7.5e-4 (SR/Rcut=2a) orthogonalises trivially: min eig =
  min sv for Hermitian PSD; the √ shows up only on the orthogonaliser `V=S^-1/2` (`cond(V)=√cond(S)`, amplifies
  ≤ 1/√min_eig ≈ 36×). You'd need min eig ~1e-16 to matter. Confirm with the residual `‖V·S·Vᵀ − I‖`. A
  slightly-more-SR basis is cheap insurance, not the fix. (Earlier "near-singular metric → instability" and
  "unoccupied → redundant → instability" were BOTH wrong — user corrections.)
- **Our density "fit" basis IS PLANE WAVES** (`GPW_IBS.C:41` `PlaneWaveFit_IBS`), the SAME family as CP2K
  (whose "no fit" is really a hidden PW/grid fit). So we are NOT in a different regime from CP2K — the Hartree
  (W-tensor `FT[χ_iχ_j]` × `4π/G²`) is exact given the grid; XC is grid quadrature. **Any residual
  non-variationality is PROCEDURAL** (the `FittedVee`/`FittedVxc` projection/consistency), not the basis. And
  the DIIS≠GDM "proof" of non-variationality was WEAK (both runs hit the 60-iter cap → not two minima, two
  unfinished trajectories). **⇒ whole-density collocation (ditch the `Fitted*` wrappers, collocate ρ directly
  on the grid) should reproduce CP2K under the hood.** This is the deep fix.
- **SCF instability = ionic charge-transfer oscillation, NOT conditioning.** NaF wants Na⁺F⁻ but the test
  seeds `SeedStrategy::Uniform` (line 365) — the SCF must move a whole electron Na→F from a flat start →
  oscillation that mixing-factor throttling can't damp (user saw: reduced DIIS EMax, throttled relax, no
  success). The relax auto-tune KEYS OFF [F,D] (`SCFIterator.C:207-216`), an unreliable signal here, so it
  misfires. `IonicSAD` is NOT a ready fix: documented "Phase 3, not implemented" and the dcmplx/GPW path falls
  back to Uniform (`Seed.C:29-31`); even in PW it's WORSE for NaF (crude too-compact ionic ρ → high-G noise,
  `PlaneWaveDFTUT.C:1473`). Real fixes: a properly-DIFFUSE ionic seed (real F⁻ is diffuse), Kerker/preconditioned
  mixing (damps charge sloshing — linear mixing amplifies it), electronic smearing, or the variational-energy +
  direct-min path once collocation lands.

**Get CP2K converging (the oracle — `UnitTests/CP2K/naf_gpw.inp`, currently diverges to +400 Ha under OT):**
1. **Isolate the variable:** single Na q1 atom-in-box, then F, then NaF. Atoms converge but NaF doesn't ⇒ the
   ionic charge transfer, not the basis/PP.
2. **Kill the overshoot:** `MINIMIZER CG` (OT-CG doesn't extrapolate → no +400).
3. **Robust preconditioner:** `PRECONDITIONER FULL_ALL`.
4. **Fix the guess** (the ATOMIC guess is far from Na⁺F⁻; CP2K printed `electrons 11→9→rescale 8`): traditional
   diagonalisation + Broyden mixing + a little electronic smearing at LOW CUTOFF → `SCF_GUESS RESTART` into OT.
5. **Converge CUTOFF upward** (100→200→400 Ry) to separate grid effects from SCF stability.

**Fit-quality metrics (for ρ and Vxc) — and what CP2K reports:**
- **ρ (Hartree side):** the rigorous metric is the **Coulomb-metric residual** `‖ρ−ρ̃‖_C = √(∬ Δρ(r)Δρ(r′)/|r−r′|)`
  (the RI-V norm; Hartree-energy error is 2nd-order in it → near-variational). On a grid = the Fourier tail
  beyond G_max. **Practical scalar: `∫ρ_grid − N`** (grid charge conservation). **Never ΔE_total** (non-var).
- **Vxc:** nonlinear → its quality is grid resolution where ρ is sharp (tight F); watch **Exc vs CUTOFF and vs
  REL_CUTOFF** (the denser-Vxc-grid knob for ∇ρ).
- **CP2K reports it directly:** `Electronic density on regular grids: -7.9963  0.0037` — integrated grid ρ and
  its **error (0.0037 e⁻ lost to truncation)**; the `Re-scaling ... Number of electrons: 8` step corrects it and
  the rescale magnitude IS that error. Plus the multigrid (4 levels + `REL_CUTOFF 30`) = per-exponent grid
  mapping. **DONE 2026-07-12: `∫ρ_grid − N` readout added** — `qchem::Hamiltonian::ReportGridCharge()` (opt-in),
  printed by `PW_XC::RefreshRhoGrid` per iteration (see the DIAGNOSTICS block at the top of §0). Turns "is our
  grid good enough" into CP2K's controlled number.

**Iteration-output refactor (user-requested — diagnostic infrastructure):** the per-iteration columns in
`SCFIterator.C:148` are hardcoded (`Etotal  ε+V/K  Δ[F,D]  Δρ  …`). Atoms / Molecules / Solids (and HF vs DFT)
want DIFFERENT ideal columns. **Refactor the header + `DisplayEnergies` through VIRTUAL DISPATCH on
`tSCFIterator<T>`** so derived (per-system) classes choose the columns and their order. For solid/GPW-DFT the
useful columns are **`‖ρ−ρ̃‖_C`, `∫ρ_grid − N`, `ΔE`** — and DROP **Δ[F,D]** (non-variational, useless here) and
the **virial `2+V/T`** (meaningless under a PP / periodic). While there, fix the relax auto-tune keying off [F,D].

**Ordered experiment plan for next session:** (a) a properly-diffuse ionic seed OR Kerker mixing OR smearing to
kill the NaF charge-transfer oscillation (biggest immediate win); (b) get CP2K converging (isolate → CG →
warm-start) for the real reference; (c) **whole-density collocation** (the deep fix: match CP2K, remove the
procedural fit noise, make the energy variational so GDM/OT can win); (d) the iteration-output virtual-dispatch
refactor + the `∫ρ_grid−N` readout (diagnostics); (e) slightly-more-SR basis as conditioning insurance;
(f) magnitude-screen the overlap (correctness+speed, drops the arbitrary Rcut).

## (archived from TODO §2) valence-basis generator + NaF cross-validation, full records

**PROGRESS (2026-07-11): a valence-basis GENERATOR, not hand-rolled files.** `qchem.ValenceBasisGen`
(`src/Calculation/ValenceBasisGen.C`) generates a low-q valence Gaussian basis straight from an **atomic
pseudo-atom SCF**: `GenerateValenceBasis(recipe)` runs the spherical solver (correct l-occupation, no molecular
open-shell degeneracy) in a candidate even-tempered window to VALIDATE it, then emits the per-l shells as a
Gaussian94 element block; `AssembleBasisFile` combines blocks into one file. Enabled by `AtomCalcOptions.exponents`
(the "bring your own exponents" atom path). Output so far: **`BasisSetData/valence_lowq.bsd`** (organised by TYPE,
all elements in one file, per the BasisSetData convention) with **F** (F⁻ window, 8s+6p, E=−21.10) and **Na**
(neutral 3s¹, 5s+2p, E=−0.144). Wired as `BasisSetData::VALENCE_LOWQ` / `"valence_lowq"`. Tests: `UnitTests/
ValenceBasisGen_UT.C` (energies + round-trip load). KEY LESSONS: (a) canned bases are F⁻-optimised → don't copy;
the atom calc is the generator/validator. (b) Validate against the physically-relevant CHARGE STATE (F⁻ for NaF).
(c) Oracle GS-energy matching is the WRONG objective (user) — N≈8 windows, move on; refine later from a NaF-GPW
**orbital-coefficient heat-map**. (d) Keep per-l exponents DISJOINT: the molecular Gaussian94 reader has a
flagged inverted-condition bug (`PG_Cart/Imp/IrrepBasisSet.C`) that drops a shared-exponent p shell; fixing it
shifts every density-fit DFT anchor 10–70 mHa → its own re-pin task. NEXT: Cs/I blocks; then multi-species GPW
NaF/CsI (thread the species→q map through `RunGPW`/`GPWFactory`; `Ham_PW_DFT` multi-species ctor already exists).


**NaF cross-validation findings (2026-07-11):**
- **GDM vs DIIS = non-variational confirmed.** Our GPW-SR/Rcut=2a gives DIFFERENT iter-capped totals under
  DIIS (−23.556, Ekin 12.1) vs GDM (−23.936, Ekin 29.3). A variational energy would give the SAME minimum
  under both minimisers; different answers ⇒ the fitted GPW Etot is non-variational (fit noise), so the
  limiter is the ENERGY FUNCTIONAL, not the solver. GDM (our OT analog, now dcmplx via `89f210f0`) does NOT
  rescue it — matches the plan's prior note.
- **CP2K (`UnitTests/CP2K/naf_gpw.inp` + `VALENCE-LOWQ-BASIS`, GTH-PADE-q1/q7, LDA_X+LDA_C_VWN, Γ, CUTOFF 400):**
  the FULL diffuse basis DIVERGES the SCF under both P_Mix/Diag AND OT (energies → +200..+400 Ha); the SR basis
  also diverges under OT, but **transiently passes −23.64** — right next to our GPW-SR −23.556. So both codes
  agree the answer FOR THIS GAUSSIAN BASIS is ≈ −23.6, and the ~3.3 Ha gap to PW's complete-basis −20.3293 is
  **Gaussian-basis incompleteness** (the "GPW vs PW = basis quality" leg). Neither converges cleanly because
  the cause is a **near-singular overlap METRIC, not occupation** (an earlier note wrongly said "unoccupied Na
  functions → redundant → instability"; unoccupied functions just get small well-defined coefficients — user
  correction). Our sweep measured it: SR/Rcut=2a has **min eig(S)=7.5e-4, cond≈8000** — barely PSD. Every SCF
  step (OT geodesic, DIIS Fock inversion) goes through S^-1/S^-1/2, so a near-singular S makes the steps
  ill-conditioned: CP2K's OT gradient stays ~23 and the energy overshoots to **+400 Ha** (the minimiser
  overshooting through a broken metric, NOT variational collapse); our sharp-Rcut GPW instead makes the
  truncated S *indefinite*. Same tiny-min-eig root, two symptoms. So **magnitude-screening (fixes the
  TRUNCATION) is necessary but NOT sufficient**: if the complete-Bloch S is itself near-singular from
  over-diffuse functions, the minimiser is still ill-conditioned. Deeper fix = a **better-conditioned (less
  over-complete) basis** for ionic NaF; plus, for our GPW, the separate fit-noise floor.
- NEXT (user-directed): (1) **magnitude-screen the overlap** `(i,j,R)` by `|⟨χ_i|χ_j^R⟩|>eps` (CP2K's trick —
  PSD + fast, drops the SR/Rcut crutch); (2) reduce the fit noise that makes Etot non-variational; (3) an
  ionic-appropriate Na basis for a clean CP2K reference.

# OPEN INVESTIGATION — LARGELY RESOLVED (2026-07-13)

**Magnitude screening IMPLEMENTED + COMMITTED (`05e44fab`).** `NR_Evaluator::LatticeSum` now screens each
`(i,j,R)` term by a per-component reach `r_i=√(−ln ε/α_min,i)` (ε=1e-10), shared identically across S/T/V_nuc
(consistency is a CORRECTNESS requirement — S and H must sit on the same support for `HΨ=εSΨ`). Effect: the 1E
lattice sums are SPARSE (**~4×**: 0.37 s vs 1.46 s to Rcut=4a), so a **generous Rcut is now free** — the "pinned
at 2a for tuning" pain is gone. But screening only *removes* sub-ε terms; the caller must still ENUMERATE far
enough (screening cannot add a term never enumerated).

**KEY FINDING — the full-basis indefiniteness has TWO causes, and screening only fixes one.** Extending the
sweep to 3a/4a (now cheap): full-basis min eig converges to 0 **from below** (−0.42→−0.11→−4e-4→−4.8e-8). The
large-negative *truncation/Gibbs* part IS cured by enumerating far (screening makes it affordable), but the
residual ~0⁻ is **intrinsic OVER-COMPLETENESS** of the diffuse Bloch-summed basis — a BASIS problem, not a
cutoff one. SR is cleanly PSD (+7.5e-4→+9.6e-7, from above). So the plan's old "screening → PSD full basis" was
HALF right (kills Gibbs, exposes over-completeness).

**(1) tune basis (SR) vs (2) tune ortho (truncate eigen/SVD) — RESOLVED for GPW: (1)/SR stays, (2) is BLOCKED
at the SCF stack.** The full basis's null directions cluster at ~1e-6 in a **clean ~1000× spectral gap** below
the physical ~1e-3 spectrum, so canonical Eigen/SVD ortho with tol in the gap gives a clean transform
(‖VᴴSV−I‖=6.6e-11 vs SR+Cholesky 4e-14 — bounded but ~1000× noisier, vindicating the user's atomic-HF
truncation-noise caution). BUT the SCF validation (`DISABLED_NaFFullBasisEigenTol`) hit an **integration wall**:
truncation reduces the working dim 37→33, and the periodic stack (`Crystal_EC`/`cDM_CD`/collocation) assumes
the full `n` → `"Matrix sizes do not match"` before iter 1 (the molecular path handles rectangular V; the
periodic path does not). So **dropping SR needs rank-reduction plumbed through the periodic stack** — a future
increment. Until then SR (dimension-preserving, cleanly PD) is the GPW conditioning answer.

**AGREED DESIGN (for when the rank-reduction stack work is done):**
- **Auto-Rcut via `MaxReach(ε)`** (basis exposes one scalar, mirroring `MaxExponent`; the lattice enumerates
  `CellsInSphere(MaxReach+cell-span)` — wall (B): exponents stay behind the molecular-basis wall, k-convention
  stays lattice-side). Removes the Rcut parameter; ε (a tolerance) replaces it, exactly like CP2K's
  `EPS_PGF_ORB` (CP2K sets NO user Rcut).
- **Auto-tol via GAP DETECTION** in `LASolver` (separation of concerns — pure LA): sort eig ascending,
  force-drop `d[i]≤0`, scan the LOW region (`d[i] < √ε·d_max`) for the largest consecutive ratio `ρ=d[i+1]/d[i]`;
  if `ρ > R_threshold` (**default 30**, exposed at the Calculation facade — visible but rarely touched) it's a
  CLEAN gap → cut there; else fall back to the ε-tol and WARN (ambiguous, noise-prone — the continuum case).
  `orthoTol<0`=auto, `=0`=none, `>0`=explicit (mirrors `densityEcut`). **Auto-cut is allowed but NEVER silent** —
  always `cerr` WARN with count + gap ratio + clean/ambiguous, so the user knows what the basis was truncated by.
- **Vision:** collapse knobs to ~one physically-meaningful ε (drives auto-Rcut, and could drive grid + ortho
  tol), CP2K-like. `densityEcut` already auto; `collRcut` is the later patch/collocation axis.

---

## (superseded) original 2026-07-11 diagnosis — why is the truncated Bloch overlap S indefinite?
User's intuition (from the earlier Si session): S(k) should be PSD for **any** Rcut, and in Si an
indefinite-overlap symptom was traced to a BUG — a separable-KB projector on a **corner atom** (τ=0) whose
image/tail "outside the unit cell" was dropped; after fixing it, S was PSD at any Rcut. Asked to look for the
same bug in the NaF path. **Findings so far (uncommitted, my analysis — cross-check against that old session):**

- **New diagnostic makes this cheap:** `qchem::ReportOverlapConditioning()` (LASolver, opt-in) prints min
  eig / min sv / cond of S at `SetBasisOverlap`; `GPW_SCF.DISABLED_NaFOverlapConditioningSweep` builds ONLY
  the analytic Bloch overlap (no SCF) across Rcut in ~0.2 s. NaF full basis: min eig **−0.42** at Rcut=a,
  −0.60 at 1.5a, −0.11 at 2a; SR basis: −0.035 / −0.046 / **+7.5e-4 (PSD)** at 2a.
- **Image enumeration is CLEAN — no obvious corner-atom drop bug in the OVERLAP.** `BuildImages` uses
  `UnitCell::CellsInSphere(Rcut)` = a symmetric (`n`&`−n`), COMPLETE origin-centred sphere on `|R|≤Rcut`, with
  NO cell-membership filtering. S is Hermitian (real eigenvalues at Γ). So the overlap does not drop
  images-outside-the-cell the way the Si KB projector did.
- **The KB corner-atom bug WAS real but is a DIFFERENT term, already fixed (`95e8f4a8`):** `MakeSeparablePP`
  used the raw home orbital as the projector bra, losing the corner atom's wrapped tail (16 Ha
  translation-variance). Fixed by using the Bloch-summed orbital. That fix does NOT touch the overlap's PSD-ness.
- **The real reason S is indefinite = the analytic SINGLE lattice sum is a Dirichlet-windowed autocorrelation.**
  GPW builds `S_ij(k)=Σ_{|R|≤Rcut} e^{ik·R}⟨χ_i⁰|χ_j^R⟩` (bra home, ket imaged). The FULL sum (Rcut→∞) is the
  Gram matrix of Bloch orbitals ⇒ PSD; a SHARP `|R|≤Rcut` cutoff is the rectangular-window (Dirichlet) partial
  sum of that autocorrelation ⇒ **can go negative** (Gibbs), and does so once the dropped tail exceeds the
  basis' smallest eigenvalue — hence worse for the diffuse (ill-conditioned) full basis, cured by SR + Rcut=2a.
  This matches the code's own note ("a truncated single sum can be indefinite; a generous Rcut is the fix") and
  the Si record (PSD only at Rcut≥3a). So for the single-sum scheme, "PSD at any Rcut" does NOT hold in general.
- **Corner-atom RESONANCE that's worth a second look:** the image sphere is centred on the LATTICE ORIGIN and
  the SAME set is used for every atom pair, but the physical decay of `⟨χ_i⁰|χ_j^R⟩` is centred on the pair
  SEPARATION `τ_j−τ_i+R`. For the DIAGONAL blocks (τ_i=τ_j) the cutoff is atom-centred (symmetric); for
  OFF-DIAGONAL blocks of an offset atom (F at ¼¼¼ vs Na at the corner 0) the origin-centred `|R|` cutoff
  truncates the pair tail asymmetrically → plausibly worsens the indefiniteness for multi-atom cells. A
  **pair-separation-centred** cutoff (include images where the pair overlap is actually significant, per pair)
  would be the more symmetric truncation and is the closest thing to a "corner atom handled specially" fix.
- **The rigorous "PSD for ANY Rcut" route = the Fejér/Gram scheme (plan's "scheme B", done consistently).**
  Build S as the Gram of the TRUNCATED Bloch orbitals `⟨φ_i^k|φ_j^k⟩`, `φ_i^k=Σ_{R∈Rs}e^{ik·R}χ_i^R` — a
  double lattice sum whose image terms carry Fejér (triangular) weights `c(ΔR)=|Rs∩(Rs+ΔR)|` ⇒ PSD by
  construction, any Rcut. The plan rejected this ONLY because a scheme-B overlap was mixed with a scheme-A
  single-sum kinetic (Ekin=−300); doing ALL 1E matrices (S, ⟨p²⟩, V) in the SAME tapered Gram scheme is
  self-consistent and PSD, at the cost of a tapered (approaches-exact-as-Rcut→∞) metric and O(images²) work.
- **RESOLUTION (user insight): CP2K is fast AND PSD with "no truncation" because it screens by MAGNITUDE, not
  geometry.** CP2K's neighbour lists (`EPS_PGF_ORB`/`EPS_DEFAULT`) include an image pair `(i,j,R)` only if the
  Gaussian product `⟨χ_i⁰|χ_j^R⟩` is non-negligible — a PER-PAIR, PER-FUNCTION adaptive reach: a diffuse
  Gaussian reaches far (until its tail < eps), a tight one reaches ~nothing. This is (a) FAST (sparse — cost
  scales with real overlaps, not `Rcut³`), and (b) PSD at any Rcut (drops only sub-threshold terms, so the
  error stays below `λ_min(S)` → S ≈ the exact complete-Bloch PSD overlap; a *significant* tail is never
  dropped). **Our `|R|≤Rcut` sphere is wrong on BOTH axes:** it drags tight functions out to 2a for nothing
  (slow) AND chops diffuse tails while still significant (indefinite). SR helped because it's a crude manual
  version of magnitude screening (removes the diffuse tails by hand).
- **THE FIX (do this next): replace the fixed geometric `Rcut` with per-(i,j,R) magnitude screening** — include
  an image term only if `|⟨χ_i⁰|χ_j^R⟩| > eps` (or size each Gaussian's reach from its exponent + eps, the
  CP2K `EPS_PGF_ORB` way). Then diffuse functions get their needed reach (PSD, any effective Rcut) and tight
  functions cost nothing (fast) — CP2K's trick, and it removes the SR crutch. `BuildImages`
  (`GPW/Imp/Evaluator.C`) currently uses `UnitCell::CellsInSphere(Rcut)`; the screen belongs in
  `Molecule::LatticeSum1E` (which knows the actual pair integrals) or as a per-shell reach handed to it.
- **Short term (done, works):** SR + Rcut=2a. The Fejér/Gram scheme is an alternative but magnitude screening
  is what CP2K proves out. Cross-check the corner-atom claim against the old Si session if useful; the 0.2 s
  sweep makes any hypothesis a trivial check.

---


---

# ARCHIVED FROM GPWPlan 2026-07-23 — the full §0-campaign records (analytic rewrite → §0f)

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

## §0e NaF PRODUCTION GRID — RESOLVED (2026-07-19..22).  The full records:

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

### §0e-PP `MakeLocalPP` SETUP WALL — the CP2K local-PP split (analysis 2026-07-19; steps (a)+(b) DONE 2026-07-22)
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
    - **LADDER-SATURATION corollary (measured on the OOM-killed first rerun, 2026-07-22):** the κ bound
      holds only up to the LADDER TOP — pairs with κ·p above the finest level saturate there, carrying
      error ≈ (pair tail at top) × (FIELD tail at the top's ball).  The NaF COARSE seed stage (explicit
      Ecut=40 ladder; F pairs demand κ·p up to 2400 Ha!) shifted −24.099→−23.749 once the analytic
      short stopped cancelling the grid long's saturation error — an artifact of an absurdly cheap
      ladder, harmless in a SEED.  At the FINE 160-Ha top the long's field tail is e^{−G²r_loc²/2} ≈
      e^{−7.6} (F) — and CP2K's own V_long-on-grid saturates IDENTICALLY there (its sharp pairs also
      fall back to grid 1), so production accuracy is unaffected (the rerun verifies).  Moral: κ
      standalone-exactness is a statement about ladders that REACH κ·p OR fields that are r_loc-soft
      at the top — the sharp SHORT piece satisfies neither on cheap ladders, which is exactly why it
      had to go analytic.
    - **NaF SR2 VERIFICATION — PASSED (2026-07-22, solo rerun; OOM lesson: ONE 9-GB run at a time on
      this 14-GB box).**  Grid-matched aufbau fine stage with the analytic short: **−24.4314023**
      (20 iters, clean aufbau frontier) — 0.26 mHa from the all-grid-V_loc −24.4316608 (the fine-ladder
      saturation cost, the predicted e^{−7.6} class) and **0.19 mHa from CP2K's −24.4312134** (was
      0.45 — CLOSER, as expected: the arrangement now mirrors CP2K's own analytic-short +
      saturated-smooth-long).  Fine-stage wall ~73→~48 min (short sweep deleted, one fewer iter; even
      with the stream budget mostly consumed by the 200-iter coarse seed).  **Steps (a)+(b) CLOSED.**
    - **MIXED-RADIX FLIP LANDED SAME DAY (2026-07-22; capability `1837b21e` = PocketFFT submodule +
      `qchem.FFT` dispatch [pow2→radix-2 verbatim, else PocketFFT c2c] + `Next5Smooth`; policy flip
      `aaaf1ea0` = `FFTGrid()` pads AutoGrid to 5-smooth instead of pow2).**  MEASURED: 199/199 with
      **ZERO anchor re-pins** (every energy raster-converged — the ball is the resolution object, N is
      quadrature); suite 598→442 s, Si Γ SCF 92→48 s (32³→30³), atom-in-box 78→9.5 s (8×).
      **NaF grid-matched verification config: 5669 s → 190 s (30×)** — fine raster 128³→72³, coarse
      64³→36³, streams fully cached (0 dropped), **Etot −24.4314027 vs −24.4314023 (0.4 µHa)**, same
      22 iters, charge exact.  vs CP2K's 5.8 s the gap is now ~33× (was ~900×); the remaining raster
      factor is our alias-free DIFFERENCE-SET policy (4m+1 → 72³) vs CP2K's ball-only 36³ (~8× points)
      — revisiting that means tolerating product aliasing on the raster (the CP2K trade), a separate
      deliberate increment if ever; the rest is setup/streams machinery.
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

## §0f LEAD 1 — XC on the RAW collocated ρ → FALSIFIED; the REAL gap decomposed (2026-07-21).  The full record:
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


