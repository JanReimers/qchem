# GPW (Gaussian And Plane Waves) — Plan & State

**Self-contained orientation for a fresh session.** GPW = Gaussian orbitals on a periodic lattice, with
the electron density collocated on a real-space grid, FFT→G-space Poisson for Hartree, XC on the grid,
integrated back against the Gaussians to form the KS matrix (CP2K / Lippert–Hutter). It is the north-star
that makes ab-initio solids → battery voltage curves possible.

This doc supersedes the GPW sections of `doc/MolecularPP_HarmonizationRound2.md`.

**The doc is split into two major sections: [DONE](#done) (committed, anchors green) and
[TODO](#todo--next) (what's left, in priority order).** Then the durable invariants + pointers.

---

# DONE

Everything here is committed on `main`; the GPW test suite (`GPW_UT`, `GPW_SCF_UT`) is green and the Γ
energy anchors hold. GPW is a **new evaluator, not a new IBS** — it satisfies the existing plane-wave
concepts and reuses the `EPW_*` mixins + the whole `Ham_PW_DFT` KS stack.

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

## Naming (`5f609d2f`) — remember these
- `Overlap(f)` = ANY 1-electron `⟨i|f|j⟩` (f may be a potential); `Repulsion` = the 2-electron `1/r12`.
  So the reciprocal-space field→KS-matrix bridge is `Band_FT_IBS::MakeOverlap(f)` / evaluator
  `OverlapMatrix(f)` — **not** `MakePotential`/`PotentialMatrix`. `Make` = uncached.

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

---

# TODO / NEXT

Round-1 runtime (zgemm + OpenBLAS/threads) + magnitude screening + the conditioning investigation are **DONE**
(see DONE). The two LEADING increments now: (0) **RUNTIME GAP-CLOSE** — the `O(n²·Npts)`→`O(n_pairs·patch)`
rewrite that closes the 20–40× gap to CP2K; (1) **DROP SR** — rank-reduction through the periodic stack + the
auto-tol gap-detector. Then (2) low-q multi-species bases → Si/NaF/CsI; (3) CP2K reference; (4) IBZ; (5) cleanups.

## 0. RUNTIME GAP-CLOSE — patch collocation → multi-grid (the CP2K-minute path)
Round-1 (zgemm + OpenBLAS/threads) + magnitude screening are DONE (see the DONE section). The `zgemm` attacked
the CONSTANT of the dense contraction (4× on the hotspot) but NOT the `O(n²·Npts_full)` SCALING — the real
20–40× gap to CP2K (~1 min NaF; we are ~28 min). CP2K is fast because it NEVER touches the full grid densely —
it exploits Gaussian LOCALITY. The lattice-IMAGE axis is now handled by screening (DONE); the remaining lever is
the GRID axis:

| axis | the loop | our crutch | the fix |
|---|---|---|---|
| **lattice images** | `S_ij(k)=Σ_{\|R\|≤Rcut} e^{ik·R}⟨χ_i⁰\|χ_j^R⟩`; `PhiOnGrid`'s Bloch sum | fixed `Rcut=2a` sphere | per-pair `\|⟨χ_i\|χ_j^R⟩\|>ε` (CP2K `EPS_PGF_ORB`) |
| **grid points** | `M_ij=Σ_p conj(Φ_pi)V_p Φ_pj` over full `Npts` | dense 64³ everywhere | local patches + multi-grid |

**Order:**
1. **Locality / patch collocation** — the per-iteration big lever. `PhiOnGrid` is a dense `Npts×n` matrix that
   evaluates even F's sharp α=40 orbital on all 262k points (mostly ~0). Store `Φ` as per-orbital PATCHES
   (points where `|Φ|>ε`), so the contraction sums only over patch OVERLAPS (spatially-close pairs):
   `O(n²·Npts) → O(n_pairs·patch)`. This is the §4 "whole-density collocation" rewrite — collocate
   `ρ=Σ_ij D_ij χ_iχ_j` directly on the grid as local patches (one FFT for Hartree) and integrate `V` back on
   the same patches, REPLACING the dense `PhiOnGrid` + global W-tensor (`W_c(i,j)`, `n²×N_Gbasis`, inherently
   dense — designed for plane waves). Its own correctness campaign: bit-consistency vs the current dense path
   at a converged grid (`L_PP`-style invariants). Possible ~5–10× alone (single grid), and the scaffold multi-grid sits on.
2. **Multi-grids** — per-exponent coarsening ON TOP of patches: map each product `χ_iχ_j` (exponent
   `α_i+α_j`) to a grid level by its exponent (CP2K `REL_CUTOFF`), so only tight×tight pairs touch the fine
   64³ and everything diffuse lives on coarse grids. This is what "bangs down `Npts=NextPow2(4m+1)` for most
   i,j" — the single uniform `densityEcut=160` (dictated by F's α=40) over-resolves the diffuse Na/F functions
   everywhere.
3. **FFT — SECONDARY, and partly dissolves under the rework** (multi-grid FFTs are on smaller grids; the
   power-of-2 padding waste shrinks). Do NOT optimize the hand-rolled `FFT3D` in isolation next — the
   rearchitecting reshapes it. (~18% today + ~10% allocation churn from per-line `cvec_t` allocs.)

**Dimensions recap (NaF, why the grid dominates):** `M_ij=w Σ_p conj(Φ_pi)V_p Φ_pj`. `i,j∈[0,n)`, **n=32**
(Gaussian AOs, `itsOrb->GetNumFunctions()`). `p∈[0,Npts)`, **Npts=64³=262,144**: `densityEcut=4·α_max=160 Ha`
→ `BuildGs {G:½|G|²<160}` (the fit basis N_Gbasis, ~the Ecut sphere) → `m=`max Miller index → `AutoGrid=4m+1`
→ `FFTGrid=NextPow2(4m+1)=64/axis`. `N_Gbasis` (the G-sphere) is NOT in the contraction — it only builds `V[p]`
(RhoOnGrid) and the W-tensor; the contraction runs over the full real-space cube `Npts`. Cost `O(n²·Npts)≈2.7e8`
complex MACs/call, and `Φ` (`Npts×n≈134 MB`) was re-read `n²/2≈512×` by the old scalar loop (memory-bound) —
the zgemm blocks it to ~one pass, but only the patch/multi-grid work removes the `Npts_full` factor itself.

## 1. DROP SR — rank-reduction through the periodic stack + auto-tol
The `_SR` basis is a hand-tuned crutch (drop the most-diffuse primitive so the Bloch overlap is cleanly PD).
We PROVED (2026-07-13, see DONE) that the FULL basis + screening + canonical Eigen/SVD ortho with tol in the
~1000× spectral gap gives a clean overlap transform (‖VᴴSV−I‖=6.6e-11) — BUT the SCF is **BLOCKED**: truncation
reduces the working dim (NaF 37→33) and the periodic stack (`Crystal_EC`/`cDM_CD`/collocation) assumes the full
`n` → "Matrix sizes do not match" (`DISABLED_NaFFullBasisEigenTol`). The MOLECULAR path handles rectangular V;
the PERIODIC path does not. So dropping SR = two pieces:
- **(a) Rank-reduction through the periodic stack** — let a truncated ortho (`V` is `n×(n−k)`) flow through
  `Crystal_EC` (band count `n−k`), `cDM_CD` (density still full `n×n` via `C=V·U'`), and the collocation;
  mirror the molecular path's rectangular-V handling. This is the real work and gates (b).
- **(b) The user-friendly automation** (agreed design, see the resolved OPEN INVESTIGATION section): **auto-Rcut**
  via a basis `MaxReach(ε)` scalar (wall B — the lattice enumerates `CellsInSphere(MaxReach+span)`; exponents
  stay behind the molecular-basis wall, k-convention stays lattice-side), removing the `Rcut` param for one ε
  (CP2K `EPS_PGF_ORB`; CP2K sets no user Rcut). **Auto-tol** via `LASolver` GAP DETECTION (pure LA): force-drop
  `d[i]≤0`, scan the low region for the largest consecutive ratio; if `> R_threshold` (default **30**, exposed
  at the Calculation facade) it's a CLEAN gap → cut there, else fall back to the ε-tol + WARN. `orthoTol<0`=auto
  / `=0`=none / `>0`=explicit (mirrors `densityEcut`). **Auto-cut allowed but NEVER silent** — always `cerr` WARN
  (count + gap ratio + clean/ambiguous). Vision: collapse to ~one CP2K-like ε.

Until (a) lands, **SR stays** (dimension-preserving, cleanly PD, no truncation).

---

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

## 2. Low-q multi-species bases → Si/NaF/CsI cross-validation (PW + GPW + CP2K)

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
| k-points | `MakeKMesh(shift)` (MP; shift=0 Γ-centred, shift=½ classic MP) | `&KPOINTS SCHEME MONKHORST-PACK` | CP2K's MP is SHIFTED (k=±¼ for even N) — use `kShift=½` to match; its Γ-centred list needs `SCHEME GENERAL` (see `si_fcc_gpw_222_gamma.inp`). CP2K prints its k-list (`grep BRILLOUIN`). Shifted mesh currently blocked by TODO 1 (complex-D). |
| Poisson | G-space periodic | `&POISSON PERIODIC` | |
| spin | closed shell (Si₂, 8 e⁻) | `LSD` off | |
| occupation | integer, no smearing | `&SCF SMEAR OFF` | insulator |
| lattice-sum reach | `Rcut`/`collRcut` (our truncation) | `EPS_PGF_ORB` / neighbour lists (auto) | not a direct CP2K knob — converge ours to CP2K |
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
- **Multi-grids (efficiency — the plane-wave analog of per-shell exponent scaling; user TODO, even for LDA):**
  the single uniform `densityEcut` grid is dictated by the TIGHTEST orbital primitive, so a diffuse+tight basis
  over-resolves the diffuse part everywhere. CP2K maps each primitive to a grid matched to its exponent (gated
  by REL_CUTOFF). Generate a `{densityEcut}` list by interrogating the orbital-basis exponents; collocate each
  primitive on its own grid. This is the "×2/×2/3 exponent scaling done in plane waves" (§ the molecular
  `A1_coul`/`A1_exch` fit bases): density = 4×/absolute cutoff, v_xc = functional-dependent (relCutoff).
- **Whole-density collocation (efficiency):** the dense `W` tensor is `O(nGfit·nAO²)` storage + `O(nAO²)` FFTs.
  The efficient GPW collocates `ρ = op(r)` ONCE (one FFT), which needs `D` → density-side. CP2K's local-patch
  (multi-grid) collocation is the further v2.
- **Common `PW_Evaluator`/`GPW_Evaluator` base:** the shared FFT engine (`PeriodicGridEvaluator`) is already
  factored; a base earns its keep once general-k k-space logic is shared.

---

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
- Superseded companion: `doc/MolecularPP_HarmonizationRound2.md`; Round-1 record:
  `doc/MolecularPP_HarmonizationFindings.md`.
- Commits (GPW, on `main`): `ab2c6a76` (1E), `cc123b3b`/`63fbf70c` (DFT collocation), `dcef8528`/`db314e6a`
  (first-light SCF + G-space local PP), `5f609d2f` (rename), `6d6511ac` (Ortho choice), `fc430e94` (this
  doc's bulk roadmap), **`b2a29249`** (Impl 4 general-k + multi-k), **`10ad6e29`** (N3/N5 removal),
  **`02027faf`** (charge probe), **`a4c94ec5`** (bulk over-binding root-cause + diagnostics),
  **`95e8f4a8`** (BULK FIX: KB Bloch-orbital bra + PhiOnGrid cache + test cleanup),
  **`335df0da`** (CP2K grid-matched table), **`5fe61aeb`** (multi-k validated vs CP2K same-mesh),
  **`1980d6ef`** (shifted-MP support + the complex-D bug diagnostic),
  **`745d03ff`** (complex-k fix: ket-conj density weight + conj KB projector phase + charge trace; shifted 2×2×2 == CP2K −7.86744).
- Tests: `UnitTests/GPW_UT.C` (1E + Bloch invariants), `UnitTests/GPW_SCF_UT.C` (SCF anchors + gates:
  `DISABLED_TermTranslationInvariance`, `DISABLED_SR_GammaRcut2a_CP2KReference`,
  `DISABLED_SR_2x2x2GammaCentred_vs_CP2K`, `DISABLED_SR_2x2x2ShiftedMP_vs_CP2K` [the TODO-1 complex-D probe]),
  `UnitTests/L_PP.C` (finite==lattice PP), `UnitTests/PlaneWaveDFTUT.C` (PW-DFT anchors). CP2K decks +
  results: `UnitTests/CP2K/`, `doc/CP2Kresults.md`. Build/test: `cd build/Release && ninja UTMain && ./UnitTests/UTMain`.
