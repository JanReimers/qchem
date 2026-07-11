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

---

# TODO / NEXT

Bulk energy (Γ), Γ-centred multi-k dispersion, AND the CP2K-default shifted-MP mesh (complex k) are all
**DONE and CP2K-validated** (see DONE — the complex-k fix landed 2026-07-10). Full-BZ GPW works at any k.
Remaining: (1) low-q multi-species bases → Si/NaF/CsI cross-validation (the active NEXT work); (2) the CP2K
reference library (the oracle for §1); (3) IBZ; (4) cleanups.

## 1. Low-q multi-species bases → Si/NaF/CsI cross-validation (PW + GPW + CP2K) — THE NEXT WORK

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
40-a.u. exponent forces a fine density grid). NEXT: converge `densityEcut`, add Rcut images, then CP2K.

**Gates / deliverables.** `doc/CP2Kresults.md` rows Si/NaF/CsI × {PW, GPW, CP2K} (Etot + runtime); `GPW_SCF`
NaF/CsI converge (charge, Etot) == CP2K same-basis; the GPW−PW gap documented (basis quality). **Pitfalls:**
iodine is the first GTH Gaussian basis for the element (validate its pseudo-atom carefully); F's tight 2p is
the hardest (needs the highest cutoff, per the PW NaF vs CsI experience — F set the cutoff, not the heavy I).

## 2. CP2K reference library (the oracle for §1) — BUILT; growing it
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

## 3. Symmorphic space groups + auto-IBZ (deferred qcSymmetry track)
Only the POINT-GROUP part matters for k-reduction. Cubic Si (O_h, centrosymmetric) → IBZ = 1/48 of the BZ,
no time-reversal factor. Validate bit-level: reduced-mesh-with-weights == full-mesh (against Step 2). The IBZ
is an *efficiency* layer, not a correctness requirement — hence it comes AFTER a working full-BZ reference.

## 4. Deferred cleanups (do once bulk works — "the working code is the definitive declaration")
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

# OPEN INVESTIGATION (2026-07-11): why is the truncated Bloch overlap S indefinite? (for the next session)
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
