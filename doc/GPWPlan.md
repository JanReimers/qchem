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

## Bulk over-binding FIXED — GPW bulk matches CP2K to 1e-5 (was TODO 1)
The root cause was **one thing wearing two costumes: an incompletely-wrapped Bloch orbital.**
- **The real bug (KB nonlocal, the 16 Ha term):** `MakeSeparablePP` used the **raw home orbital `*itsOrb`** as
  the projector bra on the single-cell mesh. A boundary-straddling corner atom lost its wrapped tail → `b_i`
  ≈ half (corner trace 21 vs interior 37) → the nonlocal PP was translation-variant by ~16 Ha. **Fix: use the
  Bloch-summed orbital (`Eval`, precomputed on the mesh) as the bra.**
- **The FFT-raster `Vloc`/Hartree/XC term was a RED HERRING:** once the orbital is fully wrapped (`Rcut ≥ 2a`)
  its translation-variance also vanishes (the on-node over-weighting self-corrects when the full periodic
  density is present). So **the voxel-grid-shift (old TODO 1b, Option A) was reverted entirely** — simpler.
  Both terms go to Δ = 0.0000 at `Rcut ≥ 2a`.
- **Validation vs CP2K (Γ, SIPP_SR, Rcut=2a):** Etot **−7.11505** (CP2K −7.11506), charge 8, Exc −2.544
  (CP2K −2.544). Nonlocal-PP term hits CP2K's +0.9406. **Committed anchors safe:** at `Rcut=0`, `Eval` = the
  raw orbital, so `SiliconGammaConverges` (−8.24758) and the atom-in-box are unchanged.

## Naming (`5f609d2f`) — remember these
- `Overlap(f)` = ANY 1-electron `⟨i|f|j⟩` (f may be a potential); `Repulsion` = the 2-electron `1/r12`.
  So the reciprocal-space field→KS-matrix bridge is `Band_FT_IBS::MakeOverlap(f)` / evaluator
  `OverlapMatrix(f)` — **not** `MakePotential`/`PotentialMatrix`. `Make` = uncached.

---

# TODO / NEXT

The bulk-energy blocker is **DONE** (GPW bulk matches CP2K to 1e-5 — see "Bulk over-binding FIXED" in DONE).
What remains, roughly in order: (1) full-BZ multi-k GPW (the next real feature); (2) grow the CP2K reference
library; (3) IBZ; (4) deferred cleanups.

## 1. Full-BZ multi-k GPW (the next feature)
General-k GPW + the multi-k `GPW_BasisSet` plumbing are DONE (Impl 4); the single-k dispersive total is now
correct (bulk fix). Remaining: run the actual BZ dispersion — one `GPW_IBS` per k over `lat.MakeKMesh()`
(already wired), converge `Rcut ≥ 2a` (needed for the full orbital wrap that the bulk fix relies on), and
gate the BZ-summed total against a CP2K Monkhorst-Pack run (TODO 2). NB CP2K validation only needs a k-mesh
once qchem has multi-k GPW dispersion — a single-k CP2K Γ reference (DONE) already validated the fix.

## 2. CP2K reference library (the oracle) — BUILT; growing it
CP2K's Quickstep **is** the reference GPW implementation (Lippert–Hutter); its per-term breakdown points
straight at a bug (as this session's did).
CP2K's Quickstep **is** the reference GPW implementation (Lippert–Hutter); its per-term breakdown points
straight at a bug (as this session's hand-rolled breakdown did: Een ×15.7 → local PP → the raster).
- **DONE — CP2K 2026.1 built** (serial ssmp, gcc 15.2) at `~/Code/cp2k` (sibling to qchem6, outside the git
  tree). Toolchain: OpenBLAS+FFTW+libxc+libxsmm+DBCSR, no MPI/libint. Build: `tools/toolchain/build_cp2k.sh`
  (CMake, NOT the old arch-file `make`). Run needs `source install/setup` +
  `LD_LIBRARY_PATH=install/lib`.
- **DONE — FCC-Si Γ reference (SIPP_SR, GTH-PADE-q4, LDA_X+VWN5):** **Etot = −7.11506 Ha, charge 8**,
  converged by `CUTOFF` 80 Ry (≈40 Ha). Breakdown: Core-H (kin+PP) +5.565, Hartree +10.380, XC −2.544;
  PP total −7.548 (local −8.489, nonlocal +0.941); core self-energy −20.516. (CP2K's GPW electrostatic split
  differs from ours — compare the TOTAL + the cleaner sub-terms kin/XC/nonlocal-PP.) **This is TODO 1's gate.**
  Input files: **`UnitTests/CP2K/`** (`si_fcc_gpw.inp` + `SIPP-SR-BASIS` + README) — a growing library of
  reference decks.
- **PP already aligned:** our `src/Pseudopotential/Data/gth_potentials.json` IS the CP2K GTH-PADE database
  (Si GTH-PADE-q4 params match ours exactly — verified). **Basis: same exponents, transcribed to CP2K
  `BASIS_SET` format** (uncontracted → one set per primitive; see `UnitTests/CP2K/SIPP-SR-BASIS`).
- **NEXT with CP2K:** multi-k (uncomment `&KPOINTS MONKHORST-PACK`) as the full-BZ reference (TODO 3); NaF/CsI
  decks; term-by-term breakdown once our fixed GPW runs the same geometry.

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
| k-points | `MakeKMesh` (MP, Γ-centred) | `&KPOINTS SCHEME MONKHORST-PACK` | match mesh + Γ-centring; single-point = `GAMMA` |
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
- **Whole-density collocation (efficiency):** the dense `W` tensor is `O(nGfit·nAO²)` storage + `O(nAO²)` FFTs.
  The efficient GPW collocates `ρ = op(r)` ONCE (one FFT), which needs `D` → density-side. CP2K's local-patch
  (multi-grid) collocation is the further v2.
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
- Superseded companion: `doc/MolecularPP_HarmonizationRound2.md`; Round-1 record:
  `doc/MolecularPP_HarmonizationFindings.md`.
- Commits (GPW, on `main`): `ab2c6a76` (1E), `cc123b3b`/`63fbf70c` (DFT collocation), `dcef8528`/`db314e6a`
  (first-light SCF + G-space local PP), `5f609d2f` (rename), `6d6511ac` (Ortho choice), `fc430e94` (this
  doc's bulk roadmap), **`b2a29249`** (Impl 4 general-k + multi-k), **`10ad6e29`** (N3/N5 removal),
  **`02027faf`** (charge probe), **`a4c94ec5`** (bulk over-binding root-cause + diagnostics).
- Tests: `UnitTests/GPW_UT.C` (1E + Bloch invariants), `UnitTests/GPW_SCF_UT.C` (SCF + the DISABLED bulk
  diagnostics), `UnitTests/L_PP.C` (finite==lattice PP), `UnitTests/PlaneWaveDFTUT.C` (PW-DFT anchors incl.
  Si/NaF/CsI + the 2×2×2 multi-k reference). Build/test: `cd build/Release && ninja UTMain && ./UnitTests/UTMain`.
