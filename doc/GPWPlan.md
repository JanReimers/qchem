# GPW (Gaussian And Plane Waves) — Plan & State

**Self-contained orientation for a fresh session.** GPW = Gaussian orbitals on a periodic lattice, with
the electron density collocated on a real-space grid, FFT→G-space Poisson for Hartree, XC on the grid,
integrated back against the Gaussians to form the KS matrix (CP2K / Lippert–Hutter). It is the north-star
that makes ab-initio solids → battery voltage curves possible.

This doc supersedes the GPW sections of `doc/MolecularPP_HarmonizationRound2.md` (now misnamed — that doc's
job was harmonizing the molecular ↔ plane-wave pseudopotential paths; GPW outgrew it). The durable
leftovers from it are folded in below (§5).

---

## 1. Where we are (committed on `main`, all green)

Two increments done, 181/181 UTMain green, `allTests` builds. GPW is a **new evaluator, not a new IBS** —
it satisfies the existing plane-wave concepts and reuses the `EPW_*` mixins.

**Increment 1 — periodic Gaussian 1-electron integrals at Γ** (`ab2c6a76`)
- `GPW_Evaluator` (`src/BasisSet/Lattice_3D/Evaluators/GPW/`) satisfies `isPW_1E_Evaluator`; `GPW_IBS`
  (`src/BasisSet/Lattice_3D/GPW_IBS.C`) = `EPW_Orbital1E_IBS<GPW_Evaluator>` + identity. Scalar type is
  **`dcmplx`** (reuse the PW stack; Γ = k=0, real values stored complex).
- The overlap/kinetic(`⟨p²⟩`)/nuclear matrices are real-space Bloch lattice sums `M_ij=Σ_R⟨χ_i|Ô|χ_j(·−R)⟩`.
  These are **delegated** to the molecular Gaussian basis via a new engine-neutral capability
  **`Molecule::LatticeSum1E`** (`src/BasisSet/Molecule/LatticeSum1E.C`): `MakeOverlap/MakeKinetic/
  MakeNuclear(Rs[,cl])`, realised by `PG_Cart::Orbital_IBS` forwarding to `PG_Cart_MnD::NR_Evaluator`,
  which lattice-sums the existing analytic McMurchie–Davidson kernels via `GaussianRF::AtCenter`. GPW reaches
  it by an **abstract→abstract cross-cast** — no Gaussian internals cross into qcLattice_BS.
  - New library edge `qcLattice_BS → qcMolecule_BS` (no cycle).
  - Engine seam: `LatticeSum1E` is integral-engine-neutral. MnD implements it today; **libCint (`PG_LibCint`)
    is the faster follow-up** — implement the face there, flip the `Molecule::Factory` `Engine::` arg, GPW
    unchanged.
- Validated (`UnitTests/GPW_UT.C`): home cell `R={0}` reproduces the finite matrices exactly (`<1e-12`);
  images give textbook large-cell convergence (a=14→1.5e-2, 22→7.9e-6, 30→5.9e-11).

**Increment 2 — DFT-tier collocation (Hartree/XC machinery)** (`cc123b3b`, cleaned `63fbf70c`)
- GPW satisfies `isPW_DFT_Evaluator` and reuses the **entire** PW-DFT stack (`PW_Hartree`, `PW_XC`,
  `IrrepCD`) by filling the `Repulsion3C`/`Overlap3C` 3-centre tensors with dense collocation weights.
- `G_ERI3` (`src/BasisSet/Internal/GMap.C`) gained `std::vector<mat_t<dcmplx>> weights` (empty = PW delta;
  filled = GPW); `ContractG_ERI3` branches.
- `GPW_Evaluator` holds a density grid (`PW_Grid_Evaluator` at `densityEcut`, Γ). It collocates
  `W_c(i,j)=(1/Ω)∫χ_iχ_j e^{-iG_c·r}=GridCoeff(ForwardFFT(χ_iχ_j),G_c)`, and `OverlapMatrix(f)` grid-
  integrates the adjoint `⟨χ_i|V|χ_j⟩`. Tensor caching is delegated to the framework (`Band_FT_IBS::
  Repulsion3C/Overlap3C` via `theCache<dcmplx>()`); `densityEcut` is in the cache key (`IDFragment`).
- Validated: G=0 collocation weight `W_0·Ω` == analytic overlap to grid accuracy (4.2% @ Ecut=30, Si), and a
  constant field integrates back to `V0·S` with the **identical** residual (forward-collocation == adjoint
  consistency).

**Naming (`5f609d2f`) — remember these:**
- `Overlap(f)` = ANY one-electron `⟨i|f|j⟩` (f may be a potential or not); `Repulsion` = the two-electron
  `1/r12` integral. So the reciprocal-space field→KS-matrix bridge is `Band_FT_IBS::MakeOverlap(f)` /
  evaluator `OverlapMatrix(f)` — **not** the old `MakePotential`/`PotentialMatrix`. `Make` = uncached (f
  changes every SCF iteration). It coexists with the no-arg `MakeOverlap()`/`OverlapMatrix()` (`⟨i|j⟩`) via
  `using` declarations (needed in any IBS combining the 1E + DFT mixins).
- **The Coulomb factorisation (an important GPW feature):** `W_c(i,j)` is a SINGLE-`r` integral (the density
  side). The second electron coordinate and `1/r12` are the diagonal Poisson kernel `4π/|G_c|²` (the `r_2`
  integral in reciprocal space). Full two-electron Coulomb = weight × kernel, factorised through G-space —
  never a 2-`r` integral.

---

## 2. THE mesh insight — the next increment builds on this (do not use the stopgap)

Increment 2's collocation quadrature (`GPW_Evaluator::BuildWeights`/`OverlapMatrix(f)`, hand-rolled
`PhiOnGrid` + `Integral` on the FFT grid) is a **stopgap**. The intended, clean design routes GPW's
real-space quadrature through the **existing `qcMesh` machinery**, which already has exactly the right
abstraction:

- `src/Mesh/Quadrature.C` exports `WeightedOverlap(const Mesh& m, const VectorFunction<T>& χ, const
  ScalarFunction<double>& V) -> hmat_t<T>` = `⟨χ_i|V|χ_j⟩`, plus `Overlap(m,χ)` (`⟨i|j⟩`), `KineticGrad2`,
  `Integrate`. These take the basis + field **directly** (no `GridPoints()` leak; the mesh owns points+weights).
- `Structure::CreateIntegrationMesh(mp)` is THE seam. `UnitCell::CreateIntegrationMesh`
  (`src/Structure/Imp/UnitCell.C`) returns a uniform `qcMesh::Mesh`: `mp.eCut>0` → `n=ceil(2·maxEdge·
  √(2 eCut)/π)` — its own comment calls this **"the GPW / Nyquist path"** and notes *"plane-wave DFT
  integrates in G-space so it never asks for this"* — **GPW is the intended caller.** Weights `Ω/n³`, points
  at fractional midpoints `A·((i+½)/n)` (a proper quadrature mesh, offset ½-cell from the FFT raster `A·(i/n)`
  — so it is NOT literally the FFT grid).

**Target design (decouple "FFT for Poisson" from "mesh for the KS matrix"):**
- KS-matrix quadrature (`⟨χ_i|V|χ_j⟩`, `∫ρ`, `∫ε_xc ρ`) → `qcMesh::WeightedOverlap`/`Integrate` on
  `cell.CreateIntegrationMesh({.eCut = densityEcut})`, with the G-space potential `V_H`/`v_xc` presented as a
  `ScalarFunction` (evaluated via `G_FieldEvaluator::EvalField`) — exactly the molecular `FittedVxc` path.
- The FFT grid (`PeriodicGridEvaluator`) does **only** the Poisson solve `ρ→ρ̃→V_H(G)`.

This encapsulates `{r}` (the grid-encapsulation principle: expose operations, keep the raw point list
inside; `GridPoints()` raw exposure is the soft spot), and **unifies GPW's real-space assembly with the
molecular DFT/PP path** — `L_PP` already runs `WeightedOverlap`-style terms bit-consistently on a `UnitCell`
uniform mesh. It reframes GPW cleanly: **kinetic / nuclear / XC are real-space mesh quadrature (reuse the
molecular machinery); only Hartree needs the FFT (for Poisson).** The `GridPoints()` consumers that feed the
FFT (the XC fitter's batch-sample, GPW's `PhiOnGrid`) legitimately keep a "sample on grid" op — but on the
grid engine, `{r}` internal.

---

## 3. Next increment — the full periodic SCF total energy

Build it on §2 (route the KS quadrature through `qcMesh::WeightedOverlap` on `CreateIntegrationMesh`), since
the SCF needs that real-space quadrature anyway. Remaining pieces:
1. A `dcmplx` **external** term for the GPW orbitals. Kinetic reuses `PW_Kinetic` (calls `MakeKinetic`).
   For the external potential: either an all-electron nuclear term, or (recommended, GPW's whole point) the
   **molecular pseudopotential** on the lattice (`PP_Local`/`PP_NonLocal` already run on a `UnitCell` uniform
   mesh — `L_PP` proves it — but they are `<double>`; need the `dcmplx`/Bloch path, or the G-space `PW_Pseudo`
   if the orbital basis can supply the form factors, which Gaussians cannot → use the real-space PP terms).
2. The **G=0 neutralising background / alignment** (periodic only; `PW_Pseudo` already carries the pattern).
3. **Seeding** (SAD / uniform) + the `cSCFIterator<dcmplx>` drive + `Ham_GPW_DFT` (or reuse `Ham_PW_DFT`
   term-assembly with the GPW basis; the terms are basis-agnostic via the abstract faces).
4. **Validation gate:** Γ molecule-in-a-box GPW total energy == density-fit DFT on the same system to
   grid-cutoff tolerance, and variational-stable as the cutoff rises. Mirror `L_PP`: converge to the finite
   molecular DFT energy as the box grows.

**Open structural decision (settle first):** `double` vs `dcmplx`. The old spec argued GPW is all-`double`
with the FFT a term-private Poisson technique (a `<double>` `GPW_Hartree`); the built evaluator network reuses
the `dcmplx` PW stack (Γ = k=0). At Γ/real orbitals `double` is enough and cleaner; general-k Bloch orbitals
force `dcmplx`. Increment 2 committed to `dcmplx` (reuse); revisit if a `<double>` Γ path proves worth it.

### First SCF validation system — recommendation: **Si**, then NaF, then CsI
- **Si (do first).** Single-species, covalent, 8 valence e⁻ (2 × Zion=4), `l≤1` GTH-LDA projectors. Uses the
  **SIPP** valence Gaussian basis already exercised by `L_PP` and the GPW tests — well-conditioned, moderate
  cutoff. Covalent → no Madelung/charge-transfer convergence trouble. There is a trusted reference: the same
  SIPP basis + GTH PP as a finite molecule via the molecular density-fit DFT path (large-box limit), and the
  existing PW-DFT Si anchors (`PlaneWaveDFTUT`: Si-Γ, Si 2×2×2). This is the cleanest first light.
- **NaF (follow-up).** Multi-species (Na Zion=1 + F Zion=7) and ionic. Exercises the multi-species PP routing
  and ionic-SCF convergence (needed a raised DIIS `EMax` in the PW path — ionic charge-transfer). F's tight 2p
  (`r_loc=0.219`) is the hard atom → higher `densityEcut`.
- **CsI (follow-up).** The **d-projector** (`l=2`) validation (both species carry `l=2` KB channels — Si/NaF
  are `l≤1`). Heavy, soft → a LOWER cutoff than NaF. Use Cs-q1 (the semicore q9 has an `l=3`/f projector the
  analytic HGH `Qli` table omits).

Rationale: get first light on the simplest, most-anchored covalent single-species case; add multi-species +
ionic (NaF) and d-projectors (CsI) as separate, targeted follow-ups. (All three GTH PPs + the multi-species
PW machinery already exist and are green in `PlaneWaveDFTUT`.)

---

## 4. Deferred cleanups (do once the SCF works — "the working code is the definitive declaration")
- **Route quadrature through `qcMesh` (§2)** — do it WITH the SCF increment (it needs the quadrature).
- **Whole-density collocation** — the dense `W` tensor is `O(nGfit·nAO²)` storage + `O(nAO²)` FFTs. The
  efficient GPW collocates `ρ = op(r)` **once** (one FFT), which needs `D` → density-side. CP2K's local-patch
  (multi-grid) collocation is a further v2.
- **Design-pure framing** (debated, parked in favour of the working reuse): GPW Hartree/XC as a **distinct
  Coulomb-strategy term** (`GPW_Hartree`/`GPW_XC`, peer of `Vee`/`FittedVee`, grid term-private — the Factory
  picks it by structure type) rather than reuse-via-`G_ERI3`. Related: a `ProjectedDensity_R` fitter variant
  (also parked).
- **Common `PW_Evaluator`/`GPW_Evaluator` base** — DEFERRED to general-k. The shared FFT engine
  (`PeriodicGridEvaluator`) is already factored (both hold it); the residual orbital-evaluator overlap is just
  the k-point. It earns a base only when general-k adds real shared k-space logic (`e^{ik·R}`, `|k+G|`).

---

## 5. Durable pins / invariants (carry into all GPW work)
- **PP-smoothness is GPW's enabler; GAPW is out of scope (first pass).** Plain GPW is clean only when the
  density is smooth enough to collocate on a modest grid — pseudopotentials give that. All-electron cores are
  too sharp (would need PAW-like augmentation = GAPW, deferred). Validate with a well-conditioned GTH valence
  basis, never all-electron.
- **GPW is a Coulomb/Hartree STRATEGY orthogonal to the orbital basis** — a third one beside exact-4-centre
  (`Vee`) and density-fitting (`FittedVee`). Same `⟨χ|V_H|χ⟩` out, different internals; it does not restructure
  the SCF. Mirror `FittedVee` (Factory picks by structure type: molecules → density-fit; solids/supercells →
  GPW's O(N log N) FFT Coulomb).
- **Never assume `orbital == fit`.** Any fit/aux basis comes from the orbital basis via `Create{CD,Vxc}
  FitBasisSet(...)` — the factory is the seam even when the answer is trivial.
- **Fit quality is measured by grid-convergence of ρ, NEVER by ΔE_total** (the fit is non-variational).
- **Spin-polarized is the native formulation**; unpolarized is the ζ=0 collapse. New GPW/periodic terms are
  spin-native (`FittedVxcPol`/`FittedVcorrPol`), unpolarized as the efficiency corner.
- **Use well-conditioned bases for SCF** (Slater/High for atoms; a cleanly-converted GTH valence basis for
  PP). "LASolver" symptoms are basis conditioning, not a solver gap. `N3`/`N5` are test-only pools.
- **Regression style:** periodic/GPW energies are "did-E-move" anchors (pin the converged value, no
  `Converged()` guard). Where a real-space-on-lattice quantity must equal its finite counterpart, assert
  bit-consistency (`L_PP`-style) rather than an absolute oracle.

### Symmetry comes AFTER a working GPW (independent optimisation layer, does not gate GPW)
GPW runs at Γ / a small explicit k-mesh with no space-group machinery. Once GPW works, the symmetry track
(from the Round-2 doc) bolts on without rework: symmorphic space groups in `qcSymmetry` → BZ reduction
(irreducible wedge, the highest-value solid speedup, `reduced+weighted == full`) → SALC with plane waves
(star-of-G blocking). None of these gate GPW.

---

## 6. Pointers
- Superseded companion: `doc/MolecularPP_HarmonizationRound2.md` (§2.4/§2.5 GPW discussion; §2.1–2.2 the
  deferred symmetry track). Round-1 record: `doc/MolecularPP_HarmonizationFindings.md`.
- Commits (this GPW work, on `main`): `09cf4170` facade-preserves-UnitCell, `17fd8f0c` MakePotential→orbital
  face (pre-GPW), evaluator-network refactor (`4eb26bff`/`b45532f7`/`85765dfb`/`65128554`), then
  **`ab2c6a76`** (1E), **`cc123b3b`** (DFT collocation), **`63fbf70c`** (cache/doc), **`5f609d2f`** (rename).
- Tests: `UnitTests/GPW_UT.C` (GPW), `UnitTests/L_PP.C` (finite==lattice PP), `UnitTests/PlaneWaveDFTUT.C`
  (PW-DFT anchors incl. Si/NaF/CsI). Build/test: `cd build/Release && ninja UTMain && ./UnitTests/UTMain`.
