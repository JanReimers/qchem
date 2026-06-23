# Plane-Wave Calculation — Implementation Plan

Prep notes for implementing a plane-wave (PW) electronic-structure calculation,
written 2026-06-22 at the end of the `src/Structure` groundwork session. Symbol
and units conventions are documented in `src/Structure/Lattice_3D/Lattice_3D.C`
(Doxygen).

## 1. Geometry groundwork — DONE

`src/Structure` is ready for periodic work:
- `UnitCell` stores the cell matrix `A` (columns = lattice vectors `aᵢ`); `ToCartesian(f)=A·f`.
- `Lattice_3D` (was `Lattice`; module `qchem.Lattice_3D`, folder `src/Structure/Lattice_3D/`) — renamed in anticipation of `Lattice_1D` (polymers) and `Lattice_2D` (graphene) as sibling folders.
- `ReciprocalLattice` (`Lattice_3D::Reciprocal()`): `B = 2π·A⁻ᵀ`, `GetGVectors(Gmax)`, `GetGLength(m)`.
- `KMesh` (`Lattice_3D::MakeKMesh()`): Monkhorst–Pack `(k,weight)` list (Γ-centred or shifted), ready for later IBZ symmetry reduction.
- `UnitCell::CellsInSphere` for real-space lattice sums (Ewald real part later).

## 2. Validation targets (staged, simplest first)

Do them in this order — each isolates a different piece:

1. **Empty lattice (free electron), V = 0.** Eigenvalues are exactly `½|k+G|²`.
   Validates G-vector generation, k-points, the kinetic operator, and the
   eigensolver with zero potential ambiguity. The standard first PW check; exact.
   **DONE (2026-06-22):** `BasisSet::Lattice_3D::PlaneWave_IBS` (a concrete
   `Orbital_1E_IBS<dcmplx>`) in `src/BasisSet/Lattice_3D/`, lib `qcLattice_BS`,
   test `UnitTests/PlaneWaveUT.C` (Γ + off-Γ) matches `½|k+G|²` to 1e-10.
2. **Separable cosine potential** `V₀·Σᵢcos(2π xᵢ/a)`. Its PW matrix elements are
   nonzero only between `G` and `G ± bᵢ` (sparse, exactly representable). 1D case
   = Mathieu equation with known bands. Validates potential-matrix assembly with
   no singularity.
   **DONE (2026-06-22):** reusable `PlaneWave_IBS::MakePotential(Vtilde)` assembles
   `⟨G|V|G'⟩ = Ṽ(m−m')` (the cosine and the milestone-2.3 nuclear structure factor
   are both just `Ṽ` suppliers). Tests: exact sparse V structure, traceless ⇒
   `ΣE = Σ½|k+G|²` exactly, `V₀→0` recovers milestone 1, and the Γ ground state
   matches 2nd-order PT `−3V₀²a²/(4π²)` to 1%.
3. **Hydrogen, bare Coulomb, large `a` → −0.5 Ha/cell** (the headline physics).
   Caveats — bake into the tolerance:
   - The 1s **cusp converges slowly** in PW (`FT(−1/r)=−4π/G²`); approach −0.5
     **from above**, needs high `E_cut`. Won't hit −0.5 to 1e-4 at modest cutoffs
     without a smooth/pseudized nucleus.
   - The **G=0 Coulomb term diverges** and is conventionally dropped (uniform
     neutralizing background) → a finite-cell shift that → 0 as `a → ∞`.
   - Large `a` ⇒ flat bands ⇒ **Γ-point only suffices** (so hydrogen does NOT
     stress k-points; the empty-lattice test does that better).
   **DONE (2026-06-22):** `PlaneWave_IBS::MakeNuclear` assembles the bare-Coulomb
   structure factor `V(ΔG) = −(4π/Ω)·Σ_a Z_a e^{−iΔG·τ_a}/|ΔG|²`, `ΔG=0` dropped.
   Tests: exact nearest-G coupling `−1/(πa)` (single H at origin), and end-to-end
   variational convergence (E0 monotone-decreasing in `E_cut`, bound, above −0.5).
   Confirmed the warned slow cusp convergence empirically — at `a=8`, `E_cut`=4→6
   gives only ≈ −0.136→−0.149, so the test asserts convergence *behaviour*, not a
   −0.5 value (which needs far higher `E_cut` / a pseudized nucleus).
   The 1E integral building blocks now return `hmat_t<T>` (`chmat_t` for complex):
   Hermitian for `dcmplx`, identical `SymmetricMatrix` for real `T` (so atom/molecule
   bases are unaffected). `MakeNuclear` fills a `HermitianMatrix`, so off-origin /
   multi-atom cells (complex phases) assemble correctly — the earlier `csmat_t`
   symmetric-storage seam is closed.

## 2b. Beyond bare Coulomb — two lineages, and the pseudopotential ladder

Two distinct families improve on a plain plane-wave (PW) calculation. They have
**different destinations**, and the project aims to demonstrate the shared OOD
framework spans *both* within the solids sector:

- **Lineage A — smooth-potential / pseudopotential.** The basis stays plane waves;
  the *potential* is made smooth so the PW basis converges fast.
  Ladder: local pseudopotential → Kleinman–Bylander separable nonlocal projectors
  (norm-conserving) → ultrasoft → PAW.
- **Lineage B — augmented all-electron bases.** Keep the bare Coulomb, enrich the
  *basis* with atom-centred functions inside muffin-tin spheres: APW → LAPW →
  FLAPW (and LMTO/NMTO). Heavier; the all-electron "gold standard". A *sequence* of
  IBS types, each its own class (and eventually its own Evaluator) under
  `src/BasisSet/Lattice_3D/`.
  **Lineage B, first IBS (DONE) — `APW_IBS`** (`src/BasisSet/Lattice_3D/APW_IBS.C`):
  plane waves in the interstitial, value-matched to free-particle radial solutions
  inside a single origin muffin-tin sphere. Fixed-energy demonstrator: the secular
  matrix `Γ(E)=H(E)−E·O(E)` (`MakeSecular(E)`) is built per energy parameter and is
  singular exactly at the free-electron energies `½|k+G|²` (empty-lattice limit). The
  radial integral collapses to a surface log-derivative `u_l'(R)/u_l(R)` (no radial
  quadrature). Test: `Γ(E)` singular (min|eig|→0, machine-zero by l_max=8) at the
  free-electron energies, non-singular between them, converging with l_max. `APW_IBS`
  shares the `IrrepBasisSet<dcmplx>` base with `PlaneWave_IBS` but **not** the
  energy-independent `Orbital_1E_IBS` (APW's overlap is E-dependent) — confirming each
  lineage-B IBS carries its own interface needs.
  **Lineage B, second IBS (DONE) — `LAPW_IBS`** (`src/BasisSet/Lattice_3D/LAPW_IBS.C`):
  the linearised method. Inside the sphere the radial function is `a_l·u_l + b_l·u̇_l`
  at a fixed linearization energy `E_l`, with value+derivative matching at `R` (2×2
  Wronskian → `a_l,b_l`). `H` and `O` are then energy-independent (`MakeHamiltonian()`,
  `MakeOverlap()`), so a single generalised eigenproblem `Hc=εOc` (`LASolver`) yields
  the whole band — a true band solver, unlike APW's per-energy determinant. Built from
  the symmetric-gradient form (Hermitian; sidesteps the `u/u̇` Hamiltonian asymmetry);
  radial integrals by Simpson (integrands vanish at r=0). Test: the empty-lattice bands
  reproduce `½|k+G|²`, with the linearization error smallest near `E_l` (≈1e-9 on the
  level at `E_l`, growing to ≈1e-4 far from it). Reuses the `qchem::Math` Bessel/Legendre.
  Next: a real (l>0) potential where −0.5 becomes an achievable check.

**Rung 1 (DONE, lineage A) — local potential abstraction.** `LocalPotential`
(`src/BasisSet/Lattice_3D/LocalPotential.C`) is the open/closed extension point:
it supplies only a per-species reciprocal form factor `v(Z,|G|²)`, while
`PlaneWave_IBS::MakeLocalPotential` folds in `1/Ω`, the structure factor, and the
`G=0` drop. Implementations: `BareCoulomb` (`−4πZ/G²`) and `GaussianSmearedNucleus`
(`−4πZ·e^{−σ²G²/2}/G²`, the rung-1 local pseudopotential). `MakeNuclear` is now just
`MakeLocalPotential(cl, BareCoulomb())`. Tests confirm the cusp is tamed: at `a=8,
σ=0.5` the smeared energy is converged by `E_cut=6` (drift ~1e-4 to `E_cut=12`)
while bare Coulomb keeps crawling; `σ→0` reproduces `BareCoulomb` element-by-element.
**Rung 2 (DONE, lineage A) — separable nonlocal.** `SeparablePotential`
(`src/BasisSet/Lattice_3D/SeparablePotential.C`) is the nonlocal sibling of
`LocalPotential`: per-species projector form factors `β̃_p(|q|)` + KB coefficients
`D_p`. `PlaneWave_IBS::MakeSeparablePotential` assembles the Kleinman–Bylander
form `V_NL = Σ_a Σ_p |β^a_p⟩ D_p ⟨β^a_p|` = `(1/Ω) Σ_a e^{−iΔG·τ_a} Σ_p
β̃_p(|k+G|) D_p β̃_p(|k+G'|)`. The external one-body potential is now `V = V_loc +
V_NL` — both model-parameterized contributions summed into `H(k)`. Demonstrator
implementation: `GaussianProjector` (one s-channel projector, `β̃=e^{−σ²q²/2}`);
l>0 needs `Y_lm(q̂)` + spherical-Bessel transforms (the production-PP extension).
Tests: exact rank-1 matrix element `(D/Ω)β̃ᵢβ̃ⱼ`; a single projector yields exactly
one nonzero eigenvalue equal to the trace (separability); and a repulsive projector
added to `V_loc` raises the ground state. The same separable structure underlies
USPP and PAW.

## 3. Scope insight — the 1-electron target needs much less than full SCF

`H = T + V_ext` is a **single diagonalisation** — no self-consistency, no
Hartree/XC, no ERIs. So the first milestone needs only:
- a complex **Hermitian H(k)** in the PW basis (kinetic `½|k+G|²` on the diagonal;
  `V_ext` via the structure factor / G-space potential), and
- a (complex Hermitian) **eigensolver**.

It does **not** need the complex `ChargeDensity` / `WaveFunction` / `SCFIterator`
stack. Those are for the *later* self-consistent calc.

## 4. Complex (`T=dcmplx`) readiness — smoke-test findings

PW Bloch states are complex, so the templated cores must build for `dcmplx`.
Current state (two camps):

- **Templated cores** (`Orbitals`, `ChargeDensity`): templated on `T`, but the
  `dcmplx` explicit instantiations are commented out:
  - `src/Orbitals/Internal/Imp/TOrbital.C:76`  `// template class TOrbitalImp<std::complex<double>>;`
  - `src/Orbitals/Internal/Imp/TOrbitals.C:161` (add `TOrbitalsImp<dcmplx>`)
  - `src/ChargeDensity/Internal/Imp/IrrepCD.C:150` (add `IrrepCD<dcmplx>`)
- **Hardwired modules** (`Hamiltonian`, `WaveFunction`, `SCFIterator`): NOT
  templated — their `Types.C` pin `BasisSet<double>` etc. Making these `dcmplx`
  is a real refactor (templating them), part of the PW work, not a smoke test.

**Smoke test = uncomment those 3 instantiation lines and build qcOrbitals +
qcChargeDensity.** When green, leave them uncommented as the permanent build
check. Uncommenting today surfaces these issues (all bounded, no algorithmic
blockers):

| Issue | Location | Needed for 1e target? |
|---|---|---|
| `Iterate<TOrbital<double>>` + `smat_t<double>` hardcoded — should be `<T>` | `TOrbitals.C:103` | yes (clear bug, fix now) |
| `IrrepCD_Factory(…, const obs_t*)` hardwires `obs_t`=…`<double>`; itsBasisSet is `tobs_t<T>*` | `TOrbitals.C:109`, factory sig | yes |
| `Static_CC/Dynamic_CC::GetMatrix(const obs_t*)` hardwires double | `IrrepCD.C:77,84` | self-consistent stage |
| `blazem::norm(complex matrix)` → needs `std::real(...)` for `double` return | `IrrepCD.C:118` | self-consistent stage |
| `Vector3D<double> = Vector3D<complex>` in density gradient path | `vector3d.C:53` via `IrrepCD::Gradient`, `TOrbitals::Gradient` | no (no UT cov) |
| complex `AccumulateDirect/Exchange` (HF/ERI), `GradientContraction`, `Repulsion3C` | `IrrepCD.C` (double-only specializations / commented) | no (ERI/DFT, later) |

Note IrrepCD already half-handles complex: `DM_Contract`/`GetTotalCharge`/
`operator()` use `std::real(...)` and `blazem::conj(phir)`. The remaining core
work is mostly **templating the hardwired `obs_t`/`smat_t<double>` cross-module
signatures** (the `Types.C` double seam) + the `TOrbitals.C:103` `<double>`→`<T>`
fix. Watch for `blazem::trans` vs `ctrans` (conjugate transpose) once it compiles
— e.g. `TOrbital.C:48 blazem::trans(itsCoeff)*gr`.

## 5. Evaluator framework — the PW basis set should fit it

Integral/value evaluation in this codebase is converging on a shared **Evaluator**
concept that spans Atom and Molecule basis sets. The PW basis set should slot into
the same framework rather than grow a one-off path. See
[MolecularBasisSetPlan.md](MolecularBasisSetPlan.md), "Longer-term direction":

- **B. Lift the Evaluator idea to PG basis sets** — move
  `src/BasisSet/Atom/Evaluators/Evaluator.C` out to
  `src/BasisSet/Internal/Evaluator.C` so any basis set (incl. plane waves) can reuse it.
- **C. Hoist the 1E integral loops** (`Integrals_Overlap/Kinetic/Nuclear`) to a
  basis-set-agnostic level.

That plan **already names plane waves as the risk case**: the Lattice tree has its
own `IBS_Evaluator` (complex-valued, with k-point phase factors), and item C carries
a ⚠ spike — *confirm a single 1E loop signature genuinely unifies Atom / Molecule /
Lattice, because k-point phases may break the "just loop i,j" assumption.* So when
building the PW `IrrepBasisSet` / its evaluator (§6.1), design it to fit that shared
interface, and treat the **k-point-phase unification as the open risk to validate
early** — it constrains the shared Evaluator interface for all three sectors. See
also `[[project_evaluator_framework_friction]]` and `[[project_pg_refactor_plan]]`.

## 6. Fresh-session checklist

1. PW `IrrepBasisSet`: a k-point's `{G : ½|k+G|² < E_cut}` set (from
   `ReciprocalLattice::GetGVectors`), each a normalised plane wave `e^{i(k+G)·r}/√V`.
   New code goes under `src/BasisSet/Lattice_3D/` (placeholder folder reserved;
   `# add_subdirectory(Lattice_3D)` in `src/BasisSet/CMakeLists.txt` to wire it in).
2. Kinetic `MakeGrad2`/`MakeKinetic`: diagonal `½|k+G|²` (see the `Lattice_3D.C` conventions; mind the ½ — `[[project_kinetic_p2_conventions]]`).
3. `V_ext` matrix: `V(G−G')` from the structure factor `Σ_atoms Z·e^{-i(G−G')·τ}` × `FT(potential)`; drop/define the `G=0` term.
4. Complex Hermitian eigensolver per k-point → bands; energy/cell.
5. Targets 2.1 → 2.2 → 2.3 with the tolerances above.
6. For self-consistency later: do the `dcmplx` core fixes in §4, then template
   `Hamiltonian`/`WaveFunction`/`SCFIterator` off their `Types.C` double seam.

See `[[project_structure_groundwork]]` for what was built this session.
