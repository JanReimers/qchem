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
2. **Separable cosine potential** `V₀·Σᵢcos(2π xᵢ/a)`. Its PW matrix elements are
   nonzero only between `G` and `G ± bᵢ` (sparse, exactly representable). 1D case
   = Mathieu equation with known bands. Validates potential-matrix assembly with
   no singularity.
3. **Hydrogen, bare Coulomb, large `a` → −0.5 Ha/cell** (the headline physics).
   Caveats — bake into the tolerance:
   - The 1s **cusp converges slowly** in PW (`FT(−1/r)=−4π/G²`); approach −0.5
     **from above**, needs high `E_cut`. Won't hit −0.5 to 1e-4 at modest cutoffs
     without a smooth/pseudized nucleus.
   - The **G=0 Coulomb term diverges** and is conventionally dropped (uniform
     neutralizing background) → a finite-cell shift that → 0 as `a → ∞`.
   - Large `a` ⇒ flat bands ⇒ **Γ-point only suffices** (so hydrogen does NOT
     stress k-points; the empty-lattice test does that better).

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

## 5. Fresh-session checklist

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
