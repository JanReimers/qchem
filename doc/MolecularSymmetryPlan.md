# Molecular Symmetry / SALC Plan

Design notes for adding point-group symmetry (and, eventually, space / magnetic /
Brillouin-zone symmetry) to the molecular SCF path. Captured from a design discussion;
this is a living document.

## ⟢ RESUME HERE (current status; reindexed 2026-06-21 after the unit-test migration)

**Done and committed** — the whole group-theory engine + the bridge to the real basis:

| Stage | What | Where | Tests |
|---|---|---|---|
| 1 detect | `DetectPointGroup` → Schoenflies + abelian subgroup | `src/Symmetry/PointGroup.C` | `M_PointGroup.C` |
| 2 represent | `CartesianShellRep`, `BuildOperationRep` (AO rep `M(g)`) | `src/Symmetry/CartesianRep.C` | `M_CartesianRep.C` |
| 3a tables | `AbelianCharacterTable` (8 abelian, Mulliken labels) | `src/Symmetry/CharacterTable.C` | `M_CharacterTable.C` |
| 3a ops | `BuildAbelianGroup` (concrete ops in molecule frame) | `src/Symmetry/AbelianGroup.C` | `M_PointGroup.C` |
| 3b SALC | `BuildSALCs` → transform **O** (irrep-blocked, labelled) | `src/Symmetry/SALC.C` | `M_SALC.C` |
| 5a bridge | `ExtractAoShells` + `ClusterToSymPoints` (PG basis → symmetry) | `src/BasisSet/Molecule/PolarizedGaussian/Symmetry.C` | `M_PGSymmetry.C` |

`qcSymmetry` stays LAPACK-free. ~40 tests, all green. **O is produced and validated on a real
24-function H₂O basis** (every SALC column an irrep eigenvector `M(g)v = χ(g)v`).
**Build/test (per-library targets, post-migration — `UTMain` is gone):** symmetry =
`cmake --build build/Debug --target UTSymmetry -j4` then
`./build/Debug/src/Symmetry/tests/UTSymmetry`; the PG bridge =
`--target UTMolecule_BS` then `./build/Debug/src/BasisSet/Molecule/tests/UTMolecule_BS`. Test
sources live in `src/<Lib>/tests/`. (`src/BasisSet1` was renamed to `src/BasisSet`.)

**Decision locked (the 2-electron fork): Option 1 — the decorator.** Build J and K in the AO
basis (existing engine untouched), then transform the 2-index density/Fock in/out:
`D_AO = Σ_Γ O_Γ D_Γ O_Γᵀ`, build `F_AO`, slice `F_Γ = O_Γᵀ F_AO O_Γ`. RAM/CPU identical to the
simpler "symmetry-adapted diagonalization" option but it packages each irrep as an IBS, so it
reuses the whole per-irrep SCF/accelerator stack AND the Mulliken labels reach the orbital
tables for free. (4-index ERI transform rejected on RAM — transformed ERIs ≈ raw ERI tables.)
Caching: use the existing global cache; the decorator's `MakeOverlap()` returns
`O_Γᵀ·raw->Overlap()·O_Γ` (raw is public+cached) and memoizes the transformed block under an
irrep-specific `(RadialID, AngularID)` key — both raw and transformed cached, no new cache code.
The density is block-diagonal by irrep (one isolated `D_Γ` per irrep); blocks couple only
through the totally-symmetric mean field, which is why "build `F_AO`, slice" is natural.

**STAGE 4 IS DONE (1-e + 2-e, all tested with real integrals).**  `SymmetryAdapted_IBS`
(`src/BasisSet/SymmetryAdapted_IBS.C`) wraps the raw whole-molecule IBS + an `O_Gamma` block:
1-electron = `O^T M_raw O` (cached per irrep); 2-electron = build the AO Coulomb/exchange from
each cd-irrep density block via the raw IBS (linear in D, NO 4-index transform) and slice
`O^T J_AO O` -- `Orbital_HF_IBS::AccumulateDirect/Exchange` made virtual for this.
`SymmetryAdaptedBasisSet` (`.../PolarizedGaussian/SymmetryAdaptedBasisSet.C`) presents one
labelled `SymmetryAdapted_IBS` per irrep; `MolecularIrrep` carries the Mulliken label.  Tests:
`src/BasisSet/Molecule/tests/M_PGSymmetry.C` (UTMolecule_BS) -- 1-e blocks vs real overlap,
SAB iteration + labels, and the 2-e Coulomb matches `O^T J_AO O` on a real H2O s+p basis.

**END-TO-END SYMMETRIC HF WORKS (water).**  `UnitTests/M_PG_Sym.C` (in UTMain) runs HF for
water two ways -- plain single-IBS vs `SymmetryAdaptedBasisSet` with a fixed C2v occupation
`{A1:6,B1:2,B2:2}` (`MolecularSym_EC`, `src/ElectronConfigurations/MolecularSym_EC.C`) -- and
the energies agree to 1e-5.  So the whole stack is validated; molecules reuse the atomic
per-irrep WaveFunction/ChargeDensity/accelerator machinery unchanged.  (Cache note: the
decorator's `MakeOverlap/Kinetic/Nuclear` call the raw COMPUTE `raw->MakeX()`, NOT the cached
`raw->X()` -- the global cache's stateful Has/Set "last key" gets clobbered by a nested cached
access mid-Set.  So raw 1-e is recomputed per irrep (cheap); the 2-e AO ERIs stay cached.)

**NEXT STEPS:**
1. **General global aufbau** -- a `CompositeAufbauWF` (subclass of CompositeWF) that overrides
   `FillOrbitals` to collect all irreps' orbital energies, fill the globally-lowest N/2, and
   set per-irrep occupations each iteration.  Replaces the hand-coded `MolecularSym_EC`.
2. Molecular Factory hook (bridge -> detect -> BuildSALCs -> wrap in SAB) + a polarized path.
3. Optimisation: the 2-e decorator rebuilds the AO Coulomb/exchange once per (irrep,cd-irrep)
   pair (N^2 per iteration); sum the back-transformed densities first, build once, slice all.

(Below: the original blocker analysis, now resolved.)
**(resolved) the OCCUPATION/AUFBAU blocker:**
- `CompositeWF::FillOrbitals` does `w->FillOrbitals(itsEC)` per irrep, i.e. each irrep takes a
  FIXED electron count `ec->GetN(irrep)`.  Right for atoms; wrong for a symmetric molecule,
  where the occupied MOs fill the *globally-lowest* across irrep blocks (the per-irrep
  occupation is an SCF result, not known a priori).  Need a **global aufbau across irreps**:
  collect all irreps' orbital energies, fill the lowest N/2 globally, set per-irrep occupations.
  (First-validation shortcut: hand-code a known molecule's per-irrep occupation -- e.g. water
  a1:3,b1:1,b2:1 -- into the EC, run HF with fixed occupations, compare energy; then generalise.)
- Then a molecular Factory hook that wraps a basis in a `SymmetryAdaptedBasisSet` (extract
  shells+points via the bridge -> detect -> BuildSALCs -> wrap), and an H2O HF run vs the
  non-symmetric energy, confirming irrep-labelled orbitals.

(Historical: the original stage-4 plan -- the 1-electron decorator -- below; now done.)
**(done) stage 4, the 1-electron decorator:**
1. `SymmetryAdapted_IBS` (per irrep): IS-A `Orbital_1E_IBS`, holds `{raw Orbital_1E_IBS*, O_Γ
   block (nAO×dΓ), irrep sym_t}`. `MakeOverlap()=O_Γᵀ·raw->Overlap()·O_Γ` (same for Kinetic/
   Nuclear), `GetNumFunctions()=dΓ`, `GetSymmetry()`=the irrep, `RadialID()/AngularID()`=raw's
   + irrep suffix (cache key). Lives in BasisSet (PG area or a generic SymmetryAdapted module).
2. `SymmetryAdaptedBasisSet` (IS-A `BasisSet<double>`): holds raw molecular basis + O (from
   `BuildSALCs` via the stage-5a bridge), creates one `SymmetryAdapted_IBS` per irrep block.
3. Test: build the H₂O/N₂ basis, make the irrep IBSs, assert `O_Γᵀ S_raw O_Γ` block-diagonal
   and the per-irrep overlap blocks are correct (block dims sum to nAO; cross-irrep ≈ 0).
Then: 2-electron (`F_AO` once, slice per irrep) and the molecular Factory + SCF wiring.

Key types: `Symmetry::AoShell{shellType,center,monomials,norm,offset}`,
`Symmetry::SALCs{rmat_t O, vector<string> irrep, vector<size_t> blockStart}`. The molecular
orbital basis is a single PG `Orbital_IBS` (IS-A `PGData`, public `radials[]/pols[]/ns[]`).
Build a basis without a file: `Orbital_IBS(rvec_t exponents, size_t LMax, const Cluster*)`.

## Goal

Let the existing SCF machinery solve **molecules** with point-group symmetry, the same way
it already solves atoms (`l`/`κ` channels) and solids (Bloch `k`-irreps): with a
**block-diagonal Hamiltonian, one block per irrep**. The molecular wrinkle is that raw
Cartesian Gaussian AOs are *not* eigenfunctions of the point group, so we need a transform
**O** (Symmetry-Adapted Linear Combinations, SALCs) to rotate the AO basis into irrep blocks.

## The key architectural insight: `IrrepBasisSet` is already the seam

The whole SCF/Hamiltonian/Orbitals/Accelerator stack is organized around
`IrrepBasisSet<T>` (IBS) — one per irrep, with H block-diagonal across IBSs. Each IBS
already carries a polymorphic `Symmetry` and exposes every matrix the SCF needs as a mixin:

- 1-electron: `Integrals_Overlap::MakeOverlap()/Overlap()` (+ Kinetic / Nuclear siblings)
- 2-electron: `Orbital_HF_IBS::MakeDirect/MakeExchange` + `AccumulateDirect/AccumulateExchange`

For **atoms** the basis functions already *are* symmetry eigenfunctions (each `l`/`κ`
channel is its own IBS); for **solids** the Bloch functions already are (`BlochQN`). Only
**molecules** need O ≠ identity.

Unifying statement: **every IBS is an irrep; O = identity for atoms/solids, O = SALC blocks
for molecules.** That is the project thesis in one sentence, and the infrastructure was
designed for it from the start.

Therefore the symmetry-adapted molecular basis is a **Decorator** over the raw molecular
IBS — and *nothing downstream changes*. SCFIterator, the accelerators (DIIS/GDM/Ladder),
Hamiltonian and Orbitals already iterate IBSs and build block-diagonal H; presenting each
point-group irrep as a `SymmetryAdapted_IBS` is transparent to them.

## Decisions

1. **Build from scratch (OOD/SOLID)** rather than wrap a foreign library. Reasons:
   - Neither candidate library reaches the long-term goal. **libmsym** (C, manual
     malloc/free) does molecular point groups + SALCs only; **symmol** (Fortran) only
     *detects* the point group / symmetrizes geometry (no SALCs). Neither does space
     groups, magnetic groups, or BZ k-point symmetry — exactly the planned extensions.
   - Symmetry is this codebase's core differentiating abstraction (`Symmetry`/`Irrep`/
     `BlochQN` already unify atoms and solids); outsourcing the molecular case fragments it.
   - The math is well-bounded — a few hundred lines of textbook group theory, smaller than
     the relativistic angular integrals or GDM already in the codebase.
   - **Use libmsym as a test oracle** (offline): generate SALCs for H₂O / NH₃ / benzene and
     assert our O matches. Cheap insurance on the convention details.
2. **Abelian point groups first** — D2h and its subgroups (C1, Ci, Cs, C2, C2h, C2v, D2).
   Essentially every production QC code runs SCF in only these 8 groups; their irreps are
   all 1-D, so SALCs are just ± sign combinations — no degenerate-irrep / partner-function
   bookkeeping. Design the `Group` abstraction to *admit* non-abelian / space / magnetic
   groups, but implement abelian concretely first.
3. **Cartesian Gaussians now, spherical later.** The current `PolarizedGaussian` basis is
   Cartesian (McMurchie–Davidson recursions). Keep spherical-harmonic Gaussians in mind as
   a *second* molecular basis-set option for the future; it changes only the angular-rotation
   step (step 2 below), not the architecture.

## The 5-step pipeline (with refinements)

**1) Point group from `Cluster`.** `Cluster` gives centers + Z, and `GetAtomIndex(r, tol)`
is exactly the "does operation g map this atom onto another of the same Z" test. Detection
is a self-contained geometric algorithm, basis-agnostic ⇒ lives in **Symmetry** as
`DetectPointGroup(const Cluster&, tol)`.

**2) Classify basis functions under group ops** (the one bridge step needing both worlds).
For each operation g, its action on the AO basis is
`R(g) = (permutation of atomic centers) ⊗ (rotation of the angular shell)`.
The center permutation comes from step 1; the angular part is the **real-Cartesian** (later
real-spherical) rotation matrix `D^l(g)`. This is the bug-prone step (d/f ordering,
Cartesian vs spherical conventions). Building blocks already exist (`Wigner3j`/`RelWigner3j`
in the atomic angular code; `D^l` is adjacent). Needs basis knowledge ⇒ lives in
**BasisSet** (or a thin BasisSet↔Symmetry adapter), consuming Symmetry's `Group`.

**3) SALCs → O.** Projection operator `P^Γ = (l_Γ/h) Σ_g χ^Γ(g)* R(g)` applied to the AOs
gives irrep Γ's combinations.
*Correction to the naive picture:* O does **not** generally produce an orthonormal basis,
because Gaussian AOs overlap (S ≠ I). O's real job is to **block-diagonalize** S and every
symmetry-commuting operator into irrep blocks. Within-block orthonormalization (S^{-1/2} /
Löwdin) is **already done per block by `LASolver`** (the same `Transform`/`BackTransform`
GDM uses). So O only block-diagonalizes; orthonormalization stays where it already lives.

**4) Where O lives** — split by responsibility:
- **Symmetry library**: abstract `Group` (ops, multiplication, character table) +
  `PointGroup` + detection. Basis-agnostic; the home for future `SpaceGroup` / `MagneticGroup`.
- **BasisSet**: the rep-builder (step 2) and the `SymmetryAdapted_IBS` decorator that
  stores and applies O. O is *produced* by the bridge (needs the basis) and *consumed* by
  the decorator. Symmetry never learns what a Gaussian is; BasisSet never learns what a
  character table is.

**5) `SymmetryAdaptedBasisSet` + `SymmetryAdapted_IBS`** (the decorator). Performance
refinement that makes it cheap — **do not 4-index-transform the ERIs**:
- 1-electron (S, T, V): transform once, `Oᵀ M O`. Trivial.
- 2-electron: keep the ERI engine + cache in the AO basis **untouched**. Each SCF
  iteration the decorator back-transforms the density to AO (`D_AO = O D_Γ Oᵀ`), calls the
  existing `AccumulateDirect/Exchange` in AO, then transforms the resulting Fock block to
  the irrep basis (`F_Γ = Oᵀ F_AO O`).

This delivers the **block-diagonal eigenproblem** benefit (smaller diagonalizations,
per-irrep DIIS/GDM) immediately, with *zero* changes to the integral engine. Skipping
symmetry-zero ERIs (the "skeleton / petite list" speedup) is a separate, later optimization
on the engine — do **not** couple it to this.

## Suggested staging (each stage validates independently)

1. **[DONE]** `PointGroup` + `DetectPointGroup(points)` in the Symmetry library
   (`src/Symmetry/PointGroup.C`, `Imp/PointGroup.C`; tests `UnitTests/M_PointGroup.C`,
   21 tests in UTMain).  Operates on `SymPoint{species, position}` (decoupled from
   `Cluster`; a `Cluster->SymPoint` adapter is deferred to the wiring step).  Built bottom
   up: `SymOp` + `IsSymmetryOf` + `Centroid` (1a); `InertiaTensor` + dependency-free Jacobi
   3x3 eigensolver + `ClassifyTop` (1b-1); `FindRotationAxes` (1b-2a); `FindMirrorPlanes` /
   `HasInversion` / `FindImproperAxes` (1b-2b); `DetectPointGroup` -> Schoenflies symbol +
   abelian-subgroup descent (1b-2c).  Validated: H₂O C2v, NH₃ C3v, benzene D6h, CH₄ Td,
   CO₂ D∞h, rectangle D2h, scalene-planar Cs (symbols, orders, abelian subgroups).
2. **[DONE]** Rep-builder, in the Symmetry library (`CartesianRep.C`, tests
   `M_CartesianRep.C`).  `CartesianShellRep(R, exps)` (one shell's monomial transform) and
   `BuildOperationRep(shells, R, origin, tol)` (full AO basis = center permutation ⊗ shell
   rep ⊗ normalization), on a generic `AoShell` layout.  Validated: p-shell rep == R,
   identity, the faithful-representation law `D(R1)D(R2)=D(R1 R2)`, and the full-basis law on
   H₂O.  (PG `AoShell` extraction from the real molecular basis is the wiring step, stage 5.)
3. **[DONE]** SALCs → O, **abelian first**, in the Symmetry library.  `AbelianCharacterTable`
   (8 abelian tables, Mulliken labels — `CharacterTable.C`); `BuildAbelianGroup` (concrete
   ops aligned to the molecule's axes, incl. cubic→C2v and linear→D2h descents —
   `AbelianGroup.C`); `BuildSALCs` -> the transform O, block-structured by irrep, each column
   carrying its Mulliken label (`SALC.C`).  Validated on H₂O C2v: block dims a1=3,a2=0,
   {b1,b2}={1,2}, and every column is an irrep eigenvector `M(g)v = chi^Gamma(g) v`.
   (Mulliken labels flow to the orbital eigenvalue/occupation tables.)
4. `SymmetryAdapted_IBS` decorator (1-e transform + 2-e density/Fock wrapping) +
   `SymmetryAdaptedBasisSet`.
5. Wire into the molecular `Factory`; run H₂O HF, compare total energy to a reference;
   confirm per-irrep blocks shrink the diagonalizations.

## Future extensions (design `Group` to admit these)

- Non-abelian point groups (2-D / 3-D irreps; degenerate-partner bookkeeping).
- Spherical-harmonic Gaussian molecular basis (second basis option; changes only step 2).
- Space groups (solids) and magnetic space groups.
- Brillouin-zone k-point symmetry (likely a facet of space groups).
- ERI skeleton / petite-list speedup (engine-level, exploits symmetry-zero integrals).

## Relevant existing code (entry points)

- `src/Cluster/Cluster.C` — `Cluster` interface (centers, Z, `GetAtomIndex`).
- `src/Symmetry/` — `Symmetry`, `Irrep`, `BlochQN` (where `Group`/`PointGroup` go).
- `src/BasisSet/IrrepBasisSet.C` — the IBS interface + 1-e integral mixins.
- `src/BasisSet/Orbital_HF_IBS.C` — Direct/Exchange + `Accumulate*` (2-e Fock build).
- `src/BasisSet/Molecule/PolarizedGaussian/` — the existing Cartesian molecular basis
  (`Block`, `Polarization`, `IntegralEngine`, McMurchie–Davidson).
