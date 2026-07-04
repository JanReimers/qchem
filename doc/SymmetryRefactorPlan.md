# qcSymmetry refactor — folder/namespace reorg + Angular-math extraction

**Status: IN PROGRESS (branch `symmetry-refactor`).** This is the single consolidated plan; it supersedes:

- `doc/SymmetryConsolidationPlan.md` (SALC pipeline notes + the `GetAoShells` DIP — that part is DONE)
- `doc/AngularMathPlan.md` (the `Monomial`/`SolidHarmonics` → qcMath extraction — folded in as Stages 4–5)
- memory `project_symmetry_naming_cleanup` (the `Symmetry::Symmetry`→`Irrep`/`SpinIrrep` semantic rename —
  still a **separate later pass**; only the `CarriesSpin()` contract fix is pulled forward here)

Goal: get `src/Symmetry` (project-module `qcSymmetry`) as organized as we can *before* molecular
pseudopotentials. Two independent axes:

1. **Structural** — nest the flat `src/Symmetry` into `Atom/ Molecule/ Lattice_3D/` sub-areas, with matching
   nested namespaces (`qchem::Symmetry::Atom::` …) and C++20 module names (`qchem.Symmetry.Atom.*`), and an
   `Imp/` impl folder in every sub-area.
2. **Layering** — push the shared *angular math* (`Monomial`, real solid harmonics) DOWN into `qcMath`
   (`qchem.Math.Angular`), where both qcSymmetry and the basis evaluators already link — zero new lib edges.

The user has no fear of the rename ripple ("the compiler and unit tests are good guards"), so stages are
sequenced for *review clarity*, not to minimise rebuilds.

## Target layout

```
Symmetry/                     namespace qchem::Symmetry::        module qchem.Symmetry
  Symmetry.C  Irrep.C  Orbital_QNs.C  Spin.C  Factory.C  Unit.C   ← shared sym_t spine (root, all public)
  Imp/
  Atom/                       qchem::Symmetry::Atom::            qchem.Symmetry.Atom
    Spherical.C               ← AtomicSymmetry / Spherical / SphericalSpinor QNs (public)
    Imp/
    Internal/                 qchem::Symmetry::Atom::Internal::  qchem.Symmetry.Atom.Internal
      RealSpinorHarmonics.C   ← was Internal/Spherical.C (Ylm + Ωκ specifics); only Factory uses it
      Imp/
  Molecule/                   qchem::Symmetry::Molecule::        qchem.Symmetry.Molecule
    Irrep.C                   ← was MolecularIrrep.C
    SALC.C  OperationRep.C  ShellRep.C  CartesianRep.C  SphericalRep.C
    PointGroup.C  AbelianGroup.C  CharacterTable.C
    Imp/
  Lattice_3D/                 qchem::Symmetry::Lattice_3D::      qchem.Symmetry.Lattice_3D
    BlochQN.C                 ← promoted out of Internal/ (it is cross-lib public)
    Imp/
```

### Why no `Internal/` under `Molecule/` (deviation from the first sketch)

The first sketch filed the reps + `PointGroup`/`AbelianGroup`/`CharacterTable` under `Molecule/Internal/`.
Verified against the code, that is not viable: the SALC pipeline **spans the library boundary**. In
`qcMolecule_BS`, `PG_Cart/Symmetry.C` declares `std::vector<Symmetry::SymPoint> StructureToSymPoints(...)`
and calls `Symmetry::BuildAbelianGroup` / `BuildSALCs`; `IrrepBasisSet.C` uses `AoShell`; the PG evaluators
use `CartesianShellRep`/`SphericalShellRep`. All of it reaches the molecule side through the
`SALC → AbelianGroup → PointGroup` (and `OperationRep → ShellRep`) `export import` chain. Per CLAUDE.md,
`.Internal.` modules must not be re-exported across a library boundary, so **none of these can be
`Internal`** without first narrowing the public API — the "hide them behind a narrower face" design track,
which was explicitly declined for this pass. So everything molecule-side is public and flat in `Molecule/`.

The **only** genuinely-internal module is the atom `Internal/Spherical` (imported solely by `Factory`), which
becomes `Atom/Internal/RealSpinorHarmonics`. `Unit` looks internal (it lived in `Internal/`) but is imported
by four `qcMolecule_BS` files → it is public and moves to the root. `BlochQN` likewise → `Lattice_3D/`.

### `Internal`-folder ↔ `.Internal.`-module rule

A file under any `Internal/` folder gets `.Internal.` in its module name and `::Internal::` in its
namespace; nothing outside qcSymmetry may `import` (or `export import`) it. Tests may cheat.

## Angular-math extraction (was AngularMathPlan.md)

The same angular math is duplicated: `IVec3 = std::array<int,3>` in `CartesianRep.C` vs `Polarization`
`(n,l,m)` in `PG_Cart_MnD`; and the real-solid-harmonic c2s map as `HarmonicC2S` (type, in `SphericalRep.C`)
vs `SphericalShell(l)` (data, in `PG_Spherical_MnD/SolidHarmonics.C`). Both go DOWN to a new qcMath module
`qchem.Math.Angular` — **not** qcSymmetry (that would force the integral evaluators to depend on group
theory, the wrong direction). qcMath is already a `PUBLIC` dep of both consumers → zero new lib edges.

- `Monomial` — the `(n,l,m)` exponents + `operator[]` / lexicographic `operator<` (≡ the old `LMax`-radix
  order, so `map<Polarization,…>` iteration is preserved) / `==`. `IVec3` becomes `Monomial`;
  `Polarization : Math::Monomial` (keeps arithmetic returning `Polarization`, `operator Index3()`, real-space
  eval — those are evaluator-coupled and must NOT move).
- `CartTerm = {Monomial p; double c;}`, `HarmonicC2S = vector<vector<CartTerm>>`, `SphericalShell(int l)`.
  `SphericalShellRep` consumes the qcMath c2s type; the qcSymmetry tests stop hard-coding d-harmonics.

Regression bar is the **whole UTMain energy suite** byte-identical (the integral core is a `Polarization`
consumer), not just the symmetry tests.

## Stages (build + test after each; commit at green)

0. **[DONE]** Consolidated doc (this file). (Delete the three superseded docs/notes once AngularMath lands.)
1. **[DONE]** `CarriesSpin()` contract fix — false "no spin" claim removed from `Symmetry.C`;
   `virtual bool CarriesSpin()` added (false on base, true on `SphericalSpinor`).
2. **[DONE]** `Lattice_3D/` — `BlochQN` moved (+ `Getk` nested to `::Lattice_3D`).
3. **[DONE]** Root tidy — `Unit.C` lifted from `Internal/` to root.
4. **[DONE]** `Atom/` — `Spherical`→`Atom/`; `Internal/Spherical`→`Atom/Internal/SphericalQNs` (renamed for
   clarity). The pry-out helpers `Getl/Getmls/Getκ/Getmjs` deliberately STAY at the `qchem::Symmetry`
   root (ADL from many bare call sites + the `Evaluator::Getl()` member name clash).
5. **[DONE]** `Molecule/` — `MolecularIrrep`→`Molecule::Irrep` (class renamed), `SALC`, the reps,
   `PointGroup`/`AbelianGroup`/`CharacterTable`, all public and flat in `Molecule/`.

   *Reorg complete: 167/167 UTMain green (full run incl. `A_*`). Branch `symmetry-refactor`, 4 commits.*

6. **[DONE]** AngularMath Stage 1 — `Monomial` in new `qchem.Math.Angular`; `IVec3` is now an alias for it;
   `Polarization : qchem::Math::Monomial`. Lexicographic `Monomial::operator<` verified byte-compatible with
   both the old `std::array` order and `Polarization`'s `LMax`-radix order (the sole Polarization-keyed map,
   Hermite `indexCache`, has non-negative keys). **167/167 UTMain green — byte-identical** (incl. the
   A_HF_dfPin 1e-8 self-pins). Commit on branch `symmetry-refactor`.

7. **AngularMath Stage 2 (REMAINING, lower value)** — move `SphericalShell(l)` + `CartTerm` + the c2s type
   into `qchem.Math.Angular`; unify `SphericalRep`'s `HarmonicC2S` (`pair<IVec3,double>`) with the evaluator's
   `vector<CartTerm>`. Note: `CartTerm`'s member becomes `Monomial` (not `Polarization`), so the
   PG_Spherical_MnD site that builds `SphData` from `SphericalShell` needs a `Monomial`→`Polarization` step
   (add a `Polarization(const Monomial&)` ctor). Byte-identical bar = spherical SCF energies. The plan calls
   this the softer duplication — fine to leave for when convenient.

Class names are preserved verbatim in Stages 1–5 (pure move); the `MolecularIrrep`→`Irrep` file/class
rename is the one exception (the sub-namespace makes the prefix redundant). Semantic renames
(`Symmetry::Symmetry`→`Irrep`, `Irrep`→`SpinIrrep`, `Orbital_QNs`→`OrbitalQNs`) stay in the deferred
`project_symmetry_naming_cleanup` pass.
