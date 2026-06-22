# qchem

**One quantum-chemistry engine for atoms, molecules, and solids — for any Hamiltonian,
any basis set, with or without relativity — built to show that object-oriented and SOLID
design, (hopefully) *done well*, scales to hard computational physics instead of getting in the way.**

The thesis of this project is that with the right abstractions, a *single* code base
can span the whole matrix of

```
{ atoms, molecules, solids } × { any Hamiltonian } × { any basis set } × { non-rel, relativistic }
```

without combinatorial duplication — each new basis set, Hamiltonian, or symmetry is added
*once* and composes with everything already there.  Obviously these attributes do not perfectly decouple, atoms use very 
different basis sets than do molecules and solids.  The Hamiltonian requires the basis set to deliver tables of specific types of integrals.  The required integrals are different for HF/DFT/Relatistic/non-relativistic etc.  The goal is to handle this gracefully in the software architecture by replacing a combinatorial explosion of procedural "if" statements with modern virtual dispatch and template meta programming.

Another question the author is interested in addressing: Is AI now ready to develop reliable computational physics software?

The implementation language selected is c++.  julia would also be a very good language for demonstrating these ideas.

---

## What works today

| Axis | Implemented | In progress / planned |
|------|-------------|-----------------------|
| **Systems** | Atoms (full); molecules (Gaussian, partial) | Periodic solids |
| **Hamiltonians** | Hartree–Fock (RHF/UHF), Kohn–Sham DFT (via libxc), Dirac–Hartree–Fock (Dirac–Coulomb) | Post-HF, Breit/Gaunt |
| **Basis sets** | Slater (STO), Gaussian (GTO), B-splines; relativistic RKB large/small-component variants | plane waves (solids) |
| **Relativity** | 4-component Dirac–Coulomb with Restricted Kinetic Balance; one-electron spin–orbit fine structure | `(LL\|SS)` / `(SS\|SS)` ERIs for heavy-atom accuracy |
| **SCF** | DIIS and GDM accelerators, configuration-averaged & spin-polarized | |

Atomic results are validated against published references — non-relativistic
Hartree–Fock–Roothaan totals (Saito B-spline data) to ~10 significant figures, and
relativistic Dirac–Fock totals, spinor energies, and 2p₁/₂–2p₃/₂ fine-structure
splittings across the periodic table.

## The design idea (why it composes)

The combinatorial explosion is tamed by keeping three concerns **orthogonal**, each behind
a small interface, and composing them with C++20 concepts and mixins:

1. **Basis function ("Evaluator")** — *how* a radial/spatial function and its integrals are
   computed (Slater, Gaussian, B-spline, RKB small component, …). Defined by C++20
   `concept`s (`is1E_Evaluator`, `isHF_Evaluator`, `isRKBL_Evaluator`, …) so a basis type
   only has to provide the integral primitives it can support — and the type system enforces
   the rest.

2. **Integral capabilities** — kinetic, nuclear, overlap, two-electron (Coulomb/exchange),
   DFT fitting — are each a **separate mixin** templated on the Evaluator
   (`Integrals_Kinetic<E>`, `Orbital_HF_IBS<E>`, `Orbital_RKBL_IBS<E>`, …). A basis set is
   *assembled* from exactly the capabilities its Hamiltonian needs. Add a capability once;
   every basis that satisfies the concept gets it.

3. **Symmetry** — spherical harmonics for non-relativistic atoms, spherical *spinors*
   (`SphericalSpinor`, labelled by κ) for the Dirac case, molecular point groups, Bloch function wave vectors — all abstracted so the same orbital and SCF machinery drives both. 

`Factory` functions select the concrete combination at run time from JSON, so a calculation
is "pick a system, a Hamiltonian, a basis, and a relativity flag." The relativistic
`Slater_RKB`/`Gaussian_RKB` path, for example, reuses the *same* non-relativistic radial
evaluators, two-electron engines, and SCF — it only adds the small-component evaluator and
the jj-coupled angular coefficients.

This is the SOLID payoff in practice: single-responsibility pieces, open for extension
(new Evaluator or Hamiltonian) but closed for modification, substitutable behind concepts.

## Current Layout

```
src/
  Common/        constants, periodic table, math
  Symmetry/      spherical harmonics & spherical spinors, electron configurations
  Structure/       atoms and molecules (nuclear framework)
  Mesh/          numerical integration grids
  BasisSet/     Evaluators (Slater/Gaussian/BSpline ± RKB), IrrepBasisSets, factories
  Hamiltonian/   HF, DFT (libxc), Dirac terms; energy breakdown
  ChargeDensity/ density matrices, direct/exchange accumulation
  Orbitals/ WaveFunction/ SCFIterator/ SCFAccelerator/   the SCF loop
  Fitting/ LASolver/ Factory/
UnitTests/       GoogleTest suite (build target: UTMain)
gtkapp/          GTK front-end for atoms
doc/             reference data (HF, DFT, DHF energies)
```

## Building

Requires a C++20 compiler with modules support, CMake ≥ 3.31, a Fortran compiler
(for wignerSymbols/libxc/GaussLegendreIntegration), and BLAS/LAPACK. Third-party dependencies are vendored as git
submodules (blaze, libxc, wignerSymbols, json, googletest, tabulate, BSplinebasis).

```sh
git clone --recursive <repo-url> qchem
cd qchem
cmake -S . -B build/Debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build/Debug --target UTMain -j4
# run the tests (from the UnitTests build dir, so the doc/ data resolves):
( cd build/Debug/UnitTests && ./UTMain )
```

## Status

This is an active research/learning code base — an experiment in how far disciplined
design can carry a from-scratch quantum-chemistry engine. Atoms (non-relativistic and
relativistic) are the most mature; molecules are partial; solids are the long-term target.

## License

GPL-3.0 (see [LICENSE](LICENSE)).
