# Structure members by caller (V-STRUCTURE-ISP step 1)

RECORD (lane S, 2026-11-03, on main b1c8d632).  What the 120-odd dependents of `qchem.Structure` actually CALL, so the step-2
faces are cut from use, not from the declarations.  Open work stays in the V-STRUCTURE-ISP row (doc/OOD-SOLID-Cleanup.md).

## Method
- `scripts/structureprobe` runs `clang-query` on each translation unit's REAL compile command: the `.modmap` compiles from the
  `ninja -t compdb` of a built tree.  So overloads and inherited calls resolve as the compiler resolves them.
- It covers the 282 TUs outside `src/Structure` that name the family.
- It records three kinds of use:
  - every call of a member function DECLARED in `src/Structure`, credited to the declaring class (a UnitCell calling an inherited Structure member is a Structure use);
  - every read of a data member declared there;
  - every call of a free function declared there.
- Only call sites in the probed file itself count.
- One TU was not probed (`src/BasisSet/Gaussian/bench/MnDBench.C`: not built).
- `scripts/structuretable` renders the rows.  The raw rows are `StructureUses.tsv`.
- `src/Structure` also holds Lattice_3D, the imposition, the k mesh and Ewald.  Their rows are in the tables, but the analysis below is about the SITE-SET family (Structure / Atom / Molecule / UnitCell).
- Regenerate:
  `ninja -C build/Release -t compdb > cc/all.json` (keep the `.modmap` compiles as `cc/compile_commands.json`), then
  `scripts/structureprobe cc src uses.tsv -j4` (about 25 min), then
  `scripts/structuretable uses.tsv --commit SHA`.

## What the measurement says
**Production use is small.**
- PRODUCTION code (not tests, ITs or CLI apps) uses the family at 308 sites in 39 files.
- Of the 131 files that `import qchem.Structure` directly, 63 call nothing on it: they only pass a `Structure*` through a signature.
- 13 modules re-export it.  Two are tier-1 hubs: `qchem.BasisSet` and `qchem.BasisSet.Orbital_1E_IBS`, where `MakeNuclear(const Structure*)` is the signature that justifies it.

**Production site access (the step-2 site-set face):**
| route | files | sites |
|---|---:|---:|
| `ForEachSite(Z, R, spinFlip)` (the V1.19 face) | 10 | 21 |
| the CONCRETE-Atom route: range-for `begin`/`end` + `operator[]` handing out `Atom*` | 10 + 8 | 32 + 14 |
| ...then reading its public fields `itsR` / `itsZ` / `itsCharge` | 14 / 12 / 1 | 31 / 30 / 1 |
| `GetNumAtoms` | 12 (+3 on Molecule) | 21 (+8) |
| `GetNumElectrons` / `GetNuclearCharge` | 3 / 2 | 5 / 2 |
| `CreateIntegrationMesh` (the mesh type is the geometry's) | 5 | 7 |
| `Clone` | 3 | 8 |

(`GetAtom()` itself is protected and has no outside caller.  Its public doors `operator[]` and the iterators are the concrete-`Atom*` route.)

**Production cell geometry (UnitCell only, all in BasisSet and ChargeDensity):**
- `ToCartesian` 11 files / 50 sites; `GetCellVolume` 5 / 9; `CellsInSphere` 3 / 14; `GetCellMatrix` 3 / 13; `MakeReciprocalCell` 3 / 4; `ToFractional` 2 / 10; `GetMaximumCellEdge` 4 / 5; `GetDistance` 2 / 3.
- All of these are periodic questions.  In the user's model they belong with periodicity (the Lattice), not on the site set.

**`isFinite()`: 4 production sites, as traced in the row.**
- `IonIon.C:87` only words its report line.
- `PWTerms_PP.C:53` is the G = 0 alignment.
- `GPW_Evaluator.C:1110` subtracts the V_short cell mean.
- `Hamiltonian/Imp/Factory.C:23` is admission.
- Plus `NuclearRepulsion` inside `src/Structure`, and 3 test sites (`L_PP.C` ×2, `PW_Terms.C`).
- CORRECTION to the row: the Hartree and XC terms do NOT ask it.  `PWTerms_Hartree.C:31` and `PWTerms_XC.C:30` only import `qchem.Structure`, under a stale comment; only `PWTerms_PP` makes the G = 0 choice (`SumFormFactors`: 1 file).

**Test-only builders dominate the raw counts:** `UnitCell::AddAtom` (155 sites), `Molecule::Insert` (121) and `StructureData::GetMolecule` (50).

---

### Structure members by caller (V-STRUCTURE-ISP step 1, clang-query on b1c8d632)

1018 call sites in 115 files across 14 libraries/dirs.  A member DECLARED in a base counts as the base's (a UnitCell calling an inherited Structure member is a Structure use).

### Per member

| class | member | kind | files | libraries | sites | libraries using it |
|---|---|---|---:|---:|---:|---|
| Molecule | `Insert` | method | 30 | 8 | 121 | BasisSet/tests, CLIapps, Calculation/tests, ChargeDensity, Hamiltonian/tests, IntegrationTests, Response/tests, SCFIterator/tests |
| UnitCell | `AddAtom` | method | 26 | 7 | 155 | BasisSet/tests, Calculation/tests, ChargeDensity/tests, Hamiltonian/tests, IntegrationTests, Response/tests, SCFIterator/tests |
| UnitCell | `ToCartesian` | method | 20 | 5 | 112 | BasisSet, BasisSet/tests, ChargeDensity, Hamiltonian/tests, IntegrationTests |
| Lattice_3D | `GetStructure` | method | 19 | 7 | 39 | BasisSet/tests, CLIapps, Calculation, ChargeDensity/tests, Hamiltonian/tests, IntegrationTests, Response/tests |
| Atom | `itsR` | field | 17 | 5 | 34 | BasisSet, Calculation, Calculation/tests, Hamiltonian, IntegrationTests |
| Structure | `CreateIntegrationMesh` | method | 14 | 7 | 16 | BasisSet, BasisSet/tests, Calculation, Hamiltonian, IntegrationTests, Response/tests, SCFIterator |
| Structure | `ForEachSite` | method | 13 | 7 | 29 | BasisSet, Calculation, Calculation/tests, ChargeDensity, Hamiltonian, IntegrationTests, Response/tests |
| Structure | `GetNumAtoms` | method | 13 | 6 | 22 | BasisSet, Calculation, ChargeDensity, Hamiltonian, IntegrationTests, Pseudopotential |
| Atom | `itsZ` | field | 12 | 4 | 30 | BasisSet, Calculation, Hamiltonian, Pseudopotential |
| (free) | `GetMolecule` | free | 11 | 4 | 50 | BasisSet/tests, Calculation/tests, IntegrationTests, SCFIterator/tests |
| Lattice_3D | `Reciprocal` | method | 11 | 4 | 48 | BasisSet, BasisSet/tests, Hamiltonian/tests, IntegrationTests |
| Structure | `begin` | method | 11 | 2 | 17 | BasisSet, IntegrationTests |
| Structure | `end` | method | 11 | 2 | 17 | BasisSet, IntegrationTests |
| ReciprocalLattice | `GetCell` | method | 10 | 3 | 36 | BasisSet, BasisSet/tests, IntegrationTests |
| Structure | `operator[]` | method | 10 | 6 | 16 | BasisSet, Calculation, Calculation/tests, ChargeDensity, Hamiltonian, Pseudopotential |
| UnitCell | `GetCellVolume` | method | 10 | 4 | 24 | BasisSet, BasisSet/tests, Calculation/tests, IntegrationTests |
| ReciprocalLattice | `CoulombKernel` | method | 6 | 3 | 8 | BasisSet, BasisSet/tests, ChargeDensity |
| Structure | `GetNumElectrons` | method | 6 | 5 | 9 | Calculation, Calculation/tests, ChargeDensity, Hamiltonian, IntegrationTests |
| Structure | `isFinite` | method | 6 | 4 | 7 | BasisSet, Hamiltonian, Hamiltonian/tests, IntegrationTests |
| UnitCell | `GetCellMatrix` | method | 6 | 6 | 19 | BasisSet, BasisSet/tests, Calculation, Calculation/tests, ChargeDensity, Response/tests |
| UnitCell | `MakeReciprocalCell` | method | 6 | 3 | 9 | BasisSet, BasisSet/tests, ChargeDensity/tests |
| Molecule | `GetNumAtoms` | method | 5 | 2 | 13 | Calculation, Calculation/tests |
| Structure | `ID` | method | 5 | 4 | 19 | BasisSet, CLIapps, Calculation/tests, IntegrationTests |
| UnitCell | `ToFractional` | method | 5 | 4 | 16 | BasisSet, BasisSet/tests, ChargeDensity, IntegrationTests |
| (free) | `ResolveImposition` | free | 4 | 3 | 6 | BasisSet, BasisSet/tests, Calculation |
| Structure | `GetNuclearCharge` | method | 4 | 3 | 5 | BasisSet, BasisSet/tests, Calculation |
| UnitCell | `GetMaximumCellEdge` | method | 4 | 2 | 5 | BasisSet, Calculation |
| (free) | `GetReciprocalLattice` | free | 3 | 2 | 6 | ChargeDensity, ChargeDensity/tests |
| Lattice_3D | `GetUnitCell` | method | 3 | 2 | 18 | BasisSet, Calculation |
| Structure | `Clone` | method | 3 | 2 | 8 | Calculation, SCFIterator |
| UnitCell | `CellsInSphere` | method | 3 | 1 | 14 | BasisSet |
| (free) | `BravaisCell` | free | 2 | 2 | 3 | BasisSet/tests, IntegrationTests |
| (free) | `EwaldEnergy` | free | 2 | 2 | 2 | Hamiltonian/tests, IntegrationTests |
| Atom | `ID` | method | 2 | 1 | 3 | Calculation/tests |
| Lattice_3D | `GetLimits` | method | 2 | 2 | 5 | BasisSet, Calculation |
| Lattice_3D | `GetSpaceGroup` | method | 2 | 2 | 2 | BasisSet/tests, Calculation |
| Lattice_3D | `MakeKMesh` | method | 2 | 2 | 3 | BasisSet, Calculation |
| Lattice_3D | `ShubnikovOps` | method | 2 | 2 | 5 | Calculation, Response/tests |
| Molecule | `CreateIntegrationMesh` | method | 2 | 2 | 3 | BasisSet/tests, IntegrationTests |
| UnitCell | `Clone` | method | 2 | 1 | 3 | Hamiltonian/tests |
| UnitCell | `CreateIntegrationMesh` | method | 2 | 1 | 3 | IntegrationTests |
| UnitCell | `GetDistance` | method | 2 | 1 | 3 | BasisSet |
| (free) | `AtomInBox` | free | 1 | 1 | 1 | Calculation |
| (free) | `CellFromEntry` | free | 1 | 1 | 1 | Calculation |
| (free) | `CellNames` | free | 1 | 1 | 1 | Calculation |
| (free) | `DimerInBox` | free | 1 | 1 | 1 | Calculation |
| (free) | `GetCell` | free | 1 | 1 | 1 | Calculation |
| (free) | `GetUnitCell` | free | 1 | 1 | 2 | ChargeDensity |
| (free) | `KindOf` | free | 1 | 1 | 1 | Calculation |
| (free) | `MoleculeNames` | free | 1 | 1 | 1 | Calculation/tests |
| (free) | `NuclearRepulsion` | free | 1 | 1 | 1 | Hamiltonian |
| (free) | `isPeriodicCell` | free | 1 | 1 | 1 | ChargeDensity |
| Atom | `CreateIntegrationMesh` | method | 1 | 1 | 1 | IntegrationTests |
| Atom | `itsCharge` | field | 1 | 1 | 1 | Calculation |
| CrystalImposition | `average` | field | 1 | 1 | 8 | BasisSet |
| CrystalImposition | `cell` | field | 1 | 1 | 1 | BasisSet |
| CrystalImposition | `decision` | field | 1 | 1 | 2 | BasisSet |
| CrystalImposition | `from` | field | 1 | 1 | 1 | BasisSet |
| CrystalImposition | `imposition` | field | 1 | 1 | 3 | BasisSet |
| CrystalImposition | `ops` | field | 1 | 1 | 1 | BasisSet |
| ImpositionOptions | `operator==` | method | 1 | 1 | 1 | BasisSet |
| KMesh | `begin` | method | 1 | 1 | 2 | BasisSet |
| KMesh | `end` | method | 1 | 1 | 2 | BasisSet |
| KMesh | `size` | method | 1 | 1 | 1 | Calculation |
| KPoint | `coordinates` | field | 1 | 1 | 9 | BasisSet |
| KPoint | `weight` | field | 1 | 1 | 3 | BasisSet |
| Lattice_3D | `SiteEnvironmentRotations` | method | 1 | 1 | 1 | Calculation |
| Lattice_3D | `SiteRotations` | method | 1 | 1 | 1 | Calculation |
| ReciprocalLattice | `GetGLength` | method | 1 | 1 | 2 | BasisSet |
| ReciprocalLattice | `GetGVectors` | method | 1 | 1 | 1 | BasisSet |
| Structure | `SumFormFactors` | method | 1 | 1 | 1 | Hamiltonian |
| UnitCell | `EmitUniformGridReport` | method | 1 | 1 | 1 | BasisSet |
| XCMeshSource | `MakeXCMesh` | method | 1 | 1 | 1 | BasisSet |

### Per dependent library

| library | files | members used (sites) |
|---|---:|---|
| BasisSet/tests | 30 | `UnitCell::AddAtom` 91, `Molecule::Insert` 55, `UnitCell::ToCartesian` 47, `Lattice_3D::Reciprocal` 29, `UnitCell::GetCellVolume` 9, `Lattice_3D::GetStructure` 9, `Structure::CreateIntegrationMesh` 6, `UnitCell::ToFractional` 5, `ReciprocalLattice::GetCell` 4, `(free)::ResolveImposition` 3, `UnitCell::GetCellMatrix` 3, `ReciprocalLattice::CoulombKernel` 3, `UnitCell::MakeReciprocalCell` 3, `Structure::GetNuclearCharge` 3, `(free)::GetMolecule` 2, `(free)::BravaisCell` 2, `Molecule::CreateIntegrationMesh` 1, `Lattice_3D::GetSpaceGroup` 1 |
| IntegrationTests | 21 | `Molecule::Insert` 36, `(free)::GetMolecule` 29, `UnitCell::AddAtom` 21, `Lattice_3D::Reciprocal` 14, `Lattice_3D::GetStructure` 12, `UnitCell::ToCartesian` 8, `ReciprocalLattice::GetCell` 4, `UnitCell::CreateIntegrationMesh` 3, `UnitCell::GetCellVolume` 3, `Structure::CreateIntegrationMesh` 2, `Molecule::CreateIntegrationMesh` 2, `Structure::isFinite` 2, `Structure::GetNumElectrons` 2, `Atom::CreateIntegrationMesh` 1, `UnitCell::ToFractional` 1, `Structure::begin` 1, `Structure::end` 1, `Structure::ForEachSite` 1, `Atom::itsR` 1, `(free)::BravaisCell` 1, `Structure::ID` 1, `Structure::GetNumAtoms` 1, `(free)::EwaldEnergy` 1 |
| BasisSet | 21 | `UnitCell::ToCartesian` 46, `ReciprocalLattice::GetCell` 28, `Atom::itsZ` 24, `Atom::itsR` 24, `Structure::begin` 16, `Structure::end` 16, `UnitCell::CellsInSphere` 14, `Structure::ID` 12, `UnitCell::GetCellVolume` 9, `KPoint::coordinates` 9, `CrystalImposition::average` 8, `Lattice_3D::GetUnitCell` 7, `UnitCell::ToFractional` 6, `UnitCell::GetMaximumCellEdge` 4, `ReciprocalLattice::CoulombKernel` 4, `UnitCell::GetCellMatrix` 4, `UnitCell::MakeReciprocalCell` 4, `Structure::GetNumAtoms` 3, `Structure::operator[]` 3, `Lattice_3D::GetLimits` 3, `KPoint::weight` 3, `CrystalImposition::imposition` 3, `UnitCell::GetDistance` 3, `Structure::CreateIntegrationMesh` 2, `Lattice_3D::MakeKMesh` 2, `KMesh::begin` 2, `KMesh::end` 2, `CrystalImposition::decision` 2, `(free)::ResolveImposition` 2, `ReciprocalLattice::GetGLength` 2, `Structure::isFinite` 1, `XCMeshSource::MakeXCMesh` 1, `Lattice_3D::Reciprocal` 1, `ImpositionOptions::operator==` 1, `Structure::ForEachSite` 1, `CrystalImposition::ops` 1, `CrystalImposition::from` 1, `CrystalImposition::cell` 1, `UnitCell::EmitUniformGridReport` 1, `ReciprocalLattice::GetGVectors` 1, `Structure::GetNuclearCharge` 1 |
| Calculation/tests | 10 | `(free)::GetMolecule` 18, `Structure::ForEachSite` 6, `UnitCell::AddAtom` 6, `Molecule::GetNumAtoms` 5, `Structure::ID` 4, `Molecule::Insert` 4, `Atom::ID` 3, `UnitCell::GetCellVolume` 3, `UnitCell::GetCellMatrix` 2, `Structure::operator[]` 2, `Atom::itsR` 2, `Structure::GetNumElectrons` 2, `(free)::MoleculeNames` 1 |
| Calculation | 8 | `Structure::ForEachSite` 14, `Lattice_3D::GetUnitCell` 11, `Molecule::GetNumAtoms` 8, `Lattice_3D::GetStructure` 8, `Structure::GetNumAtoms` 6, `Structure::operator[]` 6, `Atom::itsR` 5, `UnitCell::GetCellMatrix` 5, `Structure::Clone` 3, `Lattice_3D::ShubnikovOps` 3, `Atom::itsZ` 2, `Lattice_3D::GetLimits` 2, `(free)::KindOf` 1, `(free)::GetCell` 1, `(free)::CellFromEntry` 1, `(free)::CellNames` 1, `(free)::AtomInBox` 1, `(free)::DimerInBox` 1, `Structure::CreateIntegrationMesh` 1, `Structure::GetNumElectrons` 1, `Structure::GetNuclearCharge` 1, `Atom::itsCharge` 1, `Lattice_3D::SiteRotations` 1, `Lattice_3D::SiteEnvironmentRotations` 1, `Lattice_3D::GetSpaceGroup` 1, `(free)::ResolveImposition` 1, `UnitCell::GetMaximumCellEdge` 1, `KMesh::size` 1, `Lattice_3D::MakeKMesh` 1 |
| Hamiltonian | 8 | `Structure::GetNumAtoms` 9, `Structure::isFinite` 3, `Structure::operator[]` 3, `Atom::itsZ` 3, `Structure::GetNumElectrons` 2, `Atom::itsR` 2, `Structure::CreateIntegrationMesh` 1, `Structure::ForEachSite` 1, `Structure::SumFormFactors` 1, `(free)::NuclearRepulsion` 1 |
| ChargeDensity/tests | 4 | `UnitCell::AddAtom` 6, `Lattice_3D::GetStructure` 4, `(free)::GetReciprocalLattice` 4, `UnitCell::MakeReciprocalCell` 2 |
| Hamiltonian/tests | 4 | `UnitCell::AddAtom` 11, `Molecule::Insert` 8, `UnitCell::ToCartesian` 7, `Lattice_3D::Reciprocal` 4, `UnitCell::Clone` 3, `Structure::isFinite` 1, `(free)::EwaldEnergy` 1, `Lattice_3D::GetStructure` 1 |
| ChargeDensity | 3 | `Structure::ForEachSite` 5, `UnitCell::GetCellMatrix` 4, `UnitCell::ToFractional` 4, `UnitCell::ToCartesian` 4, `(free)::GetReciprocalLattice` 2, `Structure::GetNumElectrons` 2, `Structure::GetNumAtoms` 2, `(free)::GetUnitCell` 2, `(free)::isPeriodicCell` 1, `Molecule::Insert` 1, `Structure::operator[]` 1, `ReciprocalLattice::CoulombKernel` 1 |
| CLIapps | 2 | `Lattice_3D::GetStructure` 4, `Structure::ID` 2, `Molecule::Insert` 1 |
| Pseudopotential | 1 | `Structure::GetNumAtoms` 1, `Structure::operator[]` 1, `Atom::itsZ` 1 |
| Response/tests | 1 | `UnitCell::AddAtom` 14, `Molecule::Insert` 3, `Lattice_3D::ShubnikovOps` 2, `Structure::CreateIntegrationMesh` 1, `UnitCell::GetCellMatrix` 1, `Structure::ForEachSite` 1, `Lattice_3D::GetStructure` 1 |
| SCFIterator | 1 | `Structure::Clone` 5, `Structure::CreateIntegrationMesh` 3 |
| SCFIterator/tests | 1 | `Molecule::Insert` 13, `UnitCell::AddAtom` 6, `(free)::GetMolecule` 1 |

### isFinite() call sites: 7 sites in 6 files

- `IntegrationTests/L_PP.C:153`
- `IntegrationTests/L_PP.C:162`
- `src/BasisSet/Gaussian/Lattice/Imp/GPW_Evaluator.C:1110`
- `src/Hamiltonian/Imp/Factory.C:23`
- `src/Hamiltonian/Internal/Imp/PWTerms_PP.C:53`
- `src/Hamiltonian/Internal/IonIon.C:87`
- `src/Hamiltonian/tests/PW_Terms.C:92`

### GetAtom() call sites (a concrete Atom* handed out): 0 sites in 0 files


### Direct data-member reads (Atom/Molecule fields read from outside): 93 sites in 21 files

- `Atom::itsR`: 34 reads in 17 files
- `Atom::itsZ`: 30 reads in 12 files
- `KPoint::coordinates`: 9 reads in 1 files
- `CrystalImposition::average`: 8 reads in 1 files
- `KPoint::weight`: 3 reads in 1 files
- `CrystalImposition::imposition`: 3 reads in 1 files
- `CrystalImposition::decision`: 2 reads in 1 files
- `CrystalImposition::ops`: 1 reads in 1 files
- `CrystalImposition::from`: 1 reads in 1 files
- `CrystalImposition::cell`: 1 reads in 1 files
- `Atom::itsCharge`: 1 reads in 1 files
