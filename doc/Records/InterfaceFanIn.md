# Interface fan-in: the compile-checked import table (V-INTERFACE-TIERS step 1)

RECORD (lane S, 2026-11-03).  Measured on main 86eb4895.  Produced by `scripts/importprobe` (greedy per-import removal, each trial
re-checked by the TU's own compile command at `-fsyntax-only` against a frozen build snapshot) and rendered by `scripts/importtable`.
Raw probe rows: `InterfaceFanIn-imports.tsv` (file, unit, module, import, verdict, seconds).  To regenerate on a later tree:

    cp -al build/Release <snap>
    scripts/memsafe -H 4G scripts/importprobe <snap> imports.tsv -j 2
    scripts/importtable <snap> imports.tsv InterfaceFanIn.md --stale stale.tsv

Coverage: 667 TUs (593 with imports).  Six baseline failures are snapshot artefacts, not code: the CLIapps (rundeck, scfrun, valgen,
gth2upf, gpwprobe) and the MnD bench had stale or missing modmaps in the snapshot.  Every library and test TU compiled.

## Rulings on it (controller + user, 2026-11-03)

- **STALE imports in implementation units: remove** (413: 142 in module implementation units, 271 in plain TUs -- tests and CLIapps).
  Interface-unit stale removals were a separate step behind a whole-build check (a downstream template can need them).  DONE 2026-11-05 (78743353): of 86 candidates in 63 interface units, 59 went, 22 moved to the implementation unit that uses them, and 5 stayed because the build needed them in the interface (qchem.Iterators in BasisSet.C, qchem.Fitting.FitOperations and qchem.Fitting.Types in DeltaFunctionFitter.C, qchem.Blaze in Radial/IrrepBasisSet.C and Krylov.C).  The single-TU check could not see those 5 and the 22: it does not compile the module's implementation units or instantiate downstream templates.
- **REDUNDANT imports STAY.**  Under the re-export rule (V-INTERFACE-TIERS: CRP + include-what-you-use) a consumer imports what it
  uses directly, so a redundant direct import is the correct style and the RE-EXPORT it duplicates is the suspect.  The follow-up
  is to shrink the re-exports that fail the signature test (the list after the table), not to delete direct imports.
  DONE 2026-11-05 (c0cb59d2) for the non-hub rows.  The list below is a TEXT heuristic and overstates: re-tested name by name, 43 of its 94 rows are real signature uses or the `:Fit_IBS` partition and stay.  Of the rest, 20 changed: 11 re-exports went, and 9 became the narrower module that carries the type (`qchem.BasisSet.Projector3` for `GMap`, `qchem.Response.Reference` and `qchem.LASolver.Krylov` for `qchem.Response`, `OperationRep` and `PointGroup` for `SALC`); `SolidCalculation` gained `qchem.Response.Solver`, which `qchem.Response` no longer carries.  6 are on tier-1 hubs and need a ruling (`qchem.BasisSet.Orbital_DFT_IBS` to `Fit_Types`, `qchem.BasisSet` to `ElectronConfiguration` and `Matrix3D`, `qchem.ChargeDensity` to `Symmetry.Irrep`, `qchem.Fitting.FunctionFitter` to `GMap`, `qchem.WaveFunction` to `ElectronConfiguration`).  9 are deliberate umbrellas (`qchem`, `qchem.Math`, `qchem.CMath`) and 16 are test harnesses (`qchem.Tests.*`): both stay.
- **The cast column is the side-door list.**  `qchem.BasisSet.Lattice.BasisSet`'s 6x `BasisSetImp<dcmplx>` + 3x `GPW_Evaluator`
  casts (abstract to CONCRETE) go with V-BS-RASTER step 1.  The `TOrbitals<double|dcmplx>` casts in `SolidCalculation.C` (6),
  `SolidState.C` (4), `MoleculeCalculation.C` (1) and `PseudoAtomOrbital.C` (1) are NOT side doors: `TOrbitals<T>`
  (`src/Orbitals/Orbitals.C:122`) is an abstract face (pure virtuals, virtual bases; `TOrbitalsImp` is the concrete class), so
  each is an abstract-to-abstract SCALAR SWITCH -- a type-dispatch smell, not a reach into an implementation.
- **Martin's A, as counted here** (user, 2026-11-03): a class is abstract iff it declares a pure virtual FUNCTION (binary per
  class); VALUE STRUCTS (no virtual member, and a `struct` or no member function beyond operators / defaulted members /
  constructors) leave the denominator, and stay in their files.  `pv` beside it is a gradation, NOT Martin's A.
- **Tier 1** (user): the abstract hubs to protect and the volatile concrete hubs to fix toward abstraction -- see the
  V-INTERFACE-TIERS row of `doc/OOD-SOLID-Cleanup.md`.

Heuristic caveat on the no-sig list: a re-exported module PARTITION (`:Fit_IBS`) is how a module assembles its own interface,
not a side door; the text test cannot tell.


Import verdicts (compiler, -fsyntax-only): {'redundant': 1200, 'stale': 497, 'keep': 1229, 'export': 529}.  `stale` = used by no route; `redundant` = reaches the TU only through another import's re-export (still a dependency).
A, the sig mark and cast targets are TEXT heuristics.

A = Martin's abstractness (binary per class: abstract iff it declares a pure virtual function; value structs excluded).  pv = the mean pure-virtual fraction of the non-value classes -- NOT Martin's A.  — = the module declares no non-value class (A and D undefined).

| module | lib | files | ext libs | Ca | Ce | I | A | D | pv | casts | cast targets | re-exports | Get* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `qchem.Types` | qcMath | 213 | 33 | 204 | 3 | 0.01 | — | — | — | 0 |  | qchem.Matrix2D; qchem.Matrix3D; qchem.Vector3D |  |
| `qchem.Structure` | qcStructure | 136 | 20 | 120 | 8 | 0.06 | 0.33 | 0.60 | 0.11 | 0 |  |  | GetAtom, GetAtomIndex, GetNetCharge, GetNuclearCharge, GetNumAtoms, GetNumElectrons |
| `qchem.Blaze` | qcMath | 130 | 30 | 119 | 1 | 0.01 | 0.00 | 0.99 | 0.00 | 0 |  | qchem.Types |  |
| `qchem.BasisSet.Orbital_DFT_IBS` | qcBasisSet | 61 | 16 | 53 | 8 | 0.13 | 1.00 | 0.13 | 0.67 | 0 |  | :Fit_IBS (no-sig?); qchem.BasisSet.Fit_Types (no-sig?); qchem.BasisSet.IrrepBasisSet; qchem.BasisSet.Orbital_1E_IBS; qchem.BasisSet.Projector3 |  |
| `qchem.BasisSet` | qcBasisSet | 60 | 16 | 58 | 14 | 0.19 | 1.00 | 0.19 | 0.56 | 6 | GridResolution×2, Orbital_DFT_IBS<double,dcmplx>×1, Orbital_DFT_IBS<dcmplx>×1, Orbital_DFT_IBS<double>×1, XCMeshSource×1 | qchem.BasisSet.Fit_Types; qchem.BasisSet.GridResolution (no-sig?); qchem.BasisSet.Orbital_1E_IBS; qchem.BasisSet.Orbital_DFT_IBS; qchem.ElectronConfiguration (no-sig?); qchem.Matrix3D (no-sig?); qchem.Streamable (no-sig?); qchem.Structure; qchem.Symmetry; qchem.Symmetry.Imposition | GetIBS, GetImposition, GetIrreps, GetNumFunctions, GetNumIBS, GetRealIBS |
| `qchem.Mesh` | qcMesh | 60 | 19 | 60 | 1 | 0.02 | 0.00 | 0.98 | 0.00 | 0 |  | qchem.Types |  |
| `qchem.Symmetry.Irrep` | qcSymmetry | 51 | 18 | 50 | 4 | 0.07 | 0.00 | 0.93 | 0.00 | 0 |  | qchem.Symmetry; qchem.Symmetry.Spin | GetDegeneracy |
| `qchem.Lattice_3D` | qcStructure | 46 | 13 | 46 | 6 | 0.12 | 0.00 | 0.88 | 0.00 | 0 |  | qchem.KMesh; qchem.ReciprocalLattice; qchem.Structure; qchem.Symmetry.Lattice_3D.SpaceGroup; qchem.UnitCell | GetBasisNumber, GetBasisVector, GetBonds, GetBondsInSphere, GetCellCoord, GetCellVolume, GetCellsInSphere, GetCoordinate, GetDistances, GetGVectors, GetLatticeVolume, GetLimits … |
| `qchem.ChargeDensity` | qcChargeDensity | 44 | 7 | 32 | 6 | 0.16 | 0.81 | 0.03 | 0.55 | 2 | tSpinResolved_CD<T>×1, tDM_CD<T>×1 | qchem.ChargeDensity.FourierDensity; qchem.Fitting.FunctionFitter; qchem.Symmetry.Irrep (no-sig?); qchem.Symmetry.Spin | GetChangeFrom, GetChannel, GetDensityMatrix, GetEMatrix, GetEMatrixR, GetMatrix, GetMaxChangeFrom, GetTotalCharge, GetTotalSpin |
| `qchem.UnitCell` | qcStructure | 44 | 11 | 42 | 12 | 0.22 | 0.00 | 0.78 | 0.00 | 0 |  | qchem.Types | GetCellMatrix, GetCellVolume, GetDistance, GetMaximumCellEdge, GetMinimumCellEdge, GetNumCells |
| `qchem.BasisSet.Orbital_1E_IBS` | qcBasisSet | 43 | 12 | 43 | 3 | 0.07 | 1.00 | 0.07 | 1.00 | 0 |  | qchem.BasisSet.IrrepBasisSet; qchem.Structure |  |
| `qchem.LASolver` | qcLASolver | 38 | 14 | 36 | 2 | 0.05 | 1.00 | 0.05 | 0.88 | 0 |  | qchem.Types | GetOrthoDim |
| `qchem.Matrix3D` | qcMath | 37 | 16 | 36 | 2 | 0.05 | 0.00 | 0.95 | 0.00 | 0 |  | qchem.Vector3D | GetCol, GetRow |
| `qchem.BasisSet.Gaussian.Point.Factory` | qcGaussian_BS | 32 | 10 | 32 | 8 | 0.20 | — | — | — | 0 |  | qchem.BasisSet; qchem.BasisSet.Gaussian.Point.ShellTrim; qchem.Structure |  |
| `qchem.Energy` | qcHamiltonian | 32 | 5 | 22 | 1 | 0.04 | 0.00 | 0.96 | 0.00 | 0 |  |  | GetBandEnergy, GetConstantEnergy, GetElectronicEnergy, GetEnergy, GetKineticEnergy, GetPotentialEnergy, GetTotalEnergy, GetVirial |
| `qchem.SCFIterator` | qchem | 30 | 6 | 30 | 21 | 0.41 | 0.00 | 0.59 | 0.00 | 4 | tDM_CD<T>×2, tAdaptiveMixer<T>×1, OrbitalView<T>×1 | qchem.ChargeDensity.DensityMixer; qchem.ChargeDensity.Seed; qchem.SCFAccelerator; qchem.SCFParams; qchem.WaveFunction | GetEnergy, GetError, GetIterationCount, GetOrbitals, GetQNs, GetTotalEnergy, GetVirial, GetWaveFunction |
| `qchem.BasisSet.Lattice.BasisSet` | qcLattice_BS | 30 | 10 | 30 | 26 | 0.46 | — | — | — | 12 | BasisSetImp<dcmplx>×6, GPW_Evaluator×3, StreamFoldable×1, GaussianSharpness×1, AoShellSource×1 | qchem.BasisSet; qchem.BasisSet.Gaussian.Lattice.GPWTolerances (no-sig?); qchem.BasisSet.Gaussian.Lattice.GPW_IBS; qchem.BasisSet.Gaussian.Lattice.RasterLadder; qchem.BasisSet.Gaussian.Point.ShellTrim; qchem.BasisSet.PlaneWave.Evaluators; qchem.BasisSet.PlaneWave.PlaneWave_IBS; qchem.Lattice_3D; qchem.Lattice_3D.Imposition |  |
| `qchem.BasisSet.G_FieldEvaluator` | qcBasisSet | 29 | 9 | 27 | 4 | 0.13 | 1.00 | 0.13 | 0.55 | 0 |  | qchem.ScalarFunction |  |
| `qchem.SCFParams` | qchem | 29 | 6 | 28 | 0 | 0.00 | — | — | — | 0 |  |  |  |
| `qchem.Fitting.FunctionFitter` | qcFitting | 28 | 7 | 22 | 10 | 0.31 | 0.85 | 0.16 | 0.68 | 2 | G_RasterTransform×1, G_FieldEvaluator×1 | qchem.BasisSet.GMap (no-sig?); qchem.BasisSet.Orbital_DFT_IBS; qchem.Mesh.Integrator; qchem.ScalarFunction | GetRepulsion3C, GetScalarFunction, GetUnconstrainedFit |
| `qchem.ScalarFunction` | qcMath | 28 | 16 | 28 | 1 | 0.03 | 1.00 | 0.03 | 0.50 | 0 |  | qchem.Types |  |
| `qchem.Mesh.Quadrature` | qcMesh | 25 | 13 | 21 | 5 | 0.19 | — | — | — | 0 |  | qchem.Mesh; qchem.Mesh.Radial; qchem.ScalarFunction; qchem.VectorFunction | GetVectorSize |
| `qchem.BasisSet.GMap` | qcBasisSet | 24 | 9 | 24 | 5 | 0.17 | — | — | — | 0 |  | qchem.BasisSet.Projector3; qchem.Matrix3D (no-sig?); qchem.Symmetry.Lattice_3D.SpaceGroup |  |
| `qchem.Pseudopotential.GTH_Potentials` | qcPseudopotential | 23 | 8 | 22 | 2 | 0.08 | — | — | — | 0 |  | qchem.Pseudopotential.LocalPotential; qchem.Pseudopotential.SeparablePotential | GetGTH, GetGTHValences |
| `qchem.Symmetry.Factory` | qcSymmetry | 23 | 13 | 23 | 4 | 0.15 | — | — | — | 0 |  | qchem.Symmetry |  |
| `qchem.Symmetry` | qcSymmetry | 22 | 11 | 22 | 2 | 0.08 | 1.00 | 0.08 | 0.38 | 0 |  | qchem.Types | GetDegeneracy, GetLabel, GetPrincipleOffset, GetWeight |
| `qchem.Hamiltonian` | qcHamiltonian | 21 | 7 | 21 | 9 | 0.30 | 0.75 | 0.05 | 0.36 | 1 | tDM_CD<T>×1 | qchem.ChargeDensity; qchem.ChargeDensity.TransitionDensity; qchem.Energy; qchem.Hamiltonian.HubbardChannels; qchem.Hamiltonian.HubbardEstimator; qchem.Hamiltonian.HubbardUTarget (no-sig?); qchem.Hamiltonian.TransitionFock; qchem.Hamiltonian.Types | GetEMatrix, GetEMatrixR, GetEnergy, GetHubbardChannels, GetHubbardUTarget, GetMatrix, GetSpinGroup, GetTotalEnergy |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF` | qcGaussian_BS | 21 | 1 | 15 | 9 | 0.38 | 0.00 | 0.62 | 0.00 | 0 |  |  | GetCenter, GetCharge, GetCoeffs, GetExponent, GetExponents, GetGData, GetH1, GetH3, GetID, GetL, GetNMLs, GetNormalization |
| `qchem.Symmetry.Spin` | qcSymmetry | 21 | 7 | 20 | 1 | 0.05 | — | — | — | 0 |  | qchem.Types | GetDegeneracy |
| `qchem.Symmetry.Atom.Spherical` | qcSymmetry | 21 | 7 | 20 | 1 | 0.05 | 1.00 | 0.05 | 0.67 | 0 |  | qchem.Symmetry | GetPrincipleOffset |
| `qchem.ChargeDensity.Seed` | qcChargeDensity | 20 | 5 | 19 | 14 | 0.42 | — | — | — | 0 |  | qchem.ChargeDensity; qchem.ChargeDensity.AtomicDensity |  |
| `qchem.BasisSet.PlaneWave.PlaneWave_IBS` | qcPlaneWave_BS | 19 | 6 | 19 | 20 | 0.51 | 0.00 | 0.49 | 0.00 | 1 | AtomCentredSet_Gaussian×1 | qchem.BasisSet.AtomCentredSet; qchem.BasisSet.Orbital_DFT_IBS; qchem.BasisSet.Orbital_PP_IBS; qchem.ReciprocalLattice |  |
| `qchem.WaveFunction` | qcWaveFunction | 19 | 4 | 18 | 7 | 0.28 | 1.00 | 0.28 | 0.88 | 0 |  | qchem.ChargeDensity; qchem.ElectronConfiguration (no-sig?); qchem.EnergyLevel; qchem.Hamiltonian; qchem.Orbitals; qchem.Symmetry.Irrep | GetChargeDensity, GetEnergyLevels, GetOrbitals, GetQNs, GetSpinDensity, GetSpinGroup |
| `qchem.ReciprocalLattice` | qcStructure | 19 | 5 | 18 | 4 | 0.18 | 0.00 | 0.82 | 0.00 | 2 | UnitCell×2 | qchem.UnitCell | GetCell, GetDistance, GetGLength, GetGVectors, GetReciprocalCell, GetReciprocalLattice, GetUnitCell |
| `qchem.ElectronConfiguration` | qcElConfig | 19 | 6 | 19 | 1 | 0.05 | 1.00 | 0.05 | 0.50 | 0 |  | qchem.Symmetry.Irrep | GetIrreps, GetN, GetPartition |
| `qchem.ChargeDensity.DensitySampler` | qcChargeDensity | 18 | 5 | 18 | 9 | 0.33 | 1.00 | 0.33 | 0.81 | 1 | G_RasterTransform×1 | qchem.ChargeDensity.TransitionDensity; qchem.Mesh.Integrator |  |
| `qchem.Symmetry.Lattice_3D.SpaceGroup` | qcSymmetry | 18 | 10 | 17 | 3 | 0.15 | 0.00 | 0.85 | 0.00 | 0 |  | qchem.Matrix3D; qchem.Symmetry.Lattice_3D.Fold; qchem.Types |  |
| `qchem.SCFAccelerator.Factory` | qcSCFAccelerator | 18 | 5 | 17 | 6 | 0.26 | — | — | — | 0 |  | qchem.SCFAccelerator |  |
| `qchem.SolidCalculation` | qchem | 17 | 5 | 16 | 35 | 0.69 | 0.00 | 0.31 | 0.00 | 7 | TOrbitals<dcmplx>×3, Orbital_DFT_IBS<U,dcmplx>×2, TOrbitals<double>×2 | qchem.Calculation (no-sig?); qchem.Hamiltonian.HubbardEstimator; qchem.PotentialEnergySurface; qchem.Response; qchem.SCFFailure (no-sig?); qchem.SolidSetup; qchem.SolidState |  |
| `qchem.MoleculeCalculation` | qchem | 17 | 5 | 17 | 27 | 0.61 | 0.00 | 0.39 | 0.00 | 3 | sf_t×2, TOrbitals<double>×1 | qchem.Calculation (no-sig?); qchem.MoleculeSetup; qchem.PotentialEnergySurface; qchem.Response (no-sig?); qchem.SCFFailure (no-sig?) | GetLevels, GetSpinDensity, GetStructure |
| `qchem.ChargeDensity.FourierDensity` | qcChargeDensity | 17 | 2 | 14 | 3 | 0.18 | 1.00 | 0.18 | 0.33 | 0 |  | qchem.BasisSet.GMap (no-sig?) | GetFourierDensity, GetRepulsion3C, GetRepulsion3C_Raw, GetRhoOnGrid |
| `qchem.Vector3D` | qcMath | 17 | 7 | 16 | 1 | 0.06 | — | — | — | 0 |  | qchem.CMath |  |
| `qchem.Orbitals` | qcOrbitals | 16 | 5 | 15 | 8 | 0.35 | 1.00 | 0.35 | 0.91 | 0 |  | qchem.ChargeDensity; qchem.ElectronConfiguration.OrbitalView; qchem.Symmetry.Orbital | GetBasisPopulations, GetBasisSet, GetChargeDensity, GetChemicalPotential, GetCoeff, GetCoeffPrime, GetDegeneracy, GetDensityMatrix, GetEigenEnergy, GetEigenValueChange, GetLabel, GetNumOccOrbitals … |
| `qchem.Mesh.XCPolicy` | qcMesh | 16 | 7 | 15 | 2 | 0.12 | — | — | — | 0 |  | qchem.Mesh; qchem.Mesh.Angular |  |
| `qchem.Hamiltonian.Factory` | qcHamiltonian | 15 | 5 | 14 | 13 | 0.48 | — | — | — | 0 |  | qchem.Hamiltonian; qchem.Math.Gaussian (no-sig?); qchem.Outcome (no-sig?); qchem.Pseudopotential.Set; qchem.Symmetry.Spin |  |
| `qchem.BasisSet.PlaneWave.Evaluators` | qcPlaneWave_BS | 15 | 7 | 15 | 10 | 0.40 | 0.00 | 0.60 | 0.00 | 0 |  | qchem.BasisSet.GMap (no-sig?); qchem.BasisSet.G_FieldEvaluator; qchem.ReciprocalLattice | GetGCartesian, GetGIndex |
| `qchem.Hamiltonian.Types` | qcHamiltonian | 15 | 1 | 9 | 6 | 0.40 | — | — | — | 0 |  | qchem.BasisSet; qchem.BasisSet.Orbital_1E_IBS; qchem.BasisSet.Orbital_DFT_IBS; qchem.BasisSet.Orbital_DHF_IBS; qchem.BasisSet.Orbital_HF_IBS |  |
| `qchem.StructureData` | qcStructure | 15 | 7 | 14 | 5 | 0.26 | — | — | — | 0 |  | qchem.Structure; qchem.UnitCell | GetCell, GetMolecule |
| `qchem.VectorFunction` | qcMath | 15 | 8 | 14 | 1 | 0.07 | 1.00 | 0.07 | 0.60 | 0 |  | qchem.Types | GetVectorSize |
| `qchem.BasisSet.Radial.Evaluators` | qcRadial_BS | 14 | 1 | 13 | 5 | 0.28 | 1.00 | 0.28 | 0.92 | 0 |  | qchem.Streamable (no-sig?); qchem.VectorFunction | GetVectorSize |
| `qchem.Materials` | qchem | 14 | 4 | 13 | 7 | 0.35 | — | — | — | 0 |  | qchem.Structure; qchem.StructureData (no-sig?); qchem.UnitCell |  |
| `qchem.SCFAccelerator` | qcSCFAccelerator | 14 | 3 | 10 | 4 | 0.29 | 1.00 | 0.29 | 0.36 | 0 |  | qchem.HandOff (no-sig?); qchem.LASolver; qchem.Symmetry.Irrep | GetError |
| `qchem.BasisSet.Projector3` | qcBasisSet | 14 | 5 | 13 | 3 | 0.19 | 0.00 | 0.81 | 0.00 | 0 |  | qchem.Mesh.Integrator; qchem.Types |  |
| `qchem.BasisSet.Gaussian.Lattice.GPW_IBS` | qcGaussian_BS | 13 | 4 | 12 | 24 | 0.67 | 0.00 | 0.33 | 0.00 | 12 | PW_Grid_Evaluator×6, SpeciesProjectorSet_R×2, AtomCentredSet_Gaussian×2, AoShellSource×1, BareCoulombSource×1 | qchem.BasisSet; qchem.BasisSet.AoShellSource; qchem.BasisSet.AtomCentredSet; qchem.BasisSet.BareCoulombSource; qchem.BasisSet.OrbitalImageSource; qchem.BasisSet.Orbital_DFT_IBS; qchem.BasisSet.Orbital_PP_IBS; qchem.BasisSet.Transition_DFT_IBS; qchem.UnitCell | GetAoShells |
| `qchem.AtomCalculation` | qchem | 13 | 3 | 13 | 24 | 0.65 | 0.00 | 0.35 | 0.00 | 1 | sf_t×1 | qchem.Calculation (no-sig?); qchem.Symmetry.Spin | GetIrreps, GetOrbital, GetStructure, GetZ |
| `qchem.BasisSet.Internal.DB_Cache` | qcBasisSet | 13 | 2 | 13 | 8 | 0.38 | 0.50 | 0.12 | 0.47 | 0 |  | qchem.BasisSet.DBCacheClient (no-sig?) | GetCache2, GetCache3, GetCache4, GetXXX |
| `qchem.Symmetry.Lattice_3D.Fold` | qcSymmetry | 13 | 9 | 13 | 3 | 0.19 | — | — | — | 0 |  | qchem.Matrix3D; qchem.Symmetry.SymOp; qchem.Types |  |
| `qchem.Math` | qcMath | 13 | 7 | 13 | 2 | 0.13 | — | — | — | 0 |  | qchem.CMath (no-sig?); qchem.Math.Gaussian (no-sig?) |  |
| `qchem.BasisSet.Radial.Factory` | qcRadial_BS | 13 | 5 | 13 | 10 | 0.43 | — | — | — | 0 |  | qchem.BasisSet; qchem.ElectronConfiguration |  |
| `qchem.BasisSet.Internal.BasisSetImp` | qcBasisSet | 12 | 4 | 11 | 2 | 0.15 | 0.00 | 0.85 | 0.00 | 0 |  | qchem.BasisSet; qchem.Types | GetChild, GetIBS, GetImposition, GetIrrep, GetIrreps, GetNumFunctions, GetNumIBS, GetRealIBS |
| `qchem.BasisSet.Internal.IrrepBasisSetImp` | qcBasisSet | 12 | 4 | 12 | 2 | 0.14 | 0.00 | 0.86 | 0.00 | 0 |  |  | GetIrrep, GetSymmetry, GetSymt |
| `qchem.BasisSet.Gaussian.PG_Cart` | qcGaussian_BS | 11 | 2 | 11 | 22 | 0.67 | 0.00 | 0.33 | 0.00 | 5 | PGData×3, ToleranceBinding×1, Orbital_HF_IBS<double>×1 |  | GetAoShells, GetNumFunctions |
| `qchem.BasisSet.Gaussian.Lattice.LatticeSum1E` | qcGaussian_BS | 11 | 5 | 11 | 8 | 0.42 | 0.71 | 0.14 | 0.49 | 0 |  | qchem.BasisSet.Gaussian.Lattice.GPWTolerances (no-sig?); qchem.BasisSet.Gaussian.Lattice.LatticeScreener; qchem.BasisSet.GridResolution (no-sig?); qchem.Symmetry.Lattice_3D.SpaceGroup | GetNumFunctions |
| `qchem.ChargeDensity.SeedCD` | qcChargeDensity | 11 | 3 | 11 | 11 | 0.50 | 0.00 | 0.50 | 0.00 | 1 | G_StructureFactor×1 | qchem.ChargeDensity; qchem.ChargeDensity.AtomicDensity; qchem.ChargeDensity.FourierDensity | GetChannel, GetFourierDensity, GetRepulsion3C, GetTotalCharge, GetTotalSpin |
| `qchem.Ewald` | qcStructure | 11 | 7 | 11 | 3 | 0.21 | — | — | — | 1 | UnitCell×1 |  | GetCellVolume, GetNumAtoms |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization` | qcGaussian_BS | 11 | 1 | 10 | 4 | 0.29 | 0.00 | 0.71 | 0.00 | 0 |  |  | GetMaximumL, GetSign, GetTotalL |
| `qchem.BasisSet.Internal.ERI4` | qcBasisSet | 11 | 3 | 10 | 3 | 0.23 | 0.00 | 0.77 | 0.00 | 0 |  | qchem.Types |  |
| `qchem.BasisSet.Gaussian.Point.Reader` | qcGaussian_BS | 11 | 1 | 6 | 2 | 0.25 | 1.00 | 0.25 | 0.75 | 0 |  | qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF; qchem.Structure | GetLs |
| `qchem.CompositeCD` | qcChargeDensity | 10 | 5 | 10 | 8 | 0.44 | 0.00 | 0.56 | 0.00 | 9 | tHF_Pair_CD<T>×1, tStatic_CC<double>×1, Dynamic_CC_RealBlock×1, RealBlockEnergy_CD×1, CoulombMetric_ProjectedDensity×1, tComposite_CD<T>×1 | qchem.ChargeDensity; qchem.ChargeDensity.FourierDensity; qchem.Symmetry.Irrep; qchem.Symmetry.Lattice_3D.SpaceGroup | GetChangeFrom, GetChannel, GetFourierDensity, GetMaxChangeFrom, GetRepulsion3C, GetRepulsion3C_Raw, GetRhoOnGrid, GetTotalCharge |
| `qchem.Symmetry.Lattice_3D.BlochQN` | qcSymmetry | 10 | 7 | 10 | 4 | 0.29 | 0.00 | 0.71 | 0.00 | 4 | BlochQN×4 | qchem.Outcome (no-sig?); qchem.Symmetry; qchem.Symmetry.SelectionRule; qchem.Types | GetDegeneracy, GetPrincipleOffset, GetWeight |
| `qchem.ChargeDensity.Types` | qcChargeDensity | 10 | 1 | 10 | 4 | 0.29 | — | — | — | 0 |  | qchem.BasisSet.Orbital_1E_IBS; qchem.BasisSet.Orbital_DFT_IBS; qchem.BasisSet.Orbital_HF_IBS |  |
| `qchem.BasisSet.Orbital_PP_IBS` | qcBasisSet | 10 | 4 | 7 | 3 | 0.30 | 1.00 | 0.30 | 0.62 | 0 |  |  |  |
| `qchem.Hamiltonian.Internal.ExFunctional` | qcHamiltonian | 10 | 1 | 8 | 2 | 0.20 | 0.50 | 0.30 | 0.08 | 0 |  | qchem.Symmetry.Spin | GetEpsXc, GetExcDensity, GetFxc, GetVxc |
| `qchem.Mesh.Angular` | qcMesh | 10 | 2 | 10 | 3 | 0.23 | 0.00 | 0.77 | 0.00 | 0 |  | qchem.Mesh; qchem.Types |  |
| `qchem.Hamiltonian.Internal.SlaterExchange` | qcHamiltonian | 10 | 3 | 10 | 1 | 0.09 | 0.00 | 0.91 | 0.00 | 0 |  |  | GetFxc, GetVxc |
| `qchem.BasisSet.IrrepBasisSet` | qcBasisSet | 9 | 2 | 9 | 6 | 0.40 | 0.75 | 0.15 | 0.62 | 0 |  | qchem.BasisSet.DBCacheClient (no-sig?); qchem.Symmetry.Irrep; qchem.Symmetry.Representation; qchem.VectorFunction | GetID, GetIrrep, GetNumFunctions, GetSymmetry, GetSymt, GetVectorSize |
| `qchem.Symmetry.SymOp` | qcSymmetry | 9 | 3 | 9 | 2 | 0.18 | — | — | — | 0 |  | qchem.Matrix3D; qchem.Types |  |
| `qchem.Symmetry.Molecule.OperationRep` | qcSymmetry | 9 | 5 | 8 | 2 | 0.20 | — | — | — | 0 |  | qchem.Symmetry.Molecule.ShellRep |  |
| `qchem.BasisSet.Internal.Orbital_ERI4_IBS` | qcBasisSet | 9 | 4 | 9 | 3 | 0.25 | 1.00 | 0.25 | 0.33 | 1 | Orbital_ERI4_IBS<T>×1 | qchem.BasisSet.Internal.ERI4; qchem.BasisSet.Orbital_HF_IBS |  |
| `qchem.Pseudopotential.Set` | qcPseudopotential | 9 | 5 | 9 | 2 | 0.18 | 0.00 | 0.82 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.AoShellSource` | qcBasisSet | 9 | 5 | 9 | 1 | 0.10 | 1.00 | 0.10 | 0.50 | 0 |  | qchem.Symmetry.Molecule.OperationRep | GetAoShells |
| `qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView` | qcGaussian_BS | 9 | 5 | 9 | 13 | 0.59 | — | — | — | 4 | BareCoulombSource×1, AoShellSource×1, TransitionCollocation×1, Periodic_Gaussian_IBS×1 |  |  |
| `qchem.Hamiltonian.Internal.VWN_Correlation` | qcHamiltonian | 8 | 3 | 8 | 2 | 0.20 | 0.00 | 0.80 | 0.00 | 0 |  |  | GetEpsXc, GetVxc |
| `qchem.BasisSet.BareCoulombSource` | qcBasisSet | 8 | 3 | 7 | 1 | 0.12 | 0.50 | 0.38 | 0.25 | 0 |  |  |  |
| `qchem.Deck` | qchem | 7 | 2 | 7 | 20 | 0.74 | — | — | — | 0 |  |  |  |
| `qchem.ChargeDensity.TransitionDensity` | qcChargeDensity | 7 | 2 | 6 | 6 | 0.50 | 1.00 | 0.50 | 0.60 | 3 | WaveVectorShift×1, TransitionFourierDensity×1, ProjectableTransition×1 | qchem.ChargeDensity.Types; qchem.Symmetry.Irrep; qchem.Symmetry.SelectionRule; qchem.Symmetry.Spin | GetMatrix, GetTransitionRepulsion |
| `qchem.ChargeDensity.DensityMixer` | qcChargeDensity | 7 | 2 | 7 | 15 | 0.68 | 1.00 | 0.68 | 0.53 | 4 | FourierDensity×4 | qchem.ChargeDensity; qchem.HandOff (no-sig?) | GetRelax |
| `qchem.ChargeDensity.Imp.IrrepCD` | qcChargeDensity | 7 | 3 | 7 | 8 | 0.53 | 0.00 | 0.47 | 0.00 | 17 | rohfbs_t×6, IrrepCD_Core<T>×4, scalar_t,dcmplx>×3, Orbital_DFT_IBS<T,dcmplx>×2, WholeSystemFock_IBS<double>×1, Orbital_DFT_IBS<double,double>×1 | qchem.ChargeDensity; qchem.ChargeDensity.FourierDensity; qchem.Symmetry.Irrep | GetChangeFrom, GetFourierDensity, GetMaxChangeFrom, GetRepulsion3C, GetRhoOnGrid, GetTotalCharge |
| `qchem.BasisSet.AtomCentredSet` | qcBasisSet | 7 | 5 | 6 | 4 | 0.40 | 0.75 | 0.15 | 0.48 | 0 |  | qchem.BasisSet.BareCoulombSource; qchem.Math.Gaussian (no-sig?) |  |
| `qchem.ChargeDensity.Factory` | qcChargeDensity | 7 | 4 | 7 | 7 | 0.50 | — | — | — | 1 | Orbital_DFT_IBS<T,dcmplx>×1 | qchem.ChargeDensity; qchem.FittedCD |  |
| `qchem.ElectronConfiguration.Crystal` | qcElConfig | 7 | 3 | 7 | 1 | 0.12 | 0.00 | 0.88 | 0.00 | 0 |  | qchem.ElectronConfiguration | GetIrreps, GetN, GetPartition |
| `qchem.Tests.GPW_Harness` | qcGPW_Harness | 6 | 1 | 6 | 26 | 0.81 | — | — | — | 2 | Dynamic_HT_RealBlock×1, G_RasterTransform×1 | qchem.BasisSet; qchem.BasisSet.Gaussian.Point.Factory; qchem.BasisSet.Lattice.BasisSet; qchem.ChargeDensity.DensitySampler; qchem.Energy; qchem.Lattice_3D; qchem.Materials; qchem.Mesh; qchem.Mesh.XCPolicy; qchem.Outcome (no-sig?); qchem.Reporting (no-sig?); qchem.SCFIterator; qchem.SCFParams; qchem.SolidCalculation; qchem.Types | GetEnergy, GetMatrix, GetNumIBS, GetRealIBS, GetTotalEnergy |
| `qchem.Hamiltonian.Internal.PWTerms` | qcHamiltonian | 6 | 4 | 6 | 23 | 0.79 | 0.00 | 0.21 | 0.00 | 12 | Orbital_PP_IBS<U>×3, FourierDensity×2, Transition_DFT_IBS×2, Orbital_DFT_IBS<U,dcmplx>×1, G_RasterTransform×1, G_PoissonKernel×1 |  | GetEnergy, GetMatrix |
| `qchem.BasisSet.Gaussian.Point.IBS` | qcGaussian_BS | 6 | 1 | 6 | 11 | 0.65 | 0.00 | 0.35 | 0.00 | 17 | E×17 | qchem.BasisSet.AoShellSource; qchem.BasisSet.BareCoulombSource | GetAoShells |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData` | qcGaussian_BS | 6 | 1 | 5 | 5 | 0.50 | — | — | — | 0 |  |  |  |
| `qchem.CMath` | qcMath | 6 | 0 | 6 | 3 | 0.33 | — | — | — | 0 |  | qchem.Constants (no-sig?); qchem.Factorials (no-sig?); qchem.IntPow (no-sig?) |  |
| `qchem.Symmetry.Molecule.SALC` | qcSymmetry | 6 | 3 | 6 | 4 | 0.40 | — | — | — | 0 |  | qchem.Symmetry.Molecule.AbelianGroup; qchem.Symmetry.Molecule.OperationRep |  |
| `qchem.ElectronConfiguration.OccupationPolicy` | qcElConfig | 6 | 4 | 6 | 4 | 0.40 | 0.22 | 0.38 | 0.16 | 0 |  | qchem.ElectronConfiguration.OrbitalView; qchem.Symmetry.Irrep |  |
| `qchem.Mesh.Integrator` | qcMesh | 6 | 5 | 6 | 3 | 0.33 | 0.67 | 0.00 | 0.42 | 0 |  | qchem.Mesh; qchem.Mesh.Quadrature (no-sig?); qchem.VectorFunction | GetVectorSize |
| `qchem.Symmetry.Unit` | qcSymmetry | 6 | 3 | 6 | 2 | 0.25 | 0.00 | 0.75 | 0.00 | 0 |  | qchem.Symmetry; qchem.Types | GetDegeneracy, GetPrincipleOffset |
| `qchem.Symmetry.Molecule.PointGroup` | qcSymmetry | 6 | 2 | 6 | 2 | 0.25 | 0.00 | 0.75 | 0.00 | 0 |  | qchem.Matrix3D; qchem.Types | GetKind, GetOrder |
| `qchem.BasisSet.Transition_DFT_IBS` | qcBasisSet | 6 | 3 | 5 | 2 | 0.29 | 1.00 | 0.29 | 0.75 | 0 |  | qchem.BasisSet.Orbital_DFT_IBS |  |
| `qchem.ChargeDensity.AtomicDensity` | qcChargeDensity | 6 | 2 | 5 | 2 | 0.29 | 0.00 | 0.71 | 0.00 | 0 |  |  | GetAtomicDensity, GetAtomicSpinPair, GetN |
| `qchem.Symmetry.Molecule.SphericalRep` | qcSymmetry | 6 | 4 | 6 | 1 | 0.14 | 0.00 | 0.86 | 0.00 | 0 |  | qchem.Symmetry.Molecule.CartesianRep |  |
| `qchem.Fitting.Types` | qcFitting | 6 | 0 | 3 | 1 | 0.25 | — | — | — | 0 |  | qchem.BasisSet.Orbital_DFT_IBS |  |
| `qchem.Tests.PW_Harness` | qcPW_Harness | 5 | 1 | 5 | 31 | 0.86 | — | — | — | 0 |  | qchem.BasisSet.GMap (no-sig?); qchem.BasisSet.G_FieldEvaluator (no-sig?); qchem.BasisSet.Lattice.BasisSet; qchem.BasisSet.Orbital_DFT_IBS (no-sig?); qchem.BasisSet.PlaneWave.Evaluators; qchem.BasisSet.PlaneWave.PlaneWave_IBS; qchem.Blaze (no-sig?); qchem.ChargeDensity.DensitySampler; qchem.ChargeDensity.SeedCD (no-sig?); qchem.ElectronConfiguration (no-sig?); qchem.ElectronConfiguration.Crystal; qchem.Energy; qchem.Ewald (no-sig?); qchem.Fitting.FunctionFitter; qchem.Hamiltonian; qchem.LASolver; qchem.Lattice_3D; qchem.Math (no-sig?); qchem.Matrix3D (no-sig?); qchem.Mesh; qchem.Pseudopotential.GTH_Potentials; qchem.SCFIterator; qchem.SCFParams; qchem.ScalarFunction (no-sig?); qchem.Structure; qchem.Symmetry.Irrep; qchem.Tests.PW_Fields; qchem.Types; qchem.WaveFunction | GetChargeDensity, GetEnergy, GetGIndex, GetGTH, GetIrreps, GetIterationCount, GetNumFunctions, GetSpinGroup, GetStructure, GetTotalCharge, GetTotalEnergy, GetWaveFunction |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD` | qcGaussian_BS | 5 | 1 | 5 | 15 | 0.75 | 0.00 | 0.25 | 0.00 | 0 |  | qchem.Symmetry.Lattice_3D.SpaceGroup | GetCellVolume, GetCenter, GetCoeffs, GetExponents, GetMaximumCellEdge, GetTotalL |
| `qchem.BasisSet.DeltaFit_IBS` | qcBasisSet | 5 | 3 | 5 | 11 | 0.69 | 0.00 | 0.31 | 0.00 | 0 |  | qchem.BasisSet.Fit_Types; qchem.BasisSet.Orbital_DFT_IBS; qchem.BasisSet.Transition_DFT_IBS | GetMesh, GetNumFunctions |
| `qchem.ChargeDensity.Internal.FieldMixer` | qcChargeDensity | 5 | 1 | 5 | 8 | 0.62 | 0.67 | 0.28 | 0.56 | 2 | FourierDensity×1, G_SpectralFilter×1 | qchem.ChargeDensity.FourierDensity; qchem.ChargeDensity.FourierMixCD; qchem.ReciprocalLattice | GetFourierDensity, GetRhoOnGrid |
| `qchem.BasisSet.Gaussian.Evaluators` | qcGaussian_BS | 5 | 1 | 5 | 5 | 0.50 | 1.00 | 0.50 | 0.80 | 1 | E×1 | qchem.Streamable (no-sig?) |  |
| `qchem.BasisSet.Gaussian.PG_Cart.Symmetry` | qcGaussian_BS | 5 | 1 | 5 | 7 | 0.58 | — | — | — | 0 |  | qchem.Symmetry.Molecule.SALC (no-sig?) |  |
| `qchem.Hamiltonian.Internal.Hamiltonian` | qcHamiltonian | 5 | 2 | 5 | 10 | 0.67 | 0.00 | 0.33 | 0.00 | 6 | HubbardProjection×1, HubbardChannels×1, HubbardUTarget×1, tResponse_HT<T>×1, Static_HT_RealBlock×1, Dynamic_HT_RealBlock×1 | qchem.Hamiltonian | GetHubbardChannels, GetHubbardUTarget, GetMatrix, GetSpinGroup, GetTotalEnergy |
| `qchem.Hamiltonian.Internal.Term` | qcHamiltonian | 5 | 0 | 5 | 3 | 0.38 | 0.71 | 0.09 | 0.23 | 0 |  | qchem.Hamiltonian | GetEMatrixR, GetEnergy, GetIrrep, GetMatrix, GetNumIBS, GetRealIBS |
| `qchem.ElectronConfiguration.AtomNR` | qcElConfig | 5 | 4 | 5 | 6 | 0.55 | 0.00 | 0.45 | 0.00 | 1 | Spherical_Sym×1 | qchem.ElectronConfiguration | GetIrreps, GetLMax, GetN |
| `qchem.Pseudopotential.SeparablePotential` | qcPseudopotential | 5 | 3 | 5 | 2 | 0.29 | 0.00 | 0.71 | 0.00 | 3 | SeparablePotential_R×2, SeparablePotential_Gaussian×1 | qchem.BasisSet.Orbital_PP_IBS |  |
| `qchem.BasisSet.Radial.Evaluators.Internal.Rk` | qcRadial_BS | 5 | 1 | 5 | 3 | 0.38 | 1.00 | 0.38 | 0.67 | 1 | Evaluator×1 | qchem.BasisSet.Internal.Cache4 (no-sig?); qchem.Types |  |
| `qchem.Tests.Panel` | ITMain | 5 | 0 | 5 | 2 | 0.29 | — | — | — | 0 |  | qchem.Deck; qchem.SCFCap (no-sig?) |  |
| `qchem.BasisSet.Orbital_HF_IBS` | qcBasisSet | 5 | 3 | 5 | 1 | 0.17 | 1.00 | 0.17 | 0.90 | 0 |  | qchem.BasisSet.Orbital_1E_IBS (no-sig?) |  |
| `qchem.BasisSet.Radial.Evaluators.BSpline.Internal.GLQuadrature` | qcRadial_BS | 5 | 1 | 4 | 2 | 0.33 | 0.00 | 0.67 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.Point.Readers.Gaussian94` | qcGaussian_BS | 5 | 1 | 5 | 6 | 0.55 | 0.00 | 0.45 | 0.00 | 0 |  |  | GetLs |
| `qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator` | qcGaussian_BS | 4 | 3 | 4 | 24 | 0.86 | 0.00 | 0.14 | 0.00 | 7 | SpeciesRadialField_Gaussian×3, Periodic_Gaussian_IBS×1, AoShellSource×1, TransitionCollocation×1, SpeciesProjectorSet_Gaussian×1 | qchem.BasisSet.GMap (no-sig?); qchem.BasisSet.Gaussian.Lattice.GPWTolerances (no-sig?); qchem.BasisSet.PlaneWave.Evaluators |  |
| `qchem.BasisSet.Radial.Evaluators.Gaussian.IBS` | qcRadial_BS | 4 | 1 | 4 | 12 | 0.75 | 0.00 | 0.25 | 0.00 | 1 | HF_Evaluator×1 |  |  |
| `qchem.Hamiltonian.Internal.Hamiltonians` | qcHamiltonian | 4 | 1 | 4 | 19 | 0.83 | — | — | — | 0 |  |  | GetGTH |
| `qchem.BasisSet.Lattice_IBS` | qcBasisSet | 4 | 2 | 3 | 6 | 0.67 | 0.00 | 0.33 | 0.00 | 3 | E×3 |  | GetNumFunctions |
| `qchem.WaveFunction.SCF` | qcWaveFunction | 4 | 1 | 4 | 6 | 0.60 | 1.00 | 0.60 | 1.00 | 0 |  | qchem.ChargeDensity; qchem.ElectronConfiguration.OccupationPolicy; qchem.Hamiltonian; qchem.Orbitals; qchem.Symmetry.Irrep; qchem.WaveFunction | GetOrbitals |
| `qchem.Response.Reference` | qcResponse | 4 | 0 | 4 | 7 | 0.64 | 0.00 | 0.36 | 0.00 | 0 |  | qchem.ElectronConfiguration.OccupationPolicy; qchem.Outcome (no-sig?); qchem.Symmetry.Irrep; qchem.Symmetry.Lattice_3D.BlochQN; qchem.Symmetry.SelectionRule; qchem.Types |  |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Block` | qcGaussian_BS | 4 | 0 | 3 | 4 | 0.57 | 0.00 | 0.43 | 0.00 | 0 |  |  |  |
| `qchem.Symmetry.Representation` | qcSymmetry | 4 | 3 | 4 | 3 | 0.43 | 0.60 | 0.03 | 0.20 | 0 |  | qchem.Symmetry.SymOp; qchem.Types | GetMonomial, GetRep |
| `qchem.Symmetry.Lattice_3D.BlochOperationRep` | qcSymmetry | 4 | 4 | 4 | 4 | 0.50 | — | — | — | 0 |  | qchem.Symmetry.Molecule.OperationRep; qchem.Symmetry.SymOp; qchem.Types |  |
| `qchem.Math.DIIS` | qcMath | 4 | 3 | 4 | 1 | 0.20 | — | — | — | 0 |  |  |  |
| `qchem.Symmetry.Orbital` | qcSymmetry | 4 | 2 | 4 | 2 | 0.33 | 0.00 | 0.67 | 0.00 | 0 |  | qchem.Symmetry.Irrep; qchem.Types | GetDegeneracy |
| `qchem.Mesh.Builder` | qcMesh | 4 | 2 | 4 | 2 | 0.33 | 0.00 | 0.67 | 0.00 | 0 |  | qchem.Mesh |  |
| `qchem.Pseudopotential.LocalPotential` | qcPseudopotential | 4 | 2 | 4 | 2 | 0.33 | 0.20 | 0.47 | 0.04 | 1 | SpeciesRadialField_Gaussian×1 | qchem.BasisSet.Orbital_PP_IBS | GetGTH |
| `qchem.BasisSet.Radial.Evaluators.Internal.NR_Angular` | qcRadial_BS | 4 | 0 | 4 | 4 | 0.50 | 0.00 | 0.50 | 0.00 | 2 | NR_Angular×2 | qchem.BasisSet.Radial.Evaluators |  |
| `qchem.BasisSet.PlaneWave.Internal.GVectors` | qcPlaneWave_BS | 4 | 1 | 4 | 2 | 0.33 | — | — | — | 0 |  |  | GetCell, GetDistance, GetGVectors |
| `qchem.FittedCD` | qcChargeDensity | 4 | 1 | 3 | 2 | 0.40 | 1.00 | 0.40 | 1.00 | 0 |  |  | GetRepulsion, GetSelfRepulsion |
| `qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS` | qcSCFAccelerator | 4 | 3 | 4 | 5 | 0.56 | 0.33 | 0.11 | 0.28 | 0 |  | qchem.SCFAccelerator | GetError, GetNProj, GetNproj |
| `qchem.WaveFunction.Types` | qcWaveFunction | 4 | 0 | 3 | 2 | 0.40 | — | — | — | 0 |  | qchem.BasisSet; qchem.BasisSet.Orbital_1E_IBS |  |
| `qchem.Symmetry.SelectionRule` | qcSymmetry | 4 | 3 | 4 | 1 | 0.20 | 0.67 | 0.13 | 0.33 | 0 |  | qchem.Symmetry |  |
| `qchem.BasisSet.Radial.Evaluators.Internal.AngularIntegrals` | qcRadial_BS | 4 | 1 | 4 | 2 | 0.33 | — | — | — | 0 |  | qchem.Types (no-sig?) |  |
| `qchem.BasisSet.Radial.Evaluators.Internal.Wigner3j` | qcRadial_BS | 4 | 1 | 4 | 1 | 0.20 | 0.00 | 0.80 | 0.00 | 0 |  |  |  |
| `qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull` | qcSCFAccelerator | 4 | 1 | 4 | 2 | 0.33 | 0.33 | 0.33 | 0.17 | 0 |  | qchem.SCFAccelerator | GetError |
| `qchem.Tests.PW_Fields` | qcPW_TestFields | 3 | 3 | 3 | 19 | 0.86 | 0.00 | 0.14 | 0.00 | 0 |  | qchem.BasisSet.GMap (no-sig?); qchem.BasisSet.G_FieldEvaluator; qchem.BasisSet.Lattice.BasisSet (no-sig?); qchem.BasisSet.Orbital_DFT_IBS (no-sig?); qchem.BasisSet.PlaneWave.Evaluators; qchem.BasisSet.PlaneWave.PlaneWave_IBS; qchem.Blaze (no-sig?); qchem.Ewald (no-sig?); qchem.Fitting.FunctionFitter; qchem.LASolver; qchem.Lattice_3D; qchem.Math (no-sig?); qchem.Matrix3D (no-sig?); qchem.Mesh; qchem.Pseudopotential.GTH_Potentials; qchem.ScalarFunction; qchem.Structure (no-sig?); qchem.Symmetry.Irrep (no-sig?); qchem.Types | GetCell, GetFourierDensity, GetGIndex, GetGTH, GetNumFunctions, GetScalarFunction |
| `qchem.Hamiltonian.Internal.Hubbard` | qcHamiltonian | 3 | 1 | 3 | 24 | 0.89 | 0.25 | 0.14 | 0.22 | 8 | AoShellSource×1, AtomCentredOverlapSource<U>×1, cDM_CD×1, cDM_Sourced_CD×1, Orbital_DFT_IBS<double,dcmplx>×1, Orbital_DFT_IBS<dcmplx,dcmplx>×1 | qchem.Hamiltonian.Factory; qchem.Hamiltonian.HubbardChannels; qchem.Hamiltonian.HubbardUTarget (no-sig?) | GetEnergy, GetMatrix |
| `qchem.Hamiltonian.Internal.Terms` | qcHamiltonian | 3 | 1 | 3 | 19 | 0.86 | 0.10 | 0.04 | 0.05 | 8 | odftbs_t×3, tHF_System_CD<double>×2, FitContraction<double,double>×1, ImplicitAngular_IBS×1, orkbbs_t×1 |  | GetEMatrix, GetEnergy, GetEpsXc, GetMatrix |
| `qchem.BasisSet.Internal.DB_Cache_RAM` | qcBasisSet | 3 | 1 | 3 | 8 | 0.73 | 0.00 | 0.27 | 0.00 | 0 |  | qchem.BasisSet.Internal.DB_Cache | GetCache2, GetCache3, GetCache4 |
| `qchem.Response` | qcResponse | 3 | 2 | 3 | 15 | 0.83 | — | — | — | 9 | tobs_t<dcmplx>×2, tobs_t<T>×2, Orbital_DFT_IBS<U,dcmplx>×2, TOrbitals<dcmplx>×1, TOrbitals<double>×1, OrbitalImageSource×1 | qchem.Hamiltonian.HubbardChannels; qchem.Mesh; qchem.Response.LittleGroup; qchem.Response.OrbitalFrame; qchem.Response.Probe; qchem.Response.Reference; qchem.Response.Solver (no-sig?); qchem.WaveFunction |  |
| `qchem.SolidState` | qchem | 3 | 1 | 3 | 16 | 0.84 | — | — | — | 6 | AoShellSource×2, TOrbitals<double>×2, TOrbitals<dcmplx>×2 | qchem.Outcome (no-sig?) |  |
| `qchem.BasisSet.Internal.Orbital_DHF_IBS` | qcBasisSet | 3 | 2 | 3 | 6 | 0.67 | 0.50 | 0.17 | 0.29 | 7 | Orbital_ERI4_IBS<T>×4, Orbital_RKB_Pair<T>×2, Orbital_RKB_HF_IBS_Imp×1 | qchem.BasisSet.Internal.Orbital_ERI4_IBS; qchem.BasisSet.IrrepBasisSet; qchem.BasisSet.Orbital_1E_IBS; qchem.BasisSet.Orbital_DHF_IBS; qchem.Structure | GetNumFunctions |
| `qchem.BasisSet.Radial.Evaluators.BSpline.IBS` | qcRadial_BS | 3 | 1 | 3 | 4 | 0.57 | 0.00 | 0.43 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Spherical_MnD` | qcGaussian_BS | 3 | 1 | 3 | 6 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  |  | GetCharge |
| `qchem.SymmetrizeMesh` | qcStructure | 3 | 2 | 3 | 4 | 0.57 | 0.00 | 0.43 | 0.00 | 0 |  | qchem.Mesh; qchem.Symmetry.Lattice_3D.Fold; qchem.Symmetry.Lattice_3D.SpaceGroup |  |
| `qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet` | qcGaussian_BS | 3 | 1 | 3 | 6 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.BasisSet; qchem.Structure |  |
| `qchem.ChargeDensity.FourierMixCD` | qcChargeDensity | 3 | 1 | 3 | 6 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.ChargeDensity; qchem.ChargeDensity.FourierDensity; qchem.ReciprocalLattice | GetFourierDensity, GetRepulsion3C, GetRhoOnGrid, GetTotalCharge |
| `qchem.Lattice_3D.Imposition` | qcStructure | 3 | 3 | 3 | 6 | 0.67 | 0.50 | 0.17 | 0.25 | 0 |  | qchem.Lattice_3D; qchem.Symmetry.Imposition; qchem.Symmetry.Lattice_3D.SpaceGroup |  |
| `qchem.BasisSet.Radial.Evaluators.Internal.ExponentialEvaluator` | qcRadial_BS | 3 | 1 | 3 | 5 | 0.62 | 0.00 | 0.38 | 0.00 | 0 |  | qchem.BasisSet.Radial.Evaluators; qchem.Symmetry.Atom.Spherical (no-sig?) |  |
| `qchem.Symmetry.Imposition` | qcSymmetry | 3 | 3 | 3 | 3 | 0.50 | 1.00 | 0.50 | 0.67 | 0 |  | qchem.Symmetry.Lattice_3D.SpaceGroup; qchem.Symmetry.Representation |  |
| `qchem.Mesh.Product` | qcMesh | 3 | 2 | 3 | 3 | 0.50 | — | — | — | 0 |  | qchem.Mesh; qchem.Mesh.Angular; qchem.Mesh.Radial |  |
| `qchem.BasisSet.Orbital_DHF_IBS` | qcBasisSet | 3 | 2 | 3 | 4 | 0.57 | 1.00 | 0.57 | 0.50 | 0 |  | qchem.BasisSet.IrrepBasisSet; qchem.BasisSet.Orbital_1E_IBS; qchem.Structure (no-sig?) |  |
| `qchem.BasisSet.Fit_Types` | qcBasisSet | 3 | 0 | 3 | 3 | 0.50 | — | — | — | 0 |  | qchem.Mesh; qchem.Mesh.Folded; qchem.Symmetry.Lattice_3D.Fold |  |
| `qchem.BasisSet.Internal.Fit_IBS` | qcBasisSet | 3 | 2 | 3 | 7 | 0.70 | 1.00 | 0.70 | 0.50 | 0 |  | qchem.BasisSet.Orbital_DFT_IBS; qchem.Mesh |  |
| `qchem.BasisSet.PlaneWave.Internal.KPlusG` | qcPlaneWave_BS | 3 | 1 | 3 | 3 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.Lattice.LatticeScreener` | qcGaussian_BS | 3 | 1 | 3 | 2 | 0.40 | 0.25 | 0.35 | 0.17 | 0 |  |  |  |
| `qchem.Hamiltonian.HubbardEstimator` | qcHamiltonian | 3 | 1 | 3 | 3 | 0.50 | 1.00 | 0.50 | 0.83 | 0 |  |  |  |
| `qchem.SCFRecipeReport` | qchem | 3 | 0 | 3 | 2 | 0.40 | — | — | — | 0 |  | qchem.Reporting (no-sig?); qchem.SCFParams |  |
| `qchem.FFT` | qcMath | 3 | 2 | 3 | 2 | 0.40 | — | — | — | 0 |  |  |  |
| `qchem.Symmetry.Molecule.CartesianRep` | qcSymmetry | 3 | 2 | 3 | 2 | 0.40 | 0.00 | 0.60 | 0.00 | 0 |  | qchem.Math.Angular (no-sig?); qchem.Symmetry.Molecule.ShellRep |  |
| `qchem.Symmetry.Molecule.AbelianGroup` | qcSymmetry | 3 | 1 | 3 | 2 | 0.40 | — | — | — | 0 |  | qchem.Symmetry.Molecule.CharacterTable (no-sig?); qchem.Symmetry.Molecule.PointGroup |  |
| `qchem.Symmetry.Lattice_3D.BZReduction` | qcSymmetry | 3 | 2 | 3 | 3 | 0.50 | — | — | — | 0 |  | qchem.Matrix3D; qchem.Types |  |
| `qchem.Mesh.Radial` | qcMesh | 3 | 1 | 3 | 2 | 0.40 | 0.00 | 0.60 | 0.00 | 0 |  | qchem.Mesh; qchem.Types |  |
| `qchem.BasisSet.OrbitalImageSource` | qcBasisSet | 3 | 3 | 3 | 2 | 0.40 | 1.00 | 0.40 | 0.50 | 0 |  | qchem.Symmetry.SymOp; qchem.Types |  |
| `qchem.BasisSet.Gaussian.Point.ShellTrim` | qcGaussian_BS | 3 | 2 | 3 | 4 | 0.57 | 0.00 | 0.43 | 0.00 | 0 |  | qchem.BasisSet.Gaussian.Point.Reader | GetLs |
| `qchem.Hamiltonian.HubbardChannels` | qcHamiltonian | 3 | 1 | 3 | 2 | 0.40 | 1.00 | 0.40 | 0.80 | 0 |  |  |  |
| `qchem.Hamiltonian.Internal.Kinetic` | qcHamiltonian | 3 | 2 | 3 | 2 | 0.40 | 0.00 | 0.60 | 0.00 | 0 |  |  | GetEnergy |
| `qchem.EnergyLevel` | qcOrbitals | 3 | 1 | 3 | 2 | 0.40 | 0.00 | 0.60 | 0.00 | 0 |  | qchem.Orbitals; qchem.Symmetry.Orbital |  |
| `qchem.HDF5` | qcCommon | 3 | 3 | 3 | 1 | 0.25 | 0.00 | 0.75 | 0.00 | 0 |  | qchem.Outcome (no-sig?) |  |
| `qchem.SpecialFunctions` | qcMath | 3 | 2 | 3 | 1 | 0.25 | — | — | — | 0 |  |  |  |
| `qchem.ElectronConfiguration.AtomDirac` | qcElConfig | 3 | 3 | 3 | 4 | 0.57 | — | — | — | 0 |  |  |  |
| `qchem.BasisSet.ImplicitAngular_IBS` | qcBasisSet | 3 | 3 | 3 | 1 | 0.25 | 1.00 | 0.25 | 0.67 | 0 |  |  |  |
| `qchem.BasisSet.Radial.Evaluators.Gaussian.Internal.GaussianIntegrals` | qcRadial_BS | 3 | 0 | 1 | 1 | 0.50 | — | — | — | 0 |  |  |  |
| `qchem.BasisSet.Radial.Evaluators.Slater.Internal.Integrals` | qcRadial_BS | 3 | 0 | 1 | 1 | 0.50 | — | — | — | 0 |  |  |  |
| `qchem.Orbitals.Types` | qcOrbitals | 3 | 0 | 3 | 1 | 0.25 | — | — | — | 0 |  | qchem.BasisSet.Orbital_1E_IBS |  |
| `qchem.Response.Probe` | qcResponse | 3 | 1 | 3 | 2 | 0.40 | 0.33 | 0.27 | 0.27 | 0 |  | qchem.Response.Reference |  |
| `qchem.BasisSet.Radial.IBS` | qcRadial_BS | 2 | 1 | 2 | 13 | 0.87 | 0.09 | 0.04 | 0.06 | 15 | E×14, Evaluator×1 | qchem.BasisSet.ImplicitAngular_IBS | GetCache4, GetNuclearCharge, GetNumAtoms, GetNumFunctions, GetSymmetry |
| `qchem.BasisSet.Radial.Evaluators.Slater.IBS` | qcRadial_BS | 2 | 1 | 2 | 12 | 0.86 | 0.00 | 0.14 | 0.00 | 1 | HF_Evaluator×1 |  |  |
| `qchem.BasisSet.Gaussian.Point.PG_Spherical` | qcGaussian_BS | 2 | 1 | 2 | 17 | 0.89 | 0.00 | 0.11 | 0.00 | 1 | NR_Evaluator×1 |  | GetAoShells, GetNumFunctions |
| `qchem.BasisSet.Radial.Evaluators.BSpline.Internal.Common` | qcRadial_BS | 2 | 0 | 2 | 7 | 0.78 | 0.00 | 0.22 | 0.00 | 1 | EvaluatorCommon<K>×1 | qchem.BasisSet.Radial.Evaluators; qchem.Symmetry | GetGrid |
| `qchem.ChargeDensity.Internal.TransitionDensity` | qcChargeDensity | 2 | 1 | 2 | 6 | 0.75 | 0.00 | 0.25 | 0.00 | 2 | tHF_System_CD<double>×1, Transition_DFT_IBS×1 | qchem.ChargeDensity.TransitionDensity | GetChannel, GetNumFunctions, GetTransitionRepulsion |
| `qchem.BasisSet.Internal.SymmetryAdapted_IBS` | qcBasisSet | 2 | 1 | 2 | 7 | 0.78 | 0.00 | 0.22 | 0.00 | 3 | Orbital_HF_IBS <double>×1, Orbital_DFT_IBS<double>×1, WholeSystemFock_IBS<double>×1 | qchem.BasisSet.Orbital_1E_IBS; qchem.BasisSet.Orbital_DFT_IBS; qchem.BasisSet.Orbital_HF_IBS | GetNumFunctions, GetO |
| `qchem.BasisSet.Gaussian.Evaluators.PG_LibCint` | qcGaussian_BS | 2 | 1 | 2 | 8 | 0.80 | 0.00 | 0.20 | 0.00 | 0 |  |  |  |
| `qchem.ChargeDensity.Internal.KerkerMixer` | qcChargeDensity | 2 | 1 | 2 | 4 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.ChargeDensity.DensityMixer; qchem.ChargeDensity.Internal.FieldMixer | GetRelax |
| `qchem.ChargeDensity.Internal.PulayMixer` | qcChargeDensity | 2 | 1 | 2 | 5 | 0.71 | 0.00 | 0.29 | 0.00 | 0 |  | qchem.ChargeDensity.DensityMixer; qchem.ChargeDensity.Internal.FieldMixer | GetRelax |
| `qchem.WaveFunction.Internal.CompositeWF` | qcWaveFunction | 2 | 1 | 2 | 15 | 0.88 | 0.00 | 0.12 | 0.00 | 3 | TOrbitals<U>×2, OrbitalView<T>×1 | qchem.WaveFunction.SCF | GetChargeDensity, GetEnergyLevels, GetOrbitals, GetQNs, GetSpinDensity, GetSpinGroup |
| `qchem.BasisSet.Lattice.LAPW_IBS` | qcLattice_BS | 2 | 1 | 2 | 11 | 0.85 | 0.00 | 0.15 | 0.00 | 0 |  | qchem.ReciprocalLattice | GetNumFunctions |
| `qchem.BasisSet.Lattice.BandStructure` | qcLattice_BS | 2 | 1 | 2 | 5 | 0.71 | — | — | — | 0 |  |  |  |
| `qchem.ChargeDensity.Internal.PolarizedDensityMixer` | qcChargeDensity | 2 | 1 | 2 | 4 | 0.67 | 0.00 | 0.33 | 0.00 | 5 | FourierDensity×2, tFieldMixer×2, cDM_SourceSink×1 | qchem.ChargeDensity.DensityMixer; qchem.ChargeDensity.Internal.FieldMixer | GetChannel, GetFourierDensity, GetRelax, GetRepulsion3C, GetRhoOnGrid, GetTotalCharge |
| `qchem.Hamiltonian.TransitionFock` | qcHamiltonian | 2 | 1 | 2 | 3 | 0.60 | 0.50 | 0.10 | 0.33 | 0 |  | qchem.Symmetry.Irrep; qchem.Symmetry.SelectionRule |  |
| `qchem.ChargeDensity.NumericCD` | qcChargeDensity | 2 | 1 | 2 | 5 | 0.71 | 0.00 | 0.29 | 0.00 | 1 | FIT_SF_NonOrtho×1 | qchem.ChargeDensity | GetTotalCharge, GetUnconstrainedFit |
| `qchem.Hamiltonian.Internal.IonIon` | qcHamiltonian | 2 | 1 | 2 | 4 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  |  | GetEnergy, GetNumAtoms, GetNumFunctions |
| `qchem.Orbitals.Internal.OrbitalsImp` | qcOrbitals | 2 | 0 | 1 | 9 | 0.90 | 0.00 | 0.10 | 0.00 | 1 | TOrbital<T>×1 | qchem.Orbitals; qchem.Symmetry.Irrep; qchem.Types | GetBasisPopulations, GetBasisSet, GetChargeDensity, GetChemicalPotential, GetDegeneracy, GetDensityMatrix, GetEigenValueChange, GetNumOccOrbitals, GetNumOrbitals, GetOccupation, GetOrbital, GetQNs |
| `qchem.PotentialEnergySurface` | qchem | 2 | 0 | 2 | 3 | 0.60 | 1.00 | 0.60 | 0.50 | 0 |  | qchem.Outcome (no-sig?); qchem.SCFFailure (no-sig?); qchem.Types |  |
| `qchem.Mesh.Folded` | qcMesh | 2 | 2 | 2 | 3 | 0.60 | 0.00 | 0.40 | 0.00 | 0 |  | qchem.Mesh; qchem.Symmetry.Lattice_3D.Fold; qchem.Symmetry.Spin (no-sig?) | GetFold, GetMesh, GetSpinOps |
| `qchem.LASolver.Krylov` | qcLASolver | 2 | 2 | 2 | 2 | 0.50 | 1.00 | 0.50 | 0.67 | 0 |  | qchem.Outcome (no-sig?); qchem.Types |  |
| `qchem.BasisSet.Gaussian.Evaluators.Internal.MnD` | qcGaussian_BS | 2 | 0 | 1 | 3 | 0.75 | — | — | — | 0 |  | qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Index3 (no-sig?); qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM; qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D |  |
| `qchem.Response.OrbitalFrame` | qcResponse | 2 | 0 | 2 | 4 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.ChargeDensity.TransitionDensity; qchem.Hamiltonian.TransitionFock; qchem.Response.Reference | GetC |
| `qchem.Symmetry.Molecule.ShellRep` | qcSymmetry | 2 | 0 | 2 | 2 | 0.50 | 1.00 | 0.50 | 0.60 | 0 |  | qchem.Types |  |
| `qchem.Mesh.GaussLegendre` | qcMesh | 2 | 1 | 2 | 1 | 0.33 | — | — | — | 0 |  | qchem.Types |  |
| `qchem.BasisSet.Radial.Evaluators.Internal.RKBL_Angular` | qcRadial_BS | 2 | 0 | 2 | 4 | 0.67 | 0.00 | 0.33 | 0.00 | 2 | RKB_Angular×2 | qchem.BasisSet.Radial.Evaluators |  |
| `qchem.BasisSet.Radial.Internal.FillPower` | qcRadial_BS | 2 | 0 | 2 | 1 | 0.33 | — | — | — | 0 |  |  |  |
| `qchem.BasisSet.Radial.Evaluators.BSpline.Internal.SplineGrouper` | qcRadial_BS | 2 | 0 | 1 | 1 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Hermite` | qcGaussian_BS | 2 | 0 | 1 | 3 | 0.75 | 0.00 | 0.25 | 0.00 | 0 |  |  | GetAny, GetIndex, GetSize |
| `qchem.SCFAccelerator.Internal.SCFAcceleratorGDM` | qcSCFAccelerator | 2 | 1 | 2 | 3 | 0.60 | 0.33 | 0.07 | 0.28 | 0 |  | qchem.SCFAccelerator | GetError |
| `qchem.Response.LittleGroup` | qcResponse | 2 | 0 | 2 | 4 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.Response.Reference; qchem.Symmetry.SymOp |  |
| `qchem.PseudoAtomOrbital` | qchem | 2 | 1 | 2 | 11 | 0.85 | — | — | — | 2 | AtomicSymmetry×1, TOrbitals<double>×1 | qchem.Math.Gaussian (no-sig?) |  |
| `qchem.Symmetry.Atom.Internal.SphericalQNs` | qcSymmetry | 2 | 1 | 2 | 2 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  | qchem.Symmetry.Atom.Spherical | GetDegeneracy |
| `qchem.ElectronConfiguration.OrbitalView` | qcElConfig | 2 | 1 | 2 | 1 | 0.33 | 1.00 | 0.33 | 0.80 | 0 |  |  |  |
| `qchem.ElectronConfiguration.Molecule` | qcElConfig | 2 | 2 | 2 | 2 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  | qchem.ElectronConfiguration | GetIrreps, GetN, GetPartition |
| `qchem.BasisSet.Radial.Evaluators.Internal.RelAngularIntegrals` | qcRadial_BS | 2 | 1 | 2 | 5 | 0.71 | — | — | — | 0 |  | qchem.Types (no-sig?) |  |
| `qchem.BasisSet.Radial.Evaluators.Gaussian.Internal.Rk` | qcRadial_BS | 2 | 0 | 1 | 3 | 0.75 | 0.00 | 0.25 | 0.00 | 0 |  | qchem.BasisSet.Radial.Evaluators.Internal.Rk |  |
| `qchem.BasisSet.Radial.Evaluators.Slater.Internal.Rk` | qcRadial_BS | 2 | 0 | 1 | 3 | 0.75 | 0.00 | 0.25 | 0.00 | 0 |  | qchem.BasisSet.Radial.Evaluators.Internal.Rk |  |
| `qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D` | qcGaussian_BS | 2 | 0 | 2 | 1 | 0.33 | 0.00 | 0.67 | 0.00 | 0 |  |  |  |
| `qchem.Hamiltonian.Internal.Libxc_LDA` | qcHamiltonian | 2 | 1 | 2 | 1 | 0.33 | 0.00 | 0.67 | 0.00 | 0 |  | qchem.Hamiltonian.Internal.ExFunctional | GetEpsXc, GetVxc |
| `qchem.SCFAccelerator.Internal.SCFAcceleratorLadder` | qcSCFAccelerator | 2 | 1 | 2 | 2 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  | qchem.SCFAccelerator | GetError |
| `qchem.WaveFunction.EigenRows` | qcWaveFunction | 2 | 1 | 2 | 1 | 0.33 | — | — | — | 0 |  | qchem.EnergyLevel |  |
| `qchem.SolidSetup` | qchem | 1 | 0 | 1 | 19 | 0.95 | — | — | — | 2 | GaussianSharpness×1, ToleranceBinding×1 | qchem.Pseudopotential.Set | GetStructure |
| `qchem.ChargeDensity.Internal.DensitySampler` | qcChargeDensity | 1 | 0 | 1 | 19 | 0.95 | 0.00 | 0.05 | 0.00 | 25 | tDM_CD<T>×4, Orbital_DFT_IBS<dcmplx,dcmplx>×4, ProjectableTransition×2, cDM_CD×2, G_RasterTransform×2, ScalarProjector×1 | qchem.ChargeDensity.DensitySampler; qchem.Mesh.Integrator | GetScalarFunction |
| `qchem.BasisSet.Radial.BasisSet` | qcRadial_BS | 1 | 0 | 1 | 12 | 0.92 | 0.00 | 0.08 | 0.00 | 0 |  | qchem.BasisSet; qchem.BasisSet.Orbital_HF_IBS (no-sig?); qchem.ElectronConfiguration | GetIrreps, GetNumFunctions, GetSymt |
| `qchem.BasisSet.Gaussian.Point.PG_LibCint` | qcGaussian_BS | 1 | 0 | 1 | 12 | 0.92 | 0.00 | 0.08 | 0.00 | 0 |  |  | GetAoShells, GetNumFunctions |
| `qchem.WaveFunction.Internal.IrrepWF` | qcWaveFunction | 1 | 0 | 1 | 13 | 0.93 | 0.00 | 0.07 | 0.00 | 0 |  | qchem.ChargeDensity; qchem.ElectronConfiguration (no-sig?); qchem.ElectronConfiguration.OccupationPolicy; qchem.EnergyLevel; qchem.Hamiltonian; qchem.LASolver; qchem.Orbitals; qchem.Symmetry.Irrep | GetBasisPopulations, GetChargeDensity, GetIrrep, GetOrbitals |
| `qchem.BasisSet.PlaneWave.PlaneWaveFit_IBS` | qcPlaneWave_BS | 1 | 0 | 1 | 8 | 0.89 | 0.00 | 0.11 | 0.00 | 0 |  | qchem.BasisSet.Orbital_DFT_IBS |  |
| `qchem` | qcUmbrella | 1 | 1 | 1 | 8 | 0.89 | — | — | — | 0 |  | qchem.Materials (no-sig?); qchem.MoleculeCalculation; qchem.PeriodicTable (no-sig?); qchem.SCFIterator; qchem.ScalarFunction; qchem.Structure; qchem.StructureData (no-sig?); qchem.Types (no-sig?) |  |
| `qchem.MoleculeSetup` | qchem | 1 | 0 | 1 | 8 | 0.89 | — | — | — | 0 |  | qchem.Pseudopotential.Set; qchem.Symmetry.Spin |  |
| `qchem.BasisSet.Radial.Evaluators.BSpline.IBS_r` | qcRadial_BS | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.Lattice.RasterLadder` | qcGaussian_BS | 1 | 1 | 1 | 6 | 0.86 | 0.00 | 0.14 | 0.00 | 0 |  | qchem.BasisSet.Gaussian.Lattice.GPWTolerances (no-sig?); qchem.BasisSet.PlaneWave.Evaluators |  |
| `qchem.Fitting.Internal.FunctionFitterImp` | qcFitting | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 4 | CoulombMetric_ProjectedDensity×1, OverlapMetric_ProjectedDensity×1, ProjectedDensity_AO×1, Orbital_DFT_IBS<T,T>×1 | qchem.Fitting.FunctionFitter | GetNumFunctions |
| `qchem.Fitting.Internal.DeltaFunctionFitter` | qcFitting | 1 | 0 | 1 | 4 | 0.80 | 0.00 | 0.20 | 0.00 | 1 | Integrals_Overlap3C<U,T>×1 | qchem.Fitting.FunctionFitter | GetNumFunctions, GetScalarFunction, GetSymt |
| `qchem.WaveFunction.Factory` | qcWaveFunction | 1 | 1 | 1 | 7 | 0.88 | — | — | — | 0 |  | qchem.ElectronConfiguration; qchem.Hamiltonian; qchem.SCFAccelerator; qchem.WaveFunction.SCF |  |
| `qchem.Response.Solver` | qcResponse | 1 | 0 | 1 | 7 | 0.88 | 0.00 | 0.12 | 0.00 | 0 |  | qchem.Hamiltonian; qchem.LASolver.Krylov; qchem.Response.LittleGroup; qchem.Response.OrbitalFrame; qchem.Response.Probe; qchem.Response.QStars |  |
| `qchem.SCFIterator.LoopDriver` | qchem | 1 | 0 | 1 | 5 | 0.83 | 0.33 | 0.17 | 0.17 | 1 | tDM_SourceSink<T>×1 | qchem.ChargeDensity.DensityMixer | GetChangeFrom, GetChargeDensity, GetTotalCharge |
| `qchem.Fitting.Internal.OrthoFunctionFitter` | qcFitting | 1 | 0 | 1 | 4 | 0.80 | 0.00 | 0.20 | 0.00 | 2 | ProjectedDensity_G×1, G_FieldEvaluator×1 | qchem.Fitting.FunctionFitter |  |
| `qchem.Fitting.Internal.OrthoNormalFunctionFitter` | qcFitting | 1 | 0 | 1 | 4 | 0.80 | 0.00 | 0.20 | 0.00 | 2 | G_RasterTransform×1, G_FieldEvaluator×1 | qchem.Fitting.FunctionFitter | GetScalarFunction |
| `qchem.Hamiltonian.Internal.ACBN0` | qcHamiltonian | 1 | 0 | 1 | 4 | 0.80 | 0.00 | 0.20 | 0.00 | 0 |  | qchem.Hamiltonian.HubbardEstimator; qchem.Hamiltonian.Internal.Hubbard |  |
| `qchem.BasisSet.PlaneWave.Evaluators.PeriodicGrid` | qcPlaneWave_BS | 1 | 0 | 1 | 7 | 0.88 | 0.00 | 0.12 | 0.00 | 0 |  | qchem.BasisSet.GMap (no-sig?); qchem.Mesh; qchem.ReciprocalLattice | GetGCartesian |
| `qchem.BasisSet.Lattice.APW_IBS` | qcLattice_BS | 1 | 1 | 1 | 10 | 0.91 | 0.00 | 0.09 | 0.00 | 0 |  | qchem.ReciprocalLattice | GetNumFunctions |
| `qchem.ChargeDensity.Imp.FittedCD` | qcChargeDensity | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 1 | ProjectedDensity_AO×1 |  | GetRepulsion, GetSelfRepulsion |
| `qchem.Response.QStars` | qcResponse | 1 | 0 | 1 | 4 | 0.80 | — | — | — | 0 |  | qchem.Symmetry.Lattice_3D.BlochQN; qchem.Symmetry.SymOp |  |
| `qchem.BasisSet.Radial.Evaluators.BSpline.Internal.Rk` | qcRadial_BS | 1 | 0 | 1 | 4 | 0.80 | 0.00 | 0.20 | 0.00 | 0 |  | qchem.BasisSet.Radial.Evaluators.Internal.Rk |  |
| `qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM` | qcGaussian_BS | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.PG_Cart.SymmetryAdapt` | qcGaussian_BS | 1 | 1 | 1 | 5 | 0.83 | — | — | — | 0 |  | qchem.BasisSet; qchem.Structure |  |
| `qchem.Factory` | qcFactory | 1 | 1 | 1 | 3 | 0.75 | — | — | — | 0 |  | qchem.BasisSet; qchem.BasisSet.Gaussian.Point.Factory (no-sig?); qchem.BasisSet.Radial.Factory |  |
| `qchem.ChargeDensity.Internal.LinearMixer` | qcChargeDensity | 1 | 0 | 1 | 3 | 0.75 | 0.00 | 0.25 | 0.00 | 0 |  | qchem.ChargeDensity.DensityMixer | GetChangeFrom, GetRelax, GetTotalCharge |
| `qchem.Orbitals.Internal.OrbitalImp` | qcOrbitals | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 0 |  |  | GetCoeff, GetCoeffPrime, GetDegeneracy, GetEigenEnergy, GetLabel, GetOccupation, GetQNs |
| `qchem.Matrix2D` | qcMath | 1 | 0 | 1 | 2 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.Vector2D | GetCol, GetRow |
| `qchem.KMesh` | qcStructure | 1 | 0 | 1 | 2 | 0.67 | 0.00 | 0.33 | 0.00 | 0 |  | qchem.Types |  |
| `qchem.LASolver.Internal.Lapack` | qcLASolver | 1 | 0 | 1 | 3 | 0.75 | 0.00 | 0.25 | 0.00 | 0 |  | qchem.LASolver | GetOrthoDim |
| `qchem.BasisSet.Internal.ERI4T` | qcBasisSet | 1 | 0 | 1 | 1 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  | qchem.Types |  |
| `qchem.BasisSet.Radial.Evaluators.Gaussian.Internal.ExponentScaler` | qcRadial_BS | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Radial.Evaluators.Slater.Internal.ExponentScaler` | qcRadial_BS | 1 | 0 | 1 | 5 | 0.83 | 0.00 | 0.17 | 0.00 | 0 |  |  |  |
| `qchem.BasisSet.Gaussian.Point.PG_Spherical.Symmetry` | qcGaussian_BS | 1 | 0 | 1 | 5 | 0.83 | — | — | — | 0 |  | qchem.Symmetry.Molecule.SALC (no-sig?) |  |
| `qchem.Fitting.FitOperations` | qcFitting | 1 | 1 | 1 | 2 | 0.67 | — | — | — | 0 |  | qchem.BasisSet.Orbital_DFT_IBS |  |
| `qchem.Orbitals.Factory` | qcOrbitals | 1 | 1 | 1 | 2 | 0.67 | — | — | — | 0 |  |  |  |
| `qchem.SCFIterator.Types` | qchem | 1 | 0 | 1 | 1 | 0.50 | — | — | — | 0 |  | qchem.BasisSet |  |
| `qchem.ValenceBasisGen` | qchem | 1 | 1 | 1 | 5 | 0.83 | — | — | — | 1 | rSpinResolved_CD×1 |  |  |
| `qchem.Vector2D` | qcMath | 1 | 0 | 1 | 1 | 0.50 | — | — | — | 0 |  | qchem.CMath |  |
| `qchem.Symmetry.Molecule.Irrep` | qcSymmetry | 1 | 1 | 1 | 1 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  | qchem.Symmetry | GetDegeneracy, GetLabel, GetPrincipleOffset |
| `qchem.ElectronConfiguration.ElectronCounts` | qcElConfig | 1 | 0 | 1 | 1 | 0.50 | — | — | — | 0 |  | qchem.Symmetry.Spin | GetN, GetNv |
| `qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ` | qcGaussian_BS | 1 | 0 | 1 | 1 | 0.50 | 0.00 | 0.50 | 0.00 | 0 |  |  | GetFjAt |

### Re-exports failing the signature test (no-sig?)

A re-export whose target contributes no name to the re-exporter's own exported text (TEXT heuristic): the suspect under the re-export rule.  `importers` = the re-exporting module's dependent files (every one receives the target by the side door).

| re-exporter | re-exported | importers |
|---|---|---|
| `qchem.BasisSet.Orbital_DFT_IBS` | `:Fit_IBS` | 61 |
| `qchem.BasisSet.Orbital_DFT_IBS` | `qchem.BasisSet.Fit_Types` | 61 |
| `qchem.BasisSet` | `qchem.BasisSet.GridResolution` | 60 |
| `qchem.BasisSet` | `qchem.ElectronConfiguration` | 60 |
| `qchem.BasisSet` | `qchem.Matrix3D` | 60 |
| `qchem.BasisSet` | `qchem.Streamable` | 60 |
| `qchem.ChargeDensity` | `qchem.Symmetry.Irrep` | 44 |
| `qchem.BasisSet.Lattice.BasisSet` | `qchem.BasisSet.Gaussian.Lattice.GPWTolerances` | 30 |
| `qchem.Fitting.FunctionFitter` | `qchem.BasisSet.GMap` | 28 |
| `qchem.BasisSet.GMap` | `qchem.Matrix3D` | 24 |
| `qchem.Hamiltonian` | `qchem.Hamiltonian.HubbardUTarget` | 21 |
| `qchem.WaveFunction` | `qchem.ElectronConfiguration` | 19 |
| `qchem.ChargeDensity.FourierDensity` | `qchem.BasisSet.GMap` | 17 |
| `qchem.MoleculeCalculation` | `qchem.Calculation` | 17 |
| `qchem.MoleculeCalculation` | `qchem.Response` | 17 |
| `qchem.MoleculeCalculation` | `qchem.SCFFailure` | 17 |
| `qchem.SolidCalculation` | `qchem.Calculation` | 17 |
| `qchem.SolidCalculation` | `qchem.SCFFailure` | 17 |
| `qchem.BasisSet.PlaneWave.Evaluators` | `qchem.BasisSet.GMap` | 15 |
| `qchem.Hamiltonian.Factory` | `qchem.Math.Gaussian` | 15 |
| `qchem.Hamiltonian.Factory` | `qchem.Outcome` | 15 |
| `qchem.BasisSet.Radial.Evaluators` | `qchem.Streamable` | 14 |
| `qchem.Materials` | `qchem.StructureData` | 14 |
| `qchem.SCFAccelerator` | `qchem.HandOff` | 14 |
| `qchem.AtomCalculation` | `qchem.Calculation` | 13 |
| `qchem.BasisSet.Internal.DB_Cache` | `qchem.BasisSet.DBCacheClient` | 13 |
| `qchem.Math` | `qchem.CMath` | 13 |
| `qchem.Math` | `qchem.Math.Gaussian` | 13 |
| `qchem.BasisSet.Gaussian.Lattice.LatticeSum1E` | `qchem.BasisSet.Gaussian.Lattice.GPWTolerances` | 11 |
| `qchem.BasisSet.Gaussian.Lattice.LatticeSum1E` | `qchem.BasisSet.GridResolution` | 11 |
| `qchem.Symmetry.Lattice_3D.BlochQN` | `qchem.Outcome` | 10 |
| `qchem.BasisSet.IrrepBasisSet` | `qchem.BasisSet.DBCacheClient` | 9 |
| `qchem.BasisSet.AtomCentredSet` | `qchem.Math.Gaussian` | 7 |
| `qchem.ChargeDensity.DensityMixer` | `qchem.HandOff` | 7 |
| `qchem.CMath` | `qchem.Constants` | 6 |
| `qchem.CMath` | `qchem.Factorials` | 6 |
| `qchem.CMath` | `qchem.IntPow` | 6 |
| `qchem.Mesh.Integrator` | `qchem.Mesh.Quadrature` | 6 |
| `qchem.Tests.GPW_Harness` | `qchem.Outcome` | 6 |
| `qchem.Tests.GPW_Harness` | `qchem.Reporting` | 6 |
| `qchem.BasisSet.Gaussian.Evaluators` | `qchem.Streamable` | 5 |
| `qchem.BasisSet.Gaussian.PG_Cart.Symmetry` | `qchem.Symmetry.Molecule.SALC` | 5 |
| `qchem.BasisSet.Orbital_HF_IBS` | `qchem.BasisSet.Orbital_1E_IBS` | 5 |
| `qchem.BasisSet.Radial.Evaluators.Internal.Rk` | `qchem.BasisSet.Internal.Cache4` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.BasisSet.GMap` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.BasisSet.G_FieldEvaluator` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.BasisSet.Orbital_DFT_IBS` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.Blaze` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.ChargeDensity.SeedCD` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.ElectronConfiguration` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.Ewald` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.Math` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.Matrix3D` | 5 |
| `qchem.Tests.PW_Harness` | `qchem.ScalarFunction` | 5 |
| `qchem.Tests.Panel` | `qchem.SCFCap` | 5 |
| `qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator` | `qchem.BasisSet.GMap` | 4 |
| `qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator` | `qchem.BasisSet.Gaussian.Lattice.GPWTolerances` | 4 |
| `qchem.BasisSet.Radial.Evaluators.Internal.AngularIntegrals` | `qchem.Types` | 4 |
| `qchem.Response.Reference` | `qchem.Outcome` | 4 |
| `qchem.BasisSet.Orbital_DHF_IBS` | `qchem.Structure` | 3 |
| `qchem.BasisSet.Radial.Evaluators.Internal.ExponentialEvaluator` | `qchem.Symmetry.Atom.Spherical` | 3 |
| `qchem.HDF5` | `qchem.Outcome` | 3 |
| `qchem.Hamiltonian.Internal.Hubbard` | `qchem.Hamiltonian.HubbardUTarget` | 3 |
| `qchem.Response` | `qchem.Response.Solver` | 3 |
| `qchem.SCFRecipeReport` | `qchem.Reporting` | 3 |
| `qchem.SolidState` | `qchem.Outcome` | 3 |
| `qchem.Symmetry.Molecule.AbelianGroup` | `qchem.Symmetry.Molecule.CharacterTable` | 3 |
| `qchem.Symmetry.Molecule.CartesianRep` | `qchem.Math.Angular` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.BasisSet.GMap` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.BasisSet.Lattice.BasisSet` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.BasisSet.Orbital_DFT_IBS` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.Blaze` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.Ewald` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.Math` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.Matrix3D` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.Structure` | 3 |
| `qchem.Tests.PW_Fields` | `qchem.Symmetry.Irrep` | 3 |
| `qchem.BasisSet.Gaussian.Evaluators.Internal.MnD` | `qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Index3` | 2 |
| `qchem.BasisSet.Radial.Evaluators.Internal.RelAngularIntegrals` | `qchem.Types` | 2 |
| `qchem.LASolver.Krylov` | `qchem.Outcome` | 2 |
| `qchem.Mesh.Folded` | `qchem.Symmetry.Spin` | 2 |
| `qchem.PotentialEnergySurface` | `qchem.Outcome` | 2 |
| `qchem.PotentialEnergySurface` | `qchem.SCFFailure` | 2 |
| `qchem.PseudoAtomOrbital` | `qchem.Math.Gaussian` | 2 |
| `qchem` | `qchem.Materials` | 1 |
| `qchem` | `qchem.PeriodicTable` | 1 |
| `qchem` | `qchem.StructureData` | 1 |
| `qchem` | `qchem.Types` | 1 |
| `qchem.BasisSet.Gaussian.Lattice.RasterLadder` | `qchem.BasisSet.Gaussian.Lattice.GPWTolerances` | 1 |
| `qchem.BasisSet.Gaussian.Point.PG_Spherical.Symmetry` | `qchem.Symmetry.Molecule.SALC` | 1 |
| `qchem.BasisSet.PlaneWave.Evaluators.PeriodicGrid` | `qchem.BasisSet.GMap` | 1 |
| `qchem.BasisSet.Radial.BasisSet` | `qchem.BasisSet.Orbital_HF_IBS` | 1 |
| `qchem.Factory` | `qchem.BasisSet.Gaussian.Point.Factory` | 1 |
| `qchem.WaveFunction.Internal.IrrepWF` | `qchem.ElectronConfiguration` | 1 |
