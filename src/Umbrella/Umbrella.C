// File: Umbrella/Umbrella.C  The `qchem` umbrella module -- one import for consumer / binding code.
//
// `import qchem;` brings the whole molecule front door into scope at once: the Calculation facade and
// its options (+ qchem::Model/Pol), the Molecule/Atom builders, ScalarFunction sampling, and -- via
// qchem.SCFIterator -- SCFParams, the live SCFProgress trace, and the WaveFunction/Orbitals/
// ChargeDensity/Irrep/EnergyBreakdown query surface.  End-user and binding code should prefer this;
// library-internal code keeps importing the granular modules (e.g. `import qchem.Calculation;`) to
// minimise recompile coupling.  Internals (*.Internal.*) are deliberately NOT re-exported.
//
// Lives in its OWN library (qcUmbrella, above the aggregate) ON PURPOSE: an umbrella that re-exports a
// SAME-LIBRARY sibling trips CMake's module scanner into a bogus cycle (see the qcMath-split lesson).
// Here every re-exported module is in a different, already-linked library, so there is no back-edge.
export module qchem;

export import qchem.Calculation;     // Calculation, CalcOptions, AcceleratorOptions, qchem::Model/Pol
export import qchem.SCFIterator;     // SCFParams, SCFProgress, SCFIterator + (transitive) WaveFunction,
                                     //   Orbitals, ChargeDensity, Irrep, EnergyBreakdown, Hamiltonian
export import qchem.Structure;       // Molecule, Atom, Structure
export import qchem.ScalarFunction;  // qchem::ScalarFunction<double> (the density and every MO)
export import qchem.Types;           // qchem::Vector3D / rvec3_t, dcmplx
export import qchem.PeriodicTable;   // element symbol / Z lookups
