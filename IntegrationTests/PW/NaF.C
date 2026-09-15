// File: IntegrationTests/PW/NaF.C  NaF rocksalt on the plane-wave basis (multi-species PP, IonicSAD vs uniform seed).
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(PW_<Material>, <k>_[tokens]_<Claim>); `Prototype` marks the standalone
// loop.  `scripts/testgrid` renders the coverage table from these names.
//
//   PW_NaF.Γ_Anchor

#include <map>
#include <set>
#include <memory>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include "gtest/gtest.h"

import qchem.Tests.PW_Harness;   // the PW drivers + the direct-grid oracles + PWFixture
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
import qchem.BasisSet.PlaneWave.Evaluators;   // PW_Grid_Evaluator -- a UNIT TEST may reach the internal
import qchem.BasisSet.Lattice.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
import qchem.ScalarFunction;                 // ScalarFunction<double> -- arg of the moved-here field oracles
import qchem.Mesh;                           // qcMesh::MeshParams (Vee_Hartree's fit-basis factory arg; ignored)
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
import qchem.Blaze;          // mat_t<dcmplx>
import qchem.Math;           // Pi
import qchem.Hamiltonian.Internal.ExFunctional;     // the validated LDA functional interface
import qchem.Hamiltonian.Internal.SlaterExchange;   // Dirac exchange (alpha=2/3), eps_x = 3/4 v_x
import qchem.Hamiltonian.Internal.VWN_Correlation;  // VWN5 correlation (validated vs libxc)
import qchem.Hamiltonian.Internal.PWTerms;          // Ven_PP_Short/Long, Vee_Hartree (dcmplx Hamiltonian terms)
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
import qchem.Hamiltonian;                           // cStatic_HT / cDynamic_HT aliases (public term interfaces)
import qchem.Hamiltonian.Internal.Hamiltonian;      // cHamiltonianImp (the dcmplx Hamiltonian = sum of terms)
import qchem.Hamiltonian.Internal.Hamiltonians;     // Ham_PW_DFT (the assembled plane-wave LDA KS Hamiltonian)
import qchem.Hamiltonian.Internal.IonIon;           // IonIon<double> (ion-ion term: pair sum / Ewald via isFinite)
import qchem.Hamiltonian.Internal.Kinetic;          // Kinetic<dcmplx> (the shared kinetic term)
import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
import qchem.Energy;                                // EnergyBreakdown
import qchem.ChargeDensity.Imp.IrrepCD;             // PeriodicIrrepCD<dcmplx> (the periodic leaf -- 3c-2b)
import qchem.ChargeDensity.SeedCD;                  // SeedCD + PolarizedSeedCD (the spin-SAD staggering gate)
import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
import qchem.BasisSet.GMap;                             // ΔG_Map (the G-space coefficient map)
import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
import qchem.Symmetry.Irrep;                        // Irrep
import qchem.LASolver;                              // complex Hermitian eigensolver
import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
import qchem.SCFIterator;                           // cSCFIterator (the real framework SCF driver)
import qchem.SCFParams;                             // SCFParams
import qchem.WaveFunction;                          // cWaveFunction (read view of the converged state)
import qchem.ElectronConfiguration;                 // ElectronConfiguration base
import qchem.ElectronConfiguration.Crystal;         // Crystal_EC (single-k Bloch configuration)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;      // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

using namespace qchem;
using namespace qchem::tests::pw;

using BasisSet::PlaneWave::PlaneWave_IBS;
using BasisSet::PlaneWave::PW_Grid_Evaluator;   // internal grid evaluator (unit test may reach it directly)
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

namespace
{

// Multi-species ionic crystal NaF (rocksalt = FCC + 2-atom basis), through the full SCFIterator with the
// multi-species Ham_PW_DFT facade.  Na (Zion=1) + F (Zion=7) = 8 valence electrons; the per-Z router model
// dispatches Na's vs F's pseudopotential per atom.  F's tight 2p sets the (high) cutoff.
TEST(PW_NaF, Γ_Anchor)
{
    using namespace qchem::Hamiltonian;
    const double a=8.73;                          // NaF lattice constant ~4.62 A (a.u.)
    FCCUnitCell cell(a);                          // rocksalt: FCC lattice, 2-atom basis
    cell.AddAtom(11, {0,0,0});                    // Na  (true species Z=11; Zion=1 via the PP)
    cell.AddAtom(9,  {0.5,0.5,0.5});              // F   (true species Z=9;  Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    using SS=qchem::ChargeDensity::SeedStrategy;
    auto mkHam=[&](const BasisSet::Complex_BS* bs)   // multi-species one-call; bs -> the Hartree fit basis
        { return new Ham_PW_DFT(lat.GetStructure(), bs, {{"Na",1},{"F",7}}, "LDA"); };

    // IonicSAD seed: the electronegativity heuristic assigns Na+ / F-, and the per-species valence densities
    // are scaled to the formal-charge electron counts (0 e- on Na, 8 e- on F) -- the seed pre-bakes the
    // Na->F charge transfer and CONSERVES the cell electron count.  This asserts the seed is CORRECT (charge
    // sum rule + same converged energy as Uniform: the seed cannot change the answer).
    FwResult I=RunFrameworkGamma(lat, /*Ecut*/6.0, /*Nelec*/8, mkHam, "NaF IonicSAD", SS::IonicSAD);
    FwResult U=RunFrameworkGamma(lat, /*Ecut*/6.0, /*Nelec*/8, mkHam, "NaF Uniform",  SS::Uniform);
    std::cout << "[NaF seed compare] IonicSAD iters="<<I.iters<<"  Uniform iters="<<U.iters << std::endl;

    EXPECT_TRUE(I.converged);
    EXPECT_NEAR(I.charge, 8.0, 1e-6);             // 1 (Na) + 7 (F) valence electrons, conserved by the ionic seed
    // Regression anchor (Ecut=6, Gamma-only -> underconverged but deterministic), like the Si total.
    // Negative: the ionic Madelung (Enn~-14) + G=0 alignment dominate.
    EXPECT_NEAR(I.E.GetTotalEnergy(), -20.3293, 5e-3) << "Enn="<<I.E["Enn"]<<" E_alphaZ="<<I.E["E_alphaZ"];
    EXPECT_NEAR(I.E.GetTotalEnergy(), U.E.GetTotalEnergy(), 1e-3);   // seed-independence of the converged answer
    // IonicSAD now HALVES the iterations vs Uniform (17 vs 35 at Ecut=6) -- the DIFFUSE F- pseudo-valence
    // density does it.  History: the old seed scaled the NEUTRAL F valence x8/7 (too COMPACT -> high-G noise ->
    // 69 vs 58, WORSE).  Fix (2026-07-12): SeedCD/IonicSAD now pull the library's CHARGE-STATE density (F-,
    // Nelec=8, <r>=1.62 vs neutral 1.18) generated offline by qchem::ValenceBasisGen::GenerateSeedDensity, so
    // the anion diffuseness is physical, not an amplitude hack.  A guard that the win holds:
    EXPECT_LT(I.iters, U.iters) << "diffuse F- IonicSAD should converge in fewer iterations than Uniform";
}

} //namespace
