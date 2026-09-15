// File: IntegrationTests/PW/CsI.C  CsI (CsCl structure) on the plane-wave basis -- the d-projector test (Cs q1, I q7).
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(PW_<Material>, <k>_[tokens]_<Claim>); `Prototype` marks the standalone
// loop.  `scripts/testgrid` renders the coverage table from these names.
//
//   PW_CsI.Γ_Anchor

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

// Heavy multi-species ionic crystal CsI (CsCl = simple cubic + 2-atom basis).  THE d-PROJECTOR TEST: both
// Cs (q1) and I (q7) carry l=2 (d) Kleinman-Bylander channels, the first DFT exercise of the (2l+1)P_2
// angular path (Si/NaF are l<=1).  Cs q1 (Zion=1) deliberately avoids the semicore q9, whose l=3 (f)
// projector the analytic HGH Qli table doesn't tabulate.  And the PP promise: Cs/I are SOFTER than F
// (bigger r_loc), so this heavy salt needs a LOWER cutoff than NaF.
TEST(PW_CsI, Γ_Anchor)
{
    using namespace qchem::Hamiltonian;
    const double a=8.63;                          // CsI lattice constant ~4.567 A (a.u.)
    UnitCell cell(a);                             // CsCl: simple cubic, 2-atom basis
    cell.AddAtom(55, {0,0,0});                    // Cs  (true species Z=55; Zion=1, q1 -- s,p,d channels)
    cell.AddAtom(53, {0.5,0.5,0.5});              // I   (true species Z=53; Zion=7 -- s,p,d channels)
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    auto mkHam=[&](const BasisSet::Complex_BS* bs)   // bs -> the Hartree fit basis (created once in BuildTerms)
        { return new Ham_PW_DFT(lat.GetStructure(), bs, {{"Cs",1},{"I",7}}, "LDA"); };

    FwResult R=RunFrameworkGamma(lat, /*Ecut*/4.0, /*Nelec*/8, mkHam, "CsI SCFIterator-Gamma");

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);             // 1 (Cs) + 7 (I) valence electrons (d-projectors active)
    // Regression anchor (Ecut=4, Gamma-only).  The point is the d-channel assembly runs end-to-end;
    // both species' l=2 Kleinman-Bylander projectors contribute via the (2l+1)P_2(cos gamma) path.
    EXPECT_NEAR(R.E.GetTotalEnergy(), -11.3868, 5e-3) << "Enn="<<R.E["Enn"]<<" E_alphaZ="<<R.E["E_alphaZ"];
}

} //namespace
