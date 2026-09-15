// File: IntegrationTests/PW/Boxes.C  two Mn in a cubic box: the spin-SAD seed's AFM staggering at the plane-wave level (seed-level, no SCF).
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(PW_<Material>, <k>_[tokens]_<Claim>); `Prototype` marks the standalone
// loop.  `scripts/testgrid` renders the coverage table from these names.
//
//   PW_Mn2Box.Γ_Pol_SeedStaggered

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

// The SPIN-POLARIZED SAD seed (PolarizedSeedCD, doc/SCFSeedingPlan.md §10) -- the AFM staggering gate.
// Two identical Mn (q7 valence, library Hund pair with moment 5) on a CsCl arrangement, ONE flipped: the
// smallest staggered magnetic cell.  Everything below is exact seed algebra (no SCF): the per-channel
// structure-factor sums against the stored (majority, minority) pair, with the flip swapping the pair on
// the -m site.
TEST(PW_Mn2Box, Γ_Pol_SeedStaggered)
{
    using qchem::ChargeDensity::SeedCD;
    using qchem::ChargeDensity::PolarizedSeedCD;
    const double a=8.4;
    auto mkcell=[&](bool f0, bool f1)
    {
        auto c=std::make_shared<UnitCell>(a);
        c->AddAtom(25, {0,0,0},       f0);        // Mn sublattice +m (unless flipped)
        c->AddAtom(25, {0.5,0.5,0.5}, f1);        // Mn sublattice -m
        return c;
    };
    auto cellA=mkcell(false,true), cellB=mkcell(true,false);   // B = the globally spin-mirrored pattern

    Lattice_3D lat(*cellA, ivec3_t(1,1,1));
    PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS>  fb (pw.CreateCDFitBasisSet (cellA.get(), qcMesh::MeshParams{}));
    std::shared_ptr<const BasisSet::cFIT_SF_ABS>  fsf(pw.CreateVxcFitBasisSet(cellA.get(), qcMesh::MeshParams{}));

    PolarizedSeedCD A(fb, cellA.get());
    PolarizedSeedCD B(fb, cellB.get());
    SeedCD          T(fb, cellA.get());           // the spin-agnostic total (flip-blind)

    // Channel bookkeeping: each channel carries majority(one site) + minority(the other) = (6+1) = 7 e-;
    // the staggering cancels the NET spin exactly while both sites are fully staggered locally.
    EXPECT_NEAR(A.GetTotalCharge(), 14.0, 0.02);                  // 2 x 7 valence electrons
    EXPECT_NEAR(A.GetTotalSpin(),    0.0, 1e-10);                 // AFM: <up>-<down> = 0 by symmetry
    EXPECT_NEAR(A.GetChannel(Spin::Up)->GetTotalCharge(), 7.0, 0.01);

    // The channel rho-tilde maps.  A channel is a matrix-free FourierDensity; the seed ignores the SF fit
    // basis argument (its support is its own CD fit basis's difference set), so every map shares one key set.
    auto fd=[&](const qchem::ChargeDensity::cChargeDensity* cd)->ΔG_Map
    {
        auto* f=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
        EXPECT_TRUE(f) << "a seed channel must be a FourierDensity";
        return f->GetFourierDensity(*fsf);
    };
    ΔG_Map upA=fd(A.GetChannel(Spin::Up)), dnA=fd(A.GetChannel(Spin::Down));
    ΔG_Map dnB=fd(B.GetChannel(Spin::Down));
    ΔG_Map tot=T.GetFourierDensity(*fsf);

    // (1) The staggered ORDER PARAMETER at the AFM wavevector G=(1,1,1): identical species on both sites
    // means the TOTAL suffers the CsCl extinction there (F_maj+F_min alike, phases +1/-1), while the spin
    // density m-tilde = rho-up - rho-dn = 2(F_maj-F_min) is large -- the seed genuinely staggers.
    const ivec3_t q(1,1,1);
    ASSERT_TRUE(upA.count(q) && dnA.count(q));
    // rho-tilde carries 1/Omega, so scale back: |m-tilde(q)|*Omega/2 = F_maj-F_min at |G_q| -- the
    // staggered moment per site the seed expresses at the AFM wavevector (5 e- at G=0, ~4.4 by |G_q|).
    EXPECT_GT  (std::abs(upA.at(q)-dnA.at(q))*cellA->GetCellVolume()/2, 1.0)
        << "the staggered per-site moment at q_AFM must be electrons-scale";
    EXPECT_NEAR(std::abs(upA.at(q)+dnA.at(q)), 0.0, 1e-10) << "same-species extinction of the total";

    // (2) The global spin mirror: pattern B (both flips inverted) is A with the channels swapped.
    for (const auto& [dm,v] : upA)
    {
        ASSERT_TRUE(dnB.count(dm));
        EXPECT_NEAR(std::abs(v-dnB.at(dm)), 0.0, 1e-12) << "dm=("<<dm.x<<","<<dm.y<<","<<dm.z<<")";
    }

    // (3) Flip-blind total: up+dn reproduces the spin-agnostic SeedCD map coefficient by coefficient
    // (the two-group assembly must not perturb the total the Hartree term consumes).
    for (const auto& [dm,v] : tot)
    {
        dcmplx s = (upA.count(dm)?upA.at(dm):dcmplx(0)) + (dnA.count(dm)?dnA.at(dm):dcmplx(0));
        EXPECT_NEAR(std::abs(s-v), 0.0, 1e-11) << "dm=("<<dm.x<<","<<dm.y<<","<<dm.z<<")";
    }
}

} //namespace
