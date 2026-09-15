// File: IntegrationTests/PW/Model.C  the two model systems (jellium, a weak cosine potential) the prototype loop was built on.
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(PW_<Material>, <k>_[tokens]_<Claim>); `Prototype` marks the standalone
// loop.  `scripts/testgrid` renders the coverage table from these names.
//
//   PW_Jellium.Γ_Prototype_IsUniform
//   PW_Jellium.Γ_Prototype_XcUniformLimit
//   PW_Cosine.Γ_Prototype_Converges
//   PW_Cosine.Γ_Prototype_VxcFitConverges

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

// Jellium: no external potential, 2 electrons.  The self-consistent density stays uniform (only the
// G=0 component), so E_H=0, the kinetic energy is zero (G=0 band), and E_tot = E_xc = Omega eps_xc rho0.
// Also a sanity check that the two energy routes agree.
TEST(PW_Jellium, Γ_Prototype_IsUniform)
{
    PWFixture F;
    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};
    chmat_t Vzero=F.pw.OverlapMatrix([](const ivec3_t&){return dcmplx(0.0);});

    SCFResult R=RunSCF(F.pw, F.B(), F.Omega, Vzero, 2, ivec3_t(8,8,8), vxcOf, epsOf);
    ASSERT_TRUE(R.converged);

    const double rho0=2.0/F.Omega;
    for (const auto& kv : R.rho)                                   // density is uniform
        if (!(kv.first==ivec3_t(0,0,0))) EXPECT_NEAR(std::abs(kv.second), 0.0, 1e-8);
    EXPECT_NEAR(R.Ekin, 0.0, 1e-9);                               // G=0 band: zero kinetic
    EXPECT_NEAR(R.EH,   0.0, 1e-9);                               // uniform: no Hartree
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-8);               // two routes agree
    EXPECT_NEAR(R.Etot_direct, F.Omega*epsOf(rho0)*rho0, 1e-7);  // = N eps_xc(rho0)
}


// Uniform-density limit: rho = rho0 everywhere => Vxc(r)=v_xc(rho0) constant, so Vxc~(0)=v_xc(rho0),
// Vxc~(dm!=0)=0, and E_xc = Omega eps_xc(rho0) rho0.  Uses the validated Dirac exchange (eps_x=3/4 v_x).
TEST(PW_Jellium, Γ_Prototype_XcUniformLimit)
{
    PWFixture F;
    const double rho0=0.05;
    RhoG rho;
    rho[ivec3_t(0,0,0)] = rho0;

    qchem::Hamiltonian::SlaterExchange ex(2.0/3.0);              // Dirac exchange
    auto vxcOf =[&](double r){return ex.GetVxc (r);};
    auto epsOf =[&](double r){return ex.GetEpsXc(r);};

    double Exc=0.0, VxcDc=0.0;
    RhoG vtil=BuildXcVtilde(F.pw, rho, ivec3_t(6,6,6), F.Omega, vxcOf, epsOf, Exc, VxcDc);

    EXPECT_NEAR(std::real(RhoAt(vtil,ivec3_t(0,0,0))), ex.GetVxc(rho0), 1e-12);
    EXPECT_NEAR(std::abs (RhoAt(vtil,ivec3_t(1,0,0))), 0.0,             1e-12);
    EXPECT_NEAR(Exc,   F.Omega*ex.GetEpsXc(rho0)*rho0, 1e-10);
    EXPECT_NEAR(VxcDc, F.Omega*ex.GetVxc  (rho0)*rho0, 1e-10);   // integral rho v_xc
}


// A weak external cosine well + 2 electrons + Hartree + LDA(Dirac+VWN): a non-trivial self-consistent
// loop (Hartree and XC respond to a genuinely modulated density).  No external reference, so we check
// the strongest internal property: at the fixed point the band-sum and direct total energies agree.
TEST(PW_Cosine, Γ_Prototype_Converges)
{
    const double a=6.0, Ecut=4.0, Omega=a*a*a, V0=-0.3;
    UnitCell          cell(a);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), Ecut);

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};
    // V_ext(r) = 2 V0 (cos + cos + cos): only Fourier components are the unit reciprocal steps.
    chmat_t Vext=pw.OverlapMatrix([V0](const ivec3_t& dm)->dcmplx
    { return (dm.x*dm.x+dm.y*dm.y+dm.z*dm.z==1) ? dcmplx(V0) : dcmplx(0.0); });

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, 2, ivec3_t(12,12,12), vxcOf, epsOf, 0.5, 1e-9, 400);
    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-6);                       // stationarity / consistency
    EXPECT_GT (std::abs(RhoAt(R.rho, ivec3_t(1,0,0))), 1e-4);            // density genuinely modulated
}


// Item K EXPLORATION on a SELF-CONSISTENT density: converge the weak-cosine LDA density, then sweep the
// Vxc FIT grid (relCutoff 1->4->16, i.e. 16^3->32^3->64^3) and print how the XC energy quadrature E_xc=∫ε ρ
// and the SPATIAL fit residual ‖v_xc - v_xc,fit‖ behave.  Illustrative (prints a table); the density is real.
TEST(PW_Cosine, Γ_Prototype_VxcFitConverges)
{
    const double a=6.0, Ecut=4.0, Omega=a*a*a, V0=-0.3;
    UnitCell          cell(a);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), Ecut);

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};
    chmat_t Vext=pw.OverlapMatrix([V0](const ivec3_t& dm)->dcmplx
    { return (dm.x*dm.x+dm.y*dm.y+dm.z*dm.z==1) ? dcmplx(V0) : dcmplx(0.0); });

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, 2, ivec3_t(12,12,12), vxcOf, epsOf, 0.5, 1e-9, 400);
    ASSERT_TRUE(R.converged);

    ΔG_Map rho;                                              // the self-consistent density as a ΔG_Map
    for (const auto& kv : R.rho) rho[kv.first]=kv.second;
    PW_Grid_Evaluator grid=GridOf(pw);                       // the density grid matching this block's {G}
    FieldFnR vtrue([&](const rvec3_t& r){return vxcOf(grid.EvalField(rho, r));});   // exact v_xc(r) of the SCF density

    std::vector<rvec3_t> ref;                                // fit-grid-INDEPENDENT reference points (8^3)
    for (int i=0;i<8;i++) for (int j=0;j<8;j++) for (int k=0;k<8;k++)
        ref.push_back(cell.ToCartesian(rvec3_t((i+0.5)/8.0,(j+0.5)/8.0,(k+0.5)/8.0)));

    std::cout << "\n  SCF: converged in " << R.iters << " iters, Etot_direct=" << R.Etot_direct
              << ", Exc(scf grid)=" << R.Exc << "\n"
              << "  relCutoff | nG_fit |   FFT pts |      E_xc(∫ε ρ) | ‖v_xc−v_xc,fit‖\n"
              <<   "  ----------+--------+-----------+-----------------+----------------\n";
    std::vector<double> Excs, resids;
    for (double rc : {1.0, 2.0, 4.0})
    {
        qcMesh::MeshParams mp; mp.relCutoff=rc;
        auto fb=qchem::ChargeDensity::fitbasis_t(pw.CreateVxcFitBasisSet(nullptr, mp));
        // The fit BASIS counts its own {G} FUNCTIONS; the RASTER counts voxels and owns the uniform
        // quadrature rule over them.  Two different numbers, asked of the two different faces
        // (2026-08-23: NumPoints/Integrate came off the fit face, which is about functions).
        auto ge=dynamic_cast<const qchem::BasisSet::G_RasterTransform*>(fb.get());
        size_t nG=fb->GetNumFunctions(), Npts=ge->RasterSize();

        rvec_t rgrid=ge->RhoOnGrid(rho), exc(rgrid.size());
        for (size_t q=0;q<rgrid.size();q++) exc[q]=epsOf(rgrid[q])*rgrid[q];
        double Exc=ge->Integral(exc);

        auto fitter=qchem::Fitting::Factory(fb);
        fitter->DoFit(vtrue);
        // The fitted FIELD is a capability now, not part of the fitter face (a delta fit has no
        // value between its points), so ask for it: a plane-wave fit inverse-transforms and can.
        const ScalarFunction<double>& vfit=dynamic_cast<const ScalarFunction<double>&>(*fitter);
        double s=0.0;
        for (const rvec3_t& r : ref){double d=vfit(r)-vtrue(r); s+=d*d;}
        double resid=std::sqrt(s/ref.size());

        char row[256];
        std::snprintf(row,sizeof row,"  %8.0f  | %6zu | %9zu | %15.10f | %14.6e\n",rc,nG,Npts,Exc,resid);
        std::cout << row;
        Excs.push_back(Exc); resids.push_back(resid);
    }
    std::cout << std::endl;

    // The lesson, asserted: the SPATIAL fit residual converges FAST (spectrally, for this smooth density) as
    // the grid densifies -- but the XC ENERGY is essentially BLIND to it (flat to <1e-6).  So acceptance MUST
    // be field/density convergence vs a fine reference, NEVER dE_total (the non-variational-fit pin).
    EXPECT_LT(resids[1], resids[0]);              // residual converges as the fit grid densifies ...
    EXPECT_LT(resids[2], resids[1]);              // ... monotonically
    EXPECT_LT(std::abs(Excs[2]-Excs[0]), 1e-6);   // while E_xc barely moves -- energy is the wrong yardstick
}

} //namespace
