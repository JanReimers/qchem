// File: Hamiltonian/tests/GGA_Term.C  The GGA XC term on the molecular point route (GGA stage 2): H = dE/dD.
//
// tVxc_DFT_GGA's E is Σ_g w_g e(ρ, σ) on the density-gradient quadrature and its H the {v, w} gather.  The claim that makes the SCF (and
// GDM) right is H = ∂E/∂D: a central difference of E along a symmetric direction Δ of D equals Tr(Δ H).  That checks the whole chain at
// once -- libxc's v_ρ and v_σ (each FD-checked alone in GGA_XC.C), the w assembly 2 v_σσσ ∇ρ_σ + v_σ↑↓ ∇ρ_σ̄, the ρ/2 collapse, and the
// sampler's adjoint pair -- on H₂O with a real (DZVP) basis and the Becke mesh a molecular run builds.
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

import qchem.Hamiltonian.Internal.Vxc_DFT_GGA;      // the term under test (tests may cheat-import internals)
import qchem.Hamiltonian.Internal.GGAFunctional;
import qchem.Hamiltonian.Internal.Libxc_GGA;
import qchem.Hamiltonian;                         // tDynamic_HT, DensityID
import qchem.Hamiltonian.Types;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.DensitySampler;        // MakeDensityGradientSampler
import qchem.ChargeDensity.Factory;               // IrrepCD_Factory<double>
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.Tests.PW_Fields;                     // OneBlockBS: a one-block whole basis for the refresh
import qchem.Mesh;
import qchem.Structure;
import qchem.Symmetry.Irrep;
import qchem.Types;
import qchem.Blaze;
import qchem.Structure.MeshFactory;   // MakeMeshMaker: the quadrature a delta fit is built over (V-MESH-FACTORY step 4)
using namespace qchem;
using namespace qchem::Hamiltonian;

namespace
{
Molecule H2O()
{
    Molecule m;
    m.Insert(Site{8, rvec3_t( 0.000, 0.000, 0.000)});
    m.Insert(Site{1, rvec3_t( 1.430, 1.107, 0.000)});
    m.Insert(Site{1, rvec3_t(-1.430, 1.107, 0.000)});
    return m;
}
//! A deterministic PSD D = L Lᵀ (rank 5: about five doubly occupied orbitals), and a symmetric direction to move it along.
hmat_t<double> SomePSD_D(size_t n, double phase)
{
    hmat_t<double> D(n);
    for (size_t i=0; i<n; i++) for (size_t j=0; j<=i; j++)
    {
        double s=0.0;
        for (size_t m=0; m<5; m++) s+=std::sin(phase+1.1*i+0.7*m)*std::sin(phase+1.1*j+0.7*m);
        D(i,j)=2.0*s/double(n);
    }
    return D;
}

struct Rig
{
    Molecule                                         mol=H2O();
    std::unique_ptr<BasisSet::Real_BS>               bs;
    const BasisSet::Orbital_1E_IBS<double>*          block=nullptr;
    std::shared_ptr<const BasisSet::rFIT_SF_ABS>     fb;
    std::shared_ptr<const ChargeDensity::tDensityGradientSampler<double>> sampler;
    std::unique_ptr<tVxc_DFT_GGA<double>>            term;
    Rig()
    {
        bs.reset(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, &mol,
                                             BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        block=dynamic_cast<const BasisSet::Orbital_1E_IBS<double>*>((*bs)[0]);
        const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=40, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=17, .beckeOrder=2};
        BasisSet::FitQuadrature quad=MakeMeshMaker(mol)->Make(mp, qcMesh::NoSymmetry());
        fb.reset(bs->CreateVxcFitBasisSet(&mol, BasisSet::VxcFit::Delta, &quad));
        sampler=ChargeDensity::MakeDensityGradientSampler(fb, quad);
        auto pbe=std::make_shared<const CompositeGGAFunctional>(std::vector<std::shared_ptr<const GGAFunctional>>{
            std::make_shared<const Libxc_GGA>(101), std::make_shared<const Libxc_GGA>(130)});   // GGA_X_PBE + GGA_C_PBE
        term=std::make_unique<tVxc_DFT_GGA<double>>(pbe, sampler, SpinGroup::UnPolarized);
    }
    std::unique_ptr<ChargeDensity::tDM_CD<double>> Density(const hmat_t<double>& D) const
    {
        return std::unique_ptr<ChargeDensity::tDM_CD<double>>(ChargeDensity::IrrepCD_Factory<double>(D, block, bs->GetIrreps(Spin::None)[0]));
    }
    double Exc(const hmat_t<double>& D) const
    {
        const auto cd=Density(D);
        term->RefreshForDensity(cd.get());
        EnergyBreakdown e;
        term->GetEnergy(e, cd.get());
        return e["Exc"];
    }
    hmat_t<double> H(const hmat_t<double>& D) const
    {
        const auto cd=Density(D);
        term->RefreshForDensity(cd.get());
        qchem::tests::pw::OneBlockBS<double> one(block);
        term->BuildBlocks(&one, cd.get());
        const tDynamic_HT<double>& face=*term;   // the term face the Fock loop reads through
        return face.GetEMatrix(block, Spin::None, DensityID(cd.get()));
    }
};
} //anon

TEST(GGATerm, ItsMatrixIsTheDerivativeOfItsEnergy)
{
    Rig rig;
    const size_t n=rig.block->GetNumFunctions();
    const hmat_t<double> D=SomePSD_D(n, 0.3), Delta=SomePSD_D(n, 1.7);
    const hmat_t<double> H=rig.H(D);
    double trDH=0.0, scale=0.0;
    for (size_t i=0; i<n; i++) for (size_t j=0; j<n; j++) {trDH+=Delta(i,j)*H(i,j); scale=std::max(scale, std::fabs(H(i,j)));}
    const double e=1e-4;
    hmat_t<double> Dp=D, Dm=D, Dp2=D, Dm2=D;
    for (size_t i=0; i<n; i++) for (size_t j=0; j<=i; j++)
    {
        Dp(i,j)+=e*Delta(i,j); Dm(i,j)-=e*Delta(i,j);
        Dp2(i,j)+=2*e*Delta(i,j); Dm2(i,j)-=2*e*Delta(i,j);
    }
    const double d1=(rig.Exc(Dp)-rig.Exc(Dm))/(2*e), d2=(rig.Exc(Dp2)-rig.Exc(Dm2))/(4*e);
    const double fd=(4*d1-d2)/3;   // Richardson: O(e^4)
    std::printf("  [GGA term] E_xc %.10f, Tr(Delta H) %.12e, dE/de (Richardson) %.12e, max|H| %.3e\n", rig.Exc(D), trDH, fd, scale);
    ASSERT_GT(std::fabs(trDH), 1e-3) << "a vacuous direction would pass anything";
    EXPECT_NEAR(trDH, fd, 1e-8*std::fabs(fd)) << "H = dE/dD: the {v, w} gather is the derivative of E = sum_g w_g e(rho, sigma)";
}

TEST(GGATerm, ItsEnergyIsNotTheLDAsAndItsSigmaIsSampled)
{
    Rig rig;
    const hmat_t<double> D=SomePSD_D(rig.block->GetNumFunctions(), 0.3);
    const auto cd=rig.Density(D);
    rig.term->RefreshForDensity(cd.get());   // the forgotten-gradient detector runs here and must pass
    const rvec3vec_t& g=rig.sampler->GradRhoPol(cd.get(), Spin::Up);
    double gmax=0.0;
    for (const auto& v : g) gmax=std::max(gmax, std::sqrt(v*v));
    EXPECT_GT(gmax, 1.0) << "the term's density gradient is sampled (a molecule's density varies)";
}
