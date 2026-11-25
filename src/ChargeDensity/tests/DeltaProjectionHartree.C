// File: ChargeDensity/tests/DeltaProjectionHartree.C  The δ-projection Hartree fitter (Fitting::DeltaProjectionDensityFitter, lane S 2026-11-25)
// at a FIXED density, on water's first orbital block over a small δ mesh:
// - its Fock is the exact derivative of its energy: E_H is quadratic in D, so the central difference along any X is exact to roundoff,
//   and must equal Tr(X V) (the ruling's analytically correct expression: both legs, the q-leg included);
// - Tr(D V) = 2 E_H (the two legs have equal trace);
// - the charge is the mesh's Σ w_g ρ(r_g), and E_H agrees with the Poisson-free direct-sum δ fit on the same points to the mesh's error;
// - its projections are two per point, and DoFit refuses a vector that is not.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <iostream>
#include <string>

import qchem.ChargeDensity;
import qchem.BasisSet;                          // Real_BS, VxcFit, rFIT_SF_ABS, PairProjector
import qchem.BasisSet.Gaussian.Point.Factory;   // Gaussian::Factory
import qchem.BasisSet.CoulombDeltaFit_IBS;      // the δ Coulomb fit basis (the comparison route)
import qchem.Fitting.CoulombOperator;           // MakeDeltaProjectionDensityFitter / MakeDeltaDensityFitter
import qchem.Mesh;                              // qcMesh::MeshParams
import qchem.Mesh.CoulombOperator;              // MakeDirectSumCoulomb (the comparison route)
import qchem.Symmetry.Unit;                     // UnitQN
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;
import qchem.Structure.MeshFactory;   // MakeMeshMaker
using namespace qchem;

namespace
{
Molecule Water()
{
    Molecule m;
    m.Insert(Site{8, rvec3_t( 0.0,     0.0,     0.0)});
    m.Insert(Site{1, rvec3_t( 1.43042, 1.10715, 0.0)});
    m.Insert(Site{1, rvec3_t(-1.43042, 1.10715, 0.0)});
    return m;
}
const BasisSet::Real_DFT_OIBS* FirstBlock(const BasisSet::Real_BS& bs)
{
    for (auto b : bs.Iterate<BasisSet::Real_DFT_OIBS>()) return b;
    return nullptr;
}
hmat_t<double> SomeD(size_t n, double phase)
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
double Tr(const hmat_t<double>& A, const hmat_t<double>& B)
{ double t=0.0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.rows();j++) t+=A(i,j)*B(i,j); return t; }

struct HartreeFixture
{
    Molecule mol=Water();
    std::unique_ptr<BasisSet::Real_BS> bs;
    const BasisSet::Real_DFT_OIBS* block=nullptr;
    BasisSet::FitQuadrature quad;
    std::shared_ptr<const BasisSet::FIT_SF_ABS<double>> delta;
    HartreeFixture()
    {
        using namespace BasisSet::Gaussian;
        bs.reset(Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian));
        block=FirstBlock(*bs);
        const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=40, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=17, .beckeOrder=2};
        quad=MakeMeshMaker(mol)->Make(mp, qcMesh::NoSymmetry());
        delta.reset(bs->CreateVxcFitBasisSet(&mol, BasisSet::VxcFit::Delta, &quad));
    }
};
//! E_H and V of fitter \a f at D (one block).
std::pair<double,hmat_t<double>> Fit(Fitting::FirstOrderDensityFit<double>& f, const BasisSet::Real_DFT_OIBS& b, const hmat_t<double>& D)
{
    f.DoFit(f.Forward(b).Forward(D));
    return {f.HartreeEnergy([]{return 0.0;}), f.Repulsion(b)};
}
} // namespace

TEST(DeltaProjectionHartree, ItsFockIsTheExactDerivativeOfItsEnergyAndTracesToTwiceIt)
{
    HartreeFixture s;
    ASSERT_NE(s.block, nullptr);
    auto f=Fitting::MakeDeltaProjectionDensityFitter(s.delta, s.quad.GetMesh());
    const size_t n=s.block->GetNumFunctions();
    const hmat_t<double> D=SomeD(n, 0.3), X=SomeD(n, 1.9);
    const auto [E, V]=Fit(*f, *s.block, D);
    EXPECT_GT(E, 0.0);
    EXPECT_NEAR(Tr(D,V), 2*E, 1e-11*std::fabs(E)) << "Tr(D V) = 2 E_H: the two legs have the same trace";
    // E_H is a quadratic form in D, so the central difference is exact up to roundoff at any h.
    const double h=1e-3;
    const double Ep=Fit(*f, *s.block, hmat_t<double>(D+h*X)).first, Em=Fit(*f, *s.block, hmat_t<double>(D-h*X)).first;
    const double fd=(Ep-Em)/(2*h), an=Tr(X,V);
    EXPECT_NEAR(an, fd, 1e-9*std::max(1.0,std::fabs(fd))) << "V = dE_H/dD (both legs)";
}

TEST(DeltaProjectionHartree, ItsChargeIsTheMeshsAndItsEnergyIsTheDirectSumRoutesToThatRoutesError)
{
    HartreeFixture s;
    ASSERT_NE(s.block, nullptr);
    const size_t n=s.block->GetNumFunctions();
    const hmat_t<double> D=SomeD(n, 0.3);
    auto f=Fitting::MakeDeltaProjectionDensityFitter(s.delta, s.quad.GetMesh());
    const double E=Fit(*f, *s.block, D).first;
    // the operator route on the same δ points (same charges q_g; its potential the direct sum over the other points, no self term): a
    // sanity bound only, the direct sum's own error on this coarse mesh is the larger one
    const auto cdelta=std::make_shared<const BasisSet::CoulombDeltaFit_IBS>(std::make_shared<BasisSet::BasisCache>(), s.delta,
                                  qcMesh::MakeDirectSumCoulomb(*s.quad.GetMesh()), sym_t(new UnitQN));
    auto p=Fitting::MakeDeltaDensityFitter(cdelta, s.quad.GetMesh());
    const double Ep=Fit(*p, *s.block, D).first;
    EXPECT_NEAR(f->Integral(), p->Integral(), 1e-12*std::fabs(p->Integral())) << "the same δ charges";
    EXPECT_NEAR(E, Ep, 5e-2*std::fabs(Ep)) << "two routes on one coarse mesh agree to the direct sum's error";
    std::cout << "[δ-projection] E_H=" << E << "  direct-sum E_H=" << Ep << "  diff=" << E-Ep << "  charge=" << f->Integral() << std::endl;
}

TEST(DeltaProjectionHartree, ItsProjectionsAreTwoPerPointAndAOneLegVectorIsRefused)
{
    HartreeFixture s;
    auto f=Fitting::MakeDeltaProjectionDensityFitter(s.delta, s.quad.GetMesh());
    EXPECT_EQ(f->NumCoefficients(), 2*s.quad.GetMesh()->size()) << "two projections per point: [rho ; v]";
    EXPECT_THROW((void)f->DoFit(Projections(rvec_t(s.quad.GetMesh()->size(), 0.0))), std::logic_error) << "one projection per point is not this fitter's";
}
