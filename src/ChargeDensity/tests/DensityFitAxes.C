// File: ChargeDensity/tests/DensityFitAxes.C  The constraint axis of the Gaussian (Coulomb-metric) density fit, V-FIT-AXES step 1 (user
// 2026-11-25: "{constrained, unconstrained} and {overlap,coulomb/dunlap} should be orthogonal axes in code"), at a FIXED density on water's
// first orbital block with its default Coulomb fit basis:
// - UNCONSTRAINED, the fit is the Coulomb least-squares one, so (ρ|ρ̃) = (ρ̃|ρ̃): the robust energy 2 EeeFit − EeeFitFit equals the plain
//   EeeFit to roundoff; CONSTRAINED, the fit's charge is N exactly and the two energies part;
// - either way the Fock Σ c_k (ij|f_k) is the exact derivative of the robust energy (the fit is stationary in c, constrained or not), so the
//   central difference of the quadratic E along X is Tr(X V) to roundoff.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <memory>

import qchem.BasisSet;                          // Real_BS, rFIT_CD_ABS
import qchem.BasisSet.Gaussian.Point.Factory;   // Gaussian::Factory
import qchem.Fitting.CoulombOperator;           // MakeDunlapDensityFitter, DensityFitConstraint
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;
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

struct Energies { double eFit, eFitFit, robust; hmat_t<double> V; double charge; };
//! Fit D (one block) with charge target \a N; the robust energy as Vee_FitDunlap forms it.
Energies Fit(Fitting::DunlapDensityFit<double>& f, const BasisSet::Real_DFT_OIBS& b, const hmat_t<double>& D, double N)
{
    f.DoFit(f.Forward(b).Forward(D), N);
    hmat_t<double> V=f.Repulsion(b);
    const double eFit=0.5*Tr(D,V), eFitFit=f.SelfRepulsion();
    return {eFit, eFitFit, 2*eFit-eFitFit, V, f.Integral()};
}
struct AxesFixture
{
    Molecule mol=Water();
    std::unique_ptr<BasisSet::Real_BS> bs;
    const BasisSet::Real_DFT_OIBS* block=nullptr;
    std::shared_ptr<const BasisSet::rFIT_CD_ABS> fit;
    AxesFixture()
    {
        using namespace BasisSet::Gaussian;
        bs.reset(Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian));
        block=FirstBlock(*bs);
        fit.reset(bs->CreateCDFitBasisSet(&mol));
    }
};
} // namespace

TEST(DensityFitAxes, UnconstrainedTheRobustEnergyIsThePlainOneAndConstrainedTheChargeIsN)
{
    AxesFixture s;
    ASSERT_NE(s.block, nullptr);
    const hmat_t<double> D=SomeD(s.block->GetNumFunctions(), 0.3);
    const double N=10.0;
    auto plain=Fitting::MakeDunlapDensityFitter(s.fit, Fitting::DensityFitConstraint::None);
    auto dunlap=Fitting::MakeDunlapDensityFitter(s.fit, Fitting::DensityFitConstraint::Charge);
    const Energies u=Fit(*plain, *s.block, D, N), c=Fit(*dunlap, *s.block, D, N);
    EXPECT_NEAR(u.eFit, u.eFitFit, 1e-10*std::fabs(u.eFit)) << "unconstrained Coulomb least squares: (ρ|ρ̃) = (ρ̃|ρ̃)";
    EXPECT_NEAR(c.charge, N, 1e-10*N) << "the constraint holds the fit's charge to N";
    EXPECT_GT(std::fabs(u.charge-N), 1e-6) << "the gate is vacuous unless the unconstrained charge misses N";
    EXPECT_GT(std::fabs(c.eFit-c.eFitFit), 1e-8) << "constrained, the plain and robust energies part";
    EXPECT_GE(u.robust, c.robust-1e-12) << "the robust energy is a maximum over c at fixed D (−½‖ρ−ρ̃‖² ≤ 0): constraining c cannot raise it";
}

TEST(DensityFitAxes, EitherWayTheFockIsTheExactDerivativeOfTheRobustEnergy)
{
    AxesFixture s;
    ASSERT_NE(s.block, nullptr);
    const size_t n=s.block->GetNumFunctions();
    const hmat_t<double> D=SomeD(n, 0.3), X=SomeD(n, 1.9);
    const double N=10.0, h=1e-3;
    for (auto con : {Fitting::DensityFitConstraint::None, Fitting::DensityFitConstraint::Charge})
    {
        auto f=Fitting::MakeDunlapDensityFitter(s.fit, con);
        const Energies e=Fit(*f, *s.block, D, N);
        const double Ep=Fit(*f, *s.block, hmat_t<double>(D+h*X), N).robust, Em=Fit(*f, *s.block, hmat_t<double>(D-h*X), N).robust;
        const double fd=(Ep-Em)/(2*h), an=Tr(X,e.V);
        EXPECT_NEAR(an, fd, 1e-9*std::max(1.0,std::fabs(fd))) << (con==Fitting::DensityFitConstraint::None ? "unconstrained" : "constrained");
    }
}
