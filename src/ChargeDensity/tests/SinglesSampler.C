// File: ChargeDensity/tests/SinglesSampler.C  The molecular (real-axis) δ sampler's ρ FORWARD (lane H, 2026-11-12).
//
// Until 2026-11-12 the real δ fitter vended no forward face, so a molecule's XC sampler took the POINTWISE field projection
// for a density-matrix density: every basis function evaluated at every mesh point on every sampling (45 % of an H₂O₂ LDA
// single point, serial; OpenWork §4).  The fitter now vends \c PairProjector<double>, so ρ is the GEMM
// \f$\rho_g=\sum_{ij}\Phi_{gi}D_{ij}\Phi_{gj}\f$ over the Φ table the H_xc adjoint already holds.  Same ρ[D] (pin 28), a different
// summation order, so the gate is ELEMENTWISE to roundoff against the pointwise route it replaced -- for ρ, and for the H_xc it
// feeds -- and it also pins that the capability is THERE, since a missing face falls back to the pointwise route silently.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <memory>

import qchem.ChargeDensity;
import qchem.ChargeDensity.DensitySampler;      // MakeDensitySampler (the real overload)
import qchem.ChargeDensity.Factory;             // IrrepCD_Factory<double>
import qchem.BasisSet;                          // Real_BS, MakeXCQuadrature, VxcFit, rFIT_SF_ABS
import qchem.BasisSet.Gaussian.Point.Factory;   // Gaussian::Factory
import qchem.Fitting.FunctionFitter;            // MakeDeltaScalarFitter, FieldProjector, PairProjector
import qchem.Mesh;                              // qcMesh::MeshParams
import qchem.Structure;
import qchem.Symmetry.Irrep;
import qchem.Types;
import qchem.Blaze;
using namespace qchem;
using namespace qchem::ChargeDensity;

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

//! A deterministic PSD density matrix D = L Lᵀ (rank 5, ~ five doubly occupied orbitals), so ρ ≥ 0 as a real one is.
hmat_t<double> SomePSD_D(size_t n)
{
    const size_t r=5;
    hmat_t<double> D(n);
    for (size_t i=0; i<n; i++) for (size_t j=0; j<=i; j++)
    {
        double s=0.0;
        for (size_t m=0; m<r; m++) s+=std::sin(0.3+1.1*i+0.7*m)*std::sin(0.3+1.1*j+0.7*m);
        D(i,j)=2.0*s/double(n);
    }
    return D;
}

struct Rig
{
    Molecule                                      mol=H2O();
    std::unique_ptr<BasisSet::Real_BS>            bs;
    const BasisSet::Orbital_1E_IBS<double>*       block=nullptr;
    BasisSet::FitQuadrature                       quad;
    std::shared_ptr<const BasisSet::rFIT_SF_ABS>  fb;
    std::unique_ptr<tDM_CD<double>>               cd;
    Rig()
    {
        bs.reset(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, &mol,
                                             BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        block=dynamic_cast<const BasisSet::Orbital_1E_IBS<double>*>((*bs)[0]);
        const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=40, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=17, .beckeOrder=2};
        quad=BasisSet::MakeXCQuadrature(mol, mp, nullptr);
        fb.reset(bs->CreateVxcFitBasisSet(&mol, mp, BasisSet::VxcFit::Delta, &quad));
        cd.reset(IrrepCD_Factory<double>(SomePSD_D(block->GetNumFunctions()), block, bs->GetIrreps(Spin::None)[0]));
    }
    //! The route the sampler used before: the fitter's FIELD projection of the density, point by point.
    rvec_t PointwiseRho() const
    {
        std::shared_ptr<const BasisSet::rFIT_SF_ABS> f=fb;
        const auto fitter=Fitting::MakeDeltaScalarFitter(f);
        return static_cast<const Fitting::FieldProjector&>(*fitter).Project(*cd);
    }
};

double MaxAbs(const rvec_t& v) {double m=0.0; for (double x : v) m=std::max(m,std::fabs(x)); return m;}
} //anon

TEST(SinglesSampler, TheRealDeltaFitterVendsTheForward)
{
    Rig rig;
    std::shared_ptr<const BasisSet::rFIT_SF_ABS> f=rig.fb;
    const auto fitter=Fitting::MakeDeltaScalarFitter(f);
    EXPECT_NE(dynamic_cast<const BasisSet::PairProjector<double>*>(fitter.get()), nullptr)
        << "the molecular δ fitter must vend PairProjector<double>: without it the sampler falls back to the pointwise raster SILENTLY";
}

// ⚠ NOT at the zero-weight points (measured 2026-11-12: 330 of 13164 here, Becke weights that underflow to 0).  There the
// pointwise route's δ fit, \f$c_g=\langle\delta_g|f\rangle/\langle\delta_g|\delta_g\rangle\f$, is 0/0 and comes out 0, while the
// forward returns the true ρ(r_g).  Every consumer integrates with w (E_xc, H_xc, the Δρ residual, the site moments), so the
// two are the same functional; the gate says so exactly: equal to roundoff where w > 0, and every difference sits at w = 0.
TEST(SinglesSampler, TheForwardRhoIsThePointwiseOneElementwise)
{
    Rig rig;
    const auto sampler=MakeDensitySampler(rig.fb, rig.quad);
    const rvec_t rho=sampler->Rho(rig.cd.get());
    const rvec_t ref=rig.PointwiseRho();
    const rvec_t w=rig.fb->Charge();   // a δ basis's own integrals are its weights
    ASSERT_EQ(rho.size(), ref.size());
    ASSERT_EQ(w.size(), ref.size());
    const double scale=MaxAbs(rho);
    ASSERT_GT(scale, 1e-3) << "a vacuous ρ would pass anything";
    double worst=0.0;
    size_t nZeroW=0, nDiffAtZeroW=0;
    for (size_t g=0; g<ref.size(); g++)
    {
        if (w[g]!=0.0) {worst=std::max(worst, std::fabs(rho[g]-ref[g])); continue;}
        nZeroW++;
        if (rho[g]!=ref[g]) nDiffAtZeroW++;
        EXPECT_EQ(ref[g], 0.0) << "the pointwise δ fit at a zero-weight point is 0/0 -> 0 (point " << g << ")";
    }
    std::printf("  [singles ρ] %zu points (%zu at w = 0, %zu of those differ), max|ρ| %.3e, max|Δρ| at w > 0 %.3e\n",
                ref.size(), nZeroW, nDiffAtZeroW, scale, worst);
    EXPECT_LE(worst, 1e-13*scale) << "the forward and the pointwise ρ are the same ρ[D] to roundoff wherever the weight is not 0";
}

TEST(SinglesSampler, TheHxcFromTheForwardRhoIsThePointwiseOnesElementwise)
{
    Rig rig;
    const auto sampler=MakeDensitySampler(rig.fb, rig.quad);
    auto v=[](const rvec_t& rho)   // a Slater-shaped potential: nonlinear and pointwise, as v_xc is
    {
        rvec_t out(rho.size());
        for (size_t g=0; g<rho.size(); g++) out[g]=-std::cbrt(std::max(rho[g],0.0));
        return out;
    };
    const hmat_t<double> H   =sampler->Matrix(rig.block, v(rvec_t(sampler->Rho(rig.cd.get()))));
    const hmat_t<double> Href=sampler->Matrix(rig.block, v(rig.PointwiseRho()));
    double scale=0.0, worst=0.0;
    for (size_t i=0; i<H.rows(); i++) for (size_t j=0; j<H.columns(); j++)
    {
        scale=std::max(scale, std::fabs(Href(i,j)));
        worst=std::max(worst, std::fabs(H(i,j)-Href(i,j)));
    }
    std::printf("  [singles H_xc] n %zu, max|H| %.3e, max|ΔH| %.3e\n", H.rows(), scale, worst);
    ASSERT_GT(scale, 1e-3);
    EXPECT_LE(worst, 1e-12*scale);
}
