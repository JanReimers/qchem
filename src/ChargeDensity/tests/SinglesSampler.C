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
import qchem.BasisSet;                          // Real_BS, VxcFit, rFIT_SF_ABS
import qchem.BasisSet.Gaussian.Point.Factory;   // Gaussian::Factory
import qchem.Fitting.FunctionFitter;            // MakeDeltaScalarFitter, FieldProjector, PairProjector
import qchem.Mesh;                              // qcMesh::MeshParams
import qchem.Structure;
import qchem.Symmetry.Irrep;
import qchem.Types;
import qchem.Blaze;
import qchem.Structure.MeshFactory;   // MakeMeshMaker: the quadrature a delta fit is built over (V-MESH-FACTORY step 4)
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
        quad=MakeMeshMaker(mol)->Make(mp, qcMesh::NoSymmetry());
        fb.reset(bs->CreateVxcFitBasisSet(&mol, BasisSet::VxcFit::Delta, &quad));
        cd.reset(IrrepCD_Factory<double>(SomePSD_D(block->GetNumFunctions()), block, bs->GetIrreps(Spin::None)[0], qchem::ChargeDensity::RhoRoute::PivotedCholesky));
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

// ---- THE DENSITY-GRADIENT STRATEGY (GGA stage 2, the molecular point route) -----------------------------------------------------------
// ∇ρ at the XC points through the δ basis's gradient 3-centre overlap, and the {v, w} gather back.  Three claims: ∇ρ is the density's own
// (a central difference of its POINTWISE ρ, an independent route through the basis's values alone), the w half of H is the exact adjoint of
// that ∇ρ (a dot-product identity, which is H = dE/dD for the gradient part), and ρ itself is the LDA strategy's to the bit.

TEST(SinglesGradientSampler, TheGradientIsTheDensitysOwnByCentralDifference)
{
    Rig rig;
    const auto sampler=MakeDensityGradientSampler(rig.fb, rig.quad);
    const rvec3vec_t& gUp=sampler->GradRhoPol(rig.cd.get(), Spin::Up);   // a spin-agnostic density: ∇ρ/2 per channel
    const rvec3vec_t& gDn=sampler->GradRhoPol(rig.cd.get(), Spin::Down);
    const rvec3vec_t& pts=rig.quad.GetMesh()->Points();
    ASSERT_EQ(gUp.size(), pts.size());
    const rvec_t rho=sampler->Rho(rig.cd.get());
    const double h=1e-5;   // Richardson below: the O(h²) term of a central difference is ~1e-5 relative at h = 1e-4 on the O 1s core
    size_t n=0;
    double worst=0.0, scale=0.0;
    for (size_t g=0; g<pts.size(); g+=53)
    {
        if (rho[g]<1e-4) continue;   // a point where the density matters (the difference loses digits in the tail)
        const rvec3_t r=pts[g];
        const rvec3_t e[3]={rvec3_t(h,0,0), rvec3_t(0,h,0), rvec3_t(0,0,h)};
        const rvec3_t a=gUp[g]+gDn[g];
        const double an[3]={a.x, a.y, a.z};
        for (int k=0; k<3; k++)
        {
            const double d1=((*rig.cd)(r+e[k])-(*rig.cd)(r-e[k]))/(2.0*h);
            const double d2=((*rig.cd)(r+2.0*e[k])-(*rig.cd)(r-2.0*e[k]))/(4.0*h);
            const double fd=(4.0*d1-d2)/3.0;   // Richardson: O(h⁴)
            worst=std::max(worst, std::fabs(an[k]-fd));
            scale=std::max(scale, std::fabs(fd));
        }
        EXPECT_EQ(gUp[g].x, gDn[g].x) << "a spin-agnostic density splits its gradient evenly";
        n++;
    }
    std::printf("  [singles ∇ρ] %zu points checked, max|∇ρ| %.3e, max|analytic - central difference| %.3e\n", n, scale, worst);
    ASSERT_GT(n, 50u);
    ASSERT_GT(scale, 1e-2) << "a vacuous gradient would pass anything";
    EXPECT_LE(worst, 1e-8*scale) << "∇ρ through the gradient tensor is the density's own (Richardson central difference, h = 1e-5)";
}

TEST(SinglesGradientSampler, TheGradientHalfOfHIsTheAdjointOfTheGradient)
{
    Rig rig;
    const auto sampler=MakeDensityGradientSampler(rig.fb, rig.quad);
    const size_t N=sampler->NumPoints();
    const rvec3vec_t& gUp=sampler->GradRhoPol(rig.cd.get(), Spin::Up);
    const rvec3vec_t& gDn=sampler->GradRhoPol(rig.cd.get(), Spin::Down);
    rvec3vec_t w(N);
    for (size_t g=0; g<N; g++) w[g]=rvec3_t(std::sin(0.1+0.37*g), std::cos(0.2+0.11*g), std::sin(0.5+0.07*g));
    const rvec_t zero(N, 0.0), w8=rig.fb->Charge();
    const hmat_t<double> H=sampler->Matrix(rig.block, zero, w);
    const hmat_t<double> D=SomePSD_D(rig.block->GetNumFunctions());
    double trDH=0.0;
    for (size_t i=0; i<H.rows(); i++) for (size_t j=0; j<H.columns(); j++) trDH+=D(i,j)*H(i,j);
    double ref=0.0, mag=0.0;
    for (size_t g=0; g<N; g++)
    {
        const rvec3_t gr=gUp[g]+gDn[g];
        const double t=w8[g]*(w[g].x*gr.x+w[g].y*gr.y+w[g].z*gr.z);
        ref+=t; mag+=std::fabs(t);
    }
    std::printf("  [singles {v,w}] Tr(D H_w) %.12e, sum_g w_g w.grad rho %.12e, sum|terms| %.3e\n", trDH, ref, mag);
    ASSERT_GT(mag, 1e-2);
    EXPECT_NEAR(trDH, ref, 1e-12*mag) << "Tr(D H_w) = sum_g w_g w_g . grad rho_g: the w gather is the exact transpose of the gradient";
}

TEST(SinglesGradientSampler, ItsRhoAndItsZeroWGatherAreTheLDAStrategys)
{
    Rig rig;
    const auto lda=MakeDensitySampler(rig.fb, rig.quad);
    const auto gga=MakeDensityGradientSampler(rig.fb, rig.quad);
    const rvec_t r0=lda->Rho(rig.cd.get()), r1=gga->Rho(rig.cd.get());
    ASSERT_EQ(r0.size(), r1.size());
    double worst=0.0;
    for (size_t g=0; g<r0.size(); g++) worst=std::max(worst, std::fabs(r0[g]-r1[g]));
    EXPECT_LE(worst, 1e-15*MaxAbs(r0)) << "ρ = ρ↑ + ρ↓ off the spin cache is the LDA strategy's ρ (ρ/2 + ρ/2)";
    rvec_t v(r0.size());
    for (size_t g=0; g<v.size(); g++) v[g]=-std::cbrt(std::max(r0[g],0.0));
    const hmat_t<double> H0=lda->Matrix(rig.block, v);
    const hmat_t<double> H1=gga->Matrix(rig.block, v, rvec3vec_t(v.size(), rvec3_t(0,0,0)));
    double dh=0.0, sc=0.0;
    for (size_t i=0; i<H0.rows(); i++) for (size_t j=0; j<H0.columns(); j++) {dh=std::max(dh, std::fabs(H0(i,j)-H1(i,j))); sc=std::max(sc, std::fabs(H0(i,j)));}
    EXPECT_LE(dh, 1e-14*sc) << "a zero w adds nothing: the {v, 0} gather is the LDA one";
}
