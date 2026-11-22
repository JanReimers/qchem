// File: Hamiltonian/tests/GGA_Periodic.C  The PERIODIC point route's density gradient (GGA stage 3): the Bloch and real TRIM blocks.
//
// The δ basis's gradient 3-centre overlap on the complex fit axis, through the periodic density-gradient sampler: ∇ρ at Becke points is the
// density's own (a Richardson central difference of its POINTWISE ρ, the block's image sum evaluated directly -- an independent route), and
// the w half of H is the exact adjoint of it (Tr(D H_w) = Σ_g w_g w_g·∇ρ_g).  The Bloch block sits at k = (1/3, 0, 0), NOT a TRIM point, so
// its Φ and D are genuinely complex (a TRIM k's Bloch block has real Φ and would not test the conjugations); the real TRIM block (Γ, the 3c-3
// mixed shape) goes through the periodic face's real gather.
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>

import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;    // tGPW_IBS<T>
import qchem.BasisSet.DeltaFit_IBS;
import qchem.Symmetry.Factory;                     // BlochFactory
import qchem.ChargeDensity;
import qchem.ChargeDensity.DensitySampler;         // MakeDensityGradientSampler (the periodic overload)
import qchem.ChargeDensity.Imp.IrrepCD;            // PeriodicIrrepCD
import qchem.CompositeCD;                          // tComposite_CD<dcmplx>
import qchem.Mesh;
import qchem.Mesh.XCPolicy;                       // BeckeXCParams
import qchem.Types;
import qchem.Blaze;
import qchem.Structure.MeshFactory;   // MakeMeshMaker (the GPW block demands its cell's mesh maker, V-MESH-FACTORY §4.4)
using namespace qchem;

namespace
{
template <class U> U Entry(size_t i, size_t m)
{
    const double a=std::sin(1.0+0.7*double(i)+1.3*double(m))/(1.0+0.1*double(i));
    if constexpr (std::is_same_v<U,dcmplx>) return dcmplx(a, std::cos(0.4+0.9*double(i)-0.6*double(m))*0.5);
    else                                    return a;
}
//! A deterministic Hermitian PSD D = L L^dagger (rank 3).
template <class U> hmat_t<U> SomeD(size_t n)
{
    hmat_t<U> D(n);
    for (size_t i=0; i<n; i++) for (size_t j=i; j<n; j++)
    {
        U s{};
        for (size_t m=0; m<3; m++)
            if constexpr (std::is_same_v<U,dcmplx>) s+=Entry<U>(i,m)*std::conj(Entry<U>(j,m));
            else                                    s+=Entry<U>(i,m)*Entry<U>(j,m);
        D(i,j)=s;
    }
    return D;
}

struct Rig
{
    UnitCell cell;
    std::shared_ptr<const Structure> st;
    std::shared_ptr<const BasisSet::Real_BS> mol;
    std::unique_ptr<BasisSet::Gaussian::tGPW_IBS<dcmplx>> cx;   // Bloch, k = (1/3, 0, 0)
    std::unique_ptr<BasisSet::Gaussian::tGPW_IBS<double>> re;   // real TRIM, Gamma
    BasisSet::FitQuadrature quad;
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fb;
    std::shared_ptr<const ChargeDensity::DensityGradientSampler> sampler;
    Rig() : cell(8.0)
    {
        cell.AddAtom(14, {0.5,0.5,0.5});
        st=cell.Clone();
        namespace M = qchem::BasisSet::Gaussian;
        mol.reset(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
        cx=std::make_unique<M::tGPW_IBS<dcmplx>>(std::make_shared<BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(3,1,1), ivec3_t(1,0,0), mol, 20.0);
        re=std::make_unique<M::tGPW_IBS<double>>(std::make_shared<BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, 20.0);
        quad=BasisSet::MakeXCQuadrature(*st, qcMesh::BeckeXCParams(30, -1.0, 11), nullptr);
        fb=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<BasisSet::BasisCache>(), quad,
                                                                  Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
        sampler=ChargeDensity::MakeDensityGradientSampler(fb, quad);
    }
    template <class U> std::unique_ptr<ChargeDensity::cDM_CD> Density(const hmat_t<U>& D) const
    {
        using namespace qchem::ChargeDensity;
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        if constexpr (std::is_same_v<U,dcmplx>)
            cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D, cx.get(), cx->GetIrrep(Spin::None))), cx->GetIrrep(Spin::None));
        else
            cd->Insert(std::unique_ptr<tDM_CD<double>>(new PeriodicIrrepCD<double>(D, re.get(), re->GetIrrep(Spin::None))), re->GetIrrep(Spin::None));
        return cd;
    }
};

//! ∇ρ through the gradient tensor vs a Richardson central difference of the density's own pointwise ρ, on a sample of the mesh.
void ExpectGradientIsTheDensitysOwn(const Rig& rig, const ChargeDensity::cDM_CD& cd, const char* what)
{
    const rvec3vec_t& gu=rig.sampler->GradRhoPol(&cd, Spin::Up);
    const rvec3vec_t& gd=rig.sampler->GradRhoPol(&cd, Spin::Down);
    const rvec_t rho=rig.sampler->Rho(&cd);
    const rvec3vec_t& pts=rig.quad.GetMesh()->Points();
    const double h=1e-5;
    double worst=0.0, scale=0.0;
    size_t n=0;
    for (size_t g=0; g<pts.size(); g+=29)
    {
        if (rho[g]<1e-4) continue;
        const rvec3_t e[3]={rvec3_t(h,0,0), rvec3_t(0,h,0), rvec3_t(0,0,h)};
        const rvec3_t a=gu[g]+gd[g];
        const double an[3]={a.x, a.y, a.z};
        for (int k=0; k<3; k++)
        {
            const double d1=(cd(pts[g]+e[k])-cd(pts[g]-e[k]))/(2*h), d2=(cd(pts[g]+2.0*e[k])-cd(pts[g]-2.0*e[k]))/(4*h);
            const double fd=(4*d1-d2)/3;
            worst=std::max(worst, std::fabs(an[k]-fd)); scale=std::max(scale, std::fabs(fd));
        }
        n++;
    }
    std::printf("  [periodic grad rho, %s] %zu points, max|grad rho| %.3e, max|analytic - Richardson| %.3e\n", what, n, scale, worst);
    ASSERT_GT(n, 20u) << what;
    ASSERT_GT(scale, 1e-3) << what << ": a vacuous gradient would pass anything";
    EXPECT_LE(worst, 1e-7*scale) << what << ": grad rho through the gradient tensor is the density's own";
}

//! Tr(D H_w) = sum_g w_g w_g . grad rho_g, with H_w the w half of the {v, w} gather (v = 0).
template <class U, class H> void ExpectWGatherIsTheAdjoint(const Rig& rig, const hmat_t<U>& D, const ChargeDensity::cDM_CD& cd, H gather,
                                                             const char* what)
{
    const size_t N=rig.sampler->NumPoints();
    const rvec3vec_t& gu=rig.sampler->GradRhoPol(&cd, Spin::Up);
    const rvec3vec_t& gd=rig.sampler->GradRhoPol(&cd, Spin::Down);
    rvec3vec_t w(N);
    for (size_t g=0; g<N; g++) w[g]=rvec3_t(std::sin(0.1+0.37*g), std::cos(0.2+0.11*g), std::sin(0.5+0.07*g));
    const hmat_t<U> Hw=gather(rvec_t(N, 0.0), w);
    U tr{};
    for (size_t i=0; i<Hw.rows(); i++) for (size_t j=0; j<Hw.columns(); j++) tr+=D(i,j)*Hw(j,i);
    const vec_t<dcmplx> q=rig.fb->Charge();   // a delta basis's own integrals are its weights
    double ref=0.0, mag=0.0;
    for (size_t g=0; g<N; g++)
    {
        const rvec3_t gr=gu[g]+gd[g];
        const double t=std::real(q[g])*(w[g]*gr);
        ref+=t; mag+=std::fabs(t);
    }
    std::printf("  [periodic {v,w}, %s] Tr(D H_w) %.12e (imag %.1e), sum_g w w.grad rho %.12e\n", what, std::real(tr), std::imag(dcmplx(tr)), ref);
    ASSERT_GT(mag, 1e-3) << what;
    EXPECT_NEAR(std::real(tr), ref, 1e-11*mag) << what << ": the w gather is the exact transpose of the gradient";
    EXPECT_NEAR(std::imag(dcmplx(tr)), 0.0, 1e-11*mag) << what << ": Tr(D H) of two Hermitian matrices is real";
}
} //anon

TEST(GGAPeriodic, TheBlochBlocksGradientIsTheDensitysOwn)
{
    Rig rig;
    const hmat_t<dcmplx> D=SomeD<dcmplx>(rig.cx->GetNumFunctions());
    const auto cd=rig.template Density<dcmplx>(D);
    ExpectGradientIsTheDensitysOwn(rig, *cd, "Bloch k=(1/3,0,0)");
}

TEST(GGAPeriodic, TheBlochBlocksWGatherIsTheAdjoint)
{
    Rig rig;
    const hmat_t<dcmplx> D=SomeD<dcmplx>(rig.cx->GetNumFunctions());
    const auto cd=rig.template Density<dcmplx>(D);
    ExpectWGatherIsTheAdjoint<dcmplx>(rig, D, *cd, [&](const rvec_t& v, const rvec3vec_t& w){return rig.sampler->Matrix(rig.cx.get(), v, w);},
                                      "Bloch k=(1/3,0,0)");
}

TEST(GGAPeriodic, TheRealTRIMBlocksGradientAndWGather)
{
    Rig rig;
    const hmat_t<double> D=SomeD<double>(rig.re->GetNumFunctions());
    const auto cd=rig.template Density<double>(D);
    ExpectGradientIsTheDensitysOwn(rig, *cd, "real TRIM Gamma");
    ExpectWGatherIsTheAdjoint<double>(rig, D, *cd, [&](const rvec_t& v, const rvec3vec_t& w){return rig.sampler->Matrix(rig.re.get(), v, w);},
                                      "real TRIM Gamma");
}
