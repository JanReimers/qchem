// File: BasisSet/Gaussian/tests/M_LibCint.C
//
// Guard for the PG_LibCint matrix-delivery evaluator: it integrates the SAME Cartesian PG basis as
// PG_Cart_MnD but via the external libcint library.  Because it reproduces PGData's component ORDER and
// unit-self-overlap NORMALIZATION, every matrix / ERI it delivers must equal -- element for element -- what
// the M&D scalar IBS builds through its public accessors (the M&D path is itself oracle-validated to machine
// precision in M_PG_Oracle, so this makes libcint an independent second engine cross-checked against it).
//
// The 1E case deliberately includes d functions (LMax=2): libcint's Cartesian order (xx,xy,xz,yy,yz,zz)
// differs from PG's (xx,xy,yy,xz,yz,zz) for l>=2, so a wrong permutation would show up here (s,p alone
// would not catch it).
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <utility>
#include <vector>

import qchem.StructureData;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;     // PGData (the IBS IS-A one)
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;        // the evaluator under test
import qchem.BasisSet.Gaussian.PG_Cart;                           // Orbital_IBS (M&D reference)
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;
import qchem.BasisSet.PairProjection;
import qchem.BasisSet.Internal.ERI4;
import qchem.Structure;
import qchem.Types;
using namespace qchem;

using BasisSet::Gaussian::PG_Cart::Orbital_IBS;
using LibCintEval = BasisSet::Gaussian::Evaluators::PG_LibCint::NR_Evaluator;
using PGData      = BasisSet::Gaussian::Evaluators::PG_Cart_MnD::PGData;

// The shared water geometry (bohr) from the structure data -- D-STRUCTDATA; was a private Angstrom-valued copy.
static Molecule* MakeWater() { return new Molecule(qchem::StructureData::GetMolecule("H2O")); }

TEST(M_LibCint, matrix_1E_matches_scalar)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 2, h2o);                           // s + p + d (exercises the d reorder)
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), h2o);
    ASSERT_EQ(lc.size(), ibs.GetNumFunctions());

    const BasisSet::Orbital_1E_IBS<double>& bs1e = ibs;
    rsmat_t So=lc.OverlapMatrix(), Ko=lc.KineticMatrix(), Vo=lc.NuclearMatrix(h2o);
    for (size_t i=0;i<lc.size();++i)
        for (size_t j=i;j<lc.size();++j)
        {
            EXPECT_NEAR(So(i,j), bs1e.Overlap()(i,j),    1e-10) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(Ko(i,j), bs1e.Kinetic()(i,j),    1e-10) << "kinetic ("<<i<<","<<j<<")";
            EXPECT_NEAR(Vo(i,j), bs1e.Nuclear(h2o)(i,j), 1e-10) << "nuclear ("<<i<<","<<j<<")";
        }
}

// Cartesian g (LMax=4) through the M&D engine vs libcint (OpenWork §3 "l=4 (g) shells", 2026-10-21): the g row is
// the first one no CP2K/QE gate exercises, so libcint is its oracle before a g shell is trusted in a lattice run.
TEST(M_LibCint, matrix_1E_matches_scalar_g)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.2, 0.4};
    Orbital_IBS ibs(exps, 4, h2o);                           // s..g
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), h2o);
    ASSERT_EQ(lc.size(), ibs.GetNumFunctions());

    const BasisSet::Orbital_1E_IBS<double>& bs1e = ibs;
    rsmat_t So=lc.OverlapMatrix(), Ko=lc.KineticMatrix(), Vo=lc.NuclearMatrix(h2o);
    const auto S=bs1e.Overlap(), K=bs1e.Kinetic(), V=bs1e.Nuclear(h2o);
    for (size_t i=0;i<lc.size();++i)
        for (size_t j=i;j<lc.size();++j)
        {
            EXPECT_NEAR(So(i,j), S(i,j), 1e-10) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(Ko(i,j), K(i,j), 1e-9)  << "kinetic ("<<i<<","<<j<<")";
            EXPECT_NEAR(Vo(i,j), V(i,j), 1e-9)  << "nuclear ("<<i<<","<<j<<")";
        }
}

TEST(M_LibCint, matrix_3C_4C_match_scalar)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, h2o);                           // s + p (keeps the dense 4C tensor small)
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), h2o);

    // 4-centre (HF): partner = self.  Reference = the public, energy-validated scalar Direct/Exchange.
    const BasisSet::Orbital_ERI4_IBS<double>& hf = ibs;
    EXPECT_LT(fnorm(lc.DirectMatrix(lc),   hf.Direct(ibs)),   1e-10);
    EXPECT_LT(fnorm(lc.ExchangeMatrix(lc), hf.Exchange(ibs)), 1e-10);

    // 3-centre (DFT): a real Coulomb-fit basis from the orbital IBS; reference = public Overlap3C/Repulsion3C.
    auto  cdfit = ibs.CreateCDFitBasisSet(h2o,qcMesh::MeshParams{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3});
    auto vxcfit = ibs.CreateVxcFitBasisSet(h2o,qcMesh::MeshParams{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3});
    LibCintEval lcdfit(dynamic_cast<const PGData&>(*cdfit), h2o);
    LibCintEval lvxcfit(dynamic_cast<const PGData&>(*vxcfit), h2o);
    const BasisSet::Orbital_DFT_IBS<double>& dft = ibs;
    EXPECT_LT(fnorm(lc.OverlapThreeC_Matrix(lvxcfit),   dft.Overlap3C  (*vxcfit)),   1e-10);
    EXPECT_LT(fnorm(lc.RepulsionThreeC_Matrix(lcdfit), dft.Repulsion3C(* cdfit)), 1e-10);
    delete vxcfit;
    delete  cdfit;
}

// Spherical grid-eval (an HF run samples rho on a mesh for its D-MINDRHO residual, so the spherical libcint oracle needs
// it): the values and gradients, integrated on a mesh, must reproduce libcint's OWN spherical overlap and kinetic
// (<grad i|grad j> = <-nabla^2>) matrices -- that pins the cart->sph map, the order and the normalization together, with no
// convention-matching of our own.  ONE off-origin atom, so every product is a polynomial of degree <= 6 on one centre and
// Lebedev 35 is exact: the residual is radial quadrature only (on water the Becke partition of the cross-centre products
// leaves ~1e-6, which would hide a small convention error).  LMax=3: f is the first shell whose cart->sph map is not a
// relabelling of a few monomials.
TEST(M_LibCint, spherical_grid_eval_matches_matrices)
{
    Molecule* o = new Molecule;
    o->Insert(Site{8, rvec3_t(0.3, -0.2, 0.5)});     // off the origin: exercises the r - R_A shift
    rvec_t exps{1.0, 0.3};
    Orbital_IBS ibs(exps, 3, o);                             // s..f
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), o);
    lc.Init(o, true);                                        // re-pack as libcint-native spherical
    const size_t n=lc.size();
    ASSERT_EQ(n, (1+3+5+7)*exps.size());

    const qcMesh::Mesh m=o->CreateIntegrationMesh(qcMesh::MeshParams{.radial=qcMesh::RadialKind::MHL, .nRadial=120, .mhl_m=2,
                                     .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=35, .beckeOrder=3});
    rmat_t S(n,n,0.0), K(n,n,0.0);
    for (size_t g=0; g<m.size(); ++g)
    {
        const rvec_t     f =lc.SphericalValues  (m.Points()[g]);
        const rvec3vec_t df=lc.SphericalGradient(m.Points()[g]);
        const double w=m.Weights()[g];
        for (size_t i=0;i<n;++i)
            for (size_t j=i;j<n;++j)
            {
                S(i,j)+=w*f[i]*f[j];
                K(i,j)+=w*(df[i]*df[j]);
            }
    }
    const rsmat_t So=lc.OverlapMatrix(), Ko=lc.KineticMatrix();
    for (size_t i=0;i<n;++i)
        for (size_t j=i;j<n;++j)
        {
            EXPECT_NEAR(S(i,j), So(i,j), 1e-10) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(K(i,j), Ko(i,j), 1e-9)  << "kinetic ("<<i<<","<<j<<")";
        }
    delete o;
}

// THE 1E GEOMETRY GRADIENTS (the analytic force's Pulay, kinetic and nuclear terms; doc/GeometryOptimizationPlan.md §4-§5): libcint's
// dM/dR_A against a central finite difference of OUR engine's matrices (McMurchie-Davidson) at moved atoms -- two engines, two methods.
// The nuclear matrix moves twice: its basis centres and its nuclei.  Translation invariance: the gradients sum to zero over the atoms.
TEST(M_LibCint, matrix_1E_geometry_gradients_match_finite_differences_of_the_MnD_matrices)
{
    Molecule* h2o = MakeWater();
    const rvec_t exps{1.0, 0.25};
    const size_t L=2;                                           // s + p + d (the d reorder and its renormalization)
    Orbital_IBS ibs(exps, L, h2o);
    LibCintEval lc(dynamic_cast<const PGData&>(ibs), h2o);
    const auto dS=lc.OverlapGradient(), dK=lc.KineticGradient(), dV=lc.NuclearGradient();
    const size_t nAtoms=h2o->GetNumAtoms(), n=lc.size();
    ASSERT_EQ(dS.size(), nAtoms);
    std::vector<rvec3_t> R0; h2o->ForEachSite([&](int, const rvec3_t& r, bool){ R0.push_back(r); });
    const double h=1e-4;
    double largest=0.0;
    for (size_t A=0;A<nAtoms;++A)
        for (int c=0;c<3;++c)
        {
            auto at=[&](double s)
            {
                auto R=R0;
                (c==0 ? R[A].x : c==1 ? R[A].y : R[A].z)+=s;
                std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(R);
                Orbital_IBS b(exps, L, st.get());
                const BasisSet::Orbital_1E_IBS<double>& b1=b;
                return std::array<rsmat_t,3>{b1.Overlap(), b1.Kinetic(), b1.Nuclear(st.get())};
            };
            const auto p=at(h), m=at(-h);
            for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j)
            {
                const double fS=(p[0](i,j)-m[0](i,j))/(2*h), fK=(p[1](i,j)-m[1](i,j))/(2*h), fV=(p[2](i,j)-m[2](i,j))/(2*h);
                EXPECT_NEAR(dS[A][c](i,j), fS, 1e-7) << "dS/dR atom " << A << " axis " << c << " (" << i << "," << j << ")";
                EXPECT_NEAR(dK[A][c](i,j), fK, 2e-7) << "dK/dR atom " << A << " axis " << c << " (" << i << "," << j << ")";
                EXPECT_NEAR(dV[A][c](i,j), fV, 2e-7) << "dV/dR atom " << A << " axis " << c << " (" << i << "," << j << ")";
                largest=std::max({largest, std::fabs(fS), std::fabs(fV)});
            }
        }
    EXPECT_GT(largest, 0.1) << "the derivatives tested are not trivially zero";
    for (int c=0;c<3;++c)
        for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j)
        {
            double s=0.0, k=0.0, v=0.0;
            for (size_t A=0;A<nAtoms;++A) { s+=dS[A][c](i,j); k+=dK[A][c](i,j); v+=dV[A][c](i,j); }
            EXPECT_NEAR(s, 0.0, 1e-12); EXPECT_NEAR(k, 0.0, 1e-11); EXPECT_NEAR(v, 0.0, 1e-11);
        }
}

// THE DUNLAP FORCE'S INTEGRALS, contracted: Σ D_ij c_k ∂(ij|k)/∂R_A and Σ c_k c_l ∂(k|l)/∂R_A, against a central finite difference of our
// engine's (MnD) contractions at moved atoms, with a random symmetric D and random c (no SCF: the integrals alone).  The 2-centre
// metric itself is gated against the fit basis's own Repulsion first.
TEST(M_LibCint, coulomb_fit_geometry_gradients_match_finite_differences_of_the_MnD_contractions)
{
    Molecule* h2o = MakeWater();
    const rvec_t exps{1.0, 0.25};
    const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3};
    auto build=[&](const Structure* st)
    {
        auto ibs=std::make_unique<Orbital_IBS>(exps, 1, st);                 // s + p
        std::unique_ptr<BasisSet::rFIT_CD_ABS> fit(ibs->CreateCDFitBasisSet(st, mp));
        return std::pair{std::move(ibs), std::move(fit)};
    };
    auto [ibs, fit]=build(h2o);
    LibCintEval lc(dynamic_cast<const PGData&>(*ibs), h2o), lf(dynamic_cast<const PGData&>(*fit), h2o);
    const size_t n=lc.size(), nf=lf.size();
    const auto& metric=dynamic_cast<const BasisSet::FIT_CD_NonOrtho&>(*fit).Repulsion();
    {
        const rsmat_t M2=lf.RepulsionTwoC_Matrix();
        double worst=0.0; for (size_t k=0;k<nf;++k) for (size_t l=0;l<nf;++l) worst=std::max(worst, std::fabs(M2(k,l)-metric(k,l)));
        EXPECT_LT(worst, 1e-10) << "the 2-centre Coulomb metric";
    }

    std::mt19937 rng(7); std::uniform_real_distribution<double> u(-1.0, 1.0);
    rmat_t D(n, n); for (size_t i=0;i<n;++i) for (size_t j=0;j<=i;++j) D(i,j)=D(j,i)=u(rng);
    rvec_t c(nf);   for (size_t k=0;k<nf;++k) c[k]=u(rng);
    const auto g3=lc.RepulsionThreeC_Gradient(lf, D, c), g2=lf.RepulsionTwoC_Gradient(c);
    auto contractions=[&](const Structure* st)
    {
        auto [b, f]=build(st);
        const BasisSet::Orbital_DFT_IBS<double>& dft=*b;
        const auto P=dft.Repulsion3C(*f);
        double s3=0.0; for (size_t k=0;k<nf;++k) for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) s3+=D(i,j)*c[k]*P.dense[k](i,j);
        const auto& G=dynamic_cast<const BasisSet::FIT_CD_NonOrtho&>(*f).Repulsion();
        double s2=0.0; for (size_t k=0;k<nf;++k) for (size_t l=0;l<nf;++l) s2+=c[k]*c[l]*G(k,l);
        return std::pair{s3, s2};
    };
    std::vector<rvec3_t> R0; h2o->ForEachSite([&](int, const rvec3_t& r, bool){ R0.push_back(r); });
    const double h=1e-4;
    double largest=0.0;
    for (size_t A=0;A<R0.size();++A)
        for (int x=0;x<3;++x)
        {
            auto at=[&](double s){ auto R=R0; (x==0 ? R[A].x : x==1 ? R[A].y : R[A].z)+=s; std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(R); return contractions(st.get()); };
            const auto p=at(h), m=at(-h);
            const double f3=(p.first-m.first)/(2*h), f2=(p.second-m.second)/(2*h);
            const double a3 = x==0 ? g3[A].x : x==1 ? g3[A].y : g3[A].z, a2 = x==0 ? g2[A].x : x==1 ? g2[A].y : g2[A].z;
            EXPECT_NEAR(a3, f3, 1e-6*std::max(1.0, std::fabs(f3))) << "Σ D c ∂(ij|k)/∂R, atom " << A << " axis " << x;
            EXPECT_NEAR(a2, f2, 1e-6*std::max(1.0, std::fabs(f2))) << "Σ c c ∂(k|l)/∂R, atom " << A << " axis " << x;
            largest=std::max({largest, std::fabs(f3), std::fabs(f2)});
        }
    EXPECT_GT(largest, 0.1) << "the contractions move with the atoms";
    rvec3_t s3(0,0,0), s2(0,0,0); for (size_t A=0;A<R0.size();++A) { s3+=g3[A]; s2+=g2[A]; }
    EXPECT_LT(norm(s3), 1e-9); EXPECT_LT(norm(s2), 1e-9);
}
