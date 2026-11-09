// File: BasisSet/Lattice/tests/GPW_ForceUT.C  THE GPW HARTREE FORCE AT FIXED D, END TO END THROUGH THE TENSORS (doc/GeometryOptimizationPlan.md
// §13.2; lane B 2026-11-09).  The first-order fitted Hartree energy on GPW's plane-wave fit, E_H = ½ Ω Σ_G k(G) |ρ̃(G)|² with
// ρ̃ = Contract(Overlap3C, D) (the GSpaceDensityFitter's energy), is quadratic in D; at fixed D its R-derivative is the pullback of the
// physics round at the fixed fit, ContractAdjointPullback(Repulsion3C, ρ̃, D) -- the collocation's ladder, level restriction, Coulomb
// kernel and the gather's jet walk all on the path.  Gated against the central difference of E_H with the GPW block rebuilt on displaced
// atoms, at a NON-TRIM k (complex phases).  The kernel itself is gated at the LatticeCollocation face in M_CollocationPullback.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <vector>

import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet;                                // Real_BS
import qchem.BasisSet.Orbital_1E_IBS;                 // Complex_OIBS
import qchem.BasisSet.Gaussian.Point.Factory;         // Gaussian::Factory
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;       // GPW_IBS (the block under test)
import qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator; // its tensors (tests may cheat-import internals)
import qchem.BasisSet.GMap;                           // PairProjection<dcmplx>, ΔG_Map, Contract / ContractAdjoint(Pullback)
import qchem.Vector3D;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Complex_OIBS;
using BasisSet::Gaussian::GPW_IBS;
using BasisSet::Gaussian::GPW_Evaluator;
using qchem::BasisSet::Gaussian::BasisSetData;

namespace
{
const ivec3_t kMesh(4,4,4), kIndex(1,0,2);   // k = (1/4, 0, 1/2): not a TRIM point

struct Si2 { rvec3_t f0{0.1,0.05,0.12}, f1{0.35,0.3,0.37}; };   // off the origin: no corner atom to re-wrap under a displacement
UnitCell CellOf(size_t moved=size_t(-1), const rvec3_t& dr=rvec3_t(0,0,0))
{
    FCCUnitCell c(10.26);
    const Si2 s;
    const rvec3_t f[2]={s.f0, s.f1};
    for (size_t a=0;a<2;a++) c.AddAtom(14, a==moved ? c.ToFractional(c.ToCartesian(f[a])+dr) : f[a]);
    return c;
}
std::shared_ptr<const Real_BS> Mol(const Structure& st)
{
    return std::shared_ptr<const Real_BS>(BasisSet::Gaussian::Factory(BasisSetData::SIPP, &st, BasisSet::Gaussian::Engine::MnD,
                                                                      BasisSet::Gaussian::Angular::Cartesian));
}
std::unique_ptr<GPW_IBS> Block(const UnitCell& cell)
{
    return std::make_unique<GPW_IBS>(cell, kMesh, kIndex, Mol(cell), /*densityEcut*/6.0);
}
chmat_t RandomHermitian(size_t n)
{
    unsigned long long r=0x9E3779B97F4A7C15ULL;
    auto next=[&]{ r^=r<<13; r^=r>>7; r^=r<<17; return double(r%2000)/1000.0-1.0; };
    chmat_t D(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)= i==j ? dcmplx(0.2+0.1*next(),0.0) : 0.03*dcmplx(next(),next());
    return D;
}
//! The GSpaceDensityFitter's energy of D on this block: ½ Ω Σ k |ρ̃|², ρ̃ the overlap projection.
double EH(const GPW_IBS& gpw, const chmat_t& D)
{
    const GPW_Evaluator& ev=gpw;
    const auto& grid=ev.DensityGrid();
    const ΔG_Map rho=Contract(ev.Overlap3CTensor(), D);
    double e=0.0;
    for (const auto& [dm,r] : rho) if (const double k=grid.Recip().CoulombKernel(dm); k>0.0) e+=k*std::norm(r);
    return 0.5*grid.Volume()*e;
}
} // namespace

TEST(GPW_Force, TheHartreeForceAtFixedDIsTheCoulombTensorsPullback)
{
    const UnitCell cell=CellOf();
    const auto gpw=Block(cell);
    const GPW_Evaluator& ev=*gpw;
    const size_t n=static_cast<const Complex_OIBS&>(*gpw).GetNumFunctions();
    const chmat_t D=RandomHermitian(n);
    const ΔG_Map rho=Contract(ev.Overlap3CTensor(), D);
    const PairProjection<dcmplx> rep=ev.Repulsion3CTensor();
    ASSERT_TRUE(RealisesAdjointPullback(rep)) << "the GPW Coulomb tensor must carry its pullback";
    const rvec3vec_t G=ContractAdjointPullback(rep, [&rho](const ivec3_t& dm)->dcmplx
        { auto it=rho.find(dm); return it==rho.end() ? dcmplx(0.0) : it->second; }, D, &cell);
    ASSERT_EQ(G.size(), 2u);

    // The first-order identity the force rests on: E_H = ½ Re Tr(D Repulsion) on this orthonormal fit.
    const chmat_t V=ContractAdjoint(rep, [&rho](const ivec3_t& dm)->dcmplx
        { auto it=rho.find(dm); return it==rho.end() ? dcmplx(0.0) : it->second; });
    double tr=0.0;
    for (size_t i=0;i<n;i++)
    {
        tr+=std::real(dcmplx(D(i,i)))*std::real(dcmplx(V(i,i)));
        for (size_t j=i+1;j<n;j++) tr+=2.0*std::real(std::conj(dcmplx(D(i,j)))*dcmplx(V(i,j)));
    }
    const double e0=EH(*gpw, D);
    std::printf("  [GPW Hartree] E_H %.12e   ½Tr(DV) %.12e\n", e0, 0.5*tr);
    // To the SCREENING tier, not machine precision: the collocation screens on the signed Re[D e^{-ikR}], the gather on |D| (deliberate,
    // LatticeCollocation::IntegratePotential's note), so at a general k the two are adjoints only up to that asymmetry (4.5e-10 here).
    // The same tier bounds how far the gather's pullback can sit from the forward energy's finite difference below.
    EXPECT_NEAR(e0, 0.5*tr, 1e-8*std::fabs(e0)) << "the first-order energy and the physics round disagree";

    auto central=[&](size_t a, int c, double s)
    {
        rvec3_t e(0,0,0); (c==0 ? e.x : c==1 ? e.y : e.z)=s;
        return (EH(*Block(CellOf(a,e)), D)-EH(*Block(CellOf(a,-1.0*e)), D))/(2*s);
    };
    const double h=1e-3;
    double big=0.0, worst=0.0;
    for (size_t a=0;a<2;a++)
        for (int c=0;c<3;c++)
        {
            const double fd=(4.0*central(a,c,0.5*h)-central(a,c,h))/3.0;   // Richardson (h, h/2)
            const double an= c==0 ? G[a].x : c==1 ? G[a].y : G[a].z;
            std::printf("  [GPW Hartree pullback] atom %zu %c  analytic % .10e  FD % .10e  diff %.2e\n", a, "xyz"[c], an, fd, an-fd);
            big=std::max(big,std::fabs(fd)); worst=std::max(worst,std::fabs(an-fd));
        }
    ASSERT_GT(big, 1e-4) << "a vanishing force cannot show the chain is wired";
    EXPECT_LT(worst, 1e-6*big) << "the GPW Hartree force at fixed D is not its energy's finite difference";
}
