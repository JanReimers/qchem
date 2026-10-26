// File: M_BareCoulomb.C  The BareCoulombSource face (DFT+U increment 3, ACBN0): the raw (m1 m2|m3 m4) tensor
// over a chosen function subset, against the textbook shell averages of ONE normalised d Gaussian shell.
//
// For a single radial function R(r) and the five real d harmonics, every two-electron integral is a
// combination of three Slater integrals F^0, F^2, F^4 -- and three averages are exactly known:
//     (1/25) Sum_ij (ii|jj)        = F^0                  (the shell-averaged Coulomb U_bare; basis-invariant)
//     (1/20) Sum_{i!=j} (ij|ij)    = (5/98)(F^2 + F^4)    (the shell-averaged exchange J_bare in the REAL basis)
//     (ii|ii) = F^0 + 4F^2/49 + 36F^4/441, the same for all five real d orbitals.
// ⚠ The often-quoted J = (F^2+F^4)/14 is a DIFFERENT average (Anisimov's, over complex-harmonic pairs with
// U_mm' - J_mm' combined); the real-basis pair average is (5/7) of it -- verified 2026-09-21 by an
// independent angular quadrature of the multipole expansion (F^2 and F^4 coefficients both 5/98).  ACBN0's own
// J-bar (eq 13) is a real-PAO-basis quantity of exactly this kind, so the real-basis number is the one to gate.
// F^k comes from an INDEPENDENT 2-D radial quadrature here, so this is the whole chain in three numbers:
// M&D FourC on the Cartesian block -> the view's cart->sphere map on all four indices -> normalisation.
#include "gtest/gtest.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
import qchem.BasisSet;                                          // Real_BS / Real_OIBS iteration
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.BareCoulombSource;
import qchem.BasisSet.AoShellSource;                        // the face under test (+ ERI4Block)
import qchem.BasisSet.AtomCentredSet;                       // the pseudo-atomic set's own one-centre integrals (HubbardProjectorPlan inc 4)
import qchem.BasisSet.Gaussian.PG_Cart;                         // the Cartesian family
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;    // MakeSphericalLatticeView
import qchem.Structure;
import qchem.Types;

using namespace qchem;

namespace
{
//! Slater integral F^k of the normalised radial R(r) = N r^2 e^{-a r^2} (a d Gaussian):
//! F^k = int int R^2(r1) R^2(r2) r1^2 r2^2 r_<^k / r_>^{k+1}  = 2 int dr1 R^2(r1) r1^{1-k} int_0^{r1} R^2(r2) r2^{2+k} dr2.
double SlaterF(int k, double a)
{
    const size_t N=40000; const double R=12.0/std::sqrt(a), h=R/N;
    // N^2 from int R^2 r^2 dr = N^2 int r^6 e^{-2 a r^2} dr = N^2 * Gamma(7/2) / (2 (2a)^{7/2}) = 1
    const double N2 = 2.0*std::pow(2.0*a,3.5)/std::tgamma(3.5);
    auto R2=[&](double r){ return N2*std::pow(r,4)*std::exp(-2.0*a*r*r); };
    double inner=0.0, F=0.0, prevIn=0.0;
    for (size_t i=1;i<=N;i++)
    {
        const double r=i*h, rp=(i-1)*h;
        const double gi=R2(r)*std::pow(r,2+k), gp=R2(rp)*std::pow(rp,2+k);
        inner+=0.5*h*(gi+gp);                                  // cumulative int_0^r R^2 r2^{2+k}
        const double fi=R2(r)*std::pow(r,1-k)*inner;
        F+=0.5*h*(fi+prevIn); prevIn=fi;
    }
    return 2.0*F;
}
const BasisSet::Real_OIBS& FirstBlock(const std::shared_ptr<const BasisSet::Real_BS>& bs)
{
    for (auto b : const_cast<BasisSet::Real_BS&>(*bs).Iterate<BasisSet::Real_OIBS>()) return *b;
    throw std::logic_error("no block");
}
//! The columns of every l=\a L shell in the block (the (exponents, LMax) ctor builds EVERY shell 0..LMax).
std::vector<size_t> ShellColumns(const BasisSet::Real_OIBS& blk, int L)
{
    std::vector<size_t> c;
    for (const auto& sh : dynamic_cast<const BasisSet::AoShellSource&>(blk).GetAoShells())
        if (sh.rep->L()==L) for (size_t k=0;k<sh.nComponents();k++) c.push_back(sh.offset+k);
    return c;
}
}

TEST(BareCoulomb, OneDShellReproducesTheSlaterAveragesF0And5F2F4Over98)
{
    const double a=0.9;
    Molecule mol; mol.Insert(new Atom(25, 0, Vector3D<double>(0,0,0)));
    auto* cart=new BasisSet::Gaussian::PG_Cart::BasisSet;
    cart->Insert(new BasisSet::Gaussian::PG_Cart::Orbital_IBS(rvec_t{a}, 2, &mol));
    std::shared_ptr<const BasisSet::Real_BS> cbs(cart);
    auto view=BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(cbs);
    const auto& sph=FirstBlock(view);
    const std::vector<size_t> d=ShellColumns(sph, 2);             // the FIVE harmonics of the one d shell
    ASSERT_EQ(d.size(), 5u);
    const BasisSet::ERI4Block eri=dynamic_cast<const BasisSet::BareCoulombSource&>(sph).BareCoulomb(d);
    ASSERT_EQ(eri.Size(), 5u);
    const double F0=SlaterF(0,a), F2=SlaterF(2,a), F4=SlaterF(4,a);
    ASSERT_GT(F0, 0.0); ASSERT_GT(F2, 0.0); ASSERT_GT(F4, 0.0);
    double U=0.0, J=0.0;
    for (size_t i=0;i<5;i++) for (size_t j=0;j<5;j++)
    {
        U+=eri(i,i,j,j);
        if (i!=j) J+=eri(i,j,i,j);
        EXPECT_GE(eri(i,i,j,j), 0.0);
        EXPECT_GE(eri(i,j,i,j), -1e-12);
    }
    U/=25.0; J/=20.0;
    EXPECT_NEAR(U, F0, 2e-6*F0)            << "the shell-averaged Coulomb is F^0";
    EXPECT_NEAR(J, 5.0*(F2+F4)/98.0, 2e-6*F0)  << "the real-basis shell-averaged exchange is (5/98)(F^2+F^4)";
    for (size_t i=0;i<5;i++) EXPECT_NEAR(eri(i,i,i,i), F0+4*F2/49+36*F4/441, 2e-6*F0) << "orbital " << i;
    // The 8-fold symmetry survives the transform.
    for (size_t a1=0;a1<5;a1++) for (size_t b=0;b<5;b++) for (size_t c=0;c<5;c++) for (size_t d=0;d<5;d++)
    {
        EXPECT_NEAR(eri(a1,b,c,d), eri(b,a1,c,d), 1e-12);
        EXPECT_NEAR(eri(a1,b,c,d), eri(c,d,a1,b), 1e-12);
        EXPECT_NEAR(eri(a1,b,c,d), eri(a1,b,d,c), 1e-12);
    }
    // U_bare - J_bare for a 0.9 bohr^-2 d Gaussian: printed for the record (eV).
    std::cout << "[BareCoulomb] a=" << a << ": F0=" << F0 << " F2=" << F2 << " F4=" << F4
              << "  U_bare=" << U*27.211386245988 << " eV  J_bare=" << J*27.211386245988 << " eV" << std::endl;
}

TEST(BareCoulomb, TheCartesianBlockAnswersDirectlyAndASubsetIsASubTensor)
{
    Molecule mol; mol.Insert(new Atom(25, 0, Vector3D<double>(0,0,0)));
    auto* cart=new BasisSet::Gaussian::PG_Cart::BasisSet;
    cart->Insert(new BasisSet::Gaussian::PG_Cart::Orbital_IBS(rvec_t{0.9, 2.5}, 2, &mol));
    std::shared_ptr<const BasisSet::Real_BS> cbs(cart);
    const auto& blk=dynamic_cast<const BasisSet::BareCoulombSource&>(FirstBlock(cbs));
    const std::vector<size_t> d=ShellColumns(FirstBlock(cbs), 2);   // 2 exponents x 6 Cartesian d
    ASSERT_EQ(d.size(), 12u);
    const BasisSet::ERI4Block all=blk.BareCoulomb(d);
    const std::vector<size_t> pick{d[3],d[7],d[10]};
    const BasisSet::ERI4Block sub=blk.BareCoulomb(pick);
    ASSERT_EQ(all.Size(), 12u); ASSERT_EQ(sub.Size(), 3u);
    const size_t cols[3]={3,7,10};
    for (size_t a=0;a<3;a++) for (size_t b=0;b<3;b++) for (size_t c=0;c<3;c++) for (size_t dd=0;dd<3;dd++)
        EXPECT_DOUBLE_EQ(sub(a,b,c,dd), all(cols[a],cols[b],cols[c],cols[dd]));
    for (size_t i=0;i<12;i++) EXPECT_GT(all(i,i,i,i), 0.0);
    // Two exponents: the six Cartesian components of ONE exponent share (ii|ii) only within {xx,yy,zz} and
    // {xy,xz,yz} (a Cartesian xx is not a pure harmonic) -- the view is where the five harmonics live.
}

// --- the ATOM-CENTRED SET's own bare integrals (HubbardProjectorPlan increment 4) ------------------------------
// The set answers one-centre integrals in closed form (multipole expansion, real Gaunt factors, Slater integrals of
// Gaussian densities as a finite incomplete-beta sum).  The ORACLE is the M&D engine above: the same functions built
// as basis functions and carried to the set's functions on all four indices.
namespace
{
//! The engine's spherical-view block over \a exps with every shell 0..\a LMax at the origin.
std::shared_ptr<const BasisSet::Real_BS> EngineView(const rvec_t& exps, int LMax)
{
    static Molecule mol; static bool once=false;
    if (!once) { mol.Insert(new Atom(25, 0, Vector3D<double>(0,0,0))); once=true; }
    auto* cart=new BasisSet::Gaussian::PG_Cart::BasisSet;
    cart->Insert(new BasisSet::Gaussian::PG_Cart::Orbital_IBS(exps, size_t(LMax), &mol));
    std::shared_ptr<const BasisSet::Real_BS> cbs(cart);
    return BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(cbs);
}
//! The view's columns of the l=\a L shell with exponent \a alpha, in the SET's m order (p: y,z,x; the view keeps a
//! pass-through p shell in its Cartesian x,y,z order).
std::vector<size_t> ViewShell(const BasisSet::Real_OIBS& blk, int L, double alpha)
{
    for (const auto& sh : dynamic_cast<const BasisSet::AoShellSource&>(blk).GetAoShells())
        if (sh.rep->L()==L && sh.exponents.size()==1 && std::fabs(sh.exponents[0]-alpha)<1e-12)
        {
            std::vector<size_t> c;
            for (size_t k=0;k<sh.nComponents();k++) c.push_back(sh.offset+k);
            if (L==1) c={c[1],c[2],c[0]};
            return c;
        }
    throw std::logic_error("ViewShell: no such shell");
}
void ExpectSameTensor(const BasisSet::ERI4Block& a, const BasisSet::ERI4Block& b, double tol)
{
    ASSERT_EQ(a.Size(), b.Size());
    const size_t m=a.Size();
    double big=0.0; for (size_t n=0;n<b.Data().size();n++) big=std::max(big,std::fabs(b.Data()[n]));
    for (size_t i=0;i<m;i++) for (size_t j=0;j<m;j++) for (size_t k=0;k<m;k++) for (size_t l=0;l<m;l++)
        EXPECT_NEAR(a(i,j,k,l), b(i,j,k,l), tol*big) << i << j << k << l;
}
}

TEST(BareCoulomb, TheAtomCentredSetMatchesTheEngineOnOnePrimitiveDShell)
{
    const double a=0.9;
    BasisSet::GaussianAtomCentredSet set; set.Add(rvec3_t{0,0,0}, 2, {Math::Gaussian{1.0,0,a}});
    auto view=EngineView(rvec_t{a}, 2);
    const auto& blk=FirstBlock(view);
    const BasisSet::ERI4Block eng=dynamic_cast<const BasisSet::BareCoulombSource&>(blk).BareCoulomb(ViewShell(blk,2,a));
    ExpectSameTensor(set.BareCoulomb({0,1,2,3,4}), eng, 1e-12);
}

TEST(BareCoulomb, TheAtomCentredSetMatchesTheEngineOnAContractedDShell)
{
    // R = 0.6 r^2 e^{-0.9 r^2} + 0.4 r^2 e^{-2.5 r^2} (normalised by the set); the engine's two unit-normalised d shells
    // carry it through T(s m, m) = c_s / (N_s Nset), N_s the shell's own radial norm.
    const rvec_t al{0.9, 2.5}, c{0.6, 0.4};
    BasisSet::GaussianAtomCentredSet set; set.Add(rvec3_t{0,0,0}, 2, {Math::Gaussian{c[0],0,al[0]}, Math::Gaussian{c[1],0,al[1]}});
    const double nset=std::sqrt(BasisSet::GaussianAtomCentredSet::RadialNorm2(2, {Math::Gaussian{c[0],0,al[0]}, Math::Gaussian{c[1],0,al[1]}}));
    auto view=EngineView(al, 2);
    const auto& blk=FirstBlock(view);
    std::vector<size_t> cols; rmat_t T(10, 5, 0.0);
    for (size_t s=0;s<2;s++)
    {
        const std::vector<size_t> sh=ViewShell(blk, 2, al[s]);
        const double ns=1.0/std::sqrt(BasisSet::GaussianAtomCentredSet::RadialNorm2(2, {Math::Gaussian{1.0,0,al[s]}}));
        for (size_t mm=0;mm<5;mm++) { cols.push_back(sh[mm]); T(5*s+mm, mm)=c[s]/(ns*nset); }
    }
    const BasisSet::ERI4Block eng=dynamic_cast<const BasisSet::BareCoulombSource&>(blk).BareCoulomb(cols).Transform(T);
    ExpectSameTensor(set.BareCoulomb({0,1,2,3,4}), eng, 1e-12);
}

TEST(BareCoulomb, TheAtomCentredSetMatchesTheEngineAcrossPAndDOnOneCentre)
{
    // Mixed l on one centre: the odd-k multipoles and the p-d Gaunt factors.  Set manifolds: p (0..2), d (3..7).
    const double a=1.3;
    BasisSet::GaussianAtomCentredSet set;
    set.Add(rvec3_t{0,0,0}, 1, {Math::Gaussian{1.0,0,a}});
    set.Add(rvec3_t{0,0,0}, 2, {Math::Gaussian{1.0,0,a}});
    auto view=EngineView(rvec_t{a}, 2);
    const auto& blk=FirstBlock(view);
    std::vector<size_t> cols=ViewShell(blk,1,a); for (size_t x : ViewShell(blk,2,a)) cols.push_back(x);
    const BasisSet::ERI4Block eng=dynamic_cast<const BasisSet::BareCoulombSource&>(blk).BareCoulomb(cols);
    ExpectSameTensor(set.BareCoulomb({0,1,2,3,4,5,6,7}), eng, 1e-12);
}

TEST(BareCoulomb, TheAtomCentredSetRefusesTwoCentres)
{
    BasisSet::GaussianAtomCentredSet set;
    set.Add(rvec3_t{0,0,0}, 2, {Math::Gaussian{1.0,0,0.9}});
    set.Add(rvec3_t{0,0,2}, 2, {Math::Gaussian{1.0,0,0.9}});
    EXPECT_THROW(set.BareCoulomb({0,5}), std::invalid_argument);
    EXPECT_NO_THROW(set.BareCoulomb({5,6,7}));
}
