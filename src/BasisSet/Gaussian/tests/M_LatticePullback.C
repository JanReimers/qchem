// File: BasisSet/Gaussian/tests/M_LatticePullback.C  THE LATTICE-SUM 1E PULLBACKS AND VIRIALS (GeometryOptimizationPlan §13.2;
// doc/Records/CellHalf_Design.md §5; lane B 2026-11-10): LatticeSum1E::MakeLatticeOverlapPullback / MakeLatticeKineticPullback, libcint's
// ip integrals per shell pair and image over the energy series' own image list, against the central finite difference of what they
// differentiate, E = Re Tr(Ā† M) with M = MakeOverlap(phase,A) or MakeKinetic(phase,A) at a FIXED adjoint Ā and fixed fractional k:
//   - the GRADIENT: each atom moved along each axis, the basis rebuilt;
//   - the VIRIAL: the homogeneous strain A -> (1+ε)A with the fractional positions held (so every atom AND image rides it), all nine
//     components.  At fixed Ā the virial is NOT symmetric (a rotation turns the centres but not a p shell's axes), so all nine are gated;
//     a rotation costing nothing is a property of the converged total energy, gated where the stress is summed, not here.
//   - ATTRIBUTION (controller condition, 2026-11-10): the energy is MnD's lattice sum and the force is libcint's, two integral codes, so
//     libcint's VALUE integral is gated against MnD's image by image (an indicator phase makes MakeOverlap return one image's matrix);
//     a force/energy mismatch then points at the code that is wrong.
//   - TRANSLATION: Σ_A g_A = 0 to roundoff (the image list moves with the atoms, so nothing changes under a rigid shift).
// The same three fixtures as M_CollocationPullback: uncontracted Cartesian (Si SIPP_SR), CONTRACTED shells (H3 dzvp), and d shells through
// the SPHERICAL lattice view (MnO), each at a NON-TRIM k (complex phases: a sign or conjugation slip cannot hide).
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

import qchem.BasisSet.Gaussian.Point.Factory;               // Factory(BasisSetData, cell, engine, angular)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;        // LatticeSum1E
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;   // the block's own PGData (a test may cheat: the libcint twin is built on it)
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;           // ShellPairImageValue (the per-image attribution gate)
import qchem.BasisSet;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.UnitCell;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;

namespace {

constexpr double kTwoPi=6.283185307179586476925286766559;

struct Atom { int Z; rvec3_t f; };
struct Fixture { const char* what; Matrix3D<double> A; std::vector<Atom> atoms; BasisSetData data; bool spherical; };

//! The cell, optionally strained (A -> F A, F = 1 + s e_uv, fractional positions held) or with atom \a moved displaced by \a dr.
UnitCell CellOf(const Fixture& x, size_t moved=size_t(-1), const rvec3_t& dr=rvec3_t(0,0,0), int u=-1, int v=-1, double s=0.0)
{
    Matrix3D<double> A=x.A;
    if (u>=0)
    {
        Matrix3D<double> B=x.A;
        for (int j=1;j<=3;j++) B(u+1,j)=x.A(u+1,j)+s*x.A(v+1,j);   // (F A)_uj = A_uj + s A_vj
        A=B;
    }
    UnitCell c(A);
    for (size_t a=0;a<x.atoms.size();a++)
    {
        rvec3_t f=x.atoms[a].f;
        if (a==moved || moved==size_t(-2)) f=c.ToFractional(c.ToCartesian(f)+dr);   // -2: every atom (the translation)
        c.AddAtom(x.atoms[a].Z, f);
    }
    return c;
}
struct Block { std::shared_ptr<const BasisSet::Real_BS> bs; const LatticeSum1E* ls=nullptr; size_t n=0; };
Block Make(const Fixture& x, const UnitCell& cell)
{
    Block b;
    std::shared_ptr<const BasisSet::Real_BS> cart(Factory(x.data, &cell, Engine::MnD, Angular::Cartesian));
    // THE SCREEN AT 1e-16, not the 1e-10 default.  The default screens the KINETIC series on the overlap prefactor, so it drops kinetic terms
    // up to 4.4e-8 (Si SIPP_SR, measured 2026-11-10; doc/OpenWork.md row "the 1E lattice-sum screen drops kinetic terms").  A dropped term
    // that crosses the screen between the ±h geometries is a 4e-8 STEP in the screened energy, which a central difference reads as a
    // 2e-5 slope.  These gates test the derivative of the series, so they run it where no such step exists.
    {
        GPWTolerances tight; tight.screenEps=1e-16;
        auto* bind=dynamic_cast<ToleranceBinding*>(const_cast<BasisSet::Real_BS*>(cart.get()));
        if (!bind) throw std::logic_error("M_LatticePullback: the MnD Cartesian set offers no ToleranceBinding door");
        bind->BindTolerances(tight);
    }
    b.bs = x.spherical ? std::shared_ptr<const BasisSet::Real_BS>(PG_Spherical::MakeSphericalLatticeView(cart)) : cart;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(*b.bs).Iterate<BasisSet::Real_OIBS>())
    {
        b.ls=dynamic_cast<const LatticeSum1E*>(ibs);
        b.n =ibs->GetNumFunctions();
        break;
    }
    return b;
}
LatticeSum1E::cellphase_t PhaseOf(const rvec3_t& k)
{ return [k](const ivec3_t& n){ const double a=kTwoPi*(k.x*n.x+k.y*n.y+k.z*n.z); return dcmplx(std::cos(a),std::sin(a)); }; }
const rvec3_t kNonTRIM(0.25,0.0,0.125);

chmat_t RandomHermitian(size_t n, unsigned long long seed)
{
    unsigned long long r=seed;
    auto next=[&]{ r^=r<<13; r^=r>>7; r^=r<<17; return double(r%2000)/1000.0-1.0; };
    chmat_t D(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)= i==j ? dcmplx(0.1*next(),0.0) : 0.05*dcmplx(next(),next());
    return D;
}
//! Re Tr(Ā† M) over the full Hermitian matrices, from the upper triangles.
double ReTr(const chmat_t& D, const chmat_t& h)
{
    double e=0.0;
    for (size_t i=0;i<D.rows();i++)
    {
        e+=std::real(dcmplx(D(i,i)))*std::real(dcmplx(h(i,i)));
        for (size_t j=i+1;j<D.rows();j++) e+=2.0*std::real(std::conj(dcmplx(D(i,j)))*dcmplx(h(i,j)));
    }
    return e;
}
double Comp(const rvec3_t& v, int c) {return c==0 ? v.x : c==1 ? v.y : v.z;}

//! Both operators from ONE basis per geometry (the energies E_S, E_T).
struct Pair { double S, T; };
Pair Energies(const Fixture& x, const UnitCell& c, const chmat_t& Abar)
{
    const Block b=Make(x, c);   // a FRESH basis: no memo carries across geometries
    return {ReTr(Abar, b.ls->MakeOverlap(PhaseOf(kNonTRIM), c)), ReTr(Abar, b.ls->MakeKinetic(PhaseOf(kNonTRIM), c))};
}

void PullbacksAreTheFiniteDifference(const Fixture& x)
{
    const UnitCell cell=CellOf(x);
    const Block b=Make(x, cell);
    ASSERT_TRUE(b.ls) << x.what << ": the block realises no lattice-sum face";
    const chmat_t Abar=RandomHermitian(b.n, 0x5DEECE66DULL+b.n);
    const auto PS=b.ls->MakeLatticeOverlapPullback(Abar, PhaseOf(kNonTRIM), cell, &cell);
    const auto PT=b.ls->MakeLatticeKineticPullback(Abar, PhaseOf(kNonTRIM), cell, &cell);
    ASSERT_EQ(PS.gradient.size(), x.atoms.size());
    ASSERT_EQ(PT.gradient.size(), x.atoms.size());

    // RICHARDSON on the central difference (h, h/2), as M_CollocationPullback: the d shells and the contracted H core are narrow.
    auto rich=[](auto E, double h)
    {
        auto central=[&](double s){ const Pair p=E(s), m=E(-s); return Pair{(p.S-m.S)/(2*s), (p.T-m.T)/(2*s)}; };
        const Pair a=central(0.5*h), c=central(h);
        return Pair{(4.0*a.S-c.S)/3.0, (4.0*a.T-c.T)/3.0};
    };
    double bigS=0.0, bigT=0.0, worstS=0.0, worstT=0.0;
    auto tally=[&](const char* what, double anS, double anT, const Pair& fd)
    {
        std::printf("  [lattice pullback] %-22s %-10s S: analytic % .10e FD % .10e diff %.2e   T: analytic % .10e FD % .10e diff %.2e\n",
                    x.what, what, anS, fd.S, anS-fd.S, anT, fd.T, anT-fd.T);
        bigS=std::max(bigS,std::fabs(fd.S)); worstS=std::max(worstS,std::fabs(anS-fd.S));
        bigT=std::max(bigT,std::fabs(fd.T)); worstT=std::max(worstT,std::fabs(anT-fd.T));
    };
    const double h=1e-3;
    for (size_t a=0;a<x.atoms.size();a++)
        for (int c=0;c<3;c++)
        {
            const Pair fd=rich([&](double s){ rvec3_t e(0,0,0); (c==0 ? e.x : c==1 ? e.y : e.z)=s; return Energies(x, CellOf(x,a,e), Abar); }, h);
            char what[16]; std::snprintf(what, sizeof what, "atom %zu %c", a, "xyz"[c]);
            tally(what, Comp(PS.gradient[a],c), Comp(PT.gradient[a],c), fd);
        }
    for (int u=0;u<3;u++)
        for (int v=0;v<3;v++)
        {
            const Pair fd=rich([&](double s){ return Energies(x, CellOf(x,size_t(-1),rvec3_t(0,0,0),u,v,s), Abar); }, h);
            char what[16]; std::snprintf(what, sizeof what, "strain %c%c", "xyz"[u], "xyz"[v]);
            tally(what, PS.virial(u+1,v+1), PT.virial(u+1,v+1), fd);
        }
    ASSERT_GT(bigS, 1e-4) << x.what << ": a vanishing overlap pullback cannot show the walk is wired";
    ASSERT_GT(bigT, 1e-4) << x.what << ": a vanishing kinetic pullback cannot show the walk is wired";
    EXPECT_LT(worstS, 1e-7*bigS) << x.what << ": the overlap pullback / virial is not the finite difference of Re Tr(Ā† S)";
    EXPECT_LT(worstT, 1e-7*bigT) << x.what << ": the kinetic pullback / virial is not the finite difference of Re Tr(Ā† T)";

    rvec3_t sS(0,0,0), sT(0,0,0);
    for (size_t a=0;a<x.atoms.size();a++) { sS=sS+PS.gradient[a]; sT=sT+PT.gradient[a]; }
    const double tS=std::max({std::fabs(sS.x),std::fabs(sS.y),std::fabs(sS.z)}), tT=std::max({std::fabs(sT.x),std::fabs(sT.y),std::fabs(sT.z)});
    std::printf("  [lattice pullback] %-22s |Σ_A g_A|  S %.2e  T %.2e\n", x.what, tS, tT);
    EXPECT_LT(tS, 1e-12*std::max(1.0,bigS)) << x.what << ": a rigid translation moves nothing, so the overlap gradients sum to zero";
    EXPECT_LT(tT, 1e-12*std::max(1.0,bigT)) << x.what << ": a rigid translation moves nothing, so the kinetic gradients sum to zero";
}

Fixture Si()   // FCC primitive, a = 10.26, the atoms OFF the origin (no corner atom: the displaced cell must not re-wrap one)
{ return {"Si SIPP_SR cartesian", Matrix3D<double>(0.0,5.13,5.13, 5.13,0.0,5.13, 5.13,5.13,0.0), {{14,{0.1,0.05,0.12}},{14,{0.35,0.3,0.37}}},
          BasisSetData::SIPP_SR, false}; }
Fixture H3()   // CONTRACTED dzvp shells (H 1s: 4 primitives) in an orthorhombic box
{ return {"H3 dzvp (contracted)", Matrix3D<double>(9.0,0.0,0.0, 0.0,9.5,0.0, 0.0,0.0,10.0),
          {{1,{0.50,0.50,0.50}},{1,{0.60,0.52,0.50}},{1,{0.47,0.59,0.53}}}, BasisSetData::DZVP, false}; }
Fixture MnO()   // rocksalt, d shells through the spherical view's T
{ return {"MnO VA spherical view", Matrix3D<double>(0.0,4.20,4.20, 4.20,0.0,4.20, 4.20,4.20,0.0), {{25,{0.1,0.1,0.1}},{8,{0.6,0.58,0.62}}},
          BasisSetData::VALENCE_LOWQ_VA, true}; }

} // namespace

TEST(M_LatticePullback, OverlapAndKineticPullbacksAndVirialsAreTheirFiniteDifferences_Uncontracted) {PullbacksAreTheFiniteDifference(Si());}
TEST(M_LatticePullback, OverlapAndKineticPullbacksAndVirialsAreTheirFiniteDifferences_Contracted)   {PullbacksAreTheFiniteDifference(H3());}
TEST(M_LatticePullback, OverlapAndKineticPullbacksAndVirialsAreTheirFiniteDifferences_Spherical)    {PullbacksAreTheFiniteDifference(MnO());}

namespace {
void LibcintImagesAreTheEnergyEnginesImages(const Fixture& x)
{
    using LC=Evaluators::PG_LibCint::NR_Evaluator;
    const UnitCell cell=CellOf(x);
    const Block b=Make(x, cell);
    ASSERT_TRUE(b.ls);
    const auto* pg=dynamic_cast<const Evaluators::PG_Cart_MnD::PGData*>(b.ls);
    ASSERT_TRUE(pg) << x.what << ": the Cartesian block is no PGData";
    const LC lc(*pg, &cell);
    double big=0.0, worst=0.0;
    for (const ivec3_t& n0 : {ivec3_t(0,0,0), ivec3_t(1,0,0), ivec3_t(0,1,0), ivec3_t(-1,0,1), ivec3_t(1,1,-1)})
    {
        const LatticeSum1E::cellphase_t only=[n0](const ivec3_t& n){ return (n.x==n0.x && n.y==n0.y && n.z==n0.z) ? dcmplx(1.0) : dcmplx(0.0); };
        const rvec3_t R=cell.ToCartesian(rvec3_t(n0.x,n0.y,n0.z));
        for (int which=0;which<2;which++)
        {
            const chmat_t M = which==0 ? b.ls->MakeOverlap(only, cell) : b.ls->MakeKinetic(only, cell);   // the upper triangle: K_ij(n0), i <= j
            for (size_t si=0;si<lc.NumShells();si++)
                for (size_t sj=0;sj<lc.NumShells();sj++)
                {
                    const auto& P=lc.ShellComponents(si); const auto& Q=lc.ShellComponents(sj);
                    const rmat_t V=lc.ShellPairImageValue(which, si, sj, R);
                    for (size_t p=0;p<P.size();p++) for (size_t q=0;q<Q.size();q++)
                    {
                        if (P[p]>Q[q]) continue;   // the lower triangle holds conj(M_ji), the image -n0's
                        const double mnd=std::real(dcmplx(M(P[p],Q[q])));
                        big=std::max(big,std::fabs(mnd));
                        worst=std::max(worst,std::fabs(mnd-V(p,q)));
                    }
                }
        }
    }
    std::printf("  [lattice images] %-22s max|K| %.3e  max|libcint - MnD| %.2e\n", x.what, big, worst);
    ASSERT_GT(big, 1e-2) << x.what << ": no image integral to compare";
    EXPECT_LT(worst, 1e-10*std::max(1.0,big)) << x.what << ": libcint's lattice image integral is not the energy engine's (MnD)";
}
} // namespace
TEST(M_LatticePullback, LibcintsImageIntegralsAreTheEnergyEnginesImageByImage_Uncontracted) {LibcintImagesAreTheEnergyEnginesImages(Si());}
TEST(M_LatticePullback, LibcintsImageIntegralsAreTheEnergyEnginesImageByImage_Contracted)   {LibcintImagesAreTheEnergyEnginesImages(H3());}

TEST(M_LatticePullback, RefusesAMisSizedAdjointAndAMissingStructure)
{
    const Fixture x=Si();
    const UnitCell cell=CellOf(x);
    const Block b=Make(x, cell);
    ASSERT_TRUE(b.ls);
    EXPECT_THROW((void)b.ls->MakeLatticeOverlapPullback(chmat_t(b.n+1), PhaseOf(kNonTRIM), cell, &cell), std::invalid_argument);
    EXPECT_THROW((void)b.ls->MakeLatticeKineticPullback(chmat_t(b.n), PhaseOf(kNonTRIM), cell, nullptr), std::invalid_argument);
}
