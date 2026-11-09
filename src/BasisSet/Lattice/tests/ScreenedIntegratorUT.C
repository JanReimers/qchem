// File: BasisSet/Lattice/tests/ScreenedIntegratorUT.C  ScreenedMatrixIntegrator -- the PairProjection raw
// pair behind the mesh-level MatrixForward/MatrixAdjoint faces.
//
// Synthetic closures, deliberately: what is under test is the ADAPTER and its contract, not GPW's
// collocation.  A real screened pair needs a periodic Gaussian basis and belongs in the GPW gates; the
// plumbing is what a unit test can pin exactly (user rule: dev loop is unit, integration is acceptance).
#include "gtest/gtest.h"
#include <stdexcept>
#include <complex>
#include <functional>

import qchem.BasisSet.PairProjection;
import qchem.Types;

using namespace qchem;
using qchem::PairProjection;
using qchem::ScreenedMatrixIntegrator;

namespace {
// A 2x2 toy whose forward/adjoint pair IS exactly adjoint: rho_g = sum_ij D_ij P(g,i) P(g,j) and
// H_ij = sum_g w_g v_g P(g,i) P(g,j), with one fixed P.  Two "raster points".
const double P[2][2] = {{1.0, 2.0},{3.0,-1.0}};   // P[g][i]
const rvec_t W{0.25, 0.75};

PairProjection<double> MakeAdjointPair()
{
    PairProjection<double> g;
    g.applyRaw = [](const hmat_t<double>& D)
    {
        rvec_t rho(2,0.0);
        for (size_t p=0;p<2;p++)
            for (size_t i=0;i<2;i++) for (size_t j=0;j<2;j++) rho[p]+=D(i,j)*P[p][i]*P[p][j];
        return rho;
    };
    g.applyRawAdjoint = [](const rvec_t& v)
    {
        hmat_t<double> H(2);
        for (size_t i=0;i<2;i++) for (size_t j=i;j<2;j++)
        {
            double s=0.0;
            for (size_t p=0;p<2;p++) s+=W[p]*v[p]*P[p][i]*P[p][j];
            H(i,j)=s;
        }
        return H;
    };
    return g;
}
} // anon

//---------------------------------------------------------------------------------------------------
// ★ THE CONTRACT: an integrator is the PAIR.  A PairProjection that can collocate but not integrate back --
// or the reverse -- is not one, and saying so at CONSTRUCTION beats discovering it at the first Adjoint
// call inside an SCF iteration.
TEST(ScreenedMatrixIntegrator, AHalfPairIsRejectedAtConstruction)
{
    PairProjection<double> fwdOnly;  fwdOnly.applyRaw       =[](const hmat_t<double>&){return rvec_t(2,0.0);};
    PairProjection<double> adjOnly;  adjOnly.applyRawAdjoint=[](const rvec_t&){return hmat_t<double>(2);};
    PairProjection<double> neither;

    EXPECT_THROW((ScreenedMatrixIntegrator<double>(fwdOnly, W)), std::runtime_error);
    EXPECT_THROW((ScreenedMatrixIntegrator<double>(adjOnly, W)), std::runtime_error);
    EXPECT_THROW((ScreenedMatrixIntegrator<double>(neither, W)), std::runtime_error);
    const PairProjection<double> ok=MakeAdjointPair();
    EXPECT_NO_THROW((ScreenedMatrixIntegrator<double>(ok, W)));
}

//---------------------------------------------------------------------------------------------------
// The adapter forwards to the tensor's own closures -- no arithmetic of its own in either direction.
TEST(ScreenedMatrixIntegrator, ItForwardsBothDirectionsToTheTensor)
{
    const PairProjection<double> g=MakeAdjointPair();
    const ScreenedMatrixIntegrator<double> I(g, W);
    EXPECT_EQ(I.NumCoefficients(), 2u);

    hmat_t<double> D(2); D(0,0)=1.0; D(0,1)=0.5; D(1,1)=2.0;
    const rvec_t rho=I.Forward(D).Vector();
    // point 0: P=(1,2) -> 1*1 + 2*(0.5*1*2) + 4*2 = 1 + 2 + 8 = 11
    EXPECT_DOUBLE_EQ(rho[0], 11.0);
    // point 1: P=(3,-1) -> 9*1 + 2*(0.5*3*-1) + 1*2 = 9 - 3 + 2 = 8
    EXPECT_DOUBLE_EQ(rho[1], 8.0);

    const rvec_t ones(2,1.0);
    EXPECT_DOUBLE_EQ(I.Integrate(ones), 1.0);      // 0.25 + 0.75
}

//---------------------------------------------------------------------------------------------------
// ★ AND THE PROPERTY THE FACE EXISTS FOR, checked on the SCREENED realization exactly as it is checked on
// the dense one: <v, Forward(D)>_w == <Adjoint(v), D>.  A screened route is adjoint only if it truncates
// IDENTICALLY in both directions, so this is the gate a future screening change has to survive.
TEST(ScreenedMatrixIntegrator, ForwardAndAdjointAreExactlyAdjoint)
{
    const PairProjection<double> g=MakeAdjointPair();
    const ScreenedMatrixIntegrator<double> I(g, W);

    hmat_t<double> D(2); D(0,0)=1.3; D(0,1)=-0.7; D(1,1)=2.1;
    const rvec_t v{0.9, -1.4};

    rvec_t vrho=I.Forward(D).Vector();
    for (size_t p=0;p<2;p++) vrho[p]*=v[p];
    const double lhs=I.Integrate(vrho);

    const hmat_t<double> H=I.Adjoint(v);
    double rhs=0.0;
    for (size_t i=0;i<2;i++) for (size_t j=0;j<2;j++) rhs+=H(i,j)*D(i,j);

    EXPECT_NEAR(lhs, rhs, 1e-13) << "the screened pair is not adjoint -- H = dE/dD does not hold";
}

//===================================================================================================
// THE ONE-OF-N RECORD (V-ONE-OF-N 2/4, user rulings 2026-11-07/08): a Projections holds exactly one
// realization; the paired implementation reads it with Vector()/GMap(), a wrong read THROWS, and everyone
// else combines records with the dispatchers.
TEST(Projections, AWrongTypedReadThrows)
{
    const Projections a(rvec_t{1.0, 2.0});
    ΔG_Map m; m[ivec3_t(0,0,1)]=dcmplx(0.5,-0.25);
    const Projections g(m);
    EXPECT_NO_THROW(a.Vector());
    EXPECT_NO_THROW(g.GMap());
    EXPECT_THROW(a.GMap(),   std::logic_error);
    EXPECT_THROW(g.Vector(), std::logic_error);
    const Projections empty;
    EXPECT_THROW(empty.Vector(), std::logic_error);
    EXPECT_THROW(empty.GMap(),   std::logic_error);
}

TEST(Projections, TheEmptyRecordIsTheAdditiveIdentityAndAMismatchThrows)
{
    Projections s;                                   // a sum over blocks starts here, realization unknown
    s+=Projections(rvec_t{1.0, -2.0});
    s+=Projections(rvec_t{0.5,  4.0});
    s*=2.0;
    ASSERT_EQ(s.Vector().size(), 2u);
    EXPECT_DOUBLE_EQ(s.Vector()[0], 3.0);
    EXPECT_DOUBLE_EQ(s.Vector()[1], 4.0);
    EXPECT_DOUBLE_EQ(s.MaxAbs(), 4.0);
    Projections d; d-=Projections(rvec_t{1.0, 0.0});
    EXPECT_DOUBLE_EQ(d.Vector()[0], -1.0);

    ΔG_Map m; m[ivec3_t(1,0,0)]=dcmplx(0.0,3.0);
    Projections g(m);
    g+=Projections(m);
    EXPECT_DOUBLE_EQ(g.MaxAbs(), 6.0);
    EXPECT_THROW(s+=g, std::logic_error);                          // rvec_t += ΔG_Map
    EXPECT_THROW(s+=Projections(rvec_t{1.0}), std::logic_error); // different lengths
}

// The rvalue read MOVES the member out, so binding a reference to a temporary's member cannot dangle.
TEST(Projections, TheRvalueReadIsAValue)
{
    const rvec_t& v=Projections(rvec_t{7.0, 8.0}).Vector();   // lifetime-extended prvalue
    ASSERT_EQ(v.size(), 2u);
    EXPECT_DOUBLE_EQ(v[1], 8.0);
}

//===================================================================================================
// THE G-SPACE REALIZATION: plane-wave columns forward (Contract) and a Fourier-lookup adjoint
// (applyAdjoint), the pair Vee_Hartree now holds.  A 2-function toy with "G" indices m_0, m_1.
namespace {
const ivec3_t m0(0,0,0), m1(1,0,0);
PairProjection<dcmplx> MakeColumnPair(double volume)
{
    PairProjection<dcmplx> g;
    g.volume=volume;
    // columns: dm = m_i - m_j, the (i,j) pairs that land on it (row-major)
    g.columns.push_back({ivec3_t( 0,0,0), {{0,0},{1,1}}});
    g.columns.push_back({ivec3_t(-1,0,0), {{0,1}}});
    g.columns.push_back({ivec3_t( 1,0,0), {{1,0}}});
    // <i|f|j> = f(m_i - m_j)
    g.applyAdjoint=[](const std::function<dcmplx(const ivec3_t&)>& f)
    {
        const ivec3_t m[2]={m0,m1};
        hmat_t<dcmplx> H(2);
        for (size_t i=0;i<2;i++) for (size_t j=i;j<2;j++) H(i,j)=f(m[i]-m[j]);
        return H;
    };
    return g;
}
} // anon

TEST(GSpacePairProjectionIntegrator, AHalfPairOrADenseProjectionIsRejectedAtConstruction)
{
    PairProjection<dcmplx> noAdjoint=MakeColumnPair(2.0); noAdjoint.applyAdjoint=nullptr;
    PairProjection<dcmplx> noForward=MakeColumnPair(2.0); noForward.columns.clear();
    EXPECT_THROW(GSpacePairProjectionIntegrator<dcmplx>{noAdjoint}, std::runtime_error);
    EXPECT_THROW(GSpacePairProjectionIntegrator<dcmplx>{noForward}, std::runtime_error);
    const PairProjection<dcmplx> ok=MakeColumnPair(2.0);
    EXPECT_NO_THROW(GSpacePairProjectionIntegrator<dcmplx>{ok});
}

// Forward is Contract, returning the G-space record; the typed Adjoint is ContractAdjoint on the owner's own map.
TEST(GSpacePairProjectionIntegrator, ForwardAndAdjointCarryTheGSpaceRecord)
{
    const PairProjection<dcmplx> g=MakeColumnPair(2.0);
    const GSpacePairProjectionIntegrator<dcmplx> I(g);
    hmat_t<dcmplx> D(2); D(0,0)=1.0; D(0,1)=dcmplx(0.3,0.4); D(1,1)=3.0;

    const Projections rho=I.Forward(D);
    const ΔG_Map ref=Contract(g, D);
    ASSERT_EQ(rho.GMap().size(), ref.size());
    for (const auto& [k,v] : ref) EXPECT_EQ(rho.GMap().at(k), v);
    EXPECT_EQ(I.NumCoefficients(), 3u);
    EXPECT_THROW(rho.Vector(), std::logic_error);   // the G-space realization: a vector read throws

    const hmat_t<dcmplx> H=I.Adjoint(ref);
    EXPECT_EQ(H(0,0), ref.at(m0));
    EXPECT_EQ(H(0,1), ref.at(ivec3_t(-1,0,0)));
}
