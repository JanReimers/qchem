// File: BasisSet/Lattice_3D/tests/ScreenedIntegratorUT.C  ScreenedMatrixIntegrator -- the Projector3 raw
// pair behind the mesh-level MatrixIntegrator face.
//
// Synthetic closures, deliberately: what is under test is the ADAPTER and its contract, not GPW's
// collocation.  A real screened pair needs a periodic Gaussian basis and belongs in the GPW gates; the
// plumbing is what a unit test can pin exactly (user rule: dev loop is unit, integration is acceptance).
#include "gtest/gtest.h"
#include <stdexcept>

import qchem.BasisSet.Internal.Projector3;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using qchem::Projector3;
using qchem::ScreenedMatrixIntegrator;

namespace {
// A 2x2 toy whose forward/adjoint pair IS exactly adjoint: rho_g = sum_ij D_ij P(g,i) P(g,j) and
// H_ij = sum_g w_g v_g P(g,i) P(g,j), with one fixed P.  Two "raster points".
const double P[2][2] = {{1.0, 2.0},{3.0,-1.0}};   // P[g][i]
const rvec_t W{0.25, 0.75};

Projector3<double> MakeAdjointPair()
{
    Projector3<double> g;
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
// ★ THE CONTRACT: an integrator is the PAIR.  A Projector3 that can collocate but not integrate back --
// or the reverse -- is not one, and saying so at CONSTRUCTION beats discovering it at the first Adjoint
// call inside an SCF iteration.
TEST(ScreenedMatrixIntegrator, AHalfPairIsRejectedAtConstruction)
{
    Projector3<double> fwdOnly;  fwdOnly.applyRaw       =[](const hmat_t<double>&){return rvec_t(2,0.0);};
    Projector3<double> adjOnly;  adjOnly.applyRawAdjoint=[](const rvec_t&){return hmat_t<double>(2);};
    Projector3<double> neither;

    EXPECT_THROW((ScreenedMatrixIntegrator<double>(fwdOnly, W)), std::runtime_error);
    EXPECT_THROW((ScreenedMatrixIntegrator<double>(adjOnly, W)), std::runtime_error);
    EXPECT_THROW((ScreenedMatrixIntegrator<double>(neither, W)), std::runtime_error);
    const Projector3<double> ok=MakeAdjointPair();
    EXPECT_NO_THROW((ScreenedMatrixIntegrator<double>(ok, W)));
}

//---------------------------------------------------------------------------------------------------
// The adapter forwards to the tensor's own closures -- no arithmetic of its own in either direction.
TEST(ScreenedMatrixIntegrator, ItForwardsBothDirectionsToTheTensor)
{
    const Projector3<double> g=MakeAdjointPair();
    const ScreenedMatrixIntegrator<double> I(g, W);
    EXPECT_EQ(I.NumCoefficients(), 2u);

    hmat_t<double> D(2); D(0,0)=1.0; D(0,1)=0.5; D(1,1)=2.0;
    const rvec_t rho=I.Forward(D);
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
    const Projector3<double> g=MakeAdjointPair();
    const ScreenedMatrixIntegrator<double> I(g, W);

    hmat_t<double> D(2); D(0,0)=1.3; D(0,1)=-0.7; D(1,1)=2.1;
    const rvec_t v{0.9, -1.4};

    rvec_t vrho=I.Forward(D);
    for (size_t p=0;p<2;p++) vrho[p]*=v[p];
    const double lhs=I.Integrate(vrho);

    const hmat_t<double> H=I.Adjoint(v);
    double rhs=0.0;
    for (size_t i=0;i<2;i++) for (size_t j=0;j<2;j++) rhs+=H(i,j)*D(i,j);

    EXPECT_NEAR(lhs, rhs, 1e-13) << "the screened pair is not adjoint -- H = dE/dD does not hold";
}
