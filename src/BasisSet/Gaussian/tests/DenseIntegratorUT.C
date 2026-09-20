// File: DenseIntegratorUT.C  DenseProjector3Integrator -- the ANALYTIC forward/adjoint pair (R1.0q, 2026-09-19).
//
// The dense realization is the molecular Gaussian auxiliary-basis route: one full <ab|c> matrix per fit
// function, no grid anywhere.  These gates pin the property the face exists for on THAT realization --
// H = dE/dD holds because the adjoint is the exact transpose of the forward on the same tensor -- plus the
// two mechanics a realization can get wrong (the factored overload being hidden; a tensor with no dense part).
#include "gtest/gtest.h"
#include <stdexcept>
#include <vector>
import qchem.BasisSet.Projector3;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;

namespace
{
// Two fit functions over a 3-function orbital block: hand-picked symmetric <ab|c> matrices.
Projector3<double> MakeDense()
{
    Projector3<double> g;
    smat_t<double> A(3), B(3);
    A(0,0)=1.0; A(0,1)=0.5; A(0,2)=-0.25; A(1,1)=2.0; A(1,2)=0.75; A(2,2)=0.5;
    B(0,0)=0.3; B(0,1)=-0.4; B(0,2)=0.1; B(1,1)=1.1; B(1,2)=0.6; B(2,2)=-0.2;
    g.dense={A,B};
    return g;
}
const rvec_t Q{0.8, 1.7};   // <f_a|1> per fit function
}

//---------------------------------------------------------------------------------------------------
TEST(DenseProjector3Integrator, ATensorWithNoDensePartIsRejectedAtConstruction)
{
    Projector3<double> empty;
    EXPECT_THROW((DenseProjector3Integrator<double>(empty, Q)), std::runtime_error);
    const Projector3<double> g=MakeDense();
    EXPECT_THROW((DenseProjector3Integrator<double>(g, rvec_t(3,1.0))), std::runtime_error);   // wrong count
    EXPECT_NO_THROW((DenseProjector3Integrator<double>(g, Q)));
}

//---------------------------------------------------------------------------------------------------
// The forward is Sum_ab D_ab <ab|c>: the loop FiniteIrrepCD::GetRepulsion3C used to run, elementwise
// product over the FULL symmetric matrix (both triangles), which is what a Blaze D % A sum gives.
TEST(DenseProjector3Integrator, TheForwardIsTheFullDoubleContraction)
{
    const Projector3<double> g=MakeDense();
    const DenseProjector3Integrator<double> I(g, Q);
    EXPECT_EQ(I.NumCoefficients(), 2u);

    hmat_t<double> D(3); D(0,0)=1.0; D(0,1)=0.5; D(0,2)=0.0; D(1,1)=2.0; D(1,2)=-1.0; D(2,2)=0.5;
    const rvec_t rho=I.Forward(D);
    double a=0.0, b=0.0;
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) {a+=D(i,j)*g.dense[0](i,j); b+=D(i,j)*g.dense[1](i,j);}
    EXPECT_DOUBLE_EQ(rho[0], a);
    EXPECT_DOUBLE_EQ(rho[1], b);
    EXPECT_DOUBLE_EQ(I.Integrate(rvec_t{1.0,1.0}), 2.5);   // 0.8 + 1.7
}

//---------------------------------------------------------------------------------------------------
// ★ THE PROPERTY: <v, Forward(D)> == Tr(Adjoint(v) D) exactly -- same tensor, contracted both ways.
TEST(DenseProjector3Integrator, ForwardAndAdjointAreExactlyAdjoint)
{
    const Projector3<double> g=MakeDense();
    const DenseProjector3Integrator<double> I(g, Q);

    hmat_t<double> D(3); D(0,0)=1.3; D(0,1)=-0.7; D(0,2)=0.2; D(1,1)=2.1; D(1,2)=0.4; D(2,2)=-0.6;
    const rvec_t v{0.9, -1.4};

    const rvec_t rho=I.Forward(D);
    double lhs=0.0;
    for (size_t a=0;a<2;a++) lhs+=v[a]*rho[a];

    const hmat_t<double> H=I.Adjoint(v);
    double rhs=0.0;
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) rhs+=H(i,j)*D(i,j);

    EXPECT_NEAR(lhs, rhs, 1e-13) << "the dense pair is not adjoint -- H = dE/dD does not hold";
}

//---------------------------------------------------------------------------------------------------
// The factored overload must reach the base default (form LL^T, delegate), not be hidden by the override.
TEST(DenseProjector3Integrator, TheFactoredForwardAgreesWithTheUnfactoredOne)
{
    const Projector3<double> g=MakeDense();
    const DenseProjector3Integrator<double> I(g, Q);

    mat_t<double> L(3,2); L(0,0)=1.0; L(0,1)=0.2; L(1,0)=-0.5; L(1,1)=1.5; L(2,0)=0.3; L(2,1)=-0.4;
    hmat_t<double> D(3);
    for (size_t i=0;i<3;i++) for (size_t j=i;j<3;j++) {double s=0; for (size_t m=0;m<2;m++) s+=L(i,m)*L(j,m); D(i,j)=s;}

    const rvec_t full=I.Forward(D), fact=I.Forward(L);
    ASSERT_EQ(full.size(), fact.size());
    for (size_t a=0;a<full.size();a++) EXPECT_NEAR(full[a], fact[a], 1e-14);
}
