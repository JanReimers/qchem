// File: Hubbard.C  DFT+U unit gates -- the LowdinProjector pair on a synthetic block, no SCF.
//
// The property the term rests on is that the occupation-matrix FORWARD and the potential ADJOINT are the
// same T = S^{1/2}[:,M] contracted both ways, so H = dE/dD holds exactly; and that Loewdin occupations of
// an idempotent density (projected onto a manifold that spans it) are exactly 0 or 1.
#include "gtest/gtest.h"
#include <complex>
#include <cmath>
#include <vector>
import qchem.Hamiltonian.Internal.Hubbard;   // LowdinProjector (unit tests may import Internal)
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::Hamiltonian;

namespace
{
// A 4-function block whose overlap is a mildly non-orthogonal SPD matrix; manifold = functions {1,2}.
hmat_t<double> OverlapR()
{
    hmat_t<double> S(4);
    for (size_t i=0;i<4;i++) S(i,i)=1.0;
    S(0,1)=0.2; S(0,2)=-0.1; S(0,3)=0.05; S(1,2)=0.3; S(1,3)=-0.15; S(2,3)=0.1;
    return S;
}
hmat_t<dcmplx> OverlapC()
{
    hmat_t<dcmplx> S(4);
    for (size_t i=0;i<4;i++) S(i,i)=1.0;
    S(0,1)=dcmplx(0.2,0.05); S(0,2)=dcmplx(-0.1,0.0); S(0,3)=dcmplx(0.05,-0.02);
    S(1,2)=dcmplx(0.3,0.1);  S(1,3)=dcmplx(-0.15,0.0); S(2,3)=dcmplx(0.1,0.03);
    return S;
}
const std::vector<std::vector<size_t>> M{{1,2}};
}

//---------------------------------------------------------------------------------------------------
TEST(LowdinProjector, ForwardAndAdjointAreExactlyAdjoint_Real)
{
    const LowdinProjector<double> P(OverlapR(), M);
    EXPECT_EQ(P.NumCoefficients(), 4u);   // 2x2 flattened
    hmat_t<double> D(4); D(0,0)=1.3; D(0,1)=-0.7; D(0,2)=0.2; D(0,3)=0.1; D(1,1)=2.1; D(1,2)=0.4; D(1,3)=-0.3; D(2,2)=0.6; D(2,3)=0.05; D(3,3)=0.9;
    const rvec_t W{0.9, -0.4, -0.4, 1.7};   // a symmetric 2x2, flattened
    const rvec_t n=P.Forward(D);
    double lhs=0.0; for (size_t a=0;a<4;a++) lhs+=W[a]*n[a];
    const hmat_t<double> V=P.Adjoint(W);
    double rhs=0.0; for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++) rhs+=V(i,j)*D(i,j);
    EXPECT_NEAR(lhs, rhs, 1e-13) << "Tr(W n) != Tr(V D): the pair is not adjoint";
}

TEST(LowdinProjector, ForwardAndAdjointAreExactlyAdjoint_Complex)
{
    const LowdinProjector<dcmplx> P(OverlapC(), M);
    hmat_t<dcmplx> D(4);
    D(0,0)=1.3; D(0,1)=dcmplx(-0.7,0.2); D(0,2)=dcmplx(0.2,-0.1); D(0,3)=0.1; D(1,1)=2.1; D(1,2)=dcmplx(0.4,0.3); D(1,3)=-0.3; D(2,2)=0.6; D(2,3)=dcmplx(0.05,0.02); D(3,3)=0.9;
    const rvec_t W{0.9, -0.4, -0.4, 1.7};
    const rvec_t n=P.Forward(D);
    double lhs=0.0; for (size_t a=0;a<4;a++) lhs+=W[a]*n[a];
    const hmat_t<dcmplx> V=P.Adjoint(W);
    dcmplx rhs=0.0; for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++) rhs+=V(i,j)*std::conj(dcmplx(D(i,j)));
    // With a REAL symmetric W the forward keeps only Re(T^dag D T), so the adjoint identity holds on the
    // real part; Tr(V D) is real for Hermitian V, D.
    EXPECT_NEAR(lhs, std::real(rhs), 1e-13);
    EXPECT_NEAR(std::imag(rhs), 0.0, 1e-13);
}

//---------------------------------------------------------------------------------------------------
// An idempotent density living entirely inside the manifold: D = S^{-1/2} P S^{-1/2} with P a projector
// onto Loewdin-orthogonal functions {1}. Its Loewdin occupation matrix is P restricted -- eigenvalues {1,0}.
TEST(LowdinProjector, AnIdempotentDensityInTheManifoldHasOccupations0And1)
{
    const hmat_t<double> S=OverlapR();
    rvec_t w; rmat_t V; blazem::eigen(S, w, V);
    rmat_t Sinvh(4,4);
    for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++)
    { double s=0; for (size_t k=0;k<4;k++) s+=V(i,k)*V(j,k)/std::sqrt(w[k]); Sinvh(i,j)=s; }
    // D = S^{-1/2} e1 e1^T S^{-1/2}: one electron in Loewdin function 1 (a member of the manifold).
    hmat_t<double> D(4);
    for (size_t i=0;i<4;i++) for (size_t j=i;j<4;j++) D(i,j)=Sinvh(i,1)*Sinvh(j,1);
    const LowdinProjector<double> P(S, M);
    const rvec_t n=P.Forward(D);           // the 2x2 block over {1,2}: should be [[1,0],[0,0]]
    EXPECT_NEAR(n[0], 1.0, 1e-12);
    EXPECT_NEAR(n[1], 0.0, 1e-12);
    EXPECT_NEAR(n[2], 0.0, 1e-12);
    EXPECT_NEAR(n[3], 0.0, 1e-12);
    EXPECT_NEAR(P.Integrate(n), 1.0, 1e-12);   // the manifold charge
}

//---------------------------------------------------------------------------------------------------
TEST(LowdinProjector, TheFactoredForwardAgreesWithTheUnfactoredOne)
{
    const LowdinProjector<double> P(OverlapR(), M);
    rmat_t L(4,2); L(0,0)=1.0; L(0,1)=0.2; L(1,0)=-0.5; L(1,1)=1.5; L(2,0)=0.3; L(2,1)=-0.4; L(3,0)=0.1; L(3,1)=0.7;
    hmat_t<double> D(4);
    for (size_t i=0;i<4;i++) for (size_t j=i;j<4;j++) { double s=0; for (size_t m=0;m<2;m++) s+=L(i,m)*L(j,m); D(i,j)=s; }
    const rvec_t a=P.Forward(D), b=P.Forward(L);
    ASSERT_EQ(a.size(), b.size());
    for (size_t k=0;k<a.size();k++) EXPECT_NEAR(a[k], b[k], 1e-13);
}
