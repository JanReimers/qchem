// File: src/Hamiltonian/tests/RangeAudit.C  The quotient range audit's term body (AuditAlong), on a toy near-singular block.
//
// The audit (user rulings 2026-11-25) holds a term's Rayleigh quotient along the basis's weakest directions to its operator's rigorous
// range.  The toy: two functions chi = C e over an orthonormal pair e, C = [[1,1],[0,eps]] (nearly parallel, lambda_min(S) ~ eps^2/2),
// and an operator diagonal on e with spectrum {-1, 2}.  Its exact matrix A = C^T H C gives every quotient inside [-1, 2] however small
// lambda is; an error delta*I in the matrix moves the near-null quotient by delta/lambda, which the audit sees.  This keeps the check
// tested after the emulated Hermite3 overrun (the defect it was built for) is gone.
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
import qchem.Hamiltonian.Internal.Term;   // AuditAlong + TermRangeAudit (a unit test may import internals)
import qchem.Blaze;

using namespace qchem;

namespace
{
struct Toy { hmat_t<double> S, A; mat_t<double> v; rvec_t lambda; };
Toy MakeToy(double eps, double delta)
{
    rmat_t C(2,2); C(0,0)=1; C(0,1)=1; C(1,0)=0; C(1,1)=eps;
    rmat_t H(2,2); H(0,0)=-1; H(0,1)=0; H(1,0)=0; H(1,1)=2;
    Toy t; t.S.resize(2); t.A.resize(2);
    for (size_t i=0;i<2;i++) for (size_t j=i;j<2;j++)
    {
        double s=0, a=0;
        for (size_t k=0;k<2;k++) { s+=C(k,i)*C(k,j); a+=C(k,i)*H(k,k)*C(k,j); }
        t.S(i,j)=s; t.A(i,j)=a+(i==j ? delta : 0.0);
    }
    blazem::eigen(t.S, t.lambda, t.v);   // ascending; eigenvectors as columns
    return t;
}
} //anon

// The exact matrix: every quotient lies in the operator's range, at a lambda_min of 5e-9, and the near-null direction is still inside.
TEST(RangeAudit, AnExactMatrixStaysInItsRangeAlongANearNullDirection)
{
    const Toy t=MakeToy(1e-4, 0.0);
    ASSERT_LT(t.lambda[0], 1e-8);
    const Hamiltonian::TermRangeAudit a=Hamiltonian::AuditAlong<double>(t.A, -1.0, 2.0, t.v, t.lambda);
    EXPECT_TRUE(a.hasRange);
    EXPECT_GE(a.quotient, -1.0-1e-6);
    EXPECT_LE(a.quotient,  2.0+1e-6);
    EXPECT_LT(a.Excess(), 1e-7) << "an exact matrix needs no allowance beyond roundoff (excess " << a.Excess() << ")";
}

// A corrupted matrix (delta*I, a relative error of ~1e-6): the near-null quotient leaves the range by delta/lambda ~ 200 Ha, so the
// audit names that direction and an allowance factor of 1e-7 does not cover it.  The well-conditioned direction stays inside.
TEST(RangeAudit, ACorruptedMatrixLeavesItsRangeAlongTheNearNullDirection)
{
    const Toy t=MakeToy(1e-4, 1e-6);
    const Hamiltonian::TermRangeAudit a=Hamiltonian::AuditAlong<double>(t.A, -1.0, 2.0, t.v, t.lambda);
    EXPECT_EQ(a.direction, 0u) << "the worst direction is the near-null one";
    EXPECT_GT(a.quotient, 2.0+100.0);
    EXPECT_NEAR(a.lambda, t.lambda[0], 1e-20);
    EXPECT_GT(a.Excess(), 1e-7);
    EXPECT_NEAR(a.Excess()*a.scale, a.quotient-2.0, 1e-9*a.quotient);   // the excess is the distance outside, per unit scale
}

// An open range (Kinetic's [0, inf)) is a sign check: a positive quotient is inside whatever its size.
TEST(RangeAudit, AnOpenUpperSideIsASignCheck)
{
    const Toy t=MakeToy(1e-4, 1e-6);   // the corrupted quotient is large and POSITIVE
    const Hamiltonian::TermRangeAudit a=Hamiltonian::AuditAlong<double>(t.A, -1.0, std::numeric_limits<double>::infinity(), t.v, t.lambda);
    EXPECT_LT(a.Excess(), 0.0);
}
