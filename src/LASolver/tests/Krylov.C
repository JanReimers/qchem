// File: LASolver/tests/Krylov.C  Unit tests for the matrix-free GMRES solver (doc/LinearResponsePlan.md M1).
//
// The operators are dense matrices wrapped in the LinearOperator face, so the answer is checked against a
// direct solve.  Both scalars: the response of a real molecule is real, a q != 0 lattice response is complex
// and NON-Hermitian, which is the case GMRES (not CG/MINRES) is here for.
#include <gtest/gtest.h>
#include <complex>
#include <random>
import qchem.LASolver.Krylov;
import qchem.Blaze;
using namespace qchem;

namespace {

template <class T> class DenseOperator : public LinearOperator<T>
{
public:
    explicit DenseOperator(mat_t<T> A) : itsA(std::move(A)) {}
    virtual size_t   Dimension() const override {return itsA.rows();}
    virtual vec_t<T> Apply(const vec_t<T>& x, double) const override {itsCalls++; return itsA*x;}
    mutable size_t itsCalls=0;
private:
    mat_t<T> itsA;
};

template <class T> T Random(std::mt19937& g)
{
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    if constexpr (std::is_floating_point_v<T>) return u(g);
    else return T(u(g), u(g));
}

//! A NON-symmetric, well-conditioned test matrix: a shifted identity plus a random perturbation.
template <class T> mat_t<T> MakeMatrix(size_t n, unsigned seed)
{
    std::mt19937 g(seed);
    mat_t<T> A(n, n);
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++) A(i,j)=Random<T>(g)*(0.5/std::sqrt(double(n)));
    for (size_t i=0;i<n;i++) A(i,i)+=T(2.0);
    return A;
}
template <class T> vec_t<T> MakeVector(size_t n, unsigned seed)
{
    std::mt19937 g(seed);
    vec_t<T> v(n);
    for (size_t i=0;i<n;i++) v[i]=Random<T>(g);
    return v;
}
template <class T> double Residual(const mat_t<T>& A, const vec_t<T>& x, const vec_t<T>& b)
{
    vec_t<T> r=b-A*x;
    double num=0, den=0;
    for (size_t i=0;i<r.size();i++) {num+=std::norm(r[i]); den+=std::norm(b[i]);}
    return std::sqrt(num/den);
}

template <class T> void SolvesTo(size_t n, size_t restart)
{
    const mat_t<T> A=MakeMatrix<T>(n, 7);
    const vec_t<T> b=MakeVector<T>(n, 11);
    DenseOperator<T> op(A);
    auto o=SolveGMRES<T>(op, b, nullptr, {.tol=1e-12, .maxIter=500, .restart=restart});
    ASSERT_TRUE(o.IsOk()) << o.Error().detail;
    EXPECT_LE(o->residual, 1e-12);
    EXPECT_LE(Residual(A, o->x, b), 1e-11);                  // the TRUE residual, recomputed here
}

} // namespace

TEST(Krylov, RealNonSymmetric)          {SolvesTo<double>(60, 40);}
TEST(Krylov, ComplexNonHermitian)       {SolvesTo<dcmplx>(60, 40);}
//! A restart SMALLER than the problem must still converge (more cycles, same answer).
TEST(Krylov, RealRestarted)             {SolvesTo<double>(60, 5);}
TEST(Krylov, ComplexRestarted)          {SolvesTo<dcmplx>(60, 5);}

//! The iteration cap returns a FAILED Outcome carrying the residual it reached -- never an unconverged x.
TEST(Krylov, CapIsAFailureNotANumber)
{
    const mat_t<double> A=MakeMatrix<double>(60, 7);
    const vec_t<double> b=MakeVector<double>(60, 11);
    DenseOperator<double> op(A);
    auto o=SolveGMRES<double>(op, b, nullptr, {.tol=1e-14, .maxIter=3, .restart=40});
    ASSERT_FALSE(o.IsOk());
    EXPECT_GT(o.Error().residual, 1e-14);
    EXPECT_EQ(o.Error().iterations, 3u);
    EXPECT_THROW((void)o.Value(), std::runtime_error);
}

//! b = 0 has the solution x = 0, with no operator application at all.
TEST(Krylov, ZeroRightHandSide)
{
    DenseOperator<double> op(MakeMatrix<double>(10, 7));
    auto o=SolveGMRES<double>(op, vec_t<double>(10, 0.0), nullptr, {});
    ASSERT_TRUE(o.IsOk());
    EXPECT_EQ(o->iterations, 0u);
    EXPECT_EQ(op.itsCalls, 0u);
    for (size_t i=0;i<10;i++) EXPECT_EQ(o->x[i], 0.0);
}

//! A warm start that is already the solution costs one check (the true residual) and no Krylov step.
TEST(Krylov, ExactWarmStart)
{
    const mat_t<dcmplx> A=MakeMatrix<dcmplx>(30, 3);
    const vec_t<dcmplx> x=MakeVector<dcmplx>(30, 5);
    const vec_t<dcmplx> b=A*x;
    DenseOperator<dcmplx> op(A);
    auto o=SolveGMRES<dcmplx>(op, b, &x, {.tol=1e-12});
    ASSERT_TRUE(o.IsOk());
    EXPECT_EQ(o->iterations, 0u);
    EXPECT_EQ(op.itsCalls, 1u);
}

//! The identity converges in ONE step (the Krylov space is invariant after one vector: lucky breakdown).
TEST(Krylov, IdentityIsOneStep)
{
    mat_t<double> I(8, 8, 0.0);
    for (size_t i=0;i<8;i++) I(i,i)=1.0;
    const vec_t<double> b=MakeVector<double>(8, 2);
    DenseOperator<double> op(I);
    auto o=SolveGMRES<double>(op, b, nullptr, {.tol=1e-14});
    ASSERT_TRUE(o.IsOk());
    EXPECT_EQ(o->iterations, 1u);
    for (size_t i=0;i<8;i++) EXPECT_NEAR(o->x[i], b[i], 1e-15);
}
