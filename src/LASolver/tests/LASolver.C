// File: LASolver/tests/LASolver.C  Unit tests for LASolver (Eigen/SVD/Cholesky).
//
// Mock matrices — no BasisSet dependency:
//   S3: 3x3 SPD overlap  S = AᵀA, A=upper-tri [[2,1,0],[0,2,1],[0,0,2]]
//       → S = [[4,2,0],[2,5,2],[0,2,5]]
//   H3: 3x3 symmetric "Hamiltonian", diagonally dominant
//       → H = [[6,2,0],[2,8,2],[0,2,10]]
#include <gtest/gtest.h>
import qchem.LASolver;
import qchem.Blaze;
using namespace qchem;

static constexpr double kTol = 1e-9;

// ---------------------------------------------------------------------------
//  Matrix builders
// ---------------------------------------------------------------------------

static rsmat_t make_S3()
{
    rsmat_t S(3);
    S(0,0)=4; S(0,1)=2; S(0,2)=0;
              S(1,1)=5; S(1,2)=2;
                        S(2,2)=5;
    return S;
}

static rsmat_t make_H3()
{
    rsmat_t H(3);
    H(0,0)=6; H(0,1)=2; H(0,2)=0;
              H(1,1)=8; H(1,2)=2;
                        H(2,2)=10;
    return H;
}

static rsmat_t make_identity3()
{
    rsmat_t I(3);
    for (size_t i = 0; i < 3; i++) I(i,i) = 1.0;
    return I;
}

// ---------------------------------------------------------------------------
//  Verify the fundamental eigenvalue equation: H*U[:,j] == e[j]*S*U[:,j]
// ---------------------------------------------------------------------------

static void check_eigvec_equation(const rsmat_t& H, const rsmat_t& S,
                                   const rmat_t& U, const rvec_t& e)
{
    size_t n = U.columns();
    rmat_t HU  = rmat_t(H) * U;
    rmat_t SUe(U.rows(), n);
    for (size_t j = 0; j < n; j++)
        blazem::column(SUe, j) = e[j] * (rmat_t(S) * blazem::column(U, j));
    EXPECT_LT(blazem::norm(HU - SUe), kTol)
        << "Eigenvector equation H*U = S*U*diag(e) violated";
}

// ---------------------------------------------------------------------------
//  Parametrized fixture — same tests run for Cholesky / Eigen / SVD
// ---------------------------------------------------------------------------

class LASolverTest : public ::testing::TestWithParam<qchem::Ortho>
{
protected:
    std::unique_ptr<LASolver<double>> solver;

    void SetUp() override
    {
        solver.reset(LASolver<double>::Factory(GetParam()));
    }
};

// Transform(S) after SetBasisOverlap(S) should yield identity.
TEST_P(LASolverTest, TransformSgivesIdentity)
{
    rsmat_t S = make_S3();
    solver->SetBasisOverlap(S);
    rsmat_t Sp = solver->Transform(S);
    for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
            EXPECT_NEAR(Sp(i,j), i==j ? 1.0 : 0.0, kTol)
                << " at (" << i << "," << j << ")";
}

// Identity overlap → eigenvalues of H equal diagonal of H (trivial case).
TEST_P(LASolverTest, IdentityOverlapDiagHam)
{
    rsmat_t S = make_identity3();
    rsmat_t H(3);
    H(0,0)=1.0; H(1,1)=2.0; H(2,2)=3.0;  // diagonal

    solver->SetBasisOverlap(S);
    auto [U, e] = solver->Solve(H);

    ASSERT_EQ(e.size(), 3u);
    EXPECT_NEAR(e[0], 1.0, kTol);
    EXPECT_NEAR(e[1], 2.0, kTol);
    EXPECT_NEAR(e[2], 3.0, kTol);
    check_eigvec_equation(H, S, U, e);
}

// H = S → all generalised eigenvalues == 1.
TEST_P(LASolverTest, HequalsS_AllEigvalsOne)
{
    rsmat_t S = make_S3();
    solver->SetBasisOverlap(S);
    auto [U, e] = solver->Solve(S);   // H == S

    ASSERT_EQ(e.size(), 3u);
    for (size_t j = 0; j < 3; j++)
        EXPECT_NEAR(e[j], 1.0, kTol) << " eigenvalue " << j;
    check_eigvec_equation(S, S, U, e);
}

// General H, S — verify eigenvector equation and eigenvalue ordering.
TEST_P(LASolverTest, GeneralHSEigenvectorEquation)
{
    rsmat_t S = make_S3();
    rsmat_t H = make_H3();
    solver->SetBasisOverlap(S);
    auto [U, e] = solver->Solve(H);

    ASSERT_EQ(e.size(), 3u);
    // eigenvalues should be ascending
    for (size_t j = 1; j < 3; j++)
        EXPECT_LT(e[j-1], e[j]);
    check_eigvec_equation(H, S, U, e);
}

// BackTransform is the inverse of the orthogonalising transform:
// BackTransform(I) == V, and Transform + BackTransform ~ identity on U.
TEST_P(LASolverTest, SolveOrthoConsistentWithSolve)
{
    rsmat_t S = make_S3();
    rsmat_t H = make_H3();
    solver->SetBasisOverlap(S);

    auto [U1, e1]            = solver->Solve(H);
    rsmat_t Hp               = solver->Transform(H);
    auto [U2, Up2, e2]       = solver->SolveOrtho(Hp);

    // Both paths must give the same eigenvalues
    ASSERT_EQ(e1.size(), e2.size());
    for (size_t j = 0; j < e1.size(); j++)
        EXPECT_NEAR(e1[j], e2[j], kTol) << " eigenvalue " << j;

    // BackTransform(Up2) must equal U2
    rmat_t U2bt = solver->BackTransform(Up2);
    EXPECT_LT(blazem::norm(U2 - U2bt), kTol);

    // And both back-transformed results satisfy the eigenvalue equation
    check_eigvec_equation(H, S, U2, e2);
}

INSTANTIATE_TEST_SUITE_P(AllMethods, LASolverTest,
    ::testing::Values(qchem::Cholesky, qchem::Eigen, qchem::SVD, qchem::Auto),
    [](const ::testing::TestParamInfo<qchem::Ortho>& info) -> std::string {
        switch (info.param) {
            case qchem::Cholesky: return "Cholesky";
            case qchem::Eigen:   return "Eigen";
            case qchem::SVD:     return "SVD";
            case qchem::Auto:    return "Auto";
        }
        return "Unknown";
    });

// ---------------------------------------------------------------------------
//  Auto gap-detection tolerance + Cholesky-on-singular override (doc/GPWPlan §1)
// ---------------------------------------------------------------------------

// [[1,1-e],[1-e,1]] has eigenvalues {e, 2-e}; with the 1 on the 3rd diagonal -> eigenvalues ~{1e-8, 1, 2}:
// a PSD NEAR-NULL mode with a clean ~1e8x spectral gap above it.
static rsmat_t make_nearnull3()
{
    rsmat_t S(3);
    S(0,0)=1.0; S(0,1)=1.0-1e-8; S(1,1)=1.0; S(2,2)=1.0;
    return S;
}
// [[1,1+e],[1+e,1]] has eigenvalues {-e, 2+e} -> eigenvalues ~{-1e-4, 1, 2}: genuinely INDEFINITE, so a
// Cholesky factorisation cannot proceed.
static rsmat_t make_indefinite3()
{
    rsmat_t S(3);
    S(0,0)=1.0; S(0,1)=1.0+1e-4; S(1,1)=1.0; S(2,2)=1.0;
    return S;
}

// Explicit Eigen with AUTO tolerance (Factory tol < 0): gap detection drops the near-null mode, no user tol.
TEST(LASolverAuto, EigenAutoDropsCleanNearNull)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Eigen, -1.0));
    s->SetBasisOverlap(make_nearnull3());
    EXPECT_EQ(s->GetOrthoDim(), 2u);   // 3 -> 2: the ~1e-8 mode dropped at the clean gap
}
TEST(LASolverAuto, EigenAutoKeepsWellConditioned)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Eigen, -1.0));
    s->SetBasisOverlap(make_S3());
    EXPECT_EQ(s->GetOrthoDim(), 3u);
}
// AUTO (the default): a near-null overlap -> it picks canonical Eigen and drops the mode.
TEST(LASolverAuto, AutoPicksEigenOnNearNull)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Auto));
    s->SetBasisOverlap(make_nearnull3());   // clean gap -> Eigen(auto)
    EXPECT_EQ(s->GetOrthoDim(), 2u);
    rsmat_t Sp = s->Transform(make_nearnull3());
    EXPECT_EQ(Sp.rows(), 2u);
    for (size_t i=0;i<2;i++) for (size_t j=0;j<2;j++) EXPECT_NEAR(Sp(i,j), i==j?1.0:0.0, kTol);
}
// AUTO on an INDEFINITE overlap -> Eigen(auto) force-drops the negative mode (Cholesky would fail).
TEST(LASolverAuto, AutoHandlesIndefinite)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Auto));
    s->SetBasisOverlap(make_indefinite3());
    EXPECT_EQ(s->GetOrthoDim(), 2u);
}
// AUTO on a WELL-conditioned overlap -> it picks Cholesky, full rank kept.
TEST(LASolverAuto, AutoPicksCholeskyWellConditioned)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Auto));
    s->SetBasisOverlap(make_S3());
    EXPECT_EQ(s->GetOrthoDim(), 3u);
}
// EXPLICIT Cholesky is PLAIN -- it never truncates, even a near-null PSD overlap (the escape hatch for
// relativistic bases where mode-dropping is invalid).  It factors it (ill-conditioned but full rank).
TEST(LASolverAuto, CholeskyNeverTruncates)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Cholesky));
    s->SetBasisOverlap(make_nearnull3());   // PSD, potrf succeeds
    EXPECT_EQ(s->GetOrthoDim(), 3u);        // full rank -- NO truncation
}
