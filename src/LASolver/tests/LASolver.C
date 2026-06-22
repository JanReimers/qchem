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
        solver.reset(LASolver<double>::Factory(GetParam(), 1e-12));
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
    ::testing::Values(qchem::Cholesky, qchem::Eigen, qchem::SVD),
    [](const ::testing::TestParamInfo<qchem::Ortho>& info) -> std::string {
        switch (info.param) {
            case qchem::Cholesky: return "Cholesky";
            case qchem::Eigen:   return "Eigen";
            case qchem::SVD:     return "SVD";
        }
        return "Unknown";
    });
