// File: LASolver/tests/LASolver.C  Unit tests for LASolver (Eigen/SVD/Cholesky).
//
// Mock matrices — no BasisSet dependency:
//   S3: 3x3 SPD overlap  S = AᵀA, A=upper-tri [[2,1,0],[0,2,1],[0,0,2]]
//       → S = [[4,2,0],[2,5,2],[0,2,5]]
//   H3: 3x3 symmetric "Hamiltonian", diagonally dominant
//       → H = [[6,2,0],[2,8,2],[0,2,10]]
#include <gtest/gtest.h>
#include <memory>
#include <complex>
#include <cmath>
#include <stdexcept>
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
    ::testing::Values(qchem::Cholesky, qchem::Eigen, qchem::SVD, qchem::Auto, qchem::CholeskyPivoted),
    [](const ::testing::TestParamInfo<qchem::Ortho>& info) -> std::string {
        switch (info.param) {
            case qchem::Cholesky: return "Cholesky";
            case qchem::Eigen:   return "Eigen";
            case qchem::SVD:     return "SVD";
            case qchem::Auto:    return "Auto";
            case qchem::CholeskyPivoted: return "CholeskyPivoted";
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
// PIVOTED Cholesky SELECTS m independent raw AO functions and DROPS the redundant one (rank-revealing).
TEST(LASolverAuto, CholeskyPivotedDropsNearNull)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::CholeskyPivoted, 1e-6));
    s->SetBasisOverlap(make_nearnull3());   // the ~1e-8 near-null pair -> drop one function
    EXPECT_EQ(s->GetOrthoDim(), 2u);        // 3 -> 2: one redundant AO dropped outright
    rsmat_t Sp = s->Transform(make_nearnull3());   // V^H S V on the kept subspace == I_2
    EXPECT_EQ(Sp.rows(), 2u);
    for (size_t i=0;i<2;i++) for (size_t j=0;j<2;j++) EXPECT_NEAR(Sp(i,j), i==j?1.0:0.0, kTol);
}
TEST(LASolverAuto, CholeskyPivotedKeepsWellConditioned)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::CholeskyPivoted, 1e-6));
    s->SetBasisOverlap(make_S3());          // well-conditioned -> keep all 3
    EXPECT_EQ(s->GetOrthoDim(), 3u);
}
// The BASIS-NEUTRAL detector: gap-count -> which functions are redundant (original indices), no ortho.
TEST(LASolverAuto, PivotedCholeskyDropsDetector)
{
    auto drops = qchem::PivotedCholeskyDrops<double>(make_nearnull3());
    EXPECT_EQ(drops.size(), 1u);              // one redundant dimension (the ~1e-8 near-null pair)
    EXPECT_LT(drops[0], 3u);                  // a valid original index
    auto none = qchem::PivotedCholeskyDrops<double>(make_S3());
    EXPECT_TRUE(none.empty());                // well-conditioned -> nothing to drop
}
// The EXPLICIT-tol detector (pin 22's vet-stage trim iterates on it) drops exactly what the ortho path drops
// at the same orthoTol: the near-null pair loses one member at a floor above its residual, nothing at a floor
// below it, and a well-conditioned S keeps everything.
TEST(LASolverAuto, PivotedCholeskyDropsAtTolMatchesTheOrthoPath)
{
    auto one = qchem::PivotedCholeskyDrops<double>(make_nearnull3(), 1e-4);
    EXPECT_EQ(one.size(), 1u);
    auto none = qchem::PivotedCholeskyDrops<double>(make_nearnull3(), 1e-14);
    EXPECT_TRUE(none.empty());
    EXPECT_TRUE(qchem::PivotedCholeskyDrops<double>(make_S3(), 1e-4).empty());
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::CholeskyPivoted, 1e-4));
    s->SetBasisOverlap(make_nearnull3());
    EXPECT_EQ(s->GetOrthoDim(), 3u-one.size());      // the ortho path keeps n - (the detector's drops)
}
// AUTO pivot tol (Factory tol<0): the gap-middle rule finds the ~1e-8/1 gap and drops the near-null mode,
// and leaves a well-conditioned basis untouched -- no user tol.
TEST(LASolverAuto, CholeskyPivotedAutoGapMiddle)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::CholeskyPivoted, -1.0));
    s->SetBasisOverlap(make_nearnull3());
    EXPECT_EQ(s->GetOrthoDim(), 2u);        // gap-middle tol drops the ~1e-8 mode
    std::unique_ptr<LASolver<double>> s2(LASolver<double>::Factory(qchem::CholeskyPivoted, -1.0));
    s2->SetBasisOverlap(make_S3());
    EXPECT_EQ(s2->GetOrthoDim(), 3u);       // no gap -> keep all
}

// ---------------------------------------------------------------------------
//  Ortho::Orthonormal (S = I, the plane-wave basis) and its PARTIAL eigensolve (numBands, QE's nbnd)
// ---------------------------------------------------------------------------

// A dense complex Hermitian matrix with a spread spectrum (no degeneracy at the cut).
static chmat_t make_Hc(size_t n)
{
    chmat_t H(n);
    for (size_t i=0;i<n;i++)
    {
        H(i,i)=dcmplx(double(i)*0.37 + 0.01*double(i*i%7), 0.0);
        for (size_t j=i+1;j<n;j++) H(i,j)=dcmplx(0.05/(1.0+double(j-i)), 0.02*double((i+2*j)%5)-0.04);
    }
    return H;
}
static double cnorm(const cvec_t& v) { double s=0; for (size_t i=0;i<v.size();i++) s+=std::norm(v[i]); return std::sqrt(s); }
static chmat_t make_Ic(size_t n) { chmat_t I(n); for (size_t i=0;i<n;i++) I(i,i)=dcmplx(1.0,0.0); return I; }

TEST(LASolverOrthonormal, RefusesANonIdentityOverlap)
{
    std::unique_ptr<LASolver<double>> s(LASolver<double>::Factory(qchem::Orthonormal));
    EXPECT_THROW(s->SetBasisOverlap(make_S3()), std::invalid_argument);
}
// The full solve on S = I is the Cholesky one (V = I), eigenpair for eigenpair.
TEST(LASolverOrthonormal, FullSolveEqualsCholeskyOnIdentityOverlap)
{
    const size_t n=12;
    std::unique_ptr<LASolver<dcmplx>> o(LASolver<dcmplx>::Factory(qchem::Orthonormal)), c(LASolver<dcmplx>::Factory(qchem::Cholesky));
    o->SetBasisOverlap(make_Ic(n)); c->SetBasisOverlap(make_Ic(n));
    auto [Uo,eo]=o->Solve(make_Hc(n));
    auto [Uc,ec]=c->Solve(make_Hc(n));
    ASSERT_EQ(eo.size(), n);
    EXPECT_LT(blazem::norm(eo-ec), 1e-12);
}
// The lowest numBands eigenpairs: the SAME values as the full solve, eigenvectors satisfying H u = e u and normalised.
TEST(LASolverOrthonormal, PartialSolveIsTheLowestEigenpairsOfTheFullOne)
{
    const size_t n=40, nb=7;
    const chmat_t H=make_Hc(n);
    std::unique_ptr<LASolver<dcmplx>> full(LASolver<dcmplx>::Factory(qchem::Orthonormal)), part(LASolver<dcmplx>::Factory(qchem::Orthonormal, 0.0, nb));
    full->SetBasisOverlap(make_Ic(n)); part->SetBasisOverlap(make_Ic(n));
    auto [Uf,ef]=full->Solve(H);
    for (auto [U,Up,e] : {part->SolveOrtho(H)})
    {
        ASSERT_EQ(e.size(), nb);
        ASSERT_EQ(U.rows(), n);
        ASSERT_EQ(U.columns(), nb);
        for (size_t j=0;j<nb;j++)
        {
            EXPECT_NEAR(e[j], ef[j], 1e-12) << "band " << j;
            const cvec_t u=blazem::column(U,j);
            EXPECT_NEAR(cnorm(u), 1.0, 1e-12) << "band " << j;
            EXPECT_LT(cnorm(cvec_t(mat_t<dcmplx>(H)*u - e[j]*u)), 1e-10) << "band " << j;
            EXPECT_LT(cnorm(cvec_t(cvec_t(blazem::column(Up,j)) - u)), 1e-14) << "U' = U on S = I, band " << j;
        }
    }
}
TEST(LASolverOrthonormal, BandCountIsRefusedByEveryOtherOrtho)
{
    EXPECT_THROW(std::unique_ptr<LASolver<double>>(LASolver<double>::Factory(qchem::Cholesky, 0.0, 5)), std::invalid_argument);
}
// The real (TRIM-block) branch of the same partial solve: syevx.
TEST(LASolverOrthonormal, PartialSolveRealBlock)
{
    const size_t n=30, nb=5;
    rsmat_t H(n), I(n);
    for (size_t i=0;i<n;i++) { I(i,i)=1.0; H(i,i)=0.4*double(i)-0.01*double(i*i%5); for (size_t j=i+1;j<n;j++) H(i,j)=0.07/(1.0+double(j-i)); }
    std::unique_ptr<LASolver<double>> full(LASolver<double>::Factory(qchem::Orthonormal)), part(LASolver<double>::Factory(qchem::Orthonormal, 0.0, nb));
    full->SetBasisOverlap(I); part->SetBasisOverlap(I);
    auto [Uf,ef]=full->Solve(H);
    auto [U,e]=part->Solve(H);
    ASSERT_EQ(e.size(), nb);
    for (size_t j=0;j<nb;j++)
    {
        EXPECT_NEAR(e[j], ef[j], 1e-12) << "band " << j;
        const rvec_t u=blazem::column(U,j);
        EXPECT_NEAR(blazem::norm(u), 1.0, 1e-12) << "band " << j;
        EXPECT_LT(blazem::norm(rvec_t(rmat_t(H)*u - e[j]*u)), 1e-10) << "band " << j;
    }
}
