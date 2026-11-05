// File: LASolver/tests/OverlapConditioning.C  The one overlap analysis behind every facade's basis report.
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <limits>
import qchem.LASolver.OverlapConditioning;
import qchem.Blaze;

using namespace qchem;

namespace
{
//! S = [[1,c],[c,1]] has eigenvalues 1-c and 1+c.
template <class T> hmat_t<T> Pair(double c)
{
    hmat_t<T> S(2);
    S(0,0)=1.0; S(1,1)=1.0; S(0,1)=T(c);
    return S;
}
}

TEST(OverlapConditioning, ReadsTheSpectrumOfAKnownPair)
{
    const auto c=AnalyzeOverlap<double>(Pair<double>(0.5));
    EXPECT_NEAR(c.lambdaMin, 0.5, 1e-14);
    EXPECT_NEAR(c.lambdaMax, 1.5, 1e-14);
    EXPECT_NEAR(c.cond,      3.0, 1e-13);
    EXPECT_EQ(c.wouldDrop, 0u);
}

TEST(OverlapConditioning, ComplexBlockAgreesWithItsRealTwin)
{
    const auto r=AnalyzeOverlap<double>(Pair<double>(0.9));
    hmat_t<dcmplx> S(2);
    S(0,0)=1.0; S(1,1)=1.0; S(0,1)=dcmplx(0.0,0.9);   // a genuinely complex S, same spectrum
    const auto c=AnalyzeOverlap<dcmplx>(S);
    EXPECT_NEAR(c.lambdaMin, r.lambdaMin, 1e-14);
    EXPECT_NEAR(c.lambdaMax, r.lambdaMax, 1e-14);
}

// The pivoted Cholesky's second pivot of the pair is 1-c^2, so the floor sits either side of it.
TEST(OverlapConditioning, CountsDropsAgainstTheFloorItIsGiven)
{
    const double c=1.0-1e-8;   // 1-c^2 ~ 2e-8
    EXPECT_EQ(AnalyzeOverlap<double>(Pair<double>(c), 1e-6).wouldDrop, 1u);
    EXPECT_EQ(AnalyzeOverlap<double>(Pair<double>(c), 1e-9).wouldDrop, 0u);
    EXPECT_EQ(AnalyzeOverlap<double>(Pair<double>(c)).wouldDrop, 1u) << "the default is ConditioningFloor";
}

TEST(OverlapConditioning, AnIndefiniteOverlapShowsItsNegativeEigenvalue)
{
    const auto c=AnalyzeOverlap<double>(Pair<double>(1.0+1e-10));
    EXPECT_LT(c.lambdaMin, 0.0);
    EXPECT_GT(c.cond, 1e9) << "cond reads the smallest MAGNITUDE";
}

TEST(OverlapConditioning, AnExactlySingularOverlapIsNotAnInfiniteOrZeroCond)
{
    hmat_t<double> S(2);
    S(0,0)=1.0; S(1,1)=0.0; S(0,1)=0.0;
    EXPECT_EQ(AnalyzeOverlap<double>(S).cond, std::numeric_limits<double>::max());
}

TEST(OverlapConditioning, AnEmptyOverlapAnswersZeros)
{
    const auto c=AnalyzeOverlap<double>(hmat_t<double>(0));
    EXPECT_EQ(c.lambdaMin, 0.0);
    EXPECT_EQ(c.wouldDrop, 0u);
}
