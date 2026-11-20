// File: Calculation/tests/PhiTablesPolicy.C  The δ fit's XC orbital tables as a setup policy (lane S 2026-11-19; user via the controller,
// the FINAL SHAPE): auto takes the LOWER PEAK unless the other is faster by more than kPhiSpeedGainThreshold AND fits within
// kPhiPeakHeadroom; a stated value wins, a peculiar stated value warns naming auto's pick; nothing fits => REFUSED ON RESOURCES.
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
import qchem.PhiTablesPolicy;
import qchem.AuditFootprint;   // MemoryBudget, AvailableMemory
using namespace qchem;
using qchem::BasisSet::PhiFootprint;
using qchem::BasisSet::PhiRealization;

namespace
{
constexpr double GB=1073741824.0;
//! A footprint with \a dense and \a images GB (no transients), \a nComplex complex blocks.
PhiFootprint Foot(double dense, double images, size_t nComplex=8, size_t nReal=0)
{
    PhiFootprint f;
    f.points=1000; f.complexBlocks=nComplex; f.realBlocks=nReal; f.imageSources=1; f.maxFunctions=10;
    f.imagesAvailable=true; f.denseBytes=dense*GB; f.imagesBytes=images*GB; f.sampleStride=1;
    return f;
}
}

TEST(PhiTablesPolicy, AutoTakesTheLowerPeak)
{
    // Many k-blocks: images lower.  Few: dense lower; the calibrated images gain must beat the threshold AND fit to win.
    EXPECT_EQ(policy::Decide(policy::XCOrbitalTables::Auto, Foot(2.7, 0.4), 9*GB).realization, PhiRealization::Images);
    const auto real=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.2, 0.4, 0, 1), 9*GB);   // real blocks: no speed gain claimed
    EXPECT_EQ(real.realization, PhiRealization::Dense) << real.why;
    EXPECT_FALSE(real.warn);
}

TEST(PhiTablesPolicy, AutoTakesTheHigherPeakOnlyWhenFasterBeyondTheThresholdAndItFits)
{
    // Images (0.4 GB) above dense (0.3 GB) on 8 complex blocks: inside the calibration's regime, so its speed gain counts.
    const bool calibratedFaster = 1.0/policy::kPhiImagesWallRatio-1.0 > policy::kPhiSpeedGainThreshold;
    const size_t nk=policy::kPhiImagesWallRatioMinComplexBlocks;
    const auto roomy=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.3, 0.4, nk), 9*GB);
    EXPECT_EQ(roomy.realization, calibratedFaster ? PhiRealization::Images : PhiRealization::Dense) << roomy.why;
    // The same, with too little budget for the higher peak within the headroom: the lower peak wins whatever the speed.
    const auto tight=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.3, 0.4, nk), 0.6*GB);   // headroom 0.3 GB < the 0.4 GB images peak
    EXPECT_EQ(tight.realization, PhiRealization::Dense) << tight.why;
}

TEST(PhiTablesPolicy, NoSpeedClaimOutsideTheCalibrationsRegime)
{
    // A single complex Γ block (the MnO2 case, controller 2026-11-19): the record shares nothing across k and no run at nk = 1 was measured,
    // so the wall ratio is not claimed and the LOWER peak wins however much budget there is.
    const auto one=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.65, 2.0, 1), 64*GB);
    EXPECT_EQ(one.realization, PhiRealization::Dense) << one.why;
    EXPECT_FALSE(policy::ImagesWallClaimed(Foot(0.65, 2.0, 1)));
    EXPECT_FALSE(policy::ImagesWallClaimed(Foot(0.65, 2.0, policy::kPhiImagesWallRatioMinComplexBlocks-1)));
    EXPECT_TRUE (policy::ImagesWallClaimed(Foot(0.65, 2.0, policy::kPhiImagesWallRatioMinComplexBlocks)));
    EXPECT_FALSE(policy::ImagesWallClaimed(Foot(0.65, 2.0, policy::kPhiImagesWallRatioMinComplexBlocks, 1))) << "real blocks: unmeasured";
}

TEST(PhiTablesPolicy, TheChoiceReadsTheBudgetNotTheFreeMemory)
{
    // Same input, same budget, different free memory: the same choice (a live-headroom choice moved an anchor, 2026-11-19).  The free
    // memory decides only whether the chosen realization can run at all.
    const size_t nk=policy::kPhiImagesWallRatioMinComplexBlocks;
    const auto a=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.3, 0.4, nk), 9*GB, 8*GB);
    const auto b=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.3, 0.4, nk), 9*GB, 0.5*GB);
    EXPECT_EQ(a.realization, b.realization);
    EXPECT_FALSE(a.refuse);
    EXPECT_FALSE(b.refuse) << "0.5 GB free holds either peak (0.3 / 0.4 GB)";
    const auto c=policy::Decide(policy::XCOrbitalTables::Auto, Foot(0.3, 0.4, nk), 9*GB, 0.1*GB);
    EXPECT_EQ(c.realization, a.realization) << "still the budget's choice";
    EXPECT_TRUE(c.refuse) << "0.1 GB free cannot hold it: refused on resources";
}

TEST(PhiTablesPolicy, AStatedValueWinsAndAPeculiarOneWarns)
{
    const auto d=policy::Decide(policy::XCOrbitalTables::Dense, Foot(2.7, 0.4), 9*GB);
    EXPECT_EQ(d.realization, PhiRealization::Dense);
    EXPECT_FALSE(d.warn) << "2.7 GB of 9 is within the headroom";
    const auto dTight=policy::Decide(policy::XCOrbitalTables::Dense, Foot(2.7, 0.4), 4*GB);
    EXPECT_EQ(dTight.realization, PhiRealization::Dense) << "a stated value is never overridden";
    EXPECT_TRUE(dTight.warn);
    EXPECT_NE(dTight.why.find("auto would pick images"), std::string::npos) << dTight.why;
    const auto i=policy::Decide(policy::XCOrbitalTables::Images, Foot(0.3, 0.4, 1), 9*GB);   // one block: images LARGER than dense
    EXPECT_EQ(i.realization, PhiRealization::Images);
    EXPECT_TRUE(i.warn) << "images larger than the dense tables";
    EXPECT_NE(i.why.find("auto would pick"), std::string::npos) << i.why;
}

TEST(PhiTablesPolicy, WithoutImagesDenseIsTheOnlyRealization)
{
    PhiFootprint f=Foot(1.0, 0.0);
    f.imagesAvailable=false;
    const auto r=policy::Decide(policy::XCOrbitalTables::Images, f, 9*GB);
    EXPECT_EQ(r.realization, PhiRealization::Dense);
    EXPECT_TRUE(r.warn) << "a stated images with no images warns";
}

TEST(PhiTablesPolicy, NothingFitsIsRefusedOnResources)
{
    EXPECT_TRUE(policy::Decide(policy::XCOrbitalTables::Auto, Foot(20.0, 5.0), 4*GB).refuse);
    EXPECT_THROW(policy::XCOrbitalTablesPolicy(policy::XCOrbitalTables::Auto, Foot(20.0, 5.0), MemoryBudget{4*GB, "test"}, AvailableMemory{4*GB, "test"}),
                 std::runtime_error);
    EXPECT_EQ(policy::XCOrbitalTablesPolicy(policy::XCOrbitalTables::Auto, Foot(2.7, 0.4), MemoryBudget{9*GB, "test"}, AvailableMemory{9*GB, "test"}),
              PhiRealization::Images);
}
