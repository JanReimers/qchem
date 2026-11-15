// File: Calculation/tests/ERI4MemoryPolicy.C  The ERI4 RAM cap as a setup policy (doc/IntegralCacheOwnershipPlan.md S6; user ruling 6).
//
// The prediction is a COUNT of the doubles the ERI4 tensors store (canonical block pairs i <= j, packed triangles, Direct + Exchange), so
// it is gated against a hand count; the verdict against the two thresholds; and the door against a stated available memory: silent when
// the tensors fit, a warn row when they take more than half, a resource row and a throw when they cannot fit at all.
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
import qchem.ERI4MemoryPolicy;
import qchem.AuditFootprint;   // AvailableMemory
import qchem.Reporting;
using namespace qchem;

TEST(ERI4MemoryPolicy, TheCountIsTheStoredDoubles)
{
    // One block of 3: 6 outer pairs x a full 3x3 inner = 54 doubles per tensor, x2 (J,K), x8 bytes.
    EXPECT_EQ(policy::CountERI4Bytes({3}), 2.0*54*8);
    // Blocks 2 and 3: (0,0) 3 x 4, (0,1) 3 x 9, (1,1) 6 x 9 -> 12+27+54 = 93 doubles per tensor; only the canonical orientation is stored.
    EXPECT_EQ(policy::CountERI4Bytes({2,3}), 2.0*93*8);
    EXPECT_EQ(policy::CountERI4Bytes({}), 0.0);
}

TEST(ERI4MemoryPolicy, ThePredictionIsARangeAroundTheCount)
{
    const auto b=policy::PredictERI4Bytes({20,10});
    EXPECT_EQ(b.low, policy::CountERI4Bytes({20,10})) << "the low end is the rigorous count";
    EXPECT_GT(b.estimate, b.low);
    EXPECT_DOUBLE_EQ(b.high, b.estimate*policy::kERI4OverheadFactor);
}

TEST(ERI4MemoryPolicy, TheVerdictReadsTheRangesEnds)
{
    const policy::ERI4Bytes b{.low=40.0, .estimate=45.0, .high=50.0};
    EXPECT_EQ(policy::JudgeERI4(b, 200.0), policy::ERI4Verdict::Fits);     // high 50 <= half of 200
    EXPECT_EQ(policy::JudgeERI4(b,  90.0), policy::ERI4Verdict::Warn);     // high 50 >  half of 90, low 40 fits
    EXPECT_EQ(policy::JudgeERI4(b,  30.0), policy::ERI4Verdict::Refuse);   // even the low end 40 > 30
}

namespace
{
//! The ERI4Memory rows of the current run report.
std::vector<report::json> ERI4Rows()
{
    std::vector<report::json> rows;
    std::function<void(const report::json&)> walk=[&](const report::json& j)
    {
        if (j.is_object())
            for (auto it=j.begin(); it!=j.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.value("name",std::string())=="ERI4Memory") rows.push_back(r);
                walk(it.value());
            }
        else if (j.is_array()) for (const auto& e : j) walk(e);
    };
    walk(report::CurrentRunReport());
    return rows;
}
}

TEST(ERI4MemoryPolicy, TheDoorIsSilentWarnsOrRefuses)
{
    const std::vector<size_t> blocks{10,10};
    const auto b=policy::PredictERI4Bytes(blocks);
    report::Run run("ERI4MemoryPolicy");
    policy::ERI4Memory(blocks, "test", AvailableMemory{.bytes=4.0*b.high, .source="stated"});
    EXPECT_TRUE(ERI4Rows().empty()) << "a comfortable fit files nothing";
    policy::ERI4Memory(blocks, "test", AvailableMemory{.bytes=1.5*b.high, .source="stated"});
    auto rows=ERI4Rows();
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0]["verdict"], "warn");
    EXPECT_EQ(rows[0]["detail"]["lowBytes"].get<double>(), b.low);
    EXPECT_EQ(rows[0]["detail"]["availableBytes"].get<double>(), 1.5*b.high);
    EXPECT_THROW(policy::ERI4Memory(blocks, "test", AvailableMemory{.bytes=0.5*b.low, .source="stated"}), std::runtime_error);
    rows=ERI4Rows();
    ASSERT_EQ(rows.size(), 2u);
    EXPECT_EQ(rows[1]["verdict"], "resource");
    EXPECT_NE(rows[1]["reason"].get<std::string>().find("available"), std::string::npos) << "the refusal prints the memory available";
}
