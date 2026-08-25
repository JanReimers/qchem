// Outcome<T,E>: the vocabulary type for a call that can legitimately fail (qchem.Outcome).
// The gates that matter are the ones Barton & Nackman's Fallible would FAIL: the value of a failed
// Outcome must be unreadable in EVERY build, not just a debug one.
import qchem.Outcome;
#include <gtest/gtest.h>
#include <string>
#include <memory>
#include <stdexcept>

using qchem::Outcome;
namespace {

struct Why { int code; std::string what; };
using OD = Outcome<double, Why>;

TEST(Outcome, OkCarriesTheValue)
{
    auto o = OD::Ok(2.5);
    EXPECT_TRUE (o.IsOk());
    EXPECT_TRUE (bool(o));
    EXPECT_DOUBLE_EQ(o.Value(), 2.5);
    EXPECT_DOUBLE_EQ(*o, 2.5);
}

TEST(Outcome, FailCarriesTheReason)
{
    auto o = OD::Fail({7, "did not converge"});
    EXPECT_FALSE(o.IsOk());
    EXPECT_FALSE(bool(o));
    EXPECT_EQ(o.Error().code, 7);
    EXPECT_EQ(o.Error().what, "did not converge");
}

// THE POINT OF THE WHOLE TYPE.  Fallible ASSERTS here, so under NDEBUG it hands back garbage -- which is
// the mechanism that hid the site-blocks defect in Release (doc/OpenWork.md item 3).  This must THROW in
// every build, and this test runs in Release too.
TEST(Outcome, ReadingTheValueOfAFailureTHROWS_evenInRelease)
{
    auto o = OD::Fail({1, "nope"});
    EXPECT_THROW((void)o.Value(),     std::runtime_error);
    EXPECT_THROW((void)*o,            std::runtime_error);
    EXPECT_THROW((void)o.TakeValue(), std::runtime_error);
}

TEST(Outcome, ReadingTheReasonOfASuccessTHROWS)
{
    auto o = OD::Ok(1.0);
    EXPECT_THROW((void)o.Error(), std::runtime_error);
}

// T and E of the SAME type must stay unambiguous -- which is why there is no converting ctor.
TEST(Outcome, SameTypeForValueAndErrorIsUnambiguous)
{
    auto ok = Outcome<double,double>::Ok(1.0);
    auto no = Outcome<double,double>::Fail(2.0);
    EXPECT_TRUE (ok.IsOk());   EXPECT_DOUBLE_EQ(ok.Value(), 1.0);
    EXPECT_FALSE(no.IsOk());   EXPECT_DOUBLE_EQ(no.Error(), 2.0);
}

TEST(Outcome, MoveOnlyPayloadSurvives)
{
    auto o = Outcome<std::unique_ptr<int>, Why>::Ok(std::make_unique<int>(42));
    ASSERT_TRUE(o.IsOk());
    auto p = o.TakeValue();
    ASSERT_TRUE(p);
    EXPECT_EQ(*p, 42);
}

} //namespace
