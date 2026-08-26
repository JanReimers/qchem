// File: Common/tests/RunPolicy.C  The CP2K-deviation policy (doc/OpenWork.md N5).
//
// What is worth gating here is the RESOLUTION RULE, not the individual defaults: the defaults are
// recorded facts that will change as routes come and go, while the rule -- explicit knob beats the
// umbrella, and a run that is not actually at parity says so -- is the thing that makes the switch
// trustworthy.  A silently-overruled knob would turn CP2K_COMPAT into a trap, which is exactly the
// failure a one-switch design invites.
#include <gtest/gtest.h>
#include <cstdlib>

import qchem.RunPolicy;

using namespace qchem;

namespace
{
//! Set/clear an env var for the duration of a test, restoring whatever was there.  The policy is a
//! process-wide resolved-once object, so every test here re-resolves after arranging the environment.
class Env
{
public:
    Env(const char* n, const char* v) : itsName(n), itsHad(std::getenv(n)!=nullptr)
    {
        if (itsHad) itsOld=std::getenv(n);
        if (v) setenv(n,v,1); else unsetenv(n);
    }
    ~Env()
    {
        if (itsHad) setenv(itsName.c_str(), itsOld.c_str(), 1); else unsetenv(itsName.c_str());
    }
private:
    std::string itsName, itsOld;
    bool        itsHad;
};
//! Every test arranges the environment and then re-resolves; this restores the process default after.
struct Restore { ~Restore() {ReresolveRunPolicy();} };
}

// The DEFAULT run is NOT at parity, and it is able to say which routes make it so.  This is the whole
// point of the item: the low-rank rho route has been on by default since 07d13bf6 and no timing row
// since has said so, because nothing could.
TEST(RunPolicy, DefaultRunDeviatesAndNamesTheRoutes)
{
    Restore r;
    Env a("CP2K_COMPAT",nullptr), b("QCHEM_DM_LOWRANK",nullptr), c("GPW_STREAM_FOLD",nullptr),
        d("QCHEM_MIX_RHO_M",nullptr), e("GPW_XC_DM_SOURCE",nullptr), f("QCHEM_IMPOSE_SYMMETRY",nullptr);
    ReresolveRunPolicy();
    const RunPolicy& p=theRunPolicy();
    EXPECT_FALSE(p.CP2KCompat());
    EXPECT_FALSE(p.AtParity()) << "the default build runs accelerations CP2K does not -- if this ever "
                                  "passes, either a route was retired or one stopped being declared";
    EXPECT_TRUE (p.DMLowRank());
    EXPECT_TRUE (p.StreamFold());
    EXPECT_FALSE(p.MixRhoM());
    EXPECT_FALSE(p.XCFromDM());
    EXPECT_TRUE (p.SymmetryImposition()) << "by default the caller's imposeSymmetry is obeyed";
    EXPECT_EQ(p.Deviations().size(), 5u) << "a new accelerator is not finished until it is in this list";
    EXPECT_NE(p.Banner().find("DEVIATING"), std::string::npos);
}

// THE UMBRELLA: one switch turns every route off, and the policy then reports parity.
TEST(RunPolicy, CP2KCompatTurnsEveryRouteOff)
{
    Restore r;
    Env a("CP2K_COMPAT","1"), b("QCHEM_DM_LOWRANK",nullptr), c("GPW_STREAM_FOLD",nullptr),
        d("QCHEM_MIX_RHO_M",nullptr), e("GPW_XC_DM_SOURCE",nullptr), f("QCHEM_IMPOSE_SYMMETRY",nullptr);
    ReresolveRunPolicy();
    const RunPolicy& p=theRunPolicy();
    EXPECT_TRUE (p.CP2KCompat());
    EXPECT_TRUE (p.AtParity());
    EXPECT_FALSE(p.DMLowRank());
    EXPECT_FALSE(p.StreamFold());
    // THE ONE THAT OVERRULES THE CALLER: CP2K does no symmetry work at all, and every banked recipe asks
    // for an imposition -- so the veto has to come from the switch, not from editing each recipe.
    EXPECT_FALSE(p.SymmetryImposition());
    for (const Deviation& dev : p.Deviations()) EXPECT_FALSE(dev.Deviates()) << dev.knob;
    EXPECT_NE(p.Banner().find("AT PARITY"), std::string::npos);
}

// THE OVERRIDE RULE, and the reason it matters: "parity everywhere except this one" is a legitimate
// experiment, so an explicitly-named knob outranks the umbrella -- and the run must still declare that
// it is NOT at parity, or the switch quietly lies about what was measured.
TEST(RunPolicy, AnExplicitKnobOutranksTheUmbrellaAndSaysSo)
{
    Restore r;
    Env a("CP2K_COMPAT","1"), b("GPW_STREAM_FOLD","1"), c("QCHEM_DM_LOWRANK",nullptr),
        d("QCHEM_MIX_RHO_M",nullptr), e("GPW_XC_DM_SOURCE",nullptr), f("QCHEM_IMPOSE_SYMMETRY",nullptr);
    ReresolveRunPolicy();
    const RunPolicy& p=theRunPolicy();
    EXPECT_TRUE (p.CP2KCompat());
    EXPECT_TRUE (p.StreamFold())  << "an explicitly-stated knob must not be silently overruled";
    EXPECT_FALSE(p.DMLowRank())   << "...while the umbrella still governs the knobs nobody named";
    EXPECT_FALSE(p.AtParity())    << "asking for compat and then contradicting it is NOT parity";
    EXPECT_NE(p.Banner().find("(stated)"), std::string::npos);
}

// A knob set to 0 is SET.  Distinguishing "off" from "not mentioned" is what makes the override rule
// expressible at all -- without it, CP2K_COMPAT could not tell a deliberate off from a default off.
TEST(RunPolicy, SetToZeroCountsAsStated)
{
    Restore r;
    Env a("CP2K_COMPAT",nullptr), b("GPW_STREAM_FOLD","0");
    ReresolveRunPolicy();
    EXPECT_FALSE(theRunPolicy().StreamFold());
    bool found=false;
    for (const Deviation& d : theRunPolicy().Deviations())
        if (std::string(d.knob)=="GPW_STREAM_FOLD") {found=true; EXPECT_TRUE(d.stated); EXPECT_FALSE(d.Deviates());}
    EXPECT_TRUE(found);
}
