// File MixingMode.C  SCFParams::Mixing -- what Auto and the explicit modes RESOLVE to (OpenWork §3 "Kerker/Pulay as
// the DEFAULT mixer").  Pure parameter logic: no SCF is run.  The iterator calls ResolveMixing at the top of Iterate
// and the facade prints the resolved value, so these are the rules both of them see.
#include "gtest/gtest.h"
#include <stdexcept>
import qchem.SCFParams;
using namespace qchem;
using M = SCFParams::Mixing;

TEST(MixingMode, AutoIsKerkerPulayWithItsOwnRecipeOnAPeriodicRun)
{
    SCFParams p;                     // every default: Auto, StartingRelaxRo 1.0 (linear's step), no G-space knob
    const SCFParams r=ResolveMixing(p, /*gSpace*/true);
    EXPECT_EQ(r.Mixer, M::Pulay);
    EXPECT_EQ(r.StartingRelaxRo, theAutoMixing.relax);   // NOT the 1.0 linear default
    EXPECT_EQ(r.KerkerG0,        theAutoMixing.G0);
    EXPECT_EQ(r.PulayDepth,      theAutoMixing.depth);
    EXPECT_EQ(r.PulayStart,      theAutoMixing.start);
    // A linear relax set by the caller is LINEAR's step: Auto does not carry it into the G-space recipe.
    p.StartingRelaxRo=0.3;
    EXPECT_EQ(ResolveMixing(p, true).StartingRelaxRo, theAutoMixing.relax);
}

TEST(MixingMode, AutoIsLinearWithoutGSpace)
{
    SCFParams p; p.StartingRelaxRo=0.3;
    const SCFParams r=ResolveMixing(p, /*gSpace*/false);
    EXPECT_EQ(r.Mixer, M::Linear);
    EXPECT_EQ(r.StartingRelaxRo, 0.3);
    // Today's molecular behaviour: a stray G-space knob is ignored, not an error.
    p.KerkerG0=1.0; p.PulayDepth=8;
    EXPECT_EQ(ResolveMixing(p, false).Mixer, M::Linear);
}

TEST(MixingMode, AutoKeepsAHandSetRecipe)
{
    SCFParams k; k.KerkerG0=0.7; k.StartingRelaxRo=0.3;
    const SCFParams rk=ResolveMixing(k, true);
    EXPECT_EQ(rk.Mixer, M::Kerker);
    EXPECT_EQ(rk.KerkerG0, 0.7);
    EXPECT_EQ(rk.StartingRelaxRo, 0.3);
    EXPECT_EQ(rk.PulayDepth, 0);

    SCFParams h; h.PulayDepth=4; h.PulayStart=2; h.StartingRelaxRo=0.2;   // a history with no preconditioner
    const SCFParams rh=ResolveMixing(h, true);
    EXPECT_EQ(rh.Mixer, M::Pulay);
    EXPECT_EQ(rh.KerkerG0, 0.0);
    EXPECT_EQ(rh.PulayDepth, 4);
    EXPECT_EQ(rh.PulayStart, 2);
    EXPECT_EQ(rh.StartingRelaxRo, 0.2);
}

TEST(MixingMode, AutoWithTheCuspDeficitIsPlainKerker)
{
    SCFParams p; p.XCCuspDeficit=true;   // N4 is wired on the Kerker arm only; the Pulay factory throws on it
    const SCFParams r=ResolveMixing(p, true);
    EXPECT_EQ(r.Mixer, M::Kerker);
    EXPECT_EQ(r.KerkerG0, theAutoMixing.G0);
    EXPECT_EQ(r.PulayDepth, 0);
}

TEST(MixingMode, ExplicitLinearIgnoresTheGSpaceKnobs)
{
    SCFParams p; p.Mixer=M::Linear; p.KerkerG0=1.0; p.PulayDepth=8; p.StartingRelaxRo=0.3;
    const SCFParams r=ResolveMixing(p, true);
    EXPECT_EQ(r.Mixer, M::Linear);
    EXPECT_EQ(r.StartingRelaxRo, 0.3);
}

TEST(MixingMode, ExplicitGSpaceModesNeedTheirKnobAndAPeriodicRun)
{
    SCFParams k; k.Mixer=M::Kerker;
    EXPECT_THROW(ResolveMixing(k, true), std::invalid_argument);      // KerkerG0 = 0
    k.KerkerG0=1.0; k.PulayDepth=8;
    EXPECT_EQ(ResolveMixing(k, true).PulayDepth, 0);                   // Kerker = no history
    EXPECT_THROW(ResolveMixing(k, false), std::invalid_argument);     // a molecule cannot

    SCFParams h; h.Mixer=M::Pulay;
    EXPECT_THROW(ResolveMixing(h, true), std::invalid_argument);      // PulayDepth = 0
    h.PulayDepth=8;
    EXPECT_EQ(ResolveMixing(h, true).Mixer, M::Pulay);
}

TEST(MixingMode, ResolvingIsIdempotent)
{
    for (bool g : {true, false})
    {
        const SCFParams once=ResolveMixing(SCFParams{}, g), twice=ResolveMixing(once, g);
        EXPECT_NE(once.Mixer, M::Auto);
        EXPECT_EQ(twice.Mixer, once.Mixer);
        EXPECT_EQ(twice.StartingRelaxRo, once.StartingRelaxRo);
        EXPECT_EQ(twice.KerkerG0, once.KerkerG0);
        EXPECT_EQ(twice.PulayDepth, once.PulayDepth);
        EXPECT_EQ(twice.PulayStart, once.PulayStart);
    }
}

TEST(MixingMode, AutoBesideAFockHistoryIsLinear)
{
    SCFParams p; p.StartingRelaxRo=0.3;   // the measured fight: Pulay on rho beside DIIS on F
    const SCFParams r=ResolveMixing(p, /*gSpace*/true, /*fockExtrapolates*/true);
    EXPECT_EQ(r.Mixer, M::Linear);
    EXPECT_EQ(r.StartingRelaxRo, 0.3);
    // A hand-set recipe, or an explicit mode, is the caller's: kept beside DIIS too.
    SCFParams k; k.KerkerG0=1.0;
    EXPECT_EQ(ResolveMixing(k, true, true).Mixer, M::Kerker);
    SCFParams h; h.Mixer=M::Pulay; h.PulayDepth=8;
    EXPECT_EQ(ResolveMixing(h, true, true).Mixer, M::Pulay);
}
