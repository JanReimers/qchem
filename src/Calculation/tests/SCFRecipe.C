// File: Calculation/tests/SCFRecipe.C  The SCF recipe policy (qchem.SCFRecipe; doc/SCFRobustnessPlan.md R4.1 as ruled for B5).
//
// THE GOLDEN TABLE: `OldPolicy` below is the pre-B5 Auto policy COPIED VERBATIM from SolidCalculation::MakeAccelerator (the accelerator)
// and SolidCalculation::Converge (the restored stage's linear start <= 0.3), at 17ef37e6e.  B5 promised that Auto (and Insulator) move
// no anchor; the first test holds the new chooser to that over the whole product space the old code saw.
#include <gtest/gtest.h>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
import qchem.SCFRecipe;
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;
import qchem.Deck;
import qchem.SolidCalculation;

using namespace qchem;
using T = SCFAccelerators::Type;
using M = SCFParams::Mixing;

namespace {

struct Choice
{
    T         accel;
    bool      gdmRung;     //!< the Ladder forced onto its GDM rung
    bool      threw;       //!< the mixer's admission threw (R1.1)
    SCFParams resolved;    //!< what the iterator would build (valid when !threw)
};

// A Ladder keeps a DIIS rung, and DIIS is one: they extrapolate the Fock (SCFAccelerator::ExtrapolatesFock).
bool Fock(T t) {return t==T::DIIS || t==T::Ladder;}

SCFParams Resolve(const SCFParams& p, bool gSpace, T t) {return ResolveMixing(p, gSpace, Fock(t));}

//! ⚠ VERBATIM, the pre-B5 policy (SolidCalculation.C at 17ef37e6e): MakeAccelerator's Auto branch, then Converge's restored cap, then
//! the iterator's ResolveMixing.  Do not "tidy" it: it is the record the new chooser is held to.
Choice OldPolicy(T asked, bool restored, double kT, int numBands, const SCFParams& params, bool gSpace)
{
    Choice c{asked, false, false, {}};
    T t = asked;
    if (asked==T::Auto)
    {
        if (!restored)                    t=T::Null;
        else if (kT==0.0 && numBands==0) { t=T::Ladder; c.gdmRung=true; }
        else                              t=T::DIIS;   // smeared, or a PARTIAL eigensolve (GDM needs every virtual)
    }
    c.accel=t;
    try
    {
        SCFParams own=params;
        if (restored && asked==T::Auto && Resolve(params, gSpace, t).Mixer==M::Linear)
            own.StartingRelaxRo=std::min(own.StartingRelaxRo, 0.3);
        c.resolved=Resolve(own, gSpace, t);
    }
    catch (const std::invalid_argument&) {c.threw=true;}
    return c;
}

//! The B5 path, as SolidCalculation now runs it: ChooseAccelerator, FillRecipeKnobs + CapRestoredStep (Converge), ResolveMixing.
Choice NewPolicy(SCFRecipe r, T asked, SCFStageKind k, double kT, int numBands, const SCFParams& params, bool gSpace)
{
    const RecipeAccelerator a=ChooseAccelerator(r, k, asked, kT, numBands);
    Choice c{a.type, a.forceGDMRung, false, {}};
    try
    {
        const RecipeKnobs knobs=FillRecipeKnobs(r, k, params);
        const SCFParams own=CapRestoredStep(k, asked, knobs.params, Resolve(knobs.params, gSpace, a.type).Mixer==M::Linear);
        c.resolved=Resolve(own, gSpace, a.type);
    }
    catch (const std::invalid_argument&) {c.threw=true;}
    return c;
}

//! The pre-B5 policy plus the ONE change landed after it on purpose (OpenWork NEAR-TERM, 2026-11-12): a restored Auto stage that
//! resolves linear caps its CEILING at 0.3 too, not only its start.
Choice OldPolicyWithRestoredCeiling(T asked, bool restored, double kT, int numBands, const SCFParams& params, bool gSpace)
{
    Choice c=OldPolicy(asked, restored, kT, numBands, params, gSpace);
    if (restored && asked==T::Auto && !c.threw && c.resolved.Mixer==M::Linear)
        c.resolved.MaxRelaxRo=std::min(c.resolved.MaxRelaxRo, 0.3);
    return c;
}

void ExpectSame(const Choice& a, const Choice& b, const std::string& where)
{
    EXPECT_EQ(a.accel,   b.accel)   << where;
    EXPECT_EQ(a.gdmRung, b.gdmRung) << where;
    ASSERT_EQ(a.threw,   b.threw)   << where;
    if (a.threw) return;
    const SCFParams& p=a.resolved; const SCFParams& q=b.resolved;
    EXPECT_EQ(p.Mixer,           q.Mixer)           << where;
    EXPECT_EQ(p.StartingRelaxRo, q.StartingRelaxRo) << where;
    EXPECT_EQ(p.MaxRelaxRo,      q.MaxRelaxRo)      << where;
    EXPECT_EQ(p.KerkerG0,        q.KerkerG0)        << where;
    EXPECT_EQ(p.PulayDepth,      q.PulayDepth)      << where;
    EXPECT_EQ(p.PulayStart,      q.PulayStart)      << where;
    EXPECT_EQ(p.UseMOM,          q.UseMOM)          << where;
}

//! The stated-knob variants the old code met: nothing stated, a hand-set history, a Kerker, a stated linear step, a low start.
std::vector<std::pair<std::string,SCFParams>> Variants()
{
    std::vector<std::pair<std::string,SCFParams>> v;
    v.push_back({"defaults", SCFParams{}});
    { SCFParams p; p.PulayDepth=8; p.PulayStart=5; p.KerkerG0=1.0; p.StartingRelaxRo=0.45; v.push_back({"K+P knobs", p}); }
    { SCFParams p; p.PulayDepth=6;                                                        v.push_back({"Pulay-D knobs", p}); }
    { SCFParams p; p.KerkerG0=0.8;                                                        v.push_back({"Kerker knob", p}); }
    { SCFParams p; p.Mixer=M::Linear; p.StartingRelaxRo=0.5;                              v.push_back({"stated linear 0.5", p}); }
    { SCFParams p; p.StartingRelaxRo=0.2;                                                 v.push_back({"start 0.2", p}); }
    { SCFParams p; p.Mixer=M::Pulay; p.PulayDepth=8;                                      v.push_back({"stated Pulay", p}); }
    return v;
}

const char* Name(T t) {return NameOf(t);}

} // namespace

TEST(SCFRecipe, AutoAndInsulatorAreTheOldPolicyInEveryCell)
{
    int cells=0;
    for (SCFRecipe r : {SCFRecipe::Auto, SCFRecipe::Insulator})
    for (T asked : {T::Auto, T::DIIS, T::GDM, T::Ladder, T::Null})
    for (SCFStageKind k : {SCFStageKind::Cold, SCFStageKind::Carried, SCFStageKind::Restored})
    for (double kT : {0.0, 0.005})
    for (int numBands : {0, 8})
    for (bool gSpace : {true, false})
    for (const auto& [vname, p0] : Variants())
    {
        SCFParams p=p0; p.SmearingkT=kT;
        const std::string where=std::string(NameOf(r))+" asked "+Name(asked)+" "+NameOf(k)+" kT "+std::to_string(kT)
                               +" numBands "+std::to_string(numBands)+(gSpace ? " gSpace " : " D-space ")+vname;
        ExpectSame(NewPolicy(r, asked, k, kT, numBands, p, gSpace),
                   OldPolicyWithRestoredCeiling(asked, k==SCFStageKind::Restored, kT, numBands, p, gSpace), where);
        ++cells;
    }
    EXPECT_GT(cells, 0);
}

// §3 rows 2 and 3, cold and carried: Null + Kerker(G0 1)-Pulay (0.45, depth 8, start 5), no MOM.  Restored (row 4): DIIS beside linear
// D-mixing whose start is capped at 0.3 -- the same column as Insulator's (both regimes are always smeared).
TEST(SCFRecipe, NamedRegimesResolveTheirSection3Cells)
{
    for (SCFRecipe r : {SCFRecipe::MagneticInsulator, SCFRecipe::Metal})
    {
        SCFParams p; p.SmearingkT=0.005;
        for (SCFStageKind k : {SCFStageKind::Cold, SCFStageKind::Carried})
        {
            const Choice c=NewPolicy(r, T::Auto, k, p.SmearingkT, 0, p, /*gSpace*/true);
            EXPECT_EQ(c.accel, T::Null) << NameOf(r) << " " << NameOf(k);
            ASSERT_FALSE(c.threw);
            EXPECT_EQ(c.resolved.KerkerG0, 1.0);
            EXPECT_EQ(c.resolved.PulayDepth, 8);
            EXPECT_EQ(c.resolved.PulayStart, 5);
            EXPECT_EQ(c.resolved.StartingRelaxRo, 0.45);
            EXPECT_EQ(c.resolved.Mixer, M::Pulay);
            EXPECT_FALSE(c.resolved.UseMOM);
            EXPECT_TRUE(FillRecipeKnobs(r, k, p).deviations.empty());
        }
        const Choice c=NewPolicy(r, T::Auto, SCFStageKind::Restored, p.SmearingkT, 0, p, true);
        EXPECT_EQ(c.accel, T::DIIS) << NameOf(r) << " restored";
        ASSERT_FALSE(c.threw);
        EXPECT_EQ(c.resolved.Mixer, M::Linear);
        EXPECT_EQ(c.resolved.StartingRelaxRo, 0.3);
        EXPECT_EQ(c.resolved.MaxRelaxRo,      0.3);   // the CEILING too: "not 1.0!"
    }
}

TEST(SCFRecipe, AStatedKnobWinsAndIsADeviation)
{
    SCFParams p; p.SmearingkT=0.005; p.PulayDepth=12; p.UseMOM=true;
    const RecipeKnobs k=FillRecipeKnobs(SCFRecipe::MagneticInsulator, SCFStageKind::Cold, p);
    EXPECT_EQ(k.params.PulayDepth, 12) << "the stated depth wins";
    EXPECT_EQ(k.params.KerkerG0, 1.0)  << "the unstated knobs take the recipe's";
    ASSERT_EQ(k.deviations.size(), 2u);
    EXPECT_EQ(k.deviations[0], "pulayDepth 12");
    EXPECT_NE(k.deviations[1].find("useMOM"), std::string::npos);
    // A stated knob EQUAL to the recipe's is no deviation (the MnO deck's knobs beside its recipe).
    SCFParams same; same.PulayDepth=8; same.PulayStart=5; same.KerkerG0=1.0; same.StartingRelaxRo=0.45;
    EXPECT_TRUE(FillRecipeKnobs(SCFRecipe::MagneticInsulator, SCFStageKind::Cold, same).deviations.empty());
}

// The banked free MnO AFM-II deck (decks/MnO_AFM2_free.json) states Null and the K+P knobs under Auto; `scf.recipe: MagneticInsulator`
// with those knobs removed must build the same run, bit for bit.
TEST(SCFRecipe, TheMnOFreeDeckIsTheMagneticInsulatorRecipe)
{
    SCFParams deck; deck.SmearingkT=5e-3; deck.StartingRelaxRo=0.45; deck.PulayDepth=8; deck.PulayStart=5; deck.KerkerG0=1.0; deck.UseMOM=false;
    SCFParams recipe; recipe.SmearingkT=5e-3;
    for (SCFStageKind k : {SCFStageKind::Cold, SCFStageKind::Carried})
        ExpectSame(NewPolicy(SCFRecipe::MagneticInsulator, T::Auto, k, 5e-3, 0, recipe, true),
                   NewPolicy(SCFRecipe::Auto,              T::Null, k, 5e-3, 0, deck,   true), NameOf(k));
    EXPECT_EQ(FillRecipeKnobs(SCFRecipe::MagneticInsulator, SCFStageKind::Cold, recipe).params.Mixer, M::Auto)
        << "the Mixer stays Auto, as the deck's did: ResolveMixing keeps the hand-set history";
}

TEST(SCFRecipe, MagneticInsulatorAndMetalRefuseKT0AndAMolecule)
{
    for (SCFRecipe r : {SCFRecipe::MagneticInsulator, SCFRecipe::Metal})
    {
        EXPECT_THROW(AdmitRecipe(r, 0.0,   false), std::invalid_argument) << NameOf(r) << " at kT 0 (§3 row 2')";
        EXPECT_NO_THROW(AdmitRecipe(r, 1e-3, false));
        EXPECT_THROW(AdmitRecipe(r, 1e-3,  true),  std::invalid_argument) << NameOf(r) << " on a molecule (ruling (c))";
    }
    for (SCFRecipe r : {SCFRecipe::Auto, SCFRecipe::Insulator})
        for (bool molecule : {false, true}) EXPECT_NO_THROW(AdmitRecipe(r, 0.0, molecule));
}

TEST(SCFRecipe, TheDeckKeyReadsWritesAndRefuses)
{
    using deck::json;
    deck::RunSpec r;
    deck::FromJson(json{{"structure","Si_diamond"}, {"scf",{{"recipe","MagneticInsulator"},{"smearingkT",0.005}}}}, r);
    EXPECT_EQ(r.scfRecipe, SCFRecipe::MagneticInsulator);
    EXPECT_EQ(r.Stages().front().recipe, SCFRecipe::MagneticInsulator);
    const json j=deck::ToJson(r);
    EXPECT_EQ(j.at("scf").at("recipe"), "MagneticInsulator") << "the record restates the recipe";
    deck::RunSpec back; deck::FromJson(j, back); EXPECT_EQ(back.scfRecipe, SCFRecipe::MagneticInsulator);

    deck::RunSpec old; deck::FromJson(json{{"structure","Si_diamond"}, {"scf",{{"pulayDepth",8}}}}, old);
    EXPECT_EQ(old.scfRecipe, SCFRecipe::Auto) << "a deck without the key reads as before (backward compatible)";

    deck::RunSpec sch;
    deck::FromJson(json{{"structure","Si_diamond"}, {"schedule",json::array({
        {{"accelerator","Null"},{"scf",{{"recipe","Metal"},{"smearingkT",0.01}}}},
        {{"scf",{{"recipe","Insulator"}}}}})}}, sch);
    ASSERT_EQ(sch.schedule.size(), 2u);
    EXPECT_EQ(sch.schedule[0].recipe, SCFRecipe::Metal);
    EXPECT_EQ(sch.schedule[1].recipe, SCFRecipe::Insulator);

    deck::RunSpec bad;
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"}, {"scf",{{"recipe","Semiconductor"}}}}, bad), std::runtime_error);

    deck::RunSpec cold; deck::FromJson(json{{"structure","Si_diamond"}, {"scf",{{"recipe","Metal"}}}}, cold);
    EXPECT_THROW(deck::Resolve(cold), std::runtime_error) << "--check refuses a Metal at kT 0";
    deck::RunSpec warm; deck::FromJson(json{{"structure","Si_diamond"}, {"scf",{{"recipe","Metal"},{"smearingkT",0.01}}}}, warm);
    EXPECT_NO_THROW(deck::Resolve(warm));
}
