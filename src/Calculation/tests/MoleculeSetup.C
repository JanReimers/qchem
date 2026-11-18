// File: Calculation/tests/MoleculeSetup.C  The molecular RUN SETUP (doc/CalcNetworkReview.md §7, step 7): the decisions made
// once -- spin promotion, the valence electron count from the pseudopotential set, the parity and accelerator refusals.
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <nlohmann/json.hpp>
import qchem.MoleculeCalculation;
import qchem.Structure;
import qchem.Types;
import qchem.ChargeDensity.Seed;   // SeedStrategy, SeedLibrary (via AtomicDensity)
import qchem.SeedPolicy;           // SeedLibraryCovers: the Auto rule's coverage

using namespace qchem;

namespace
{
Molecule OneAtom(int Z)
{
    Molecule m;
    m.Insert(Site{Z, rvec3_t(0,0,0)});
    return m;
}
} // namespace

TEST(MoleculeSetup, AnOpenShellIsPromotedToPolarizedAndItsParityIsChecked)
{
    const auto s=MoleculeSetup::Resolve(OneAtom(7), {.multiplicity=4, .molecule={.basis="dzvp"}});
    EXPECT_EQ(s->opts.spin, SpinGroup::Polarized);
    EXPECT_EQ(s->Ne, 7);
    EXPECT_FALSE(s->pp) << "no pp named: all-electron";
    const auto refused=MoleculeSetup::Admit(OneAtom(7), {.multiplicity=3, .molecule={.basis="dzvp"}});
    ASSERT_FALSE(refused) << "a triplet of 7 electrons";
    EXPECT_EQ(refused.Error().policy, "SpinMultiplicity");
    EXPECT_THROW((void)MoleculeSetup::Resolve(OneAtom(7), {.multiplicity=3, .molecule={.basis="dzvp"}}), std::invalid_argument) << "the code door throws it";
}

TEST(MoleculeSetup, APseudopotentialSetGivesTheValenceCountAndTheValenceStructure)
{
    const Molecule na=OneAtom(11);
    const auto s=MoleculeSetup::Resolve(na, {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="valence_lowq_sr"}});
    ASSERT_TRUE(s->pp);
    EXPECT_EQ(s->Ne, s->pp->Valence("Na"));
    const auto point=s->PointStructure(na);
    EXPECT_EQ(int(point->GetNumElectrons()), s->Ne);
    EXPECT_EQ(int(s->structure->GetNumElectrons()), 11) << "the reference keeps the TRUE charges; only a point is valence-charged";
}

// A pseudopotential run seeds its OWN pseudo-atoms (OpenWork §3 "A pseudopotential run seeds from the VALENCE library for its q, on every
// structure", lane B 2026-11-12): the setup keys the seed library by each element's PP q, so a PP molecule reads the VALENCE densities at
// that q, as a solid does, and an all-electron one stays unkeyed (the all-electron library).  The tell: Si has no ALL-ELECTRON LDA entry,
// so a Si PP molecule used to take the core guess with a WARN; its q4 VALENCE entry exists, so it now takes SAD.
TEST(MoleculeSetup, APseudopotentialRunSeedsFromTheValenceLibraryAtItsQ)
{
    const Molecule si=OneAtom(14);
    const auto pp=MoleculeSetup::Resolve(si, {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="sipp"}});
    EXPECT_EQ(pp->seedKey.Q(14), pp->pp->Valence("Si")) << "the seed key carries the run's own q";
    EXPECT_EQ(ChargeDensity::SeedLibrary(pp->seedKey, 14), "atomic_valence_densities.json");
    EXPECT_EQ(pp->opts.seed, ChargeDensity::SeedStrategy::SAD) << "the valence library holds Si q4, so the Auto rule seeds from SAD";

    const auto ae=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    EXPECT_EQ(ae->seedKey.Q(8), -1) << "an all-electron run is unkeyed";
    EXPECT_EQ(ChargeDensity::SeedLibrary(ae->seedKey, 8), "atomic_densities.json");
}

// An ALL-ELECTRON molecule past neon seeds from SAD (OpenWork §3 "The default SAD seed stops at Ne", lane B 2026-11-13): the all-electron
// library held Z 1-10 only, so TeF6 on dzvp took the core guess with a WARN (before 3afed2228 it died on an uncaught throw).  It now holds
// every element NIST tabulates, so the Auto rule's coverage holds through uranium and the Te atom's setup picks SAD.
TEST(MoleculeSetup, AnAllElectronMoleculePastNeonSeedsFromSAD)
{
    const ChargeDensity::SeedKey unkeyed("LDA");
    for (int Z=1; Z<=92; Z++) EXPECT_TRUE(SeedLibraryCovers(unkeyed, Z)) << "Z=" << Z;
    const auto te=MoleculeSetup::Resolve(OneAtom(52), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    EXPECT_EQ(ChargeDensity::SeedLibrary(te->seedKey, 52), "atomic_densities.json") << "an all-electron run reads the all-electron library";
    EXPECT_EQ(te->opts.seed, ChargeDensity::SeedStrategy::SAD) << "the library holds Te, so the Auto rule seeds from SAD";
}

TEST(MoleculeSetup, ABasisWithoutAnElementIsRefusedAtResolveNamingTheBasesThatCoverIt)
{
    try
    {
        (void)MoleculeSetup::Resolve(OneAtom(6), {.molecule={.basis="sipp"}});
        FAIL() << "sipp has no carbon: the setup must refuse it before anything is built";
    }
    catch (const std::invalid_argument& e)
    {
        const std::string what=e.what();
        EXPECT_NE(what.find("basis sipp has no block for C"), std::string::npos) << what;
        EXPECT_NE(what.find("dzvp"), std::string::npos) << "the registered bases that cover carbon are named: " << what;
    }
}

// The accelerator is a TYPE in its one home (RunOptions::accelerator, V-FACADE-PARITY P2): an unknown NAME is refused where text becomes a
// type (AcceleratorTypeByName: a GUI form, the atom's json escape hatch; the deck's enum reader), and the type a molecule cannot run, Null,
// at Admit.
TEST(MoleculeSetup, AnAcceleratorAMoleculeCannotRunIsRefusedAtResolve)
{
    const auto refused=MoleculeSetup::Admit(OneAtom(2), {.accelerator=SCFAccelerators::Type::Null, .molecule={.basis="dzvp"}});
    ASSERT_FALSE(refused);
    EXPECT_EQ(refused.Error().policy, "AcceleratorType");
    EXPECT_NE(refused.Error().reason.find("Null"), std::string::npos) << refused.Error().reason;
    try { (void)AcceleratorTypeByName("Anderson", "the test"); ADD_FAILURE() << "an unknown name became a type"; }
    catch (const std::exception& e) { EXPECT_NE(std::string(e.what()).find("Anderson"), std::string::npos) << e.what(); }
}

// EVERY KNOB OF THE MOLECULAR CHAIN REACHES THE FACTORY'S DOOR (2026-11-02; the typed door since V-FACADE-PARITY P0b-1): a molecule that
// states nothing runs what it always ran -- the molecular defaults over the Z-scaled gate -- and every stated knob lands in its field.
TEST(MoleculeSetup, EveryAcceleratorKnobReachesTheFactoryDoor)
{
    const SCFAccelerators::AcceleratorParams def=ResolveAcceleratorParams(AcceleratorOptions{}, MolecularAcceleratorDefaults(0.3));
    EXPECT_EQ(def.diis.Nproj, 4u); EXPECT_EQ(def.diis.EngageBelowFD, 0.3); EXPECT_EQ(def.diis.SVTol, 0.0); EXPECT_EQ(def.diis.SVTolRel, 1e-4);
    EXPECT_EQ(def.gdm.EngageBelowFD, 0.3) << "the standalone GDM gate follows the DIIS gate, as the flat eMax always did";
    EXPECT_EQ(def.gdm.Trust, 0.1); EXPECT_EQ(def.ladderGdm.EngageBelowFD, 1e10) << "the Ladder's GDM rung always steps";
    EXPECT_EQ(def.ladderEThresh, 1e-8); EXPECT_EQ(def.ladderStall, 5); EXPECT_EQ(def.ladderFloor, 1e-8); EXPECT_EQ(def.ladderSwitchAt, 1e-6);
    EXPECT_EQ(def.ladderSignal, SCFAccelerators::ScheduleSignal::EnergyChange);
    AcceleratorOptions a; a.diis.nProj=6; a.ladder.gdm.trust=0.25; a.ladder.eThresh=2e-8; a.ladder.stall=7; a.ladder.floor=3e-8; a.ladder.switchAt=1e-6;
    const SCFAccelerators::AcceleratorParams p=ResolveAcceleratorParams(a, MolecularAcceleratorDefaults(0.5));
    EXPECT_EQ(p.diis.Nproj, 6u); EXPECT_EQ(p.diis.EngageBelowFD, 0.5); EXPECT_EQ(p.ladderGdm.Trust, 0.25); EXPECT_EQ(p.gdm.Trust, 0.1) << "stated for the rung only";
    EXPECT_EQ(p.ladderEThresh, 2e-8); EXPECT_EQ(p.ladderStall, 7); EXPECT_EQ(p.ladderFloor, 3e-8); EXPECT_EQ(p.ladderSwitchAt, 1e-6);
}

// "EMax" was the [F,D] engage gate, never an energy: RETIRED from the accelerator json door (user ruling 2026-11-17), refused BY NAME with
// the new spelling in the message; the door's keys (ToAcceleratorJson) spell it EngageBelowFD.
TEST(MoleculeSetup, TheRetiredEMaxKeyIsRefusedNamingEngageBelowFD)
{
    using SCFAccelerators::Type;
    try { delete SCFAccelerators::Factory(Type::DIIS, nlohmann::json{{"NProj",4},{"EMax",0.1},{"SVTol",0.0}}); FAIL() << "EMax must be refused"; }
    catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("EngageBelowFD"), std::string::npos) << e.what(); }
    const nlohmann::json j=ToAcceleratorJson(MolecularAcceleratorDefaults(0.2));
    EXPECT_FALSE(j.contains("EMax")); EXPECT_EQ(j.at("EngageBelowFD"), 0.2);
    std::unique_ptr<SCFAccelerators::SCFAccelerator> ok(SCFAccelerators::Factory(Type::DIIS, j));
    EXPECT_TRUE(ok);
}
