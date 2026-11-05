// Unit tests of the input deck's data layer (D-ENV step 6): typed options <-> JSON, --set, numbered revisions.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <set>
#include <nlohmann/json.hpp>
#include "forward.H"   // DeckTests: the restart fingerprint of a built calculation (SolidCalculation::Fingerprint is private)
import qchem.Hamiltonian.Factory;   // HubbardProjectorKind (the manifold projector-kind round trip)
import qchem.Deck;
import qchem.SolidCalculation;
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;   // SolidAcceleratorOptions (the acceleratorOptions round trip)
import qchem.Mesh;
import qchem.Materials;
import qchem.ChargeDensity.Seed;
import qchem.LASolver;
import qchem.SolidState;          // StateFingerprint, CompareFingerprints
import qchem.SCFCap;              // BoundedSCF: the vlocRouting test builds the evaluator with one iteration
import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // AccuracyCalibration: the vlocRouting test states eta (no pilot)
import qchem.Reporting;           // the filed run: its grids.audit* rows (the write-only gate)

using namespace qchem;
using qchem::deck::json;
namespace fs = std::filesystem;

namespace {
fs::path TmpDir(const char* tag)
{
    fs::path d=fs::temp_directory_path()/("qchem_deck_"+std::string(tag)+"_"+std::to_string(::getpid()));
    fs::remove_all(d); fs::create_directories(d); return d;
}
}

TEST(Deck, DefaultOptionsRoundTripExactlyAndAreComplete)
{
    SolidCalcOptions o;
    const json j=deck::ToJson(o);
    SolidCalcOptions back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back), j);
    EXPECT_TRUE(j.contains("tolerances")); EXPECT_TRUE(j.at("tolerances").contains("screenEps"));   // a resolved deck lists every field
    EXPECT_EQ(j.at("seed"), "IonicSAD");
    SCFParams p; EXPECT_EQ(deck::ToJson([&]{SCFParams q; deck::FromJson(deck::ToJson(p),q); return q;}()), deck::ToJson(p));
}

// SCFParams gates (2026-11-02): an OFF gate is spelled "off" -- in the deck and in the record it writes (the 1e30 sentinel is retired);
// a number >= 1e30, which every record wrote before, still reads as off; anything else that is not a number is refused.  The step
// floor minRelaxRo rides the same round trip.
TEST(Deck, AnOffGateIsSpelledOffAndTheRetiredSentinelStillReadsAsOff)
{
    SCFParams p;
    deck::FromJson(json{{"minDeltaE","off"},{"minFD",1e30},{"minVirial",4e-2},{"minRelaxRo",0.1}}, p);
    EXPECT_FALSE(p.MinΔE);
    EXPECT_FALSE(p.MinFD) << "the retired sentinel reads as off";
    ASSERT_TRUE(p.MinVirial); EXPECT_EQ(*p.MinVirial, 4e-2);
    EXPECT_EQ(p.MinRelaxRo, 0.1);
    const json j=deck::ToJson(p);
    EXPECT_EQ(j.at("minDeltaE"), "off");
    EXPECT_EQ(j.at("minFD"), "off") << "the record spells it, never 1e30";
    EXPECT_EQ(j.at("minDeltaFD"), "off") << "off by default";
    EXPECT_EQ(j.at("minRelaxRo"), 0.1);
    SCFParams back; deck::FromJson(j, back);
    EXPECT_EQ(deck::ToJson(back), j);
    EXPECT_THROW(deck::FromJson(json{{"minDeltaE","none"}}, p), std::runtime_error);
    EXPECT_THROW(deck::FromJson(json{{"minFD",true}}, p), std::runtime_error);
}

// THE ACCELERATOR CHAIN'S KNOBS ARE DECK KEYS (2026-11-02, OpenWork §3 "Accelerator knobs unreachable from the facades"): every field of
// SolidAcceleratorOptions but the facade's own ladderStartOnGDM, set to a NON-default value, survives FromJson -> ToJson; a deck that
// omits the section keeps the production defaults; a misspelled key is refused (the typo guard reaches the nested objects).
TEST(Deck, AcceleratorOptionsAreDeckKeysAndRoundTrip)
{
    const json a={{"diis",{{"nProj",6},{"engageBelowFD",0.05},{"svTol",2e-9},{"svTolRel",3e-4}}},
                  {"gdm",{{"engageBelowFD",0.5},{"trust",0.2},{"trustBackoff",0.3},{"trustMin",2e-4},{"pcFloor",0.05},{"aufbauTol",1e-5},{"aufbauMaxRepicks",4}}},
                  {"ladder",{{"gdm",{{"engageBelowFD",2.0},{"trust",0.4},{"trustBackoff",0.5},{"trustMin",1e-3},{"pcFloor",0.2},{"aufbauTol",2e-5},{"aufbauMaxRepicks",3}}},
                             {"eThresh",2e-8},{"stall",7},{"floor",3e-8},{"switchAt",2e-6},{"signal","Error"}}}};
    SCFAccelerators::SolidAcceleratorOptions o; deck::FromJson(a, o);
    EXPECT_EQ(o.diis.Nproj, 6u); EXPECT_EQ(o.diis.SVTolRel, 3e-4); EXPECT_EQ(o.gdm.AufbauMaxRepicks, 4); EXPECT_EQ(o.ladderGdm.Trust, 0.4);
    EXPECT_EQ(o.ladderStall, 7); EXPECT_EQ(o.ladderSwitchAt, 2e-6); EXPECT_EQ(o.ladderSignal, SCFAccelerators::ScheduleSignal::Error);
    EXPECT_EQ(deck::ToJson(o), a) << "every key written back, none dropped or renamed";
    const SCFAccelerators::SolidAcceleratorOptions def;
    deck::RunSpec r; deck::FromJson(json{{"structure","Si_diamond"}}, r);
    EXPECT_EQ(deck::ToJson(r.accel), deck::ToJson(def)) << "an omitted section is the production recipe";
    deck::FromJson(json{{"structure","Si_diamond"},{"acceleratorOptions",a}}, r);
    EXPECT_EQ(deck::ToJson(r).at("acceleratorOptions"), a) << "the record carries the run's knobs";
    EXPECT_THROW(deck::FromJson(json{{"diis",{{"nproj",6}}}}, o), std::runtime_error);
    EXPECT_THROW(deck::FromJson(json{{"ladder",{{"signal","Energy"}}}}, o), std::runtime_error);
}

TEST(Deck, NonDefaultValuesSurviveTheRoundTrip)
{
    SolidCalcOptions o;
    o.Nelec=48; o.multiplicity=1; o.species={{"Mn",15},{"O",6}}; o.kShift=rvec3_t(0.25,0,0.5);
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::CholeskyPivoted; o.imposeSymmetry=true; o.siteSpins={1,-1,0,0};
    o.siteCharges={2,2,-2,-2};
    o.tolerances.screenEps=1e-8; o.tolerances.mgridEcuts={53.3,17.8}; o.xcMesh.beckeEps=1e-7; o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;
    o.hubbard.push_back(HubbardU_OrthoAtomic(0,2,4.0));
    const json j=deck::ToJson(o);
    SolidCalcOptions back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back), j);
    EXPECT_EQ(back.species.size(),2u); EXPECT_EQ(back.species[0].first,"Mn"); EXPECT_EQ(back.tolerances.screenEps,1e-8);
    EXPECT_EQ(back.siteCharges,(std::vector<int>{2,2,-2,-2}));
    ASSERT_EQ(back.hubbard.size(),1u); EXPECT_TRUE(back.hubbard[0].orthoAtomic); EXPECT_NEAR(back.hubbard[0].U,4.0/27.211386245988,1e-15);
}

// HubbardProjectorPlan D7 + increment 5: a manifold names its projector kind; the DEFAULT is pseudo-atomic, siteShells
// is said, the name round-trips, and a misspelled kind is refused.  The retired site-shell contraction: an old record's
// inert defaults ("radial": [], "atomicRadial": false) still load; asking for a contraction is refused, naming the fix.
TEST(Deck, AHubbardManifoldNamesItsProjectorKind)
{
    using qchem::Hamiltonian::HubbardProjectorKind;
    SolidCalcOptions o;
    deck::FromJson(json{{"hubbard",{{{"site",0},{"l",2},{"U_eV",4.0}},{{"site",1},{"l",1},{"projector","siteShells"}}}}}, o);
    ASSERT_EQ(o.hubbard.size(),2u);
    EXPECT_EQ(o.hubbard[0].projector, HubbardProjectorKind::PseudoAtomic) << "the default since increment 5";
    EXPECT_EQ(o.hubbard[1].projector, HubbardProjectorKind::SiteShells);
    const json j=deck::ToJson(o).at("hubbard");
    EXPECT_EQ(j[0].at("projector"), "pseudoAtomic");
    EXPECT_EQ(j[1].at("projector"), "siteShells");
    SolidCalcOptions back; deck::FromJson(json{{"hubbard",j}}, back);
    EXPECT_EQ(back.hubbard[1].projector, HubbardProjectorKind::SiteShells);
    EXPECT_THROW(deck::FromJson(json{{"hubbard",{{{"site",0},{"projector","pseudoatomic"}}}}}, o), std::runtime_error);
    EXPECT_NO_THROW(deck::FromJson(json{{"hubbard",{{{"site",0},{"radial",json::array()},{"atomicRadial",false}}}}}, o));
    EXPECT_THROW(deck::FromJson(json{{"hubbard",{{{"site",0},{"atomicRadial",true}}}}}, o), std::runtime_error);
    EXPECT_THROW(deck::FromJson(json{{"hubbard",{{{"site",0},{"radial",{0.5,0.5}}}}}}, o), std::runtime_error);
}

TEST(Deck, AMissingKeyKeepsTheDefaultAndAnUnknownKeyThrowsNamingTheLegalOnes)
{
    SolidCalcOptions o; deck::FromJson(json{{"Nelec",8}}, o);
    EXPECT_EQ(o.Nelec,8); EXPECT_EQ(o.cutoffFactor,2.0);
    try { deck::FromJson(json{{"tolerances",{{"screenEPS",1e-8}}}}, o); FAIL() << "a typo'd key must throw"; }
    catch (const std::runtime_error& e)
    {
        const std::string m=e.what();
        EXPECT_NE(m.find("tolerances.screenEPS"),std::string::npos) << m;
        EXPECT_NE(m.find("screenEps"),std::string::npos) << "the message names the legal keys: " << m;
    }
    EXPECT_THROW(deck::FromJson(json{{"seed","IonicSAd"}}, o), std::runtime_error);      // an enum name is exact
    EXPECT_THROW(deck::FromJson(json{{"kShift",{1,2}}}, o), std::runtime_error);
    EXPECT_THROW(deck::FromJson(json{{"species",{{"Si"}}}}, o), std::runtime_error);
}

TEST(Deck, SetAppliesADottedPathAndParsesTheValue)
{
    json d=deck::ToJson(SolidCalcOptions{});
    deck::ApplySet(d,"tolerances.screenEps=1e-8");
    deck::ApplySet(d,"imposeSymmetry=true");
    deck::ApplySet(d,"seed=Uniform");                 // a bare word is a string
    deck::ApplySet(d,"kShift.1=0.5");                 // an array index
    SolidCalcOptions o; deck::FromJson(d,o);
    EXPECT_EQ(o.tolerances.screenEps,1e-8); EXPECT_TRUE(o.imposeSymmetry);
    EXPECT_EQ(o.seed,qchem::ChargeDensity::SeedStrategy::Uniform); EXPECT_EQ(o.kShift.y,0.5);
    EXPECT_THROW(deck::ApplySet(d,"noequals"),std::runtime_error);
    EXPECT_THROW(deck::ApplySet(d,"a..b=1"),std::runtime_error);
    EXPECT_THROW(deck::ApplySet(d,"kShift.7=1"),std::runtime_error);
    json bad=d; deck::ApplySet(bad,"tolerances.nope=1");
    EXPECT_THROW(deck::FromJson(bad,o),std::runtime_error) << "a --set typo is caught when the result is read";
}

TEST(Deck, RevisionsAreNumberedNeverReusedAndTheFileIsTheCompleteRecord)
{
    const fs::path dir=TmpDir("rev");
    const fs::path a=deck::ClaimRevision(dir,"MnO"), b=deck::ClaimRevision(dir,"MnO"), c=deck::ClaimRevision(dir,"NaF");
    EXPECT_EQ(a.filename(),"MnO.r001.json"); EXPECT_EQ(b.filename(),"MnO.r002.json"); EXPECT_EQ(c.filename(),"NaF.r001.json");
    EXPECT_TRUE(fs::exists(a)) << "a claim creates the file, so a parallel job sees the number as taken";

    SolidCalcOptions o; o.Nelec=96; o.tolerances.screenEps=1e-8;
    deck::Provenance pv; pv.commandLine="gpwprobe --deck in.json --set tolerances.screenEps=1e-8"; pv.overrides={"tolerances.screenEps=1e-8"};
    pv.codeVersion="abc123"; pv.ignoredEnvironment={"GPW_SCREEN_EPS"};
    deck::WriteRevision(b, json{{"solid",deck::ToJson(o)}}, pv);

    const json raw=json::parse(std::ifstream(b));
    EXPECT_EQ(raw.at("deck").at("schema"),deck::kSchemaVersion);
    EXPECT_EQ(raw.at("provenance").at("overrides")[0],"tolerances.screenEps=1e-8");
    EXPECT_EQ(raw.at("provenance").at("ignoredEnvironment")[0],"GPW_SCREEN_EPS");
    EXPECT_TRUE(raw.at("provenance").contains("activeEnvironment"));
    EXPECT_FALSE(raw.at("provenance").contains("parent")) << "no input deck, no parent";

    // a run that starts FROM a revision names it as its parent
    deck::Provenance child=pv; child.inputDeck=b;
    const fs::path c2=deck::ClaimRevision(dir,"MnO");
    deck::WriteRevision(c2, json{{"solid",deck::ToJson(o)}}, child);
    EXPECT_EQ(json::parse(std::ifstream(c2)).at("provenance").at("parent"),"MnO.r002.json");

    const json run=deck::LoadDeck(b,"abc123");          // header stripped; same code version = silent
    SolidCalcOptions back; deck::FromJson(run.at("solid"),back);
    EXPECT_EQ(back.Nelec,96); EXPECT_EQ(back.tolerances.screenEps,1e-8);
    fs::remove_all(dir);
}

TEST(Deck, ALaterSchemaIsRefusedAndAHandWrittenDeckIsThePayload)
{
    const fs::path dir=TmpDir("schema");
    { std::ofstream(dir/"new.json") << R"({"deck":{"schema":999},"run":{}})"; }
    EXPECT_THROW(deck::LoadDeck(dir/"new.json","x"),std::runtime_error);
    { std::ofstream(dir/"hand.json") << R"({"solid":{"Nelec":8}})"; }
    EXPECT_EQ(deck::LoadDeck(dir/"hand.json","x").at("solid").at("Nelec"),8);
    fs::remove_all(dir);
}

// ---- the structure section is a NAME ----
TEST(Deck, TheStructureIsANameAndTheMaterialSuppliesSpeciesAndElectrons)
{
    deck::RunSpec r; deck::FromJson(json{{"structure","Si_diamond"},{"scf",{{"NMaxIter",50}}}}, r);
    EXPECT_EQ(r.scf.NMaxIter,50u);
    EXPECT_EQ(r.solid.Nelec,0); EXPECT_TRUE(r.solid.species.empty()) << "unstated = derive";
    const auto m=deck::Resolve(r);
    EXPECT_EQ(r.solid.Nelec,8); ASSERT_EQ(r.solid.species.size(),1u); EXPECT_EQ(r.solid.species[0].first,"Si");
    EXPECT_EQ(m.cell->GetNumAtoms(),2u);
    const json resolved=deck::ToJson(r);                       // the record: complete, and it reloads to the same thing
    EXPECT_EQ(resolved.at("structure"),"Si_diamond"); EXPECT_EQ(resolved.at("solid").at("Nelec"),8);
    deck::RunSpec again; deck::FromJson(resolved,again); EXPECT_EQ(deck::ToJson(again),resolved);
}

// A stated per-element valence (a SEMICORE pseudopotential) changes the electron count the run derives: NaF with Na q9 holds 9+7
// electrons, not the material's 1+7 (2026-10-31: until then Nelec was taken from the material's species, a cell charged +8, silently).
TEST(Deck, AStatedSpeciesValenceSetsTheDerivedElectronCount)
{
    deck::RunSpec r; deck::FromJson(json{{"structure","NaF_rocksalt"},{"solid",{{"species",json::array({json::array({"Na",9}),json::array({"F",7})})}}}}, r);
    deck::Resolve(r);
    EXPECT_EQ(r.solid.Nelec,16);
    deck::RunSpec plain; deck::FromJson(json{{"structure","NaF_rocksalt"}}, plain);
    deck::Resolve(plain);
    EXPECT_EQ(plain.solid.Nelec,8) << "the material's own q1 Na";
}

TEST(Deck, AStatedValenceOverridesTheMaterialAndBadStructuresAreRefused)
{
    deck::RunSpec r; deck::FromJson(json{{"structure","Si_diamond"},{"solid",{{"Nelec",6}}}}, r);
    deck::Resolve(r); EXPECT_EQ(r.solid.Nelec,6) << "an explicit Nelec (a charged cell) is honoured";
    EXPECT_THROW(deck::FromJson(json{{"solid",json::object()}}, r), std::runtime_error);                 // structure is required
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"lattice",{{"a",10.0}}}}, r), std::runtime_error);   // no restating the structure
    deck::RunSpec bad; bad.structure="NoSuchMaterial";
    try { deck::Resolve(bad); FAIL(); } catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("Si_diamond"),std::string::npos) << "lists the known names"; }
    deck::RunSpec mol; mol.structure=StructureData::MoleculeNames().front();
    EXPECT_THROW(deck::Resolve(mol), std::runtime_error);
}

// ---- a deck RUNS, and runs exactly what the options say (the Si CP2K anchor, GPW_Si.Γ_Imp_CP2K's recipe) ----
import qchem.Lattice_3D;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
TEST(Deck, ARunFromADeckEqualsTheSameRunBuiltByHand_Si_Gamma)
{
    const json payload={
        {"structure","Si_diamond"}, {"basis",{{"data","SIPP_SR"}}},
        {"solid",{{"densityEcut",20.0},{"imposeSymmetry",true}}},
        {"scf",{{"NMaxIter",60},{"minDeltaRhoPerElectron",1e-3},{"minDeltaE",1e-6},{"minDeltaFD","off"},{"minVirial","off"},{"minFD","off"},{"startingRelaxRo",0.3},{"mergeTol",1e-4}}}};
    deck::RunSpec spec; deck::FromJson(payload,spec);

    const fs::path dir=TmpDir("run");
    deck::Provenance pv; pv.commandLine="test"; pv.codeVersion="t";
    const deck::RunOutcome out=deck::Run(spec,pv,dir);
    ASSERT_TRUE(out.converged) << out.summary;
    ASSERT_TRUE(out.energy.has_value());
    EXPECT_NEAR(*out.energy, -7.11506, 2e-3) << "the CP2K FCC-Si Gamma anchor (grid-gap tolerance)";

    // the record exists, is r001, and reloads to the SAME deck (defaults filled in, structure still a name)
    EXPECT_EQ(out.revision.filename(),"Si_diamond.r001.json");
    deck::RunSpec again; deck::FromJson(deck::LoadDeck(out.revision,"t"),again);
    EXPECT_EQ(again.structure,"Si_diamond"); EXPECT_EQ(again.solid.Nelec,8) << "the resolved deck records the derived electron count";

    // the SAME run assembled by hand: bit-identical energy, i.e. the deck path adds nothing
    const auto mat=Materials::Get("Si_diamond");
    Lattice_3D lat(*mat.cell, ivec3_t(1,1,1));
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, mat.cell.get()));
    SolidCalcOptions o; o.Nelec=8; o.species=mat.species; o.densityEcut=20.0; o.imposeSymmetry=true; o.label=out.revision.stem().string();
    SolidCalculation calc(SolidSetup::Resolve(lat, mol, o), spec.scf);
    auto R=calc.Result(); ASSERT_TRUE(R);
    EXPECT_DOUBLE_EQ(R->Energy(), *out.energy);
    fs::remove_all(dir);
}

TEST(Deck, ScfAndScheduleTogetherAreAmbiguousAndAScheduleRoundTrips)
{
    deck::RunSpec r;
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"scf",json::object()},{"schedule",json::array({json::object()})}},r),std::runtime_error);
    deck::FromJson(json{{"structure","Si_diamond"},{"kmesh",{2,2,2}},{"schedule",{{{"accelerator","GDM"},{"scf",{{"smearingkT",0.01}}}},{{"accelerator","DIIS"}}}}},r);
    ASSERT_EQ(r.schedule.size(),2u); EXPECT_EQ(r.schedule[0].scf.SmearingkT,0.01); EXPECT_EQ(r.kmesh.x,2);
    const json j=deck::ToJson(r); EXPECT_TRUE(j.contains("schedule")); EXPECT_FALSE(j.contains("scf"));
    deck::RunSpec back; deck::FromJson(j,back); EXPECT_EQ(deck::ToJson(back),j);
}

// ---- D-ENV 6b: the declared CP2K deviations are the deck's `policy` block; only STATED routes are recorded ----
import qchem.RunPolicy;
TEST(Deck, PolicyRecordsOnlyWhatWasStatedSoTheUmbrellaStillMeansWhatItSays)
{
    SolidCalcOptions o;
    EXPECT_EQ(deck::ToJson(o).at("policy"), json({{"cp2kCompat",false}})) << "nothing stated: only the umbrella";
    deck::FromJson(json{{"policy",{{"cp2kCompat",true},{"streamFold",true}}}}, o);
    EXPECT_TRUE(o.policy.cp2kCompat); ASSERT_TRUE(o.policy.streamFold.has_value()); EXPECT_TRUE(*o.policy.streamFold);
    EXPECT_FALSE(o.policy.dmLowRank.has_value());
    EXPECT_EQ(deck::ToJson(o).at("policy"), json({{"cp2kCompat",true},{"streamFold",true}}));
    qchem::RunPolicy p(o.policy);
    EXPECT_TRUE(p.StreamFold()) << "stated beats the umbrella"; EXPECT_FALSE(p.DMLowRank()) << "unstated follows the umbrella";
    EXPECT_THROW(deck::FromJson(json{{"policy",{{"cp2kcompat",true}}}}, o), std::runtime_error);   // a typo'd route is refused
    // 6c: the xcFromDM route's controls ride the policy block (were QCHEM_XC_DM_MIX / _BOOST)
    SolidCalcOptions c; deck::FromJson(json{{"policy",{{"xcFromDM",true},{"xcDMMix",1.0},{"xcDMBoost",2.0}}}}, c);
    qchem::RunPolicy pc(c.policy);
    EXPECT_EQ(pc.XCDMMixOverride(),1.0); EXPECT_EQ(pc.XCDMBoost(),2.0); EXPECT_TRUE(pc.XCFromDM());
    EXPECT_EQ(qchem::RunPolicy{}.XCDMMixOverride(),-1.0) << "unset = no override"; EXPECT_EQ(qchem::RunPolicy{}.XCDMBoost(),1.0);
    EXPECT_EQ(deck::ToJson(c).at("policy"), json({{"cp2kCompat",false},{"xcFromDM",true},{"xcDMMix",1.0},{"xcDMBoost",2.0}}));
    json d=deck::ToJson(SolidCalcOptions{}); deck::ApplySet(d,"policy.cp2kCompat=true"); deck::ApplySet(d,"policy.beckeXC=true");
    SolidCalcOptions viaSet; deck::FromJson(d,viaSet);
    EXPECT_TRUE(viaSet.policy.cp2kCompat); EXPECT_TRUE(viaSet.policy.beckeXC.value_or(false)) << "--set policy.<route> states it";
}

// ---- D-ENV 6d.1: eV in the file / a.u. in RAM, basis trim+autoTrim, saved states with lineage ----
TEST(Deck, UIsEvInTheFileAndAtomicUnitsInRamAndTheRecordReproducesTheBits)
{
    SolidCalcOptions o;
    deck::FromJson(json{{"hubbard",{{{"site",0},{"l",2},{"U_eV",4.0},{"Uirrep_eV",{3.0,4.5}},{"alpha_eV",0.1},{"alphaUp_eV",0.25},{"alphaSlot_eV",{0.0,0.05}}}}}}, o);
    ASSERT_EQ(o.hubbard.size(),1u);
    EXPECT_DOUBLE_EQ(o.hubbard[0].U, 4.0/27.211386245988) << "RAM is atomic units";
    EXPECT_DOUBLE_EQ(o.hubbard[0].Uirrep[1], 4.5/27.211386245988);
    const json j=deck::ToJson(o).at("hubbard")[0];
    EXPECT_EQ(j.at("U_eV"),4.0) << "the record writes the eV the user wrote (shortest value that converts back to the same double)";
    EXPECT_EQ(j.at("Uirrep_eV")[1],4.5);
    SolidCalcOptions back; deck::FromJson(deck::ToJson(o),back);
    EXPECT_EQ(back.hubbard[0].U, o.hubbard[0].U) << "bit-for-bit: a revision file re-runs with the SAME U";
    EXPECT_EQ(back.hubbard[0].Uirrep, o.hubbard[0].Uirrep);
    EXPECT_EQ(back.hubbard[0].alpha, o.hubbard[0].alpha);
    EXPECT_DOUBLE_EQ(o.hubbard[0].alphaUp, 0.25/27.211386245988); EXPECT_EQ(o.hubbard[0].alphaDown, 0.0);
    EXPECT_EQ(back.hubbard[0].alphaUp, o.hubbard[0].alphaUp);
    ASSERT_EQ(o.hubbard[0].alphaSlot.size(), 2u); EXPECT_DOUBLE_EQ(o.hubbard[0].alphaSlot[1], 0.05/27.211386245988);
    EXPECT_EQ(back.hubbard[0].alphaSlot, o.hubbard[0].alphaSlot);
    EXPECT_THROW(deck::FromJson(json{{"hubbard",{{{"site",0},{"U_Ha",0.1}}}}},o), std::runtime_error) << "the old spelling is refused, not guessed";
    // an arbitrary value round-trips exactly too
    for (double eV : {3.1,5.27,0.0,1e-3,7.123456789})
    {
        SolidCalcOptions a; deck::FromJson(json{{"hubbard",{{{"site",0},{"U_eV",eV}}}}},a);
        SolidCalcOptions b; deck::FromJson(deck::ToJson(a),b);
        EXPECT_EQ(a.hubbard[0].U,b.hubbard[0].U) << eV;
    }
}

TEST(Deck, BasisTrimVetAndStateRoundTripAndAreValidated)
{
    deck::RunSpec r;
    deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"data","VALENCE_LOWQ_VA"},{"spherical",true},{"trim",{{{"Z",25},{"l",0},{"alpha",0.06}}}}}},
                        {"state",{{"save","auto"},{"restartFrom","MnO_AFM2.r003"}}}}, r);
    ASSERT_EQ(r.basis.trim.size(),1u); EXPECT_EQ(r.basis.trim[0].Z,25); EXPECT_EQ(r.basis.trim[0].alpha,0.06);
    EXPECT_EQ(r.state.save,"auto"); EXPECT_EQ(r.state.restartFrom,"MnO_AFM2.r003");
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back); EXPECT_EQ(deck::ToJson(back),j);
    EXPECT_FALSE(deck::ToJson(deck::RunSpec{}).contains("state")) << "no state requested: no state block";

    deck::RunSpec both; deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"autoTrim",true},{"trim",{{{"Z",25},{"l",0},{"alpha",0.06}}}}}}}, both);
    EXPECT_THROW(deck::Resolve(both), std::runtime_error) << "autoTrim and a stated trim are exclusive";
    deck::RunSpec noTol; deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"autoTrim",true}}}}, noTol);
    EXPECT_THROW(deck::Resolve(noTol), std::runtime_error) << "an autoTrim with no orthoTol would drop nothing";
    deck::RunSpec retired;
    EXPECT_THROW(deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"vet",true}}}}, retired), std::runtime_error) << "basis.vet is retired, not ignored";
    deck::RunSpec sph; deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"data","VALENCE_LOWQ_SR"},{"spherical",true}}}}, sph);
    try { deck::Resolve(sph); FAIL() << "spherical + the SR transition-metal block must be refused"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("VALENCE_LOWQ_VA"),std::string::npos) << e.what(); }
    deck::RunSpec ok; deck::FromJson(json{{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"},{"spherical",true}}}}, ok);
    EXPECT_NO_THROW(deck::Resolve(ok)) << "the refusal is specific to a transition metal on the SR block";
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"state",{{"restore","x"}}}}, ok), std::runtime_error);
}

// basis.file: a path overrides the registered data file -- recorded in the resolved deck, honoured by the build, and a bad path fails loudly.
TEST(Deck, BasisFileOverridesTheRegisteredDataFileAndABadPathFails)
{
    const std::string va = (std::filesystem::path(DECKS_DIR).parent_path()/"BasisSetData"/"valence_lowq_va.bsd").string();
    deck::RunSpec reg, over;
    deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"data","VALENCE_LOWQ_SR"}}}}, reg);
    deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"data","VALENCE_LOWQ_SR"},{"file",va}}}}, over);
    EXPECT_FALSE(deck::ToJson(reg)["basis"].contains("file")) << "no override: no key";
    EXPECT_EQ(deck::ToJson(over)["basis"]["file"].get<std::string>(), va) << "the record carries the file actually used";
    deck::RunSpec back; deck::FromJson(deck::ToJson(over),back); EXPECT_EQ(deck::ToJson(back),deck::ToJson(over));
    const auto nfun=[](const deck::RunSpec& r){ return deck::Vet(r,std::filesystem::temp_directory_path()).summary; };
    EXPECT_NE(nfun(reg), nfun(over)) << "the override built a different basis (SR vs VA files)";
    deck::RunSpec bad; deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"file","/nonexistent/x.bsd"}}}}, bad);
    EXPECT_ANY_THROW(deck::Vet(bad,std::filesystem::temp_directory_path()));
}

// D-BASIS-REGISTRY: basis.data is a file STEM, any case -- an unregistered file is usable by name; a name that is not in the data
// directory is refused when the deck is read, listing the names; the record keeps the name as written.
TEST(Deck, TheBasisIsAFileStemAndAnUnknownOneIsRefusedListingTheNames)
{
    deck::RunSpec r;
    EXPECT_NO_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"basis",{{"data","valence_lowq_tm"}}}}, r)) << "unregistered, still a name";
    EXPECT_EQ(deck::ToJson(r)["basis"]["data"], "valence_lowq_tm");
    EXPECT_NO_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"basis",{{"data","sipp_sr"}}}}, r)) << "any case";
    try { deck::FromJson(json{{"structure","Si_diamond"},{"basis",{{"data","SIPP_SRR"}}}}, r); FAIL() << "a typo must throw"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("sipp_sr"), std::string::npos) << e.what(); }
}

TEST(Deck, ASavedStateRestartsAndTheRevisionNamesItsParent)
{
    const fs::path dir=TmpDir("restart");
    const json base={{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"}}},{"solid",{{"densityEcut",20.0},{"imposeSymmetry",true}}},
        {"scf",{{"NMaxIter",60},{"minDeltaRhoPerElectron",1e-6},{"minDeltaE",1e-10},{"minDeltaFD","off"},{"minVirial","off"},{"minFD","off"},{"startingRelaxRo",0.3}}}};
    json d1=base; d1["state"]={{"save","auto"}};
    deck::RunSpec s1; deck::FromJson(d1,s1);
    deck::Provenance pv; pv.codeVersion="t";
    const auto r1=deck::Run(s1,pv,dir);
    ASSERT_TRUE(r1.converged) << r1.summary;
    ASSERT_TRUE(fs::exists(dir/"states"/"Si_diamond.r001.h5")) << "save:auto writes states/<stem>.h5";

    json d2=base; d2["state"]={{"restartFrom","Si_diamond.r001"}}; d2.erase("scf");
    d2["schedule"]={{{"accelerator","GDM"},{"scf",base.at("scf")}}};          // a restart runs the FINAL stage; GDM starts at the answer
    deck::RunSpec s2; deck::FromJson(d2,s2);
    const auto r2=deck::Run(s2,pv,dir);
    ASSERT_TRUE(r2.converged) << r2.summary;
    EXPECT_EQ(r2.revision.filename(),"Si_diamond.r002.json");
    EXPECT_NEAR(*r2.energy,*r1.energy,1e-8) << "a restart from the converged state stays at the answer";
    EXPECT_EQ(json::parse(std::ifstream(r2.revision)).at("provenance").at("restartedFrom"),"Si_diamond.r001.json");
    // D-BASIS-REGISTRY: the record names the basis FILE and the checksum of its contents (a name is not an identity).
    const json pv2=json::parse(std::ifstream(r2.revision)).at("provenance");
    EXPECT_EQ(pv2.at("basisFile"),"sipp_sr.bsd");
    EXPECT_EQ(pv2.at("basisChecksum").get<std::string>().size(),16u);

    json d3=base; d3["state"]={{"restartFrom","Si_diamond.r099"}};
    deck::RunSpec s3; deck::FromJson(d3,s3);
    EXPECT_THROW(deck::Run(s3,pv,dir),std::runtime_error) << "a missing state is refused before any revision is claimed";
    EXPECT_FALSE(fs::exists(dir/"Si_diamond.r003.json"));
    fs::remove_all(dir);
}

// ---- D-ENV 6d.2: postSCF actions: JSON, pre-flight validation, and a real run ----
TEST(Deck, PostSCFRoundTripsAndIsStrict)
{
    deck::RunSpec r;
    deck::FromJson(json{{"structure","NiO_AFM2"},{"postSCF",{
        {{"estimateHubbardU",json::object()}},
        {{"hubbardLoop",{{"maxOuter",8},{"tolU_eV",1e-3}}}},
        {{"independentResponse",{{"nq",2}}}},
        {{"hubbardLinearResponse",{{"nq",{2,1,1}},{"perturb",{0,1}},{"maxIter",300},{"restart",30},{"tol",1e-9}}}},
        {{"hubbardFiniteDifference",{{"perturb",{0}},{"alpha_eV",0.1}}}}}}}, r);
    ASSERT_EQ(r.postSCF.size(),5u);
    EXPECT_EQ(r.postSCF[1].maxOuter,8u); EXPECT_EQ(r.postSCF[2].nq,ivec3_t(2,2,2)); EXPECT_EQ(r.postSCF[3].perturb, (std::vector<size_t>{0,1}));
    EXPECT_EQ(r.postSCF[3].nq,ivec3_t(2,1,1)) << "nq as [n1,n2,n3], like kmesh: a non-cubic q-mesh (the supercell gate's 2x1x1)";
    EXPECT_EQ(deck::ToJson(r)["postSCF"][2]["independentResponse"]["nq"],json(2)) << "a cubic q-mesh is written back as the integer shorthand";
    EXPECT_DOUBLE_EQ(r.postSCF[4].alpha, 0.1/27.211386245988) << "alpha is eV in the file, Hartree in RAM";
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_FALSE(deck::ToJson(deck::RunSpec{}).contains("postSCF"));

    deck::RunSpec bad;
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLop",json::object()}}}}},bad), std::runtime_error) << "an unknown action name";
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLoop",{{"maxOuters",8}}}}}}},bad), std::runtime_error) << "an unknown parameter";
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLoop",json::object()},{"estimateHubbardU",json::object()}}}}},bad), std::runtime_error) << "two actions in one entry";
}

// R4: the outer loop's estimator.  linearResponse reads the response keys; under ACBN0 they are UNKNOWN keys (a recipe typo cannot be silently ignored).
TEST(Deck, HubbardLoopEstimatorRoundTripsAndOwnsItsKeys)
{
    deck::RunSpec r;
    deck::FromJson(json{{"structure","NiO_AFM2"},{"postSCF",{
        {{"hubbardLoop",{{"estimator","linearResponse"},{"nq",2},{"perturb",{0,1}},{"tol",1e-9},{"maxOuter",6}}}},
        {{"hubbardLoop",json::object()}}}}}, r);
    ASSERT_EQ(r.postSCF.size(),2u);
    using Est=SolidCalculation::HubbardLoop::Estimator;
    EXPECT_EQ(r.postSCF[0].estimator, Est::LinearResponse);
    EXPECT_EQ(r.postSCF[0].nq, ivec3_t(2,2,2));
    EXPECT_EQ(r.postSCF[0].perturb, (std::vector<size_t>{0,1}));
    EXPECT_EQ(r.postSCF[1].estimator, Est::ACBN0) << "ACBN0 stays the default: no deck moves";
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_FALSE(j["postSCF"][1]["hubbardLoop"].contains("nq")) << "an ACBN0 loop writes no response keys";

    deck::RunSpec bad;
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLoop",{{"nq",2}}}}}}},bad), std::runtime_error) << "nq under ACBN0";
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLoop",{{"estimator","LR"}}}}}}},bad), std::runtime_error) << "an unknown estimator";

    auto msg=[](json d){ try { deck::RunSpec s; deck::FromJson(d,s); deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };
    const json hub={{"hubbard",{{{"site",0},{"l",2},{"U_eV",0.0}}}}};
    json cplx=hub; cplx["forceComplex"]=true;
    auto deckWith=[](json solid, json loop, json kmesh){ return json{{"structure","NiO_AFM2"},{"kmesh",kmesh},{"solid",solid},{"postSCF",{{{"hubbardLoop",loop}}}}}; };
    const json lr={{"estimator","linearResponse"}};
    EXPECT_NE(msg(deckWith(hub,lr,{1,1,1})).find("forceComplex"),std::string::npos) << "an LR loop runs the response: complex ansatz";
    EXPECT_EQ(msg(deckWith(cplx,lr,{1,1,1})),"NO THROW");
    EXPECT_NE(msg(deckWith(cplx,{{"estimator","linearResponse"},{"nq",2}},{3,3,3})).find("incommensurate"),std::string::npos);
    EXPECT_EQ(msg(deckWith(hub,json::object(),{1,1,1})),"NO THROW") << "an ACBN0 loop needs no complex ansatz";
    json smeared={{"structure","NiO_AFM2"},{"kmesh",{2,2,2}},{"solid",cplx},{"scf",{{"smearingkT",0.01}}},{"postSCF",{{{"independentResponse",{{"nq",2}}}}}}};
    EXPECT_EQ(msg(smeared),"NO THROW") << "one Fermi level across k is the default (D-FERMI-DEFAULT)";
    smeared["solid"]["perKFermiLevels"]=true;
    EXPECT_NE(msg(smeared).find("perKFermiLevels"),std::string::npos) << "a smeared q-mesh response needs one Fermi level";
    smeared["solid"].erase("perKFermiLevels"); smeared["solid"]["globalFermi"]=true;
    EXPECT_NE(msg(smeared).find("RETIRED"),std::string::npos) << "the old key is refused, never re-read with the opposite default";
}

TEST(Deck, PostSCFPreflightRejectsADeckThatCannotWorkBeforeAnySCF)
{
    auto resolve=[](json d){ deck::RunSpec s; deck::FromJson(d,s); deck::Resolve(s); };
    const json H={{"site",0},{"l",2},{"U_eV",4.0}};
    auto base=[&](json post, json solid, json kmesh=json::array({1,1,1}))
    { return json{{"structure","NiO_AFM2"},{"kmesh",kmesh},{"solid",solid},{"postSCF",post}}; };
    const json hub={{"hubbard",{H}}};
    auto msg=[&](json d){ try { resolve(d); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };

    EXPECT_NE(msg(base({{{"estimateHubbardU",json::object()}}},json::object())).find("no Hubbard manifold"),std::string::npos);
    EXPECT_EQ(msg(base({{{"estimateHubbardU",json::object()}}},hub)),"NO THROW");
    EXPECT_NE(msg(base({{{"hubbardLinearResponse",json::object()}}},hub)).find("forceComplex"),std::string::npos) << "the real-TRIM response face is not built";
    json cplx=hub; cplx["forceComplex"]=true;
    EXPECT_EQ(msg(base({{{"hubbardLinearResponse",json::object()}}},cplx)),"NO THROW");
    json imposed=cplx; imposed["imposeSymmetry"]=true;
    EXPECT_EQ(msg(base({{{"hubbardLinearResponse",json::object()}}},imposed)),"NO THROW")
        << "imposition does not always reduce the mesh (NiO AFM-II 2x2x2: 8 -> 8): the response decides at run time, Vet() notes it";
    EXPECT_NE(msg(base({{{"independentResponse",{{"nq",2}}}}},hub,json::array({3,3,3}))).find("incommensurate"),std::string::npos);
    EXPECT_EQ(msg(base({{{"independentResponse",{{"nq",2}}}}},hub,json::array({4,4,2}))),"NO THROW");
    EXPECT_NE(msg(base({{{"hubbardLinearResponse",{{"nq",2}}}}},cplx,json::array({2,1,1}))).find("incommensurate"),std::string::npos) << "the q loop shares the pre-flight";
    EXPECT_EQ(msg(base({{{"hubbardLinearResponse",{{"nq",{2,1,1}}}}}},cplx,json::array({2,1,1}))),"NO THROW");
    {
        json lg=base({{{"hubbardLinearResponse",{{"littleGroupK",true}}}}},cplx);
        EXPECT_NE(msg(lg).find("minDeltaRhoPerElectron <= 1e-9"),std::string::npos) << "littleGroupK on a reference converged only to the default gate";
        lg["scf"]={{"minDeltaRhoPerElectron",1e-9}};
        EXPECT_EQ(msg(lg),"NO THROW");
    }
    EXPECT_NE(msg(base({{{"hubbardLinearResponse",{{"nq",{2,1}}}}}},cplx)).find("must be n or [n1,n2,n3]"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardFiniteDifference",{{"perturb",{5}},{"alpha_eV",0.1}}}}},hub)).find("past the 1 manifolds"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardFiniteDifference",json::object()}}},hub)).find("alpha_eV must be nonzero"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardLoop",{{"maxOuter",0}}}}},hub)).find("maxOuter >= 1"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardLoop",json::object()}}},json::object())).find("postSCF.0 (hubbardLoop)"),std::string::npos) << "the message names the entry and the action";
}

// The grid audit (doc/GridEpsilonPlan.md §3) is the one postSCF action that is not about +U: it needs no Hubbard manifold, only a raster.
// C2 step 3: its reference is the run's own rasters plus an alias-exponent MARGIN (the named BallOnly/AliasFree rasters are gone).
TEST(Deck, GridAuditRoundTripsItsReferenceMargin)
{
    auto msg=[](json d){ try { deck::RunSpec s; deck::FromJson(d,s); deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };
    auto base=[](json post, json solid) { return json{{"structure","NaF_rocksalt"},{"solid",solid},{"postSCF",post}}; };
    deck::RunSpec r;
    deck::FromJson(json{{"structure","NaF_rocksalt"},{"postSCF",{{{"gridAudit",json::object()}}, {{"gridAudit",{{"referenceMargin",12.0}}}}}}}, r);
    ASSERT_EQ(r.postSCF.size(),2u);
    EXPECT_EQ(r.postSCF[0].referenceMargin, 8.0) << "the default margin: the reference's own tail ~1e-3 of the run's";
    EXPECT_EQ(r.postSCF[1].referenceMargin, 12.0);
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_EQ(msg(base({{{"gridAudit",json::object()}}},json::object())),"NO THROW") << "no Hubbard manifold needed";
    EXPECT_NE(msg(base({{{"gridAudit",{{"referenceMargin",0.0}}}}},json::object())).find("referenceMargin"),std::string::npos) << "the reference must be finer";
    EXPECT_NE(msg(base({{{"gridAudit",{{"reference","AliasFree"}}}}},json::object())).find("retired"),std::string::npos) << "a named raster is refused, saying why";
    EXPECT_NE(msg(base({{{"gridAudit",json::object()}}},json{{"raster","AliasFree"}})).find("retired"),std::string::npos) << "so is solid.raster";
}

// The ONE grid accuracy target (doc/GridEpsilonPlan.md §1), through its deprecated solid.accuracy alias, and the expert calibration round-trip; eps defaults to
// 1e-5 with NO calibration stated (C2 step 3: the run's pilot audit self-calibrates), eps = 0 states the rasters, a negative is refused.
TEST(Deck, GridAccuracyTargetRoundTripsAndRefusesNonPositive)
{
    auto read=[](json solid){ deck::RunSpec s; deck::FromJson(json{{"structure","Si_diamond"},{"solid",solid}},s); return s; };
    auto msg=[&](json solid){ try { read(solid); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };
    const deck::RunSpec d=read(json::object());
    EXPECT_EQ(d.solid.accuracy.eps, 1e-5);
    EXPECT_FALSE(d.solid.accuracy.eta.has_value()) << "the default self-calibrates";
    EXPECT_FALSE(deck::ToJson(d).at("solid").contains("accuracyCalibration")) << "an unstated calibration is not written as if stated";
    const deck::RunSpec r=read({{"accuracy",1e-6},{"accuracyCalibration",{{"vlocLong",0.1},{"density",0.01}}}});
    EXPECT_EQ(r.solid.accuracy.eps, 1e-6);
    ASSERT_TRUE(r.solid.accuracy.eta.has_value());
    EXPECT_EQ(r.solid.accuracy.eta->vlocLong, 0.1);
    EXPECT_EQ(r.solid.accuracy.eta->density, 0.01);
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j);
    EXPECT_TRUE(back.solid.accuracy==r.solid.accuracy);
    const deck::RunSpec z=read({{"accuracy",0.0},{"tolerances",{{"rasterX",16.0}}}});
    EXPECT_EQ(z.solid.accuracy.eps, 0.0) << "0 = no target: the rasters as stated";
    EXPECT_EQ(z.solid.tolerances.rasterX, 16.0);
    EXPECT_NE(msg({{"accuracy",-1.0}}).find("solid.accuracy"),std::string::npos);
    EXPECT_NE(msg({{"accuracyCalibration",{{"density",-1.0}}}}).find("accuracyCalibration"),std::string::npos);
    EXPECT_NE(msg({{"accuracyCalibration",{{"vlocLong",2.0}}}}).find("delivered/bound"),std::string::npos) << "eta > 1 would claim the bound is optimistic";
    EXPECT_NE(msg({{"accuracyCalibration",{{"xc",1.0}}}}).find("xc"),std::string::npos) << "the XC raster has no rigorous bound to calibrate";
}

// THE ACCURACY TARGET IS THE RUN'S (V-FACADE-ACCURACY, user ruling 2026-11-03): the top-level `accuracy` reaches the facade the deck
// builds; `solid.accuracy` is its deprecated alias (read, and warned about); stating both is refused; the record writes the top level.
TEST(Deck, TheAccuracyTargetIsTopLevelAndSolidAccuracyIsAWarnedAlias)
{
    auto read=[](json j){ j["structure"]="Si_diamond"; deck::RunSpec s; deck::FromJson(j,s); return s; };
    const deck::RunSpec top=read({{"accuracy",3e-6}});
    EXPECT_EQ(top.solid.accuracy.eps, 3e-6) << "the top-level key reaches the solid facade's options";
    const json rec=deck::ToJson(top);
    EXPECT_EQ(rec.at("accuracy").get<double>(), 3e-6) << "the record states it at the top level";
    EXPECT_FALSE(rec.at("solid").contains("accuracy")) << "and never as solid.accuracy";
    deck::RunSpec back; deck::FromJson(rec,back);
    EXPECT_EQ(deck::ToJson(back), rec) << "the record round-trips";

    testing::internal::CaptureStderr();
    const deck::RunSpec alias=read({{"solid",{{"accuracy",2e-6}}}});
    const std::string err=testing::internal::GetCapturedStderr();
    EXPECT_EQ(alias.solid.accuracy.eps, 2e-6) << "the deprecated alias still works";
    EXPECT_NE(err.find("'solid.accuracy' is deprecated"), std::string::npos) << err;
    EXPECT_FALSE(deck::ToJson(alias).at("solid").contains("accuracy")) << "a re-written record migrates to the top level";

    bool threw=false;
    try { (void)read({{"accuracy",1e-6},{"solid",{{"accuracy",1e-6}}}}); }
    catch (const std::runtime_error& e) { threw=true; EXPECT_NE(std::string(e.what()).find("both stated"), std::string::npos) << e.what(); }
    EXPECT_TRUE(threw) << "both stated is refused, even when they agree";
    threw=false;
    try { (void)read({{"accuracy",-1.0}}); }
    catch (const std::runtime_error& e) { threw=true; EXPECT_NE(std::string(e.what()).find("'accuracy'"), std::string::npos) << e.what(); }
    EXPECT_TRUE(threw) << "a negative target is refused";
}

// V_loc long's routing tolerance = the rigorous per-pair map (eps) corrected by a MEASURED calibration (expert
// `solid.accuracyCalibration.vlocRouting`, eta = delivered/bound, default 10): a non-default calibration reaches GridAccuracy::VlocEps(),
// the setup writes it into the tolerances the evaluator is built from (marked derived), the evaluator reports the value and its source,
// and a STATED absolute `tolerances.vlocEps` wins over it.
TEST(Deck, AVlocRoutingCalibrationReachesTheEvaluatorAndAStatedVlocEpsWins)
{
    auto read=[](json solid){ deck::RunSpec s; deck::FromJson(json{{"structure","Si_diamond"},{"solid",solid}},s); return s; };
    const deck::RunSpec d=read(json::object());
    EXPECT_EQ(d.solid.accuracy.vlocRoutingEta, 10.0);
    EXPECT_EQ(d.solid.accuracy.VlocEps(), 1e-5/10.0) << "the default is bit-identical to the old eps/10";
    const deck::RunSpec f=read({{"accuracyCalibration",{{"vlocRouting",2.0}}}});
    EXPECT_EQ(f.solid.accuracy.vlocRoutingEta, 2.0);
    EXPECT_FALSE(f.solid.accuracy.eta.has_value()) << "vlocRouting alone states no raster eta: the rasters stay self-calibrated";
    EXPECT_DOUBLE_EQ(f.solid.accuracy.VlocEps(), 5e-6);
    { const json j=deck::ToJson(f); deck::RunSpec back; deck::FromJson(j,back); EXPECT_TRUE(back.solid.accuracy==f.solid.accuracy); }
    EXPECT_FALSE(deck::ToJson(d).at("solid").contains("accuracyCalibration")) << "the default calibration is not written";
    deck::RunSpec bad;
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"solid",{{"accuracyCalibration",{{"vlocRouting",0.0}}}}}},bad),std::runtime_error);

    const auto mat=Materials::Get("Si_diamond");
    Lattice_3D lat(*mat.cell, ivec3_t(1,1,1));
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, mat.cell.get()));
    SolidCalcOptions o; o.Nelec=8; o.species=mat.species; o.densityEcut=20.0; o.label="vlocRouting";
    o.accuracy.eps=1e-5; o.accuracy.eta=BasisSet::Gaussian::AccuracyCalibration{}; o.accuracy.vlocRoutingEta=2.0;   // stated raster eta: no pilot
    const auto s=SolidSetup::Resolve(lat, mol, o);
    EXPECT_DOUBLE_EQ(s->opts.tolerances.vlocEps, 5e-6) << "the evaluator is built from these tolerances";
    EXPECT_TRUE(s->opts.tolerances.VlocEpsDerived());
    {   // the same options through a deck run, which files its report: the evaluator's grids section names the vlocEps that ran
        const fs::path dir=TmpDir("vlocrouting");
        json dj={{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"}}},{"accuracy",1e-5},
                 {"solid",{{"densityEcut",20.0},{"accuracyCalibration",{{"vlocLong",1.0},{"density",1.0},{"vlocRouting",2.0}}}}},
                 {"scf",{{"NMaxIter",1},{"minDeltaFD","off"},{"minVirial","off"},{"minFD","off"}}}};
        deck::RunSpec spec; deck::FromJson(dj,spec);
        deck::Provenance pv; pv.codeVersion="t";
        qchem::BoundedSCF bounded("the grids report's vlocEps row is the subject: one iteration builds the evaluator");
        const auto out=deck::Run(spec,pv,dir);
        const qchem::report::json* grids=nullptr;
        for (auto it=qchem::report::GlobalReport().begin(); it!=qchem::report::GlobalReport().end(); ++it)
            if (it.key().rfind(out.revision.stem().string()+"@",0)==0 && it.value().contains("grids")) grids=&it.value()["grids"];
        ASSERT_NE(grids, nullptr) << "no filed run named " << out.revision.stem() << " with a grids section";
        ASSERT_TRUE(grids->contains("vlocEps"));
        EXPECT_DOUBLE_EQ((*grids)["vlocEps"]["value"].get<double>(), 5e-6);
        EXPECT_EQ((*grids)["vlocEps"]["source"].get<std::string>(), "derived");
        EXPECT_EQ((*grids)["vlocEps"]["rigorous"].get<double>(), 1e-5) << "the rigorous piece";
        EXPECT_EQ((*grids)["vlocEps"]["etaRouting"].get<double>(), 2.0) << "the calibration";
        fs::remove_all(dir);
    }

    o.tolerances.vlocEps=3e-7;   // STATED: wins over the fraction
    const auto t=SolidSetup::Resolve(lat, mol, o);
    EXPECT_EQ(t->opts.tolerances.vlocEps, 3e-7);
    EXPECT_FALSE(t->opts.tolerances.VlocEpsDerived());
}

// The egg-box action (doc/GridEpsilonPlan.md §5 verdict 5): round-trips its shift, defaults to a generic 0.037 bohr, and is refused
// before the SCF wherever a rigid translation cannot be compared function for function (an imposed group, +U) or measures nothing.
TEST(Deck, EggBoxRoundTripsAndRefusesWhatCannotTranslate)
{
    auto msg=[](json d){ try { deck::RunSpec s; deck::FromJson(d,s); deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };
    auto base=[](json post, json solid) { return json{{"structure","NaF_rocksalt"},{"solid",solid},{"postSCF",post}}; };
    deck::RunSpec r;
    deck::FromJson(json{{"structure","NaF_rocksalt"},{"postSCF",{{{"eggBox",json::object()}}, {{"eggBox",{{"shift",{0.1,0,0}}}}}}}}, r);
    ASSERT_EQ(r.postSCF.size(),2u);
    EXPECT_GT(norm(r.postSCF[0].shift), 0.0) << "a default shift that measures something";
    EXPECT_EQ(r.postSCF[1].shift.x, 0.1);
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_EQ(msg(base({{{"eggBox",json::object()}}},json::object())),"NO THROW");
    EXPECT_NE(msg(base({{{"eggBox",{{"shift",{0,0,0}}}}}},json::object())).find("zero shift"),std::string::npos);
    EXPECT_NE(msg(base({{{"eggBox",json::object()}}},json{{"imposeSymmetry",true}})).find("free twin"),std::string::npos);
    EXPECT_NE(msg(base({{{"eggBox",json::object()}}},json{{"hubbard",{{{"site",0},{"l",1},{"U_eV",0.0}}}}})).find("U = 0 twin"),std::string::npos);
    EXPECT_NE(msg(base({{{"eggBox",{{"shift",{0.1,0}}}}}},json::object())).find("[x,y,z]"),std::string::npos);
}

// The XC-mesh audit action (doc/OpenWork.md §3 "A curvature needs an XC-MESH ladder"): round-trips its reference mesh, defaults to lane B's
// (80, 48), and is refused before the SCF where it cannot work (a plane-wave run has no Becke mesh) or names no mesh.
TEST(Deck, XCMeshAuditRoundTripsAndRefusesANonPositiveMesh)
{
    auto msg=[](json d){ try { deck::RunSpec s; deck::FromJson(d,s); deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };
    auto base=[](json post) { return json{{"structure","NaF_rocksalt"},{"postSCF",post}}; };
    deck::RunSpec r;
    deck::FromJson(base({{{"xcMeshAudit",json::object()}}, {{"xcMeshAudit",{{"nRadial",60},{"angularDegree",35}}}}}), r);
    ASSERT_EQ(r.postSCF.size(),2u);
    EXPECT_EQ(r.postSCF[0].xcRefNRadial, 80); EXPECT_EQ(r.postSCF[0].xcRefAngularDegree, 48);
    EXPECT_EQ(r.postSCF[1].xcRefNRadial, 60); EXPECT_EQ(r.postSCF[1].xcRefAngularDegree, 35);
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_EQ(msg(base({{{"xcMeshAudit",json::object()}}})),"NO THROW");
    EXPECT_NE(msg(base({{{"xcMeshAudit",{{"nRadial",0}}}}})).find("> 0"),std::string::npos);
    EXPECT_NE(msg(base({{{"xcMeshAudit",{{"angularDeg",40}}}}})).find("angularDeg"),std::string::npos) << "a typo is refused, not ignored";
}

// The displace action (OpenWork F): its points round-trip, and a deck is refused before the SCF where ONE SolidSession could not run
// them -- the free run's mesh still the symmetric one, an atom past the structure, a zero delta, no points at all.
TEST(Deck, DisplaceRoundTripsAndRefusesWhatASessionCannotRun)
{
    auto msg=[](json d){ try { deck::RunSpec s; deck::FromJson(d,s); deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); } return std::string("NO THROW"); };
    const json plain={{"xcMeshSymmetric",false}};
    auto base=[](json pts, json solid) { return json{{"structure","NaF_rocksalt"},{"solid",solid},{"postSCF",{{{"displace",{{"points",pts}}}}}}}; };
    const json pts={{{"atoms",{1}},{"delta",{0.1,0,0}}}, {{"atoms",json::array()},{"delta",{0.1,0,0}}}};   // one ion, then the rigid twin
    deck::RunSpec r; deck::FromJson(base(pts,plain), r);
    ASSERT_EQ(r.postSCF.size(),1u); ASSERT_EQ(r.postSCF[0].points.size(),2u);
    EXPECT_EQ(r.postSCF[0].points[0].atoms, std::vector<size_t>{1});
    EXPECT_TRUE(r.postSCF[0].points[1].atoms.empty()) << "empty = every atom";
    EXPECT_EQ(r.postSCF[0].points[0].delta.x, 0.1);
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_EQ(msg(base(pts,plain)),"NO THROW");
    EXPECT_NE(msg(base(pts,json::object())).find("xcMeshSymmetric false"),std::string::npos) << "the free default mesh switches along a path";
    EXPECT_EQ(msg(base(pts,json{{"imposeSymmetry",true}})),"NO THROW") << "an imposed run's mesh is the group's: the session refuses a point that breaks it";
    EXPECT_NE(msg(base({{{"atoms",{2}},{"delta",{0.1,0,0}}}},plain)).find("past the structure"),std::string::npos);
    EXPECT_NE(msg(base({{{"atoms",{0}},{"delta",{0,0,0}}}},plain)).find("zero delta"),std::string::npos);
    EXPECT_NE(msg(base(json::array(),plain)).find("no points"),std::string::npos);
    EXPECT_NE(msg(base({{{"atoms",{0}},{"delta",{0.1,0}}}},plain)).find("[x,y,z]"),std::string::npos);
    EXPECT_NE(msg(base({{{"atom",0},{"delta",{0.1,0,0}}}},plain)).find("unknown key"),std::string::npos);
}

TEST(Deck, APolarizedPlusUDeckMustStateItsMagneticOrdering)
{
    auto msg=[](json solid)
    {
        deck::RunSpec s; deck::FromJson(json{{"structure","NiO_AFM2"},{"solid",solid}},s);
        try { deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); }
        return std::string("NO THROW");
    };
    const json H={{"site",0},{"l",2},{"U_eV",4.0}};
    EXPECT_NE(msg({{"multiplicity",1},{"hubbard",{H}}}).find("siteSpins"),std::string::npos) << "never derived from the seed any more";
    EXPECT_NE(msg({{"multiplicity",1},{"hubbard",{H}}}).find("Ni Ni O O"),std::string::npos) << "the message names the sites";
    EXPECT_EQ(msg({{"multiplicity",1},{"hubbard",{H}},{"siteSpins",{1,-1,0,0}}}),"NO THROW");
    EXPECT_EQ(msg({{"multiplicity",1},{"hubbard",{H}},{"greyImposition",true}}),"NO THROW") << "the erasure control states 'no decoration'";
    EXPECT_EQ(msg({{"hubbard",{H}}}),"NO THROW") << "unpolarized: nothing to order";
    EXPECT_EQ(msg({{"multiplicity",1}}),"NO THROW") << "no +U: the imposed-without-U path is not (yet) covered";
    EXPECT_NE(msg({{"multiplicity",1},{"hubbard",{H}},{"siteSpins",{1,-1}}}).find("2 entries for 4 sites"),std::string::npos);
    EXPECT_NE(msg({{"siteSpins",{2,-1,0,0}}}).find("-1, 0 or +1"),std::string::npos);
}

// multiplicity 1 is the two-channel SINGLET, nUp = nDown (OpenWork §3, lane H 2026-11-01): an odd electron count cannot hold it, and it
// once ran a silent doublet (bcc Li, 53 mHa off, CONVERGED).  Refused at resolve -- unless one Fermi level sets the channels.
TEST(Deck, MultiplicityOneOnAnOddCountIsRefusedUnlessTheSpinsShareAFermiLevel)
{
    auto msg=[](json solid)
    {
        deck::RunSpec s; deck::FromJson(json{{"structure","NaF_rocksalt"},{"solid",solid}},s);   // 8 electrons (Na q1, F q7)
        try { deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); }
        return std::string("NO THROW");
    };
    EXPECT_EQ(msg({{"multiplicity",1}}),"NO THROW") << "an even count holds the singlet";
    EXPECT_NE(msg({{"multiplicity",1},{"charge",1}}).find("SINGLET"),std::string::npos) << "7 electrons: refused, not a silent doublet";
    EXPECT_NE(msg({{"multiplicity",1},{"charge",1}}).find("odd Nelec 7"),std::string::npos) << "the message names the count";
    EXPECT_EQ(msg({{"multiplicity",1},{"charge",1},{"spinsShareFermi",true}}),"NO THROW") << "one Fermi level: the counts float";
    EXPECT_EQ(msg({{"multiplicity",2},{"charge",1}}),"NO THROW") << "the doublet, stated";
    EXPECT_EQ(msg({{"charge",1}}),"NO THROW") << "minimal spin";
    EXPECT_NE(msg({{"multiplicity",2}}).find("parity"),std::string::npos) << "a doublet of 8 electrons";
    EXPECT_NE(msg({{"multiplicity",10}}).find("parity"),std::string::npos) << "2S = 9 > 8 electrons";
}

// The IonicSAD charge order is STATED or derived (OpenWork §3 "Spinel charge order is CHOSEN by atom index", 2026-10-30): a stated
// order is one charge per site, conserves the cell's charge, and seeds an IonicSAD run -- else the deck is refused at resolve.
TEST(Deck, AStatedChargeOrderIsOnePerSiteNeutralAndIonicSAD)
{
    auto msg=[](json solid)
    {
        deck::RunSpec s; deck::FromJson(json{{"structure","NiO_AFM2"},{"solid",solid}},s);
        try { deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); }
        return std::string("NO THROW");
    };
    EXPECT_EQ(msg({{"siteCharges",{2,2,-2,-2}}}),"NO THROW") << "IonicSAD is the default seed";
    EXPECT_NE(msg({{"siteCharges",{2,-2}}}).find("solid.siteCharges -- siteCharges has 2 entries for 4 sites (Ni Ni O O)"),std::string::npos);
    EXPECT_NE(msg({{"siteCharges",{3,2,-2,-2}}}).find("sum to 1, not 0"),std::string::npos);
    EXPECT_NE(msg({{"siteCharges",{2,2,-2,-2}},{"seed","SAD"}}).find("not IonicSAD"),std::string::npos) << "nothing else reads it";
}

//! A grid field set beside cellKind Auto used to be DROPPED (Auto resolves to a calibrated recipe): a β-MnO₂ mesh-convergence run came
//! back bit-identical.  Resolve refuses it, naming the field, so --check catches it before a revision is claimed.
TEST(Deck, AnXCMeshRecipeFieldUnderAutoIsRefusedAtResolve)
{
    auto msg=[](json xcMesh)
    {
        deck::RunSpec s; deck::FromJson(json{{"structure","Si_diamond"},{"solid",{{"xcMesh",xcMesh}}}},s);
        try { deck::Resolve(s); } catch (const std::runtime_error& e) { return std::string(e.what()); }
        return std::string("NO THROW");
    };
    EXPECT_NE(msg({{"nRadial",60}}).find("nRadial"),std::string::npos);
    EXPECT_EQ(msg({{"cellKind","Becke"},{"nRadial",60}}),"NO THROW") << "a stated kind takes its recipe";
    EXPECT_EQ(msg(json::object()),"NO THROW") << "Auto with nothing set";
}

TEST(Deck, PostSCFRunsOnTheConvergedCalculationAndTheRecordKeepsTheResults)
{
    const fs::path dir=TmpDir("post");
    json d={{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"}}},
        {"solid",{{"densityEcut",20.0},{"hubbard",{{{"site",0},{"l",1},{"U_eV",0.0}},{{"site",1},{"l",1},{"U_eV",0.0}}}}}},
        {"scf",{{"NMaxIter",60},{"minDeltaRhoPerElectron",1e-3},{"minDeltaE",1e-6},{"minDeltaFD","off"},{"minVirial","off"},{"minFD","off"},{"startingRelaxRo",0.3}}},
        {"postSCF",{{{"estimateHubbardU",json::object()}}}}};
    deck::RunSpec spec; deck::FromJson(d,spec);
    deck::Provenance pv; pv.codeVersion="t";
    const auto out=deck::Run(spec,pv,dir);
    ASSERT_TRUE(out.converged) << out.summary;
    ASSERT_EQ(out.postSCF.size(),1u);
    EXPECT_EQ(out.postSCF[0].action,"estimateHubbardU"); EXPECT_TRUE(out.postSCF[0].ok) << out.postSCF[0].summary;
    const json rec=json::parse(std::ifstream(out.revision));
    ASSERT_TRUE(rec.contains("results"));
    EXPECT_TRUE(rec.at("results").at("converged"));
    EXPECT_NEAR(rec.at("results").at("energy").get<double>(), *out.energy, 1e-12);
    const json est=rec.at("results").at("postSCF")[0].at("estimates");
    ASSERT_EQ(est.size(),2u) << "one estimate per manifold";
    EXPECT_TRUE(est[0].contains("Ueff_eV"));
    EXPECT_TRUE(rec.contains("run")) << "the deck part of the record is intact beside the results";
    fs::remove_all(dir);
}

TEST(Deck, KeysStartingWithAnUnderscoreAreCommentsAndShippedDecksLoad)
{
    deck::RunSpec r;
    deck::FromJson(json{{"_doc","a note"},{"structure","Si_diamond"},{"solid",{{"_why","x"},{"Nelec",8}}}}, r);
    EXPECT_EQ(r.solid.Nelec,8);
    EXPECT_FALSE(deck::ToJson(r).contains("_doc")) << "a record does not keep comments";
    // every deck shipped in decks/ parses and resolves (so a stale key in a recipe is caught by ctest, not by the next campaign)
    size_t n=0;
    for (const auto& e : fs::recursive_directory_iterator(fs::path(DECKS_DIR)))
        if (e.path().extension()==".json")
        {
            deck::RunSpec d; ASSERT_NO_THROW(deck::FromJson(deck::LoadDeck(e.path(),"x"),d)) << e.path();
            ASSERT_NO_THROW(deck::Resolve(d)) << e.path();
            ++n;
        }
    EXPECT_GE(n,6u) << "decks/ and decks/bench/";
}

// Every shipped deck VETS: it resolves (postSCF pre-flight included) and its lattice and basis build -- no SCF.  The campaign decks
// (decks/campaigns, reconstructed from the old logs) are the reason: a stale key, a basis file with no block for an element, or a
// postSCF the run could not honour is caught by ctest, not by the first hour of the re-run.
TEST(Deck, EveryShippedDeckVets)
{
    size_t n=0;
    for (const auto& e : fs::recursive_directory_iterator(fs::path(DECKS_DIR)))
        if (e.path().extension()==".json")
        {
            json j=deck::LoadDeck(e.path(),"x");
            // a panel deck vendors its basis beside itself (decks/panel/, resolved by IntegrationTests/GPW/Panel.C the same way)
            if (j.contains("basis") && j["basis"].contains("file"))
            {
                const fs::path f=j["basis"]["file"].get<std::string>(), beside=e.path().parent_path()/f;
                if (f.is_relative() && fs::exists(beside)) j["basis"]["file"]=beside.string();
            }
            deck::RunSpec d; ASSERT_NO_THROW(deck::FromJson(j,d)) << e.path();
            deck::VetReport v;
            ASSERT_NO_THROW(v=deck::Vet(d, fs::temp_directory_path())) << e.path();
            EXPECT_FALSE(v.summary.empty()) << e.path();
            ++n;
        }
    EXPECT_GE(n,1u);
}

// ---- an INLINE structure (2026-10-31, the GUI's door): the deck states the cell instead of naming it ----
//! THE FRIEND (src/forward.H): the restart fingerprint a built calculation would save -- what "the same run" means to a restart.
class DeckTests
{
public:
    static StateFingerprint Fingerprint(const SolidCalculation& c) { return c.Fingerprint(0.0); }
};

namespace {
// Si_diamond stated inline in the shape lane G's translator writes: explicit lattice vectors (CubicF a = 10.26: the a/2 columns),
// no name, no pp (solid.pp's default, gth-lda-valence, is the material's set).  Copied, not read: the test states what it means.
json SiInline()
{
    return {{"lattice",{{"vectors_bohr",{{0.0,5.13,5.13},{5.13,0.0,5.13},{5.13,5.13,0.0}}}}},
            {"atoms",json::array({{{"el","Si"},{"frac",{0.0,0.0,0.0}}}, {{"el","Si"},{"frac",{0.25,0.25,0.25}}}})}};
}
const json kSiRecipe={{"basis",{{"data","SIPP_SR"}}},{"solid",{{"densityEcut",20.0},{"imposeSymmetry",true}}},
                      {"scf",{{"NMaxIter",60},{"minDeltaRhoPerElectron",1e-3},{"minDeltaE",1e-6},{"minDeltaFD","off"},{"minVirial","off"},{"minFD","off"},{"startingRelaxRo",0.3},{"mergeTol",1e-4}}}};
deck::RunSpec SiDeck(const json& structure)
{
    json d=kSiRecipe; d["structure"]=structure;
    deck::RunSpec r; deck::FromJson(d,r); return r;
}
}

TEST(Deck, AnInlineStructureIsTheStructureKeyAsAnObjectAndRoundTrips)
{
    deck::RunSpec r;
    EXPECT_THROW(deck::FromJson(json{{"kmesh",{1,1,1}}},r),std::runtime_error) << "no structure";
    EXPECT_THROW(deck::FromJson(json{{"structure",42}},r),std::runtime_error) << "a name or an object";
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"inlineStructure",SiInline()}},r),std::runtime_error) << "no second key";
    deck::RunSpec both=SiDeck("Si_diamond"); both.inlineStructure=SiInline();
    EXPECT_THROW((void)deck::Resolve(both),std::runtime_error) << "a name AND an inline cell, built in C++";

    deck::RunSpec in=SiDeck(SiInline());
    EXPECT_TRUE(in.structure.empty());
    EXPECT_EQ(in.StructureName(),"inline") << "an inline cell with no name";
    json named=SiInline(); named["name"]="Si_mine";
    EXPECT_EQ(SiDeck(named).StructureName(),"Si_mine");
    const json j=deck::ToJson(in);
    EXPECT_TRUE(j.at("structure").is_object()) << "the record states the cell, not a name it does not have";
    deck::RunSpec back; deck::FromJson(j,back); EXPECT_EQ(deck::ToJson(back),j);

    // resolved like a named one: the species and Nelec are derived from the set (solid.pp here; the entry's own `pp` if it has one)
    deck::RunSpec ref=SiDeck("Si_diamond");
    const auto mi=deck::Resolve(in), mn=deck::Resolve(ref);
    EXPECT_EQ(in.solid.Nelec,ref.solid.Nelec); EXPECT_EQ(in.solid.species,ref.solid.species); EXPECT_EQ(in.solid.pp,ref.solid.pp);
    EXPECT_EQ(mi.cell->GetNumAtoms(),mn.cell->GetNumAtoms());

    // strict: a typo in the entry would build a DIFFERENT cell silently
    auto refused=[](json entry, const std::string& what)
    {
        deck::RunSpec x=SiDeck(entry);
        try { (void)deck::Resolve(x); ADD_FAILURE() << "no throw for " << what; }
        catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find(what),std::string::npos) << e.what(); }
    };
    json both2=SiInline(); both2["lattice"]["type"]="CubicF";      refused(both2,"unknown key 'type'");
    json bravais=SiInline(); bravais["lattice"]={{"type","CubicF"},{"alfa",60.0}}; refused(bravais,"unknown key 'alfa'");
    json short3=SiInline(); short3["lattice"]["vectors_bohr"]={{1,0,0},{0,1,0}}; refused(short3,"three 3-vectors");
    json extra=SiInline(); extra["kmesh"]={2,2,2};                  refused(extra,"unknown key 'kmesh'");
    json badAtom=SiInline(); badAtom["atoms"][0]["fracs"]={0,0,0};  refused(badAtom,"unknown key 'fracs'");
    json badName=SiInline(); badName["name"]="a/b";                 refused(badName,"file stem");
}

// The claim lane G binds to: an inline deck -- in lane G's shape, explicit lattice vectors -- and its named twin are the SAME run:
// the fingerprint a saved state would carry, built by deck::Build (which stops before the SCF), agrees with no difference at all,
// so either one's state resumes the other EXACTLY.
TEST(Deck, AnInlineDeckAndItsNamedTwinBuildTheSameFingerprint)
{
    const auto a=deck::Build(SiDeck("Si_diamond"));
    const auto b=deck::Build(SiDeck(SiInline()));
    auto cmp=CompareFingerprints(DeckTests::Fingerprint(*a), DeckTests::Fingerprint(*b));
    ASSERT_TRUE(cmp) << cmp.Error().details;
    EXPECT_TRUE(cmp->empty()) << "a warm difference: " << (cmp->empty() ? "" : cmp->front());
    // and a moved atom IS a different cell (the comparison is not vacuous)
    json moved=SiInline(); moved["atoms"][1]["frac"]={0.25,0.25,0.26};
    const auto c=deck::Build(SiDeck(moved));
    auto cmp2=CompareFingerprints(DeckTests::Fingerprint(*a), DeckTests::Fingerprint(*c));
    EXPECT_TRUE(!cmp2 || !cmp2->empty()) << "a moved atom left the fingerprint unchanged";
}

// deck::Build STOPS BEFORE THE SCF, and Run is Build + Converge: the observer has seen nothing until the caller converges, and the
// converged energy is Run's own, bit for bit (ONE path).
TEST(Deck, BuildStopsBeforeTheSCFAndRunIsBuildPlusConverge)
{
    deck::RunSpec spec=SiDeck("Si_diamond");
    size_t seen=0;
    const auto calc=deck::Build(spec, [&](const auto&){ ++seen; });
    EXPECT_EQ(seen,0u) << "Build ran an SCF";
    EXPECT_FALSE(calc->DidConverge());
    auto R=calc->Converge(deck::ScheduleToRun(spec));
    ASSERT_TRUE(R);
    EXPECT_GT(seen,0u) << "the observer was not attached before the first iteration";

    const fs::path dir=TmpDir("build");
    deck::Provenance pv; pv.commandLine="test"; pv.codeVersion="t";
    const deck::RunOutcome out=deck::Run(spec,pv,dir);
    ASSERT_TRUE(out.converged) << out.summary;
    EXPECT_DOUBLE_EQ(R->Energy(), *out.energy);
    fs::remove_all(dir);

    deck::RunSpec stem=spec; stem.state.restartFrom="Si_diamond.r001";
    EXPECT_THROW((void)deck::Build(stem),std::runtime_error) << "a revision stem needs Run's outDir";
    deck::RunSpec missing=spec; missing.state.restartFrom="/nonexistent/x.h5";
    EXPECT_THROW((void)deck::Build(missing),std::runtime_error);
    EXPECT_EQ(deck::ScheduleToRun(stem).size(),1u) << "a restart runs the final stage only";
}

//! THE REPORT IS WRITE-ONLY (V-REPORT-WRITEONLY, 2026-11-02): library code never reads the run report back to decide.  On Si k222 (eight
//! k-blocks sharing ONE raster ladder, V-BS-RASTER) the run files each raster verdict ONCE and renders it once -- the ladder reports
//! itself once, no dedup needed -- and the grid audit's recorded verdicts are what the rasters' owner (the ladder) answers, which are the
//! rows the run filed.
TEST(Deck, AK222RunFilesEachRasterVerdictOnceAndTheAuditAsksTheBasisForIt)
{
    const fs::path dir=TmpDir("writeonly");
    json d={{"structure","Si_diamond"},{"kmesh",{2,2,2}},{"basis",{{"data","SIPP_SR"}}},{"solid",{{"densityEcut",20.0}}},
        {"scf",{{"NMaxIter",60},{"minDeltaRhoPerElectron",1e-6},{"minDeltaE",1e-8},{"minDeltaFD","off"},{"minVirial","off"},{"minFD","off"}}},
        {"postSCF",{{{"gridAudit",json::object()}}}}};
    deck::RunSpec spec; deck::FromJson(d,spec);
    deck::Provenance pv; pv.codeVersion="t";
    testing::internal::CaptureStdout();
    const auto out=deck::Run(spec,pv,dir);
    const std::string console=testing::internal::GetCapturedStdout();
    ASSERT_TRUE(out.converged) << out.summary;

    auto count=[&](const std::string& what) { size_t n=0; for (size_t i=console.find(what); i!=std::string::npos; i=console.find(what,i+1)) n++; return n; };
    // The point's ONE ladder files both rows inside its grids section (V-BS-RASTER), so each verdict's term renders once.
    EXPECT_EQ(count("V_loc long (integrate-back)"), 1u) << "the V_long verdict renders once, not once per k-block";
    EXPECT_EQ(count("density collocation (Hartree + XC raster)"), 1u);

    const qchem::report::json* run=nullptr;
    for (auto it=qchem::report::GlobalReport().begin(); it!=qchem::report::GlobalReport().end(); ++it)
        if (it.key().rfind(out.revision.stem().string()+"@",0)==0) run=&it.value();
    ASSERT_NE(run, nullptr) << "no filed run named " << out.revision.stem();
    ASSERT_TRUE(run->contains("grids"));
    std::ifstream in(out.revision); const json rec=json::parse(in);
    const auto& verdicts=rec["results"]["postSCF"][0]["verdicts"];
    for (const char* k : {"auditDensity","auditVlocLong"})
    {
        ASSERT_TRUE((*run)["grids"].contains(k)) << k;
        ASSERT_TRUE(verdicts.contains(k)) << k << " missing from the audit's verdicts: " << verdicts.dump();
        EXPECT_EQ(verdicts[k]["status"], (*run)["grids"][k]["status"].get<std::string>()) << k;
        EXPECT_EQ(verdicts[k]["alias"].get<double>(), (*run)["grids"][k]["alias"].get<double>()) << k << ": the basis's answer is the filed row";
    }
    fs::remove_all(dir);
}
