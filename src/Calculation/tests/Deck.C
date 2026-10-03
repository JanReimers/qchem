// Unit tests of the input deck's data layer (D-ENV step 6): typed options <-> JSON, --set, numbered revisions.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <set>
#include <nlohmann/json.hpp>
import qchem.Deck;
import qchem.SolidCalculation;
import qchem.SCFParams;
import qchem.Mesh;
import qchem.Materials;
import qchem.ChargeDensity.Seed;
import qchem.LASolver;

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

TEST(Deck, NonDefaultValuesSurviveTheRoundTrip)
{
    SolidCalcOptions o;
    o.Nelec=48; o.multiplicity=1; o.species={{"Mn",15},{"O",6}}; o.kShift=rvec3_t(0.25,0,0.5);
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::CholeskyPivoted; o.imposeSymmetry=true; o.siteSpins={1,-1,0,0};
    o.tolerances.screenEps=1e-8; o.tolerances.mgridEcuts={53.3,17.8}; o.xcMesh.beckeEps=1e-7; o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;
    o.hubbard.push_back(HubbardU_OrthoAtomic(0,2,4.0));
    const json j=deck::ToJson(o);
    SolidCalcOptions back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back), j);
    EXPECT_EQ(back.species.size(),2u); EXPECT_EQ(back.species[0].first,"Mn"); EXPECT_EQ(back.tolerances.screenEps,1e-8);
    ASSERT_EQ(back.hubbard.size(),1u); EXPECT_TRUE(back.hubbard[0].orthoAtomic); EXPECT_NEAR(back.hubbard[0].U,4.0/27.211386245988,1e-15);
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
