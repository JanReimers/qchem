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
        {"scf",{{"NMaxIter",60},{"minDeltaRho",1e-3},{"minDeltaE",1e-6},{"minDeltaFD",1e30},{"minVirial",1e30},{"minFD",1e30},{"startingRelaxRo",0.3},{"mergeTol",1e-4}}}};
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
    SolidCalculation calc(lat,mol,o,spec.scf);
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

// ---- D-ENV 6d.1: eV in the file / a.u. in RAM, basis trim+vet, saved states with lineage ----
TEST(Deck, UIsEvInTheFileAndAtomicUnitsInRamAndTheRecordReproducesTheBits)
{
    SolidCalcOptions o;
    deck::FromJson(json{{"hubbard",{{{"site",0},{"l",2},{"U_eV",4.0},{"Uirrep_eV",{3.0,4.5}},{"alpha_eV",0.1}}}}}, o);
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

    deck::RunSpec both; deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"vet",true},{"trim",{{{"Z",25},{"l",0},{"alpha",0.06}}}}}}}, both);
    EXPECT_THROW(deck::Resolve(both), std::runtime_error) << "vet and a stated trim are exclusive";
    deck::RunSpec sph; deck::FromJson(json{{"structure","MnO_AFM2"},{"basis",{{"data","VALENCE_LOWQ_SR"},{"spherical",true}}}}, sph);
    try { deck::Resolve(sph); FAIL() << "spherical + the SR transition-metal block must be refused"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("VALENCE_LOWQ_VA"),std::string::npos) << e.what(); }
    deck::RunSpec ok; deck::FromJson(json{{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"},{"spherical",true}}}}, ok);
    EXPECT_NO_THROW(deck::Resolve(ok)) << "the refusal is specific to a transition metal on the SR block";
    EXPECT_THROW(deck::FromJson(json{{"structure","Si_diamond"},{"state",{{"restore","x"}}}}, ok), std::runtime_error);
}

TEST(Deck, ASavedStateRestartsAndTheRevisionNamesItsParent)
{
    const fs::path dir=TmpDir("restart");
    const json base={{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"}}},{"solid",{{"densityEcut",20.0},{"imposeSymmetry",true}}},
        {"scf",{{"NMaxIter",60},{"minDeltaRho",1e-6},{"minDeltaE",1e-10},{"minDeltaFD",1e30},{"minVirial",1e30},{"minFD",1e30},{"startingRelaxRo",0.3}}}};
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
        {{"hubbardLinearResponse",{{"perturb",{0,1}},{"maxIter",300},{"restart",30},{"tol",1e-9}}}},
        {{"hubbardFiniteDifference",{{"perturb",{0}},{"alpha_eV",0.1}}}}}}}, r);
    ASSERT_EQ(r.postSCF.size(),5u);
    EXPECT_EQ(r.postSCF[1].maxOuter,8u); EXPECT_EQ(r.postSCF[2].nq,2); EXPECT_EQ(r.postSCF[3].perturb, (std::vector<size_t>{0,1}));
    EXPECT_DOUBLE_EQ(r.postSCF[4].alpha, 0.1/27.211386245988) << "alpha is eV in the file, Hartree in RAM";
    const json j=deck::ToJson(r); deck::RunSpec back; deck::FromJson(j,back);
    EXPECT_EQ(deck::ToJson(back),j); EXPECT_TRUE(back.postSCF==r.postSCF);
    EXPECT_FALSE(deck::ToJson(deck::RunSpec{}).contains("postSCF"));

    deck::RunSpec bad;
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLop",json::object()}}}}},bad), std::runtime_error) << "an unknown action name";
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLoop",{{"maxOuters",8}}}}}}},bad), std::runtime_error) << "an unknown parameter";
    EXPECT_THROW(deck::FromJson(json{{"structure","x"},{"postSCF",{{{"hubbardLoop",json::object()},{"estimateHubbardU",json::object()}}}}},bad), std::runtime_error) << "two actions in one entry";
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
    EXPECT_NE(msg(base({{{"hubbardLinearResponse",json::object()}}},imposed)).find("FULL k-mesh"),std::string::npos);
    EXPECT_NE(msg(base({{{"independentResponse",{{"nq",2}}}}},hub,json::array({3,3,3}))).find("incommensurate"),std::string::npos);
    EXPECT_EQ(msg(base({{{"independentResponse",{{"nq",2}}}}},hub,json::array({4,4,2}))),"NO THROW");
    EXPECT_NE(msg(base({{{"hubbardFiniteDifference",{{"perturb",{5}},{"alpha_eV",0.1}}}}},hub)).find("past the 1 manifolds"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardFiniteDifference",json::object()}}},hub)).find("alpha_eV must be nonzero"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardLoop",{{"maxOuter",0}}}}},hub)).find("maxOuter >= 1"),std::string::npos);
    EXPECT_NE(msg(base({{{"hubbardLoop",json::object()}}},json::object())).find("postSCF.0 (hubbardLoop)"),std::string::npos) << "the message names the entry and the action";
}

TEST(Deck, PostSCFRunsOnTheConvergedCalculationAndTheRecordKeepsTheResults)
{
    const fs::path dir=TmpDir("post");
    json d={{"structure","Si_diamond"},{"basis",{{"data","SIPP_SR"}}},
        {"solid",{{"densityEcut",20.0},{"hubbard",{{{"site",0},{"l",1},{"U_eV",0.0}},{{"site",1},{"l",1},{"U_eV",0.0}}}}}},
        {"scf",{{"NMaxIter",60},{"minDeltaRho",1e-3},{"minDeltaE",1e-6},{"minDeltaFD",1e30},{"minVirial",1e30},{"minFD",1e30},{"startingRelaxRo",0.3}}},
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
