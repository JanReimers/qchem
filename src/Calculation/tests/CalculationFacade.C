// File: Calculation/tests/CalculationFacade.C  What the molecular facade hands a front end (lane G, 2026-10-30): the
// constructor's own SCF observed from its first iteration, every MO as a level, and the spin density -- a zero field when
// the run is unpolarized, never an absence.
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <fstream>
import qchem.MoleculeCalculation;
import qchem.AtomCalculation;   // AtomCalculation::Prepare (the atom's Build -> Run)
import qchem.Structure;         // Atom (the atom's report ID)
import qchem.StructureData;
import qchem.SCFParams;
import qchem.SCFIterator;
import qchem.Types;
import qchem.SCFCap;   // BoundedSCF: the failure arm runs a deliberately bounded SCF
import qchem.BasisSet.Gaussian.Point.BasisFiles;   // BasisFile / BasisFileChecksum (a stated basis file)
import qchem.Reporting;   // the run report the molecular mesh files itself into

using namespace qchem;
using StructureData::GetMolecule;

TEST(CalculationFacade, TheConstructorsObserverSeesItsSCFFromTheFirstIteration)
{
    std::vector<SCFIterator::IterationRow> steps;
    MoleculeCalculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}), {},
                     [&steps](const SCFIterator::IterationRow& p){ steps.push_back(p); });
    ASSERT_TRUE(calc.IsConverged());
    ASSERT_FALSE(steps.empty()) << "the observer must see the constructor's SCF, not only a later Converge";
    EXPECT_EQ(steps.front().iteration, 1u) << "iteration 1: the constructor's SCF from its start";
    EXPECT_EQ(steps.size(), calc.IterationCount()) << "every iteration, none skipped";
    EXPECT_NEAR(steps.back().ValueOf("Etotal"), calc.Energy(), 1e-10);
}

//! A molecular DFT run's report holds its XC mesh BEFORE the SCF iterates, at hamiltonian.xcQuadrature -- the path a solid files:
//! the facade builds the Hamiltonian in its "hamiltonian" context, the XC quadrature opens "xcQuadrature", the mesh builder writes
//! one level in (OpenWork §3 "A molecule's mesh never reports itself", user rulings 2026-11-02, 2026-11-04).
TEST(CalculationFacade, AMoleculesXCMeshIsInTheRunReportByIterationOne)
{
    bool seen=false, checked=false;
    qchem::report::Run run("H2O");   // the test is the driver: it opens the run the facade writes into (user ruling D, 2026-11-02)
    MoleculeCalculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA}), {},
                     [&](const SCFIterator::IterationRow& p)
                     {
                         if (p.iteration!=1) return;
                         checked=true;
                         const auto& doc=qchem::report::CurrentRunReport();
                         if (!doc.contains("hamiltonian") || !doc["hamiltonian"].contains("xcQuadrature")) return;
                         const auto& xc=doc["hamiltonian"]["xcQuadrature"];
                         seen = xc.contains("points") && xc["points"].get<long>()>0 && xc["natom"].get<long>()==3;
                     });
    ASSERT_TRUE(calc.IsConverged());
    ASSERT_TRUE(checked) << "the observer saw no iteration 1";
    EXPECT_TRUE(seen) << "at iteration 1 the run report must already hold hamiltonian.xcQuadrature (the molecule's XC mesh)";
}

//! A PSEUDOPOTENTIAL run's PP quadrature is built in the Hamiltonian's own build phase, so it files at hamiltonian.ppQuadrature -- not
//! under basis.*, where the terms' lazy first-matrix build inside the SCF iterator's construction used to put it (CLAUDE.md "the
//! run report": fix the run sequence, never the path; OpenWork, the mesh-report remainder, 2026-11-02).
// A Gaussian basis assembles the PP terms in closed form, so a molecular PP run builds no PP mesh at all.
TEST(CalculationFacade, APseudopotentialOnAGaussianBasisBuildsNoPPQuadrature)
{
    bool built=true, stray=true;
    qchem::report::Run run("H2O");   // the driver's run (user ruling D, 2026-11-02)
    MoleculeCalculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"}), {},
                     [&](const SCFIterator::IterationRow& p)
                     {
                         if (p.iteration!=1) return;
                         const auto& doc=qchem::report::CurrentRunReport();
                         built = !doc.contains("hamiltonian") || doc["hamiltonian"].contains("ppQuadrature");
                         stray = doc.contains("basis") && (doc["basis"].contains("ppLocal") || doc["basis"].contains("ppNonLocal")
                                                           || doc["basis"].contains("ppQuadrature"));
                     });
    ASSERT_TRUE(calc.IsConverged());
    EXPECT_FALSE(built) << "a PP quadrature was built for a basis that assembles the PP terms in closed form";
    EXPECT_FALSE(stray) << "a PP mesh filed under basis";
}

TEST(CalculationFacade, LevelsListEveryMOAscendingAndTheOccupiedOnesAreTheOrbitals)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    ASSERT_TRUE(calc.IsConverged());
    const auto& levels=calc.GetLevels();
    ASSERT_GT(levels.size(), calc.NumOccupied()) << "the virtuals are listed too";
    double ne=0; size_t nocc=0;
    for (size_t i=0;i<levels.size();i++)
    {
        ASSERT_NE(levels[i].phi, nullptr);
        if (i>0) EXPECT_LE(levels[i-1].e, levels[i].e) << "ascending";
        ne+=levels[i].occupation;
        if (levels[i].occupation>0.5) nocc++;
        EXPECT_EQ(levels[i].spin, Spin::None) << "an unpolarized run";
    }
    EXPECT_NEAR(ne, 10.0, 1e-12);
    ASSERT_EQ(nocc, calc.NumOccupied());
    EXPECT_EQ(levels[nocc-1].phi, &calc.HOMO()) << "the highest occupied level IS the facade's HOMO";
    EXPECT_GT(levels[nocc].e, levels[nocc-1].e) << "a gap to the LUMO";
}

TEST(CalculationFacade, AnUnpolarizedSpinDensityIsTheZeroFieldAndAPolarizedOneIsNot)
{
    MoleculeCalculation closed(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    ASSERT_TRUE(closed.IsConverged());
    EXPECT_EQ(closed.GetSpinDensity()(rvec3_t(0.0, 0.0, 0.3)), 0.0);

    MoleculeCalculation triplet(MoleculeSetup::Resolve(GetMolecule("O2"), {.basis="dzvp", .multiplicity=3}));
    ASSERT_TRUE(triplet.IsConverged());
    // The two unpaired π* electrons put the moment off the bond axis, around the atoms.
    EXPECT_GT(triplet.GetSpinDensity()(rvec3_t(0.8, 0.0, 1.141)), 1e-3);
    size_t up=0, down=0; double nUp=0, nDown=0;
    for (const auto& L : triplet.GetLevels())
    {
        if (L.spin==Spin::Up)   { up++;   nUp+=L.occupation; }
        if (L.spin==Spin::Down) { down++; nDown+=L.occupation; }
    }
    EXPECT_EQ(up, down) << "both channels listed";
    EXPECT_NEAR(nUp-nDown, 2.0, 1e-12) << "2S = multiplicity - 1";
}

// THE OUTCOME DOOR (2026-10-31): Converge and Result hand back the answers only on a PROOF of convergence, and a failure
// carries the SCFFailure both facades share -- never a plausible energy as "the answer".
TEST(CalculationFacade, ConvergeIsAnOutcomeTheAnswersLiveOnTheProof)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    auto R=calc.Result();
    ASSERT_TRUE(R) << R.Error().details;
    EXPECT_EQ(R->Energy(), calc.Energy());
    EXPECT_EQ(R->IterationCount(), calc.IterationCount());
    EXPECT_EQ(&R->Density(), &calc.Density());

    auto again=calc.Converge({.NMaxIter=60});
    ASSERT_TRUE(again);
    EXPECT_NEAR(again->Energy(), R->Energy(), 1e-8) << "a re-converge from the converged state stays there";

    qchem::BoundedSCF bounded("the failure arm is the subject: one iteration cannot converge water");
    MoleculeCalculation cut(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}), {.NMaxIter=1});
    auto F=cut.Result();
    ASSERT_FALSE(F);
    EXPECT_EQ(F.Error().why, SCFFailure::Why::NotConverged);
    EXPECT_EQ(F.Error().iterations, 1u);
    EXPECT_TRUE(std::isfinite(F.Error().lastEnergy)) << "the last iterate is kept as a DIAGNOSTIC";
    EXPECT_FALSE(F.Error().details.empty());
}

namespace
{
//! The newest FILED run whose key starts with \a prefix ("<ID> build@<time>"); null when none.
const qchem::report::json* FiledRun(const std::string& prefix)
{
    const auto& G=qchem::report::GlobalReport();
    const qchem::report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    return last;
}
}

//! BUILD, THEN RUN (OpenWork "A molecule has no build-only step", user 2026-11-02: "the GUI needs a trigger to read all setup
//! panels and start the actual work").  Prepare files the whole setup -- basis.perIrrep, every grids.* entry -- into the build run
//! its CALLER opened (a facade never opens one, user ruling D 2026-11-02; the GUI's Build opens "<ID> build"), BEFORE any SCF
//! iteration; the first Converge then reuses that graph and lands on the constructor's answer bit for bit.
TEST(CalculationFacade, PrepareFilesTheWholeSetupInItsCallersBuildRunBeforeAnySCF)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"});
    size_t rows=0;
    std::unique_ptr<MoleculeCalculation> calc;
    { qchem::report::Run run(GetMolecule("H2O").ID()+" build"); calc=MoleculeCalculation::Prepare(setup, [&](const SCFIterator::IterationRow&){ rows++; }); }
    EXPECT_EQ(rows, 0u) << "Prepare ran an SCF iteration";
    const auto* build=FiledRun(GetMolecule("H2O").ID()+" build@");
    ASSERT_NE(build, nullptr) << "no filed '<ID> build' run";
    ASSERT_TRUE(build->contains("basis") && (*build)["basis"].contains("perIrrep")) << build->dump(1);
    EXPECT_FALSE((*build)["basis"]["perIrrep"].empty());
    ASSERT_TRUE(build->contains("hamiltonian")) << build->dump(1);
    const auto& ham=(*build)["hamiltonian"];
    EXPECT_TRUE(ham.contains("xcQuadrature") && ham["xcQuadrature"]["points"].get<long>()>0) << "hamiltonian.xcQuadrature missing from the build run";
    ASSERT_TRUE(ham.contains("Hartree")) << build->dump(1);
    // A Coulomb fit has no integration grid (doc/Pins.md, user 2026-11-11): the Hartree term files its fit basis and NO mesh; the SAD seed
    // owns the points it projects itself on, and files them under its own section, never under a term.
    EXPECT_FALSE(ham["Hartree"].contains("auxiliaryFit")) << "the Hartree term files a mesh: " << ham["Hartree"].dump(1);
    EXPECT_FALSE(ham.contains("seedDensity")) << "the seed's grid files under the Hamiltonian";
    std::vector<std::string> seedAt;   // every path at which a "seedDensity" section sits
    using J=std::remove_cvref_t<decltype(*build)>;
    std::function<void(const J&, const std::string&)> find=[&](const J& j, const std::string& path)
    {
        if (!j.is_object()) return;
        for (auto it=j.begin(); it!=j.end(); ++it) { if (it.key()=="seedDensity") seedAt.push_back(path+it.key()); find(it.value(), path+it.key()+"."); }
    };
    find(*build, "");
    ASSERT_FALSE(seedAt.empty()) << "the SAD seed's grid is not filed: " << build->dump(1);
    for (const auto& at : seedAt) EXPECT_EQ(at.find("hamiltonian"), std::string::npos) << at;
    // No SCF ran: no scf header (Converge files it) -- an "scf" key may hold a setup fact the seed build files (xcRhoRoutes)
    EXPECT_FALSE(build->contains("scf") && (*build)["scf"].contains("standard")) << "the build run holds an SCF: " << (*build)["scf"].dump(1);

    auto r=calc->Converge();
    ASSERT_TRUE(r) << r.Error().details;
    EXPECT_GT(rows, 0u) << "the observer attached at Prepare sees the first Converge";
    MoleculeCalculation once(setup);
    EXPECT_DOUBLE_EQ(r->Energy(), once.Energy()) << "Prepare + Converge must be the constructor, bit for bit";
    EXPECT_EQ(r->IterationCount(), once.IterationCount());
}

//! THE CURRENT ITERATE, SAMPLED MID-RUN (OpenWork §2 "grab 3D plot data", user GUI pass 2026-11-02): from the observer -- the
//! compute thread between iterations, the binding's row callback -- \c MakeIterateDensity builds the density of the latest
//! diagonalisation.  Before any iteration there is none (a refusal, not garbage); mid-run it is a density of N electrons that is
//! not yet the answer; after the last iteration it IS the converged \c Density(), point for point.
TEST(CalculationFacade, TheCurrentIteratesDensityIsSampledFromTheObserverAndEndsAsTheConvergedOne)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"});
    const rvec3_t on(0.0, 0.0, 0.3), off(0.0, 1.0, 1.2);
    std::vector<double> mid;
    std::unique_ptr<MoleculeCalculation> calc;
    calc=MoleculeCalculation::Prepare(setup, [&](const SCFIterator::IterationRow& row)
    {
        const auto rho=calc->MakeIterateDensity();
        if (row.iteration==2) mid={(*rho)(on), (*rho)(off)};
        EXPECT_EQ(calc->MakeIterateSpinDensity()->operator()(off), 0.0) << "unpolarized: the zero field";
    });
    EXPECT_THROW((void)calc->MakeIterateDensity(), std::logic_error) << "no iteration has run: nothing to sample";
    auto r=calc->Converge();
    ASSERT_TRUE(r) << r.Error().details;
    ASSERT_EQ(mid.size(), 2u) << "the observer sampled at iteration 2";
    EXPECT_NE(mid[1], calc->Density()(off)) << "a mid-run iterate is not the answer";
    const auto last=calc->MakeIterateDensity();
    EXPECT_DOUBLE_EQ((*last)(on),  calc->Density()(on));
    EXPECT_DOUBLE_EQ((*last)(off), calc->Density()(off));
}

//! EVERY DEFAULT BUILDS (controller ruling 2026-11-07, after the sto-3g default that had no file): a molecule stating ONLY its basis
//! -- and, for the DFT defaults (XC mesh, Vxc fit, Hartree route), only its model besides -- builds its whole graph.  A default that
//! cannot build is a silent trap; this is the test that would have caught it.  The basis itself has no default: unstated, it is
//! refused at Resolve, naming the registered bases that cover the structure.
TEST(CalculationFacade, EveryDefaultBuildsAndAnUnstatedBasisIsRefusedNamingTheBasesThatCoverIt)
{
    for (const CalcOptions& opts : {CalcOptions{.basis="dzvp"}, CalcOptions{.basis="dzvp", .model=Model::LDA}})
    {
        std::unique_ptr<MoleculeCalculation> calc;
        EXPECT_NO_THROW(calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), opts)))
            << "model " << int(opts.model) << ": a default failed to build";
    }
    try
    {
        (void)MoleculeSetup::Resolve(GetMolecule("H2O"), {});
        FAIL() << "an unstated basis must be refused";
    }
    catch (const std::invalid_argument& e)
    {
        const std::string what=e.what();
        EXPECT_NE(what.find("no basis stated"), std::string::npos) << what;
        EXPECT_NE(what.find("H, O: "), std::string::npos) << "the elements are named: " << what;
        EXPECT_NE(what.find("dzvp"), std::string::npos) << "the bases that cover them are named: " << what;
    }
}

//! Two Gaussian auxiliary fits in one run: each files its fit basis and its private mesh in its own term, so neither overwrites the other
//! (one grids.auxiliaryFit key did, until the terms got scopes; user ruling 2026-11-04).
//! A STATED BASIS FILE (CalcOptions::basisFile, lane G 2026-11-10: the GUI's Basis Set Exchange downloads): the same file
//! contents at another path build the same calculation bit for bit; the report states the path and the checksum; a saved state
//! refuses an EDITED file (the checksum, not the path, identifies it) and a run on the shipped name alone.
TEST(CalculationFacade, AStatedBasisFileIsItsContentsAndASavedStateRefusesAnEditedOne)
{
    namespace fs=std::filesystem;
    const fs::path dir=fs::temp_directory_path()/"qchem_UT_basisFile";
    fs::create_directories(dir);
    const fs::path copy=dir/"mine.bsd", edited=dir/"edited.bsd", state=dir/"water.h5";
    fs::copy_file(BasisSet::Gaussian::BasisFile("dzvp.bsd"), copy, fs::copy_options::overwrite_existing);
    fs::copy_file(copy, edited, fs::copy_options::overwrite_existing);
    std::ofstream(edited, std::ios::app) << "\nan edit after the last block\n";

    MoleculeCalculation named(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    std::unique_ptr<MoleculeCalculation> mine;
    { qchem::report::Run run("H2O basisFile"); mine=std::make_unique<MoleculeCalculation>(
        MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .basisFile=copy.string()})); }
    ASSERT_TRUE(named.IsConverged() && mine->IsConverged());
    EXPECT_DOUBLE_EQ(mine->Energy(), named.Energy()) << "the same contents are the same basis";
    EXPECT_EQ(mine->IterationCount(), named.IterationCount());
    const auto* run=FiledRun("H2O basisFile@");
    ASSERT_NE(run, nullptr);
    EXPECT_EQ((*run)["basis"]["file"].get<std::string>(), copy.string());
    EXPECT_EQ((*run)["basis"]["checksum"].get<std::string>(), BasisSet::Gaussian::BasisFileChecksum(copy.string()));

    mine->SaveState(state.string());
    auto judged=[&](const CalcOptions& o){ return MoleculeCalculation(MoleculeSetup::Resolve(GetMolecule("H2O"), o)).JudgeState(state.string()); };
    const auto same=judged({.basis="dzvp", .basisFile=copy.string()});
    ASSERT_TRUE(same) << same.Error().details;
    EXPECT_TRUE(same->empty()) << "the same file: an exact resume";
    const auto ed=judged({.basis="dzvp", .basisFile=edited.string()});
    ASSERT_FALSE(ed) << "an edited file passed as the same basis";
    EXPECT_NE(ed.Error().details.find("basis.fileChecksum"), std::string::npos) << ed.Error().details;
    const auto shipped=judged({.basis="dzvp"});
    ASSERT_FALSE(shipped) << "a state made on a stated file restarted a run on the shipped name";
    EXPECT_NE(shipped.Error().details.find("basis.fileChecksum"), std::string::npos) << shipped.Error().details;
    fs::remove_all(dir);
}

TEST(CalculationFacade, EachAuxiliaryFitFilesInItsOwnTerm)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary});
    { qchem::report::Run run("H2O fits build"); auto calc=MoleculeCalculation::Prepare(setup); }
    const auto* build=FiledRun("H2O fits build@");
    ASSERT_NE(build, nullptr) << "no filed build run";
    ASSERT_TRUE(build->contains("hamiltonian")) << build->dump(1);
    const auto& ham=(*build)["hamiltonian"];
    for (const auto& [term,name] : {std::pair{"Hartree","A1_coul"}, std::pair{"Vxc_DFT","A1_exch"}})
    {
        ASSERT_TRUE(ham.contains(term)) << term;
        EXPECT_EQ(ham[term]["fitbasis"].front()["name"].get<std::string>(), name) << term;
        EXPECT_GT(ham[term]["fitbasis"].front()["nFunctions"].get<long>(), 0) << term;
    }
    // The Vxc fit projects v_xc on its own mesh; the Coulomb fit is analytic and files none (doc/Pins.md, user 2026-11-11).
    EXPECT_GT(ham["Vxc_DFT"]["auxiliaryFit"]["points"].get<long>(), 0);
    EXPECT_FALSE(ham["Hartree"].contains("auxiliaryFit")) << ham["Hartree"].dump(1);
}

//! The atom's Build -> Run twin.
TEST(CalculationFacade, AnAtomsPrepareFilesItsSetupBeforeAnySCF)
{
    const AtomCalcOptions opts{.type=AtomType::Gaussian, .model=Model::LDA};
    std::unique_ptr<AtomCalculation> calc;
    { qchem::report::Run run(Atom(8, 0.0, rvec3_t(0,0,0)).ID()+" build"); calc=AtomCalculation::Prepare(8, 0, opts); }
    const auto* build=FiledRun(Atom(8, 0.0, rvec3_t(0,0,0)).ID()+" build@");
    ASSERT_NE(build, nullptr) << "no filed '<ID> build' run";
    EXPECT_TRUE(build->contains("basis") && (*build)["basis"].contains("perIrrep")) << build->dump(1);
    EXPECT_TRUE(build->contains("hamiltonian") && (*build)["hamiltonian"].contains("xcQuadrature")) << build->dump(1);
    EXPECT_FALSE(build->contains("scf") && (*build)["scf"].contains("standard")) << "the build run holds an SCF";
    ASSERT_TRUE(calc->Converge());
    AtomCalculation once(8, 0, opts);
    EXPECT_DOUBLE_EQ(calc->Energy(), once.Energy());
    EXPECT_EQ(calc->IterationCount(), once.IterationCount());
}

//! THE ONE ACCURACY TARGET ON EVERY FACADE (V-FACADE-ACCURACY, user ruling 2026-11-03): an atom and a molecule carry eps in their
//! options, answer what their client consumes through the shared base (an energy), and -- since no term of theirs reads eps yet --
//! say so in the run report ("consumed by: none yet"), never silently.  The atom's exponent-pool preset is `basisPool`, not accuracy.
TEST(CalculationFacade, EveryFacadeCarriesTheAccuracyTargetAndSaysWhoConsumesIt)
{
    auto expectReported=[](const std::string& id, double eps)
    {
        const auto* run=FiledRun(id+"@");
        ASSERT_NE(run, nullptr) << "no filed '" << id << "' run";
        ASSERT_TRUE(run->contains("scf") && (*run)["scf"].contains("accuracy")) << run->dump(1);
        EXPECT_EQ((*run)["scf"]["accuracy"]["eps"].get<double>(), eps);
        EXPECT_NE((*run)["scf"]["accuracy"]["consumedBy"].get<std::string>().find("none yet"), std::string::npos);
    };
    {
        qchem::report::Run run("eps atom");
        AtomCalculation calc(2, 0, {.type=AtomType::Slater, .basisPool=BasisSetAccuracy::Low, .accuracy={.eps=3e-6}});
        const Calculation& base=calc;
        EXPECT_EQ(base.Consumes(), ClientConsumes::Energy);
    }
    expectReported("eps atom", 3e-6);
    {
        qchem::report::Run run("eps molecule");
        MoleculeCalculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .accuracy={.eps=4e-6}}), {});
        const Calculation& base=calc;
        EXPECT_EQ(base.Consumes(), ClientConsumes::Energy);
    }
    expectReported("eps molecule", 4e-6);
}

//! THE XC ρ AXIS'S POLICY (user ruling 2026-11-09, V-XC-ONE-TERM step 3): the Gaussian x δ cell (the A1 Vxc fit at xc.rho exact) RUNS
//! WITH A WARNING, and what needs H = dE/dD refuses on it by name -- GDM and the Ladder (its GDM rung).  The variational A1 cell
//! (xc.rho fit) reports its choice and takes GDM.
TEST(CalculationFacade, TheGaussianDeltaXCCellWarnsAndRefusesGDM)
{
    const auto water=GetMolecule("H2O");
    auto rowsNamed=[](const qchem::report::json& j, const std::string& name)
    {
        std::vector<qchem::report::json> out;
        std::function<void(const qchem::report::json&)> walk=[&](const qchem::report::json& x)
        {
            if (x.is_object())
                for (auto it=x.begin(); it!=x.end(); ++it)
                {
                    if (it.key()=="policy" && it.value().is_array())
                        for (const auto& r : it.value()) if (r.contains("name") && r["name"]==name) out.push_back(r);
                    walk(it.value());
                }
            else if (x.is_array()) for (const auto& y : x) walk(y);
        };
        walk(j);
        return out;
    };
    { qchem::report::Run run("H2O A1 rhoD resolve"); MoleculeSetup::Resolve(water, {.basis="dzvp", .model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary}); }
    const auto* warned=FiledRun("H2O A1 rhoD resolve@");
    ASSERT_NE(warned, nullptr);
    const auto w=rowsNamed(*warned, "XCRhoFittedPotential");
    ASSERT_FALSE(w.empty()) << warned->dump(1);
    for (const auto& r : w)
    {
        EXPECT_EQ(r["verdict"].get<std::string>(), "warn");
        EXPECT_NE(r["reason"].get<std::string>().find("pin 32"), std::string::npos) << "the warning names the pin it departs from";
    }

    for (const char* type : {"GDM", "Ladder"})
    {
        try
        {
            MoleculeSetup::Resolve(water, {.basis="dzvp", .model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary}, {.type=type});
            ADD_FAILURE() << type << " on the Gaussian x δ cell was not refused";
        }
        catch (const std::invalid_argument& e)
        {
            EXPECT_NE(std::string(e.what()).find("not dE/dD"), std::string::npos) << e.what();
        }
    }
    // The variational A1 cell takes GDM, and files its choice.
    { qchem::report::Run run("H2O A1 Prho resolve");
      EXPECT_NO_THROW(MoleculeSetup::Resolve(water, {.basis="dzvp", .model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary,
                                                     .xcRho=Hamiltonian::XCρ::Fit}, {.type="GDM"})); }
    const auto* chose=FiledRun("H2O A1 Prho resolve@");
    ASSERT_NE(chose, nullptr);
    const auto c=rowsNamed(*chose, "XCRhoProjected");
    ASSERT_FALSE(c.empty()) << chose->dump(1);
    for (const auto& r : c) EXPECT_EQ(r["verdict"].get<std::string>(), "choose");
    EXPECT_TRUE(rowsNamed(*chose, "XCRhoFittedPotential").empty()) << "the variational cell is not warned";
}
