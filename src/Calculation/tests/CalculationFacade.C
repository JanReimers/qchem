// File: Calculation/tests/CalculationFacade.C  What the molecular facade hands a front end (lane G, 2026-10-30): the
// constructor's own SCF observed from its first iteration, every MO as a level, and the spin density -- a zero field when
// the run is unpolarized, never an absence.
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>
#include <set>
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
import qchem.Mesh;        // qcMesh::AngularKind (the atom's radial-only XC set)
import qchem.SolidCalculation;   // the GPW twin of the no-PP-quadrature gate
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData::SIPP_SR

using namespace qchem;
using StructureData::GetMolecule;

TEST(CalculationFacade, ThePreparedObserverSeesTheSCFFromTheFirstIteration)
{
    std::vector<SCFIterator::IterationRow> steps;
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp"}}), [&steps](const SCFIterator::IterationRow& p){ steps.push_back(p); }).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->IsConverged());
    ASSERT_FALSE(steps.empty()) << "the observer handed to Prepare must see the first Converge";
    EXPECT_EQ(steps.front().iteration, 1u) << "iteration 1: the SCF from its start";
    EXPECT_EQ(steps.size(), calc->IterationCount()) << "every iteration, none skipped";
    EXPECT_NEAR(steps.back().ValueOf("Etotal"), calc->Energy(), 1e-10);
}

//! A molecular DFT run's report holds its XC mesh BEFORE the SCF iterates, at hamiltonian.nonlinearMesh -- the path a solid files:
//! the facade builds the Hamiltonian in its "hamiltonian" context, the nonlinear mesh opens "nonlinearMesh", the mesh builder writes
//! one level in (OpenWork §3 "A molecule's mesh never reports itself", user rulings 2026-11-02, 2026-11-04).
TEST(CalculationFacade, AMoleculesXCMeshIsInTheRunReportByIterationOne)
{
    bool seen=false, checked=false;
    qchem::report::Run run("H2O");   // the test is the driver: it opens the run the facade writes into (user ruling D, 2026-11-02)
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .molecule={.basis="dzvp"}}), [&](const SCFIterator::IterationRow& p)
                     {
                         if (p.iteration!=1) return;
                         checked=true;
                         const auto& doc=qchem::report::CurrentRunReport();
                         if (!doc.contains("hamiltonian") || !doc["hamiltonian"].contains("nonlinearMesh")) return;
                         const auto& xc=doc["hamiltonian"]["nonlinearMesh"];
                         seen = xc.contains("points") && xc["points"].get<long>()>0 && xc["natom"].get<long>()==3;
                     }).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->IsConverged());
    ASSERT_TRUE(checked) << "the observer saw no iteration 1";
    EXPECT_TRUE(seen) << "at iteration 1 the run report must already hold hamiltonian.nonlinearMesh (the molecule's XC mesh)";
}

//! A PSEUDOPOTENTIAL run's PP quadrature is built in the Hamiltonian's own build phase, so it files at hamiltonian.ppQuadrature -- not
//! under basis.*, where the terms' lazy first-matrix build inside the SCF iterator's construction used to put it (CLAUDE.md "the
//! run report": fix the run sequence, never the path; OpenWork, the mesh-report remainder, 2026-11-02).
// A Gaussian basis assembles the PP terms in closed form, so a molecular PP run builds no PP mesh at all.
TEST(CalculationFacade, APseudopotentialOnAGaussianBasisBuildsNoPPQuadrature)
{
    bool built=true, stray=true;
    qchem::report::Run run("H2O");   // the driver's run (user ruling D, 2026-11-02)
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="dzvp"}}), [&](const SCFIterator::IterationRow& p)
                     {
                         if (p.iteration!=1) return;
                         const auto& doc=qchem::report::CurrentRunReport();
                         built = !doc.contains("hamiltonian") || doc["hamiltonian"].contains("ppQuadrature");
                         stray = doc.contains("basis") && (doc["basis"].contains("ppLocal") || doc["basis"].contains("ppNonLocal")
                                                           || doc["basis"].contains("ppQuadrature"));
                     }).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->IsConverged());
    EXPECT_FALSE(built) << "a PP quadrature was built for a basis that assembles the PP terms in closed form";
    EXPECT_FALSE(stray) << "a PP mesh filed under basis";
}

namespace
{
//! Is \a key a member of any object at any depth of \a j?  (A report scope nests where the run sequence calls in, so a gate that a scope
//! is ABSENT searches the whole report rather than one path, which would rot with the sequence.)
template <class J> bool HasKeyAnywhere(const J& j, const std::string& key)
{
    if (j.is_object())
    {
        if (j.contains(key)) return true;
        for (auto it=j.begin(); it!=j.end(); ++it) if (HasKeyAnywhere(it.value(), key)) return true;
    }
    else if (j.is_array())
        for (const auto& e : j) if (HasKeyAnywhere(e, key)) return true;
    return false;
}
}

//! The GPW twin (V-MESH-FACTORY step 1, controller 2026-11-21): a GTH projector set has a closed-Gaussian face, so the GPW basis assembles
//! the nonlocal PP analytically (GPW_Evaluator's SpeciesProjectorSet_Gaussian branch) and builds NO `ppNonLocal` mesh.  The mesh is the
//! fallback for a projector model without that face; this gate keeps the analytic path from silently falling back when the PP mesh's
//! construction moves (its fallback arm is gated by GPW.AnalyticSeparablePPMatchesMesh).
TEST(CalculationFacade, AGTHPseudopotentialOnGPWBuildsNoPPNonLocalMesh)
{
    bool reached=false, built=true;
    qchem::report::Run run("Si gpw gth");
    RunOptions o; o.Nelec=8; o.solid.densityEcut=20.0; o.label="gpw gth no pp mesh"; o.accuracy.eps=0.0;   // stated raster: no pilot
    auto calc=SolidCalculation::Prepare(SolidSetup::Resolve(Lattice_3D(UnitCell(*Materials::Get("Si_diamond").cell), ivec3_t(1,1,1)),
                                                            {.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o),
                                        [&](const SCFIterator::IterationRow& p)
                                        {
                                            if (p.iteration!=1) return;
                                            const auto& doc=qchem::report::CurrentRunReport();
                                            reached=doc.contains("hamiltonian");
                                            built  =HasKeyAnywhere(doc, "ppNonLocal");
                                        }).TakeValue();
    ASSERT_TRUE(calc->Converge());
    ASSERT_TRUE(reached) << "the report never got to the Hamiltonian, so the absence below would prove nothing";
    EXPECT_FALSE(built) << "a GTH run built the ppNonLocal mesh: the analytic KB path fell back to the mesh quadrature";
}

TEST(CalculationFacade, LevelsListEveryMOAscendingAndTheOccupiedOnesAreTheOrbitals)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->IsConverged());
    const auto& levels=calc->GetLevels();
    ASSERT_GT(levels.size(), calc->NumOccupied()) << "the virtuals are listed too";
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
    ASSERT_EQ(nocc, calc->NumOccupied());
    EXPECT_EQ(levels[nocc-1].phi, &calc->HOMO()) << "the highest occupied level IS the facade's HOMO";
    EXPECT_GT(levels[nocc].e, levels[nocc-1].e) << "a gap to the LUMO";
}

TEST(CalculationFacade, AnUnpolarizedSpinDensityIsTheZeroFieldAndAPolarizedOneIsNot)
{
    auto closed=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp"}})).TakeValue(); (void)closed->Converge();
    ASSERT_TRUE(closed->IsConverged());
    EXPECT_EQ(closed->GetSpinDensity()(rvec3_t(0.0, 0.0, 0.3)), 0.0);

    auto triplet=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("O2"), {.multiplicity=3, .molecule={.basis="dzvp"}})).TakeValue(); (void)triplet->Converge();
    ASSERT_TRUE(triplet->IsConverged());
    // The two unpaired π* electrons put the moment off the bond axis, around the atoms.
    EXPECT_GT(triplet->GetSpinDensity()(rvec3_t(0.8, 0.0, 1.141)), 1e-3);
    size_t up=0, down=0; double nUp=0, nDown=0;
    for (const auto& L : triplet->GetLevels())
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
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    auto R=calc->Result();
    ASSERT_TRUE(R) << R.Error().details;
    EXPECT_EQ(R->Energy(), calc->Energy());
    EXPECT_EQ(R->IterationCount(), calc->IterationCount());
    EXPECT_EQ(&R->Density(), &calc->Density());

    auto again=calc->Converge([](SCFParams& p){ p.NMaxIter=60; });
    ASSERT_TRUE(again);
    EXPECT_NEAR(again->Energy(), R->Energy(), 1e-8) << "a re-converge from the converged state stays there";

    qchem::BoundedSCF bounded("the failure arm is the subject: one iteration cannot converge water");
    auto cut=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.scf={.NMaxIter=1},.molecule={.basis="dzvp"}})).TakeValue(); (void)cut->Converge();
    auto F=cut->Result();
    ASSERT_FALSE(F);
    EXPECT_EQ(F.Error().why, SCFFailure::Why::NotConverged);
    EXPECT_EQ(F.Error().iterations, 1u);
    EXPECT_TRUE(std::isfinite(F.Error().lastEnergy)) << "the last iterate is kept as a DIAGNOSTIC";
    EXPECT_FALSE(F.Error().details.empty());
}

//! A STATED linear ceiling beside the Fock DIIS converges (lane G's GUI case, 2026-11-11: water LDA dzvp at a stated start 0.2 and
//! ceiling 0.3 hit the 100-iteration cap on fa324a0b, [F,D] wandering 1e-3 .. 0.7 at Lin 0.30 with DIIS at depth 2-3; unstated it
//! converges in 8).  A user-reachable combination: it must converge or be refused at Resolve, never cap.
TEST(CalculationFacade, AStatedLinearCeilingBesideTheFockDIISConverges)
{
    SCFParams p; p.StartingRelaxRo=0.2; p.MaxRelaxRo=0.3;
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .scf=p, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    auto R=calc->Result();
    ASSERT_TRUE(R) << R.Error().details;
    auto free=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .molecule={.basis="dzvp"}})).TakeValue(); (void)free->Converge();
    ASSERT_TRUE(free->Result());
    EXPECT_NEAR(R->Energy(), free->Energy(), 1e-8) << "the ceiling changes the path, never the answer";
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
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="dzvp"}});
    size_t rows=0;
    std::unique_ptr<MoleculeCalculation> calc;
    { qchem::report::Run run(GetMolecule("H2O").ID()+" build"); calc=MoleculeCalculation::Prepare(setup, [&](const SCFIterator::IterationRow&){ rows++; }).TakeValue(); }
    EXPECT_EQ(rows, 0u) << "Prepare ran an SCF iteration";
    const auto* build=FiledRun(GetMolecule("H2O").ID()+" build@");
    ASSERT_NE(build, nullptr) << "no filed '<ID> build' run";
    ASSERT_TRUE(build->contains("basis") && (*build)["basis"].contains("perIrrep")) << build->dump(1);
    EXPECT_FALSE((*build)["basis"]["perIrrep"].empty());
    ASSERT_TRUE(build->contains("hamiltonian")) << build->dump(1);
    const auto& ham=(*build)["hamiltonian"];
    EXPECT_TRUE(ham.contains("nonlinearMesh") && ham["nonlinearMesh"]["points"].get<long>()>0) << "hamiltonian.nonlinearMesh missing from the build run";
    // V-NONLINEAR-MESH (user ruling 2026-11-20 1(a)): the HAMILTONIAN holds the mesh and lists what uses it; the XC term is one client.
    if (ham.contains("nonlinearMesh"))
    {
        const auto& nm=ham["nonlinearMesh"];
        ASSERT_TRUE(nm.contains("clients")) << nm.dump(1);
        std::set<std::string> clients;
        for (const auto& c : nm["clients"]) clients.insert(c.get<std::string>());
        for (const char* c : {"Vxc_DFT", "DeltaRho", "SiteMoments"}) EXPECT_TRUE(clients.count(c)) << "nonlinearMesh does not list " << c;
    }
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
    auto once=MoleculeCalculation::Prepare(setup).TakeValue(); (void)once->Converge();
    EXPECT_DOUBLE_EQ(r->Energy(), once->Energy()) << "Prepare + Converge must be the constructor, bit for bit";
    EXPECT_EQ(r->IterationCount(), once->IterationCount());
}

//! THE CURRENT ITERATE, SAMPLED MID-RUN (OpenWork §2 "grab 3D plot data", user GUI pass 2026-11-02): from the observer -- the
//! compute thread between iterations, the binding's row callback -- \c MakeIterateDensity builds the density of the latest
//! diagonalisation.  Before any iteration there is none (a refusal, not garbage); mid-run it is a density of N electrons that is
//! not yet the answer; after the last iteration it IS the converged \c Density(), point for point.
TEST(CalculationFacade, TheCurrentIteratesDensityIsSampledFromTheObserverAndEndsAsTheConvergedOne)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp"}});
    const rvec3_t on(0.0, 0.0, 0.3), off(0.0, 1.0, 1.2);
    std::vector<double> mid;
    std::unique_ptr<MoleculeCalculation> calc;
    calc=MoleculeCalculation::Prepare(setup, [&](const SCFIterator::IterationRow& row)
    {
        const auto rho=calc->MakeIterateDensity();
        if (row.iteration==2) mid={(*rho)(on), (*rho)(off)};
        EXPECT_EQ(calc->MakeIterateSpinDensity()->operator()(off), 0.0) << "unpolarized: the zero field";
    }).TakeValue();
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
    for (const RunOptions& opts : {RunOptions{.molecule={.basis="dzvp"}}, RunOptions{.model=Model::LDA, .molecule={.basis="dzvp"}}})
    {
        std::unique_ptr<MoleculeCalculation> calc;
        EXPECT_NO_THROW(calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), opts)).TakeValue())
            << "model " << int(opts.model.value_or(Model::HF)) << ": a default failed to build";
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
//! A STATED BASIS FILE (RunOptions::basisFile, lane G 2026-11-10: the GUI's Basis Set Exchange downloads): the same file
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

    auto named=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp"}})).TakeValue(); (void)named->Converge();
    std::unique_ptr<MoleculeCalculation> mine;
    { qchem::report::Run run("H2O basisFile"); mine=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.molecule={.basis="dzvp", .basisFile=copy.string()}})).TakeValue(); (void)mine->Converge(); }
    ASSERT_TRUE(named->IsConverged() && mine->IsConverged());
    EXPECT_DOUBLE_EQ(mine->Energy(), named->Energy()) << "the same contents are the same basis";
    EXPECT_EQ(mine->IterationCount(), named->IterationCount());
    const auto* run=FiledRun("H2O basisFile@");
    ASSERT_NE(run, nullptr);
    EXPECT_EQ((*run)["basis"]["file"].get<std::string>(), copy.string());
    EXPECT_EQ((*run)["basis"]["checksum"].get<std::string>(), BasisSet::Gaussian::BasisFileChecksum(copy.string()));

    mine->SaveState(state.string());
    auto judged=[&](const RunOptions& o){ return MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), o)).TakeValue()->JudgeState(state.string()); };
    const auto same=judged({.molecule={.basis="dzvp", .basisFile=copy.string()}});
    ASSERT_TRUE(same) << same.Error().details;
    EXPECT_TRUE(same->empty()) << "the same file: an exact resume";
    const auto ed=judged({.molecule={.basis="dzvp", .basisFile=edited.string()}});
    ASSERT_FALSE(ed) << "an edited file passed as the same basis";
    EXPECT_NE(ed.Error().details.find("basis.fileChecksum"), std::string::npos) << ed.Error().details;
    const auto shipped=judged({.molecule={.basis="dzvp"}});
    ASSERT_FALSE(shipped) << "a state made on a stated file restarted a run on the shipped name";
    EXPECT_NE(shipped.Error().details.find("basis.fileChecksum"), std::string::npos) << shipped.Error().details;
    fs::remove_all(dir);
}

TEST(CalculationFacade, EachAuxiliaryFitFilesInItsOwnTerm)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary, .coulomb={.fit=Hamiltonian::CoulombFit::Auxiliary}, .molecule={.basis="dzvp"}});
    { qchem::report::Run run("H2O fits build"); ASSERT_TRUE(MoleculeCalculation::Prepare(setup)); }
    const auto* build=FiledRun("H2O fits build@");
    ASSERT_NE(build, nullptr) << "no filed build run";
    ASSERT_TRUE(build->contains("hamiltonian")) << build->dump(1);
    const auto& ham=(*build)["hamiltonian"];
    for (const auto& [term,name] : {std::pair{"Hartree","A2_coul"}, std::pair{"Vxc_DFT","A1_exch"}})
    {
        ASSERT_TRUE(ham.contains(term)) << term;
        EXPECT_EQ(ham[term]["fitbasis"].front()["name"].get<std::string>(), name) << term;
        EXPECT_GT(ham[term]["fitbasis"].front()["nFunctions"].get<long>(), 0) << term;
    }
    // The Vxc fit projects v_xc on the run's nonlinear mesh, SHARED (V-MESH-FACTORY §4.2: it used to build a private copy at the same
    // recipe, filed as "auxiliaryFit"); the Coulomb fit is analytic and files none (doc/Pins.md, user 2026-11-11).  So exactly one
    // mesh is filed, and neither term builds a second.
    ASSERT_TRUE(ham.contains("nonlinearMesh")) << ham.dump(1);
    EXPECT_GT(ham["nonlinearMesh"]["points"].get<long>(), 0) << ham["nonlinearMesh"].dump(1);
    for (const char* term : {"Vxc_DFT", "Hartree"})
        EXPECT_FALSE(ham[term].contains("auxiliaryFit") || ham[term].contains("nonlinearMesh")) << term << ": a second mesh\n" << ham[term].dump(1);
}

//! A run that ASKED for a state and could not write it must not read green (doc/OpenWork.md "state.save into a missing directory",
//! 2026-11-22): the base's Checkpoint used to catch the writer's throw and only warn, so the facade answered CONVERGED with no file.
//! The SCF itself converges here; only the save fails, and the answer says which.
TEST(CalculationFacade, ARunWhoseStateCannotBeSavedIsNotGreen)
{
    const std::string dir=(std::filesystem::temp_directory_path()/"qchem_no_such_dir_StateNotSaved").string();
    std::filesystem::remove_all(dir);
    const std::string path=dir+"/o_atom.h5";
    RunOptions opts{.model=Model::LDA, .atom={.type=AtomType::Gaussian}};
    opts.saveStateTo=path;
    auto calc=AtomCalculation::Prepare(AtomSetup::Resolve(8, 0, opts)).TakeValue();
    const auto r=calc->Converge();
    ASSERT_FALSE(r) << "a run whose state was not saved read CONVERGED";
    EXPECT_EQ(r.Error().why, SCFFailure::Why::StateNotSaved) << r.Error().details;
    EXPECT_NE(r.Error().details.find(path), std::string::npos) << r.Error().details;
    EXPECT_GT(r.Error().iterations, 0u);
    EXPECT_FALSE(std::filesystem::exists(path));
    // The same run, its directory now there: the next Converge saves, and the answer is green again.
    std::filesystem::create_directories(dir);
    EXPECT_TRUE(calc->Converge());
    EXPECT_TRUE(std::filesystem::exists(path));
    std::filesystem::remove_all(dir);
}

//! The atom's Build -> Run twin.
TEST(CalculationFacade, AnAtomsPrepareFilesItsSetupBeforeAnySCF)
{
    const RunOptions opts{.model=Model::LDA, .atom={.type=AtomType::Gaussian}};
    std::unique_ptr<AtomCalculation> calc;
    { qchem::report::Run run(Atom(8, 0.0, rvec3_t(0,0,0)).ID()+" build"); calc=AtomCalculation::Prepare(AtomSetup::Resolve(8, 0, opts)).TakeValue(); }
    const auto* build=FiledRun(Atom(8, 0.0, rvec3_t(0,0,0)).ID()+" build@");
    ASSERT_NE(build, nullptr) << "no filed '<ID> build' run";
    EXPECT_TRUE(build->contains("basis") && (*build)["basis"].contains("perIrrep")) << build->dump(1);
    EXPECT_TRUE(build->contains("hamiltonian") && (*build)["hamiltonian"].contains("nonlinearMesh")) << build->dump(1);
    EXPECT_FALSE(build->contains("scf") && (*build)["scf"].contains("standard")) << "the build run holds an SCF";
    ASSERT_TRUE(calc->Converge());
    auto once=AtomCalculation::Prepare(AtomSetup::Resolve(8, 0, opts)).TakeValue(); (void)once->Converge();
    EXPECT_DOUBLE_EQ(calc->Energy(), once->Energy());
    EXPECT_EQ(calc->IterationCount(), once->IterationCount());
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
        auto calc=AtomCalculation::Prepare(AtomSetup::Resolve(2, 0, {.accuracy={.eps=3e-6}, .atom={.type=AtomType::Slater, .basisPool=BasisSetAccuracy::Low}})).TakeValue(); (void)calc->Converge();
        const Calculation& base=*calc;
        EXPECT_EQ(base.Consumes(), ClientConsumes::Energy);
    }
    expectReported("eps atom", 3e-6);
    {
        qchem::report::Run run("eps molecule");
        auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(GetMolecule("H2O"), {.accuracy={.eps=4e-6}, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
        const Calculation& base=*calc;
        EXPECT_EQ(base.Consumes(), ClientConsumes::Energy);
    }
    expectReported("eps molecule", 4e-6);
}

//! THE XC ρ AXIS'S POLICY (user ruling 2026-11-09, V-XC-ONE-TERM step 3): the Gaussian x δ cell (the A1 Vxc fit at xc.rho exact) RUNS
//! WITH A WARNING, and what needs H = dE/dD refuses on it by name -- GDM and the Ladder (its GDM rung).  The variational A1 cell
//! (xc.rho fit) is retired: refused at setup by pin 28 (user 2026-11-15).
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
    { qchem::report::Run run("H2O A1 rhoD resolve"); MoleculeSetup::Resolve(water, {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary, .molecule={.basis="dzvp"}}); }
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
            MoleculeSetup::Resolve(water, {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary, .accelerator=AcceleratorTypeByName(type, "the test"), .molecule={.basis="dzvp"}});
            ADD_FAILURE() << type << " on the Gaussian x δ cell was not refused";
        }
        catch (const std::invalid_argument& e)
        {
            EXPECT_NE(std::string(e.what()).find("not dE/dD"), std::string::npos) << e.what();
        }
    }
    // The A1 P rho cell is RETIRED (user 2026-11-15, pin 28): xc.rho fit on the incomplete A1 fit is refused at setup, naming the pin, and the
    // refusal files its row.
    { qchem::report::Run run("H2O A1 Prho resolve");
      try
      {
          MoleculeSetup::Resolve(water, {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary, .xcRho=Hamiltonian::XCρ::Fit, .molecule={.basis="dzvp"}});
          ADD_FAILURE() << "xc.rho fit on the A1 fit was not refused";
      }
      catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("pin 28"), std::string::npos) << e.what(); } }
    const auto* refused=FiledRun("H2O A1 Prho resolve@");
    ASSERT_NE(refused, nullptr);
    const auto c=rowsNamed(*refused, "XCRhoProjectedRefused");
    ASSERT_FALSE(c.empty()) << refused->dump(1);
    for (const auto& r : c) EXPECT_EQ(r["verdict"].get<std::string>(), "refuse");
}

//! THE RADIAL-ONLY SET'S OLD SPELLING (V-XC-ONE-TERM step 5, controller 2026-11-14: backward compatibility, gated).  The atom's XC mesh
//! is the radial-only point set, AngularKind::None.  Before it had a name the atom spelled it Lebedev at degree 0 (the 1-direction
//! rule, the same points), so a caller stating that spelling must read exactly as before: the same mesh, bit for bit.  And None is
//! one mesh whatever angularDegree says (the rule does not read it).
TEST(CalculationFacade, AnAtomStatingLebedevDegreeZeroIsTheRadialOnlySetBitForBit)
{
    auto run=[](qcMesh::AngularKind angular, int degree)
    {
        RunOptions o{.model=Model::LDA, .atom={.type=AtomType::Gaussian, .basisPool=BasisSetAccuracy::Low}};
        AtomXCMesh(o).angular=angular; AtomXCMesh(o).angularDegree=degree;
        auto calc=AtomCalculation::Prepare(AtomSetup::Resolve(2, 0, o)).TakeValue(); (void)calc->Converge();
        return std::pair{calc->Energy(), calc->IterationCount()};
    };
    const auto none  =run(qcMesh::AngularKind::None,    0);
    const auto legacy=run(qcMesh::AngularKind::Lebedev, 0);
    const auto loose =run(qcMesh::AngularKind::None,    7);
    EXPECT_EQ(legacy.first, none.first)   << "Lebedev at degree 0 on the atom must read as the radial-only set, bit for bit";
    EXPECT_EQ(legacy.second, none.second);
    EXPECT_EQ(loose.first, none.first)    << "the radial-only set does not read angularDegree";
    EXPECT_EQ(kAtomMeshDefault.angular, qcMesh::AngularKind::None) << "the atom's default XC mesh is the radial-only set";
}

// A STATED CONFIGURATION (OpenWork row "All-electron ION seeds", 2026-11-15).  The atom refuses, before any SCF, a stated configuration
// that cannot run: the wrong electron count, a hole the radial atom cannot hold, a core that is not the pseudopotential's, a Dirac model.
// Each refusal names the policy and its reason; a right one builds.
TEST(AtomConfiguration, AStatedConfigurationThatCannotRunIsRefusedBeforeAnySCF)
{
    auto refused=[](int Z, int charge, const RunOptions& o, const char* expect)
    {
        // a refusal at Resolve throws (the setup's), one while building is Prepare's BuildRefusal (R9): either names the reason
        try
        {
            auto built=AtomCalculation::Prepare(AtomSetup::Resolve(Z, charge, o));
            if (built) { ADD_FAILURE() << "must be refused: " << o.atom.configuration; return; }
            EXPECT_EQ(built.Error().policy, "AtomConfiguration");
            EXPECT_NE(built.Error().reason.find(expect), std::string::npos) << built.Error().reason;
        }
        catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find(expect), std::string::npos) << e.what(); }
    };
    const RunOptions lda{.spin=SpinGroup::Polarized, .model=Model::LDA, .atom={.type=AtomType::Slater, .basisPool=BasisSetAccuracy::Low}};
    RunOptions o=lda;
    o.atom.configuration="[Ar] 3d^6 4s^2"; refused(26, 1, o, "holds 26 electrons where the run carries 25");
    o.atom.configuration="[Ne] 3s^1 4s^1"; refused(12, 0, o, "open below the occupied 4s");
    o.atom.configuration="[Ar] 3d^6 4s^";  refused(26, 1, o, "has no electron count");
    RunOptions dirac=lda; dirac.model=Model::DHF; dirac.atom.configuration="[Ar] 3d^6 4s^1";
    refused(26, 1, dirac, "Dirac");
    RunOptions pp=lda; pp.atom.type=AtomType::Gaussian; pp.atom.pseudopotential=true; pp.atom.valence=7;   // Mn2+ q7: 5 valence electrons, charge 25 - 5
    pp.atom.configuration="[Ne] 3s^2 3p^6 3d^5"; refused(25, 20, pp, "the pseudopotential replaces 18");
    pp.atom.configuration="[Ar] 3d^5";
    EXPECT_TRUE(AtomCalculation::Prepare(AtomSetup::Resolve(25, 20, pp))) << "Mn2+ q7 stated as its valence over the PP's own [Ar] core";
    o.atom.configuration="[Ar] 3d^6 4s^1";
    EXPECT_TRUE(AtomCalculation::Prepare(AtomSetup::Resolve(26, 1, o))) << "Fe+ in NIST's configuration";
}

// AN UNSTATED ALL-ELECTRON CATION (row 196 step 3, controller ruling 2026-11-16): Fe+ is NIST's 3d^6 4s^1 (not Mn's 3d^5 4s^2), Fe2+ the
// removal rule's 3d^6, each at Hund's multiplicity, and the policy row names the branch that fired with that multiplicity.  An anion keeps
// the isoelectronic neutral.  Built only (Prepare): the configuration and its row are filed before any SCF.
TEST(AtomConfiguration, AnUnstatedCationTakesNISTsFirstCationOrTheRemovalRuleAndSaysWhich)
{
    auto row=[](const std::string& runName, int Z, int charge)
    {
        {
            qchem::report::Run run(runName);
            EXPECT_TRUE(AtomCalculation::Prepare(AtomSetup::Resolve(Z, charge, {.spin=SpinGroup::Polarized, .model=Model::LDA, .atom={.type=AtomType::Slater, .basisPool=BasisSetAccuracy::Low}})));
        }
        std::vector<qchem::report::json> out;
        std::function<void(const qchem::report::json&)> walk=[&](const qchem::report::json& x)
        {
            if (!x.is_object()) return;
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.value("name","")=="AtomConfiguration") out.push_back(r);
                walk(it.value());
            }
        };
        const auto* filed=FiledRun(runName+"@");
        if (filed) walk(*filed);
        return out;
    };
    const auto fe1=row("cation Fe+", 26, 1);
    ASSERT_EQ(fe1.size(), 1u);
    EXPECT_EQ(fe1[0]["decided"], "Fe 1s^2 2s^2 2p^6 3s^2 3p^6 3d^6 4s^1 (multiplicity 6)") << fe1[0].dump(1);
    EXPECT_EQ(fe1[0]["detail"]["branch"], "NISTFirstCation");
    EXPECT_NE(fe1[0]["reason"].get<std::string>().find("NIST's first cation"), std::string::npos);
    const auto fe2=row("cation Fe2+", 26, 2);
    ASSERT_EQ(fe2.size(), 1u);
    EXPECT_EQ(fe2[0]["decided"], "Fe 1s^2 2s^2 2p^6 3s^2 3p^6 3d^6 (multiplicity 5)") << fe2[0].dump(1);
    EXPECT_EQ(fe2[0]["detail"]["branch"], "RemovalRule");
    EXPECT_EQ(fe2[0]["detail"]["multiplicity"], 5);
    const auto f1=row("anion F-", 9, -1);
    ASSERT_EQ(f1.size(), 1u);
    EXPECT_EQ(f1[0]["decided"], "F 1s^2 2s^2 2p^6");
    EXPECT_NE(f1[0]["reason"].get<std::string>().find("isoelectronic neutral"), std::string::npos);
}

// THE ATOM'S ITERATION COUNT IS THE RUN'S (OpenWork row "After the atom pilot re-converges", 2026-11-14): the pilot rebuilds the iterator on
// the grids it raised and re-converges, so the iterator's own count is only the last SCF's.  The Mn GTH q7 pseudo-atom (the +U
// projector's) raises xcRadial 50 -> 75 in round 0 and re-converges: IterationCount(), the result line and the report must all count
// every iteration the observer saw, and the report files the pilot's share.
TEST(AtomPilot, TheIterationCountIsTheRunsTotalAndThePilotsShareIsFiled)
{
    size_t seen=0;
    RunOptions o{.atom={.type=AtomType::Gaussian, .pseudopotential=true}};   // PseudoAtomOrbital's own recipe (src/Calculation/Imp/PseudoAtomOrbital.C)
    o.atom.valence=7;
    o.vxcFit=Hamiltonian::VxcFit::Auxiliary;
    std::vector<double> pool;   // its 16-exponent pool, 0.05 .. 200
    for (int i=0;i<16;i++) pool.push_back(0.05*std::pow(200.0/0.05, i/15.0));
    for (int l=0;l<=2;l++) o.atom.exponentsByL.push_back({l, pool});
    o.scf.MinVirial=std::nullopt; o.scf.NMaxIter=120;
    std::unique_ptr<AtomCalculation> calc;
    {
        qchem::report::Run run("atom pilot iterations");
        calc=AtomCalculation::Prepare(AtomSetup::Resolve(25, 25-7, o), [&seen](const auto&){ ++seen; }).TakeValue();
        ASSERT_TRUE(calc->Converge());
    }
    const auto* filed=FiledRun("atom pilot iterations@");
    ASSERT_TRUE(filed);
    const qchem::report::json* result=nullptr;
    std::function<void(const qchem::report::json&)> find=[&](const qchem::report::json& x)
    {
        if (!x.is_object()) return;
        for (auto it=x.begin(); it!=x.end(); ++it) { if (it.key()=="result" && it.value().is_object()) result=&it.value(); find(it.value()); }
    };
    find(*filed);
    ASSERT_TRUE(result) << "the atom files a result section";
    ASSERT_TRUE(result->contains("pilotIterations")) << "the case must raise a grid, or it does not test the pilot: " << result->dump(1);
    const size_t pilot=(*result)["pilotIterations"].get<long>();
    EXPECT_GT(pilot, 0u);
    EXPECT_LT(pilot, seen) << "the entry SCF ran too";
    EXPECT_EQ(calc->IterationCount(), seen) << "every iteration the observer saw, the entry SCF's and the pilot's";
    EXPECT_EQ((*result)["iterations"].get<long>(), (long)seen);
}
