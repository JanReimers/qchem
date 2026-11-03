// File: Calculation/tests/CalculationFacade.C  What the molecular facade hands a front end (lane G, 2026-10-30): the
// constructor's own SCF observed from its first iteration, every MO as a level, and the spin density -- a zero field when
// the run is unpolarized, never an absence.
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <memory>
import qchem.Calculation;
import qchem.AtomCalculation;   // AtomCalculation::Prepare (the atom's Build -> Run)
import qchem.Structure;         // Atom (the atom's report ID)
import qchem.StructureData;
import qchem.SCFParams;
import qchem.SCFIterator;
import qchem.Types;
import qchem.SCFCap;   // BoundedSCF: the failure arm runs a deliberately bounded SCF
import qchem.Reporting;   // the run report the molecular mesh files itself into

using namespace qchem;
using StructureData::GetMolecule;

TEST(CalculationFacade, TheConstructorsObserverSeesItsSCFFromTheFirstIteration)
{
    std::vector<SCFIterator::IterationRow> steps;
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}), {},
                     [&steps](const SCFIterator::IterationRow& p){ steps.push_back(p); });
    ASSERT_TRUE(calc.IsConverged());
    ASSERT_FALSE(steps.empty()) << "the observer must see the constructor's SCF, not only a later Converge";
    EXPECT_EQ(steps.front().iteration, 1u) << "iteration 1: the constructor's SCF from its start";
    EXPECT_EQ(steps.size(), calc.IterationCount()) << "every iteration, none skipped";
    EXPECT_NEAR(steps.back().ValueOf("Etotal"), calc.Energy(), 1e-10);
}

//! A molecular DFT run's report holds its XC mesh BEFORE the SCF iterates, at grids.xcQuadrature -- the path a solid files: the
//! facade builds the Hamiltonian in its "grids" context, the XC quadrature opens "xcQuadrature", the mesh builder writes one
//! level in (OpenWork §3 "A molecule's mesh never reports itself", user ruling 2026-11-02).
TEST(CalculationFacade, AMoleculesXCMeshIsInTheRunReportByIterationOne)
{
    bool seen=false, checked=false;
    qchem::report::Run run("H2O");   // the test is the driver: it opens the run the facade writes into (user ruling D, 2026-11-02)
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA}), {},
                     [&](const SCFIterator::IterationRow& p)
                     {
                         if (p.iteration!=1) return;
                         checked=true;
                         const auto& doc=qchem::report::CurrentRunReport();
                         if (!doc.contains("grids") || !doc["grids"].contains("xcQuadrature")) return;
                         const auto& xc=doc["grids"]["xcQuadrature"];
                         seen = xc.contains("points") && xc["points"].get<long>()>0 && xc["natom"].get<long>()==3;
                     });
    ASSERT_TRUE(calc.IsConverged());
    ASSERT_TRUE(checked) << "the observer saw no iteration 1";
    EXPECT_TRUE(seen) << "at iteration 1 the run report must already hold grids.xcQuadrature (the molecule's XC mesh)";
}

//! A PSEUDOPOTENTIAL run's PP quadrature is built in the Hamiltonian's own build phase, so it files at grids.ppQuadrature -- not
//! under basis.*, where the terms' lazy first-matrix build inside the SCF iterator's construction used to put it (CLAUDE.md "the
//! run report": fix the run sequence, never the path; OpenWork, the mesh-report remainder, 2026-11-02).
TEST(CalculationFacade, APseudopotentialRunsPPQuadratureIsInGridsNotBasis)
{
    bool seen=false, stray=true;
    qchem::report::Run run("H2O");   // the driver's run (user ruling D, 2026-11-02)
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"}), {},
                     [&](const SCFIterator::IterationRow& p)
                     {
                         if (p.iteration!=1) return;
                         const auto& doc=qchem::report::CurrentRunReport();
                         seen = doc.contains("grids") && doc["grids"].contains("ppQuadrature")
                             && doc["grids"]["ppQuadrature"]["points"].get<long>()>0;
                         stray = doc.contains("basis") && (doc["basis"].contains("ppLocal") || doc["basis"].contains("ppNonLocal")
                                                           || doc["basis"].contains("ppQuadrature"));
                     });
    ASSERT_TRUE(calc.IsConverged());
    EXPECT_TRUE(seen)   << "the PP quadrature must be in grids.ppQuadrature by iteration 1";
    EXPECT_FALSE(stray) << "a PP mesh filed under basis: it was built inside the SCF iterator's construction again";
}

TEST(CalculationFacade, LevelsListEveryMOAscendingAndTheOccupiedOnesAreTheOrbitals)
{
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
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
    Calculation closed(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    ASSERT_TRUE(closed.IsConverged());
    EXPECT_EQ(closed.GetSpinDensity()(rvec3_t(0.0, 0.0, 0.3)), 0.0);

    Calculation triplet(MoleculeSetup::Resolve(GetMolecule("O2"), {.basis="dzvp", .multiplicity=3}));
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
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    auto R=calc.Result();
    ASSERT_TRUE(R) << R.Error().details;
    EXPECT_EQ(R->Energy(), calc.Energy());
    EXPECT_EQ(R->IterationCount(), calc.IterationCount());
    EXPECT_EQ(&R->Density(), &calc.Density());

    auto again=calc.Converge({.NMaxIter=60});
    ASSERT_TRUE(again);
    EXPECT_NEAR(again->Energy(), R->Energy(), 1e-8) << "a re-converge from the converged state stays there";

    qchem::BoundedSCF bounded("the failure arm is the subject: one iteration cannot converge water");
    Calculation cut(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}), {.NMaxIter=1});
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
    std::unique_ptr<Calculation> calc;
    { qchem::report::Run run(GetMolecule("H2O").ID()+" build"); calc=Calculation::Prepare(setup, [&](const SCFIterator::IterationRow&){ rows++; }); }
    EXPECT_EQ(rows, 0u) << "Prepare ran an SCF iteration";
    const auto* build=FiledRun(GetMolecule("H2O").ID()+" build@");
    ASSERT_NE(build, nullptr) << "no filed '<ID> build' run";
    ASSERT_TRUE(build->contains("basis") && (*build)["basis"].contains("perIrrep")) << build->dump(1);
    EXPECT_FALSE((*build)["basis"]["perIrrep"].empty());
    ASSERT_TRUE(build->contains("grids")) << build->dump(1);
    for (const char* g : {"xcQuadrature", "ppQuadrature", "auxiliaryFit"})
        EXPECT_TRUE((*build)["grids"].contains(g) && (*build)["grids"][g]["points"].get<long>()>0) << "grids." << g << " missing from the build run";
    // No SCF ran: no scf header (Converge files it) -- an "scf" key may hold a setup fact the seed build files (xcRhoRoutes)
    EXPECT_FALSE(build->contains("scf") && (*build)["scf"].contains("standard")) << "the build run holds an SCF: " << (*build)["scf"].dump(1);

    auto r=calc->Converge();
    ASSERT_TRUE(r) << r.Error().details;
    EXPECT_GT(rows, 0u) << "the observer attached at Prepare sees the first Converge";
    Calculation once(setup);
    EXPECT_DOUBLE_EQ(r->Energy(), once.Energy()) << "Prepare + Converge must be the constructor, bit for bit";
    EXPECT_EQ(r->IterationCount(), once.IterationCount());
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
    EXPECT_TRUE(build->contains("grids") && (*build)["grids"].contains("xcQuadrature")) << build->dump(1);
    EXPECT_FALSE(build->contains("scf") && (*build)["scf"].contains("standard")) << "the build run holds an SCF";
    ASSERT_TRUE(calc->Converge());
    AtomCalculation once(8, 0, opts);
    EXPECT_DOUBLE_EQ(calc->Energy(), once.Energy());
    EXPECT_EQ(calc->IterationCount(), once.IterationCount());
}
