// File SCFTrace.C  The SCF iteration TRACE LAYOUT -- which columns each kind of run prints.
//
// WHY THIS FILE EXISTS (V1.27).  The per-iteration trace was, until now, entirely unasserted: no test
// looked at a single column header, so dropping a column, mis-selecting a layout or breaking the UTF-8
// padding sailed through the whole suite.  That mattered most for the one layout NO enabled test even
// emitted -- the solid/GPW one, whose only Verbose runs are DISABLED_ -- which is exactly the layout the
// lattice work reads every run.
//
// It matters more now because the trace is about to be REFACTORED into column mixins, so that all four
// combos of {Molecular,Solid} x {PP,Non-PP} become expressible (user ruling 2026-08-10).  These tests are
// the specification of those four column sets, written BEFORE the refactor so it has something to be
// correct against.
//
// The assertions are deliberately on COLUMN PRESENCE, not on exact spacing: the widths are cosmetic and
// re-tuned freely, whereas "does a PP run still print a virial it has no right to?" is physics.
#include "gtest/gtest.h"
#include <sstream>
#include <ostream>
#include <iostream>
#include <string>
#include <functional>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, Pol
import qchem.SCFIterator;            // SCFParams
using namespace qchem;
using enum BasisSetAccuracy;

namespace {

// Run \a body with std::cout redirected, and hand back everything it printed.  The trace goes to cout
// (SCFIterator writes it directly), so this is the only way to see the layout from a test.
std::string CapturedTrace(const std::function<void()>& body)
{
    std::ostringstream sink;
    std::streambuf* saved = std::cout.rdbuf(sink.rdbuf());
    try { body(); } catch (...) { std::cout.rdbuf(saved); throw; }
    std::cout.rdbuf(saved);
    return sink.str();
}

// The HEADER line of the trace -- the one carrying the column labels ("#  Etotal  [F,D] ...").
// Returned empty when the run printed no trace at all, which is itself a failure worth seeing.
std::string HeaderLine(const std::string& trace)
{
    std::istringstream is(trace);
    for (std::string line; std::getline(is, line); )
        if (line.find("Etotal") != std::string::npos) return line;
    return {};
}

bool Has(const std::string& line, const std::string& col) {return line.find(col) != std::string::npos;}

} //anon namespace

// ---------------------------------------------------------------------------------------------------
// Combo 1 of 4: MOLECULAR/ATOMIC x NON-PP.  The virial theorem holds, so the run must show the 2+V/K
// column.  (He, Slater/Low: the cheapest all-electron closed shell that converges.)
TEST(SCFTrace, MolecularNonPP_ShowsTheVirialColumn)
{
    std::string trace = CapturedTrace([]{
        AtomCalculation calc(2, 0, {.type = AtomType::Slater, .accuracy = Low},
            {.NMaxIter = 40, .MinΔρ = 1e-6, .MinΔFD = 1e-6, .MinVirial = 1e30,
             .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});
        calc.Energy();
    });
    const std::string h = HeaderLine(trace);
    ASSERT_FALSE(h.empty()) << "no trace was printed at all:\n" << trace;
    EXPECT_TRUE(Has(h,"Etotal")) << h;
    EXPECT_TRUE(Has(h,"[F,D]")) << h;
    EXPECT_TRUE(Has(h,"V/K"))   << "an all-electron run must show the virial column: " << h;
    EXPECT_FALSE(Has(h,"ρ_lost")) << "the grid-charge leak is a collocation column: " << h;
}

// Combo 2 of 4: MOLECULAR/ATOMIC x PP.  A pseudopotential is not Coulombic, so the virial theorem does
// NOT hold and the column must be ABSENT -- this is V1.27's landed half (IsVirialValid), pinned.  Before
// it, a molecular PP run printed a virial the tree's own docs call invalid for it.
TEST(SCFTrace, MolecularPP_HidesTheVirialColumn)
{
    std::string trace = CapturedTrace([]{
        AtomCalculation calc(14, 10, {.type = AtomType::Slater, .accuracy = Medium, .pseudopotential = true},
            {.NMaxIter = 60, .MinΔρ = 1e-6, .MinΔFD = 1e-6, .MinVirial = 1e30,
             .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});
        calc.Energy();
    });
    const std::string h = HeaderLine(trace);
    ASSERT_FALSE(h.empty()) << "no trace was printed at all:\n" << trace;
    EXPECT_TRUE (Has(h,"Etotal")) << h;
    EXPECT_FALSE(Has(h,"V/K")) << "a PP run must NOT show a virial column (V1.27): " << h;
}

// The threshold sub-line tracks the columns: drop a column and its (tolerance) must go with it, or the
// two lines stop lining up.  Pinned separately because it is the half most easily forgotten.
TEST(SCFTrace, ThresholdSubLineMatchesTheColumns)
{
    std::string trace = CapturedTrace([]{
        AtomCalculation calc(14, 10, {.type = AtomType::Slater, .accuracy = Medium, .pseudopotential = true},
            {.NMaxIter = 60, .MinΔρ = 1e-6, .MinΔFD = 1e-6, .MinVirial = 1e30,
             .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});
        calc.Energy();
    });
    std::istringstream is(trace);
    std::string head, thresh;
    for (std::string line; std::getline(is, line); )
        if (line.find("Etotal") != std::string::npos) { head=line; std::getline(is,thresh); break; }
    ASSERT_FALSE(head.empty());
    // A PP run shows 3 gated columns ([F,D], Δ[F,D], Δρ) and no virial, so 3 (tolerance) cells.
    size_t n=0; for (size_t p=thresh.find('('); p!=std::string::npos; p=thresh.find('(',p+1)) ++n;
    EXPECT_EQ(n, 3u) << "threshold cells do not match the un-gated column set\n" << head << "\n" << thresh;
}

// ---------------------------------------------------------------------------------------------------
// Combos 3 and 4 of 4: SOLID x {PP, Non-PP}.  NOT COVERED YET, deliberately, and this is the marker.
//
// The solid layout ({..., ΔE/E, Δρ, ρ_lost/N, ..., gap} -- no virial, gap permanent) is emitted by NO
// enabled test: the only two GPW runs that set Verbose=true are both DISABLED_ (NaFGridContinuation,
// NaFFullBasisEigenTol).  So the layout the lattice work reads on every run is the one with no coverage.
//
// It is not covered HERE because driving a periodic SCF needs RunGPW/MakeBasisSR/FCCUnitCell, which are
// file-local to IntegrationTests/GPW_SCF_UT.C -- the parallel MnO campaign's working file.  Pinning these
// two needs one of:
//   (a) hoist the GPW test driver into the shared TestUtils module, then write them here; or
//   (b) add them to GPW_SCF_UT.C directly, once that file is quiet.
// (a) is the better shape -- a shared driver is useful well beyond this file -- but it edits their file
// either way, so it is a coordination call, not a technical one.
//
// What combo 4 (Solid x Non-PP, i.e. an all-electron periodic run -- the parked APW/LAPW work) would
// assert: that the virial column COMES BACK.  It currently cannot: SolidSCFIterator hard-codes the
// virial's absence instead of deriving it from IsVirialValid() the way the base layout now does.  That is
// precisely the gap the {Molecular,Solid}x{PP,Non-PP} mixin decomposition exists to close, so this test
// should be written to FAIL first and then made to pass by the refactor.
TEST(SCFTrace, DISABLED_SolidLayoutColumns) {}
TEST(SCFTrace, DISABLED_SolidNonPP_ShouldBringTheVirialColumnBack) {}
