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
import qchem.SolidCalculation;       // SolidCalculation -- the NAMED periodic front door (Step 4)
import qchem.Structure;              // FCCUnitCell
import qchem.Lattice_3D;             // Lattice_3D
import qchem.BasisSet;               // Real_BS
import qchem.BasisSet.Molecule.Factory;  // Molecule::Factory, BasisSetData/Engine/Angular
import qchem.Types;                  // ivec3_t
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
// Combo 3 of 4: SOLID x PP.  The layout that, until this test, NO enabled test emitted -- the only two
// GPW runs with Verbose=true are both DISABLED_ -- even though it is what the lattice work reads on every
// run.  Driven through SolidCalculation (the named periodic front door) rather than GPW_SCF_UT.C's
// file-local RunGPW: the facade is the thing tests should dogfood, and it keeps this file independent of
// the parallel campaign's working file.  NMaxIter is tiny on purpose -- the trace HEADER is the subject,
// so there is no reason to pay for convergence.
TEST(SCFTrace, SolidPP_ShowsGridColumnsAndNoVirial)
{
    std::string trace = CapturedTrace([]{
        FCCUnitCell cell(10.26);                       // Si diamond, 2-atom basis (the SiGamma anchor cell)
        cell.AddAtom(14, {0,0,0});
        cell.AddAtom(14, {0.25,0.25,0.25});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        auto mol = std::shared_ptr<const BasisSet::Real_BS>(
            BasisSet::Molecule::Factory(BasisSet::Molecule::BasisSetData::SIPP_SR, &cell,
                                        BasisSet::Molecule::Engine::MnD,
                                        BasisSet::Molecule::Angular::Cartesian));
        SCFParams par;
        par.NMaxIter=2; par.MinΔρ=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30;
        par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true;
        SolidCalculation calc(lat, mol, {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);
        (void)calc.Result();   // NMaxIter=2: this deliberately does NOT converge -- the trace is the subject
    });
    const std::string h = HeaderLine(trace);
    ASSERT_FALSE(h.empty()) << "the solid layout printed no trace at all:\n" << trace;
    EXPECT_TRUE (Has(h,"ΔE/E"))   << "a collocation SCF is gated on ΔE, not Δ[F,D]: " << h;
    EXPECT_TRUE (Has(h,"ρ_lost")) << "the grid-charge leak is a permanent solid health column: " << h;
    EXPECT_TRUE (Has(h,"gap"))    << "solids show the frontier gap permanently, not behind ReportBandGap: " << h;
    EXPECT_FALSE(Has(h,"V/K"))    << "a pseudopotential run must NOT show a virial column: " << h;
}

// Combo 4 of 4: SOLID x NON-PP -- an ALL-ELECTRON periodic run, which should bring the virial column BACK.
//
// DISABLED for a PHYSICS reason, not a plumbing one: no all-electron periodic path exists yet (the APW/LAPW
// tests are themselves DISABLED_), so there is no way to make a solid Hamiltonian whose IsVirialValid() is
// true.  When that path lands, enable this -- and expect it to FAIL first, because SolidSCFIterator
// hard-codes the virial's absence instead of deriving it from IsVirialValid() the way the base layout now
// does.  Closing that gap is the point of the {Molecular,Solid}x{PP,Non-PP} mixin decomposition (V1.27).
TEST(SCFTrace, DISABLED_SolidNonPP_ShouldBringTheVirialColumnBack) {}
