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
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, SpinGroup
import qchem.SCFCap;                 // BoundedSCF (the trace test runs a deliberately bounded SCF)
import qchem.SolidCalculation;       // SolidCalculation -- the NAMED periodic front door (Step 4)
import qchem.Structure;              // FCCUnitCell
import qchem.Lattice_3D;             // Lattice_3D
import qchem.BasisSet;               // Real_BS
import qchem.BasisSet.Gaussian.Point.Factory;  // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.Types;                  // ivec3_t
import qchem.SCFIterator;            // SCFParams
import qchem.Calculation;            // the molecular facade (the golden molecular trace)
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
            {.NMaxIter = 40, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = 1e30,
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
            {.NMaxIter = 60, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = 1e30,
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
            {.NMaxIter = 60, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = 1e30,
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
            BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                        BasisSet::Gaussian::Engine::MnD,
                                        BasisSet::Gaussian::Angular::Cartesian));
        SCFParams par;
        par.NMaxIter=2; par.MinΔρPerElectron=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30;
        par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true;
        qchem::BoundedSCF bounded("the trace's columns are the subject: 2 iterations print them");
        SolidCalculation calc(SolidSetup::Resolve(lat, mol, {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}), par);
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

// R2.7: a live [F,D] gate reading exactly 0 at convergence throws (the G3 shape: the Null accelerator's [F,D] = 0); a gate that is
// off (1e30) may read 0, and a tiny but COMPUTED value is a verdict, not an absence.
TEST(SCFTrace, R2_7_ADeadGateThrows)
{
    using qchem::SCFIterator::RequireLiveGates;
    EXPECT_NO_THROW(RequireLiveGates(3e-7, 1e-5, 12));
    EXPECT_NO_THROW(RequireLiveGates(1e-300, 1e-5, 12));   // computed, just small
    EXPECT_NO_THROW(RequireLiveGates(0.0, 1e30, 12));      // [F,D] gate off
    EXPECT_THROW(RequireLiveGates(0.0, 1e-5, 12), std::logic_error);   // G3
}

// ---------------------------------------------------------------------------------------------------
// THE CONSOLE IS BYTE-IDENTICAL (2026-10-31, the one-home IterationRow): the trace block -- header, tolerance line, every row,
// and the breakdown / level tables printed after it (they inherit the stream state the rows leave) -- of one molecular and one
// solid run, captured from the code BEFORE the row refactor.  Wall-clock lines ("[t=...]") are the only ones dropped.  A
// legitimate numeric change re-captures these; a layout change must not happen silently.
namespace {
std::string TraceBlock(const std::string& trace)
{
    std::istringstream is(trace); std::string out; bool in=false;
    for (std::string line; std::getline(is, line); )
    {
        if (!in && line.find("Etotal")!=std::string::npos) in=true;
        if (in && line.find("[t=")==std::string::npos) out+=line+"\n";
    }
    return out;
}
const char* const kMolecularGolden = R"G(  #               Etotal      [F,D]     Δ[F,D]        Δρ      2+V/K  ρ_mix         accel   svMin    cfg       m_site
                            (1e-06)    (1e-08)   (1e-06)    (1e+02)
----------------------------------------------------------------------------------------------------
  1    -145.500070494997   3.70e+00  -1.00e+10  7.05e-01   8.03e-02  Lin 0.50      DIIS    ----                 ----
  2    -148.746398919215   1.21e+00  -2.49e+00  1.62e-01   6.29e-03  Lin 0.81      DIIS    ----       *         ----
  3    -149.560509951632   2.96e-01  -9.11e-01  4.65e-02  -1.44e-03  Lin 0.93      DIIS:1  ----       *         ----
  4    -149.627225750406   7.86e-02  -2.17e-01  7.78e-03  -4.14e-03  Lin 0.97      DIIS:2  8.1e-02              ----
  5    -149.629512467503   2.30e-02  -5.56e-02  2.47e-03  -4.00e-03  Lin 0.99      DIIS:3  3.6e-03              ----
  6    -149.629545937427   4.15e-03  -1.88e-02  4.31e-04  -4.10e-03  Lin 1.00      DIIS:4  8.4e-05              ----
  7    -149.629546950987   1.37e-03  -2.77e-03  1.04e-04  -4.08e-03  Lin 1.00      DIIS:4  1.2e-05              ----
  8    -149.629547019861   2.67e-04  -1.11e-03  8.36e-06  -4.08e-03  Lin 1.00      DIIS:4  3.5e-07              ----
  9    -149.629547027574   8.69e-05  -1.80e-04  3.95e-06  -4.08e-03  Lin 1.00      DIIS:4  2.2e-08              ----
 10    -149.629547027643   1.02e-05  -7.67e-05  4.64e-07  -4.08e-03  Lin 1.00      DIIS:2  5.6e-09              ----
 11    -149.629547027644   1.95e-06  -8.25e-06  6.46e-08  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 12    -149.629547027645   3.08e-07  -1.64e-06  3.73e-08  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 13    -149.629547027645   2.15e-07  -9.26e-08  2.30e-08  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 14    -149.629547027645   1.55e-07  -5.99e-08  1.47e-08  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 15    -149.629547027645   1.12e-07  -4.30e-08  9.45e-09  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 16    -149.629547027645   8.15e-08  -3.10e-08  6.17e-09  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 17    -149.629547027645   5.92e-08  -2.23e-08  4.07e-09  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 18    -149.629547027645   4.31e-08  -1.61e-08  2.71e-09  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 19    -149.629547027645   3.14e-08  -1.17e-08  1.84e-09  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
 20    -149.629547027645   2.29e-08  -8.47e-09  1.26e-09  -4.08e-03  Lin 0.50      DIIS:1  1.0e+00              ----
----------------------------------------------------------------------------------------------------------
Energy    Breakdown  Total:    -149.629547028  Kinetic:     149.022225338  Potential:    -298.651772366
Potential Breakdown  Een  :    -410.569585377  Eee    :     100.155759294  Eex      :     -16.283520341
Virial               V/K  :    -2.00407537660  Eee/Exc:      -6.150743647  
+-------------+--------------+--------------+-------------+--------------+-------------+
| Occ/Degen ↑ |      ϵ↑ (au) |  n,Symmetry  | Occ/Degen ↓ |      ϵ↓ (au) |  ϵ↑-ϵ↓ (au) |
+-------------+--------------+--------------+-------------+--------------+-------------+
|     1/1     | -20.74077336 | 1No Symmetry |     1/1     | -20.68244050 | -0.05833287 |
|     1/1     | -20.74026625 | 2No Symmetry |     1/1     | -20.68139778 | -0.05886847 |
|     1/1     |  -1.72418593 | 3No Symmetry |     1/1     |  -1.59095403 | -0.13323190 |
|     1/1     |  -1.20560393 | 4No Symmetry |     1/1     |  -0.99571703 | -0.20988690 |
|     1/1     |  -0.84776689 | 5No Symmetry |     1/1     |  -0.70332220 | -0.14444470 |
|     1/1     |  -0.76809809 | 7No Symmetry |      --     |           -- |          -- |
|      --     |           -- | 6No Symmetry |     1/1     |  -0.57794358 |          -- |
|     1/1     |  -0.55935104 | 8No Symmetry |     0/1     |   0.09995070 | -0.65930175 |
+-------------+--------------+--------------+-------------+--------------+-------------+
)G";
const char* const kSolidGolden = R"G(  #               Etotal      [F,D]       ΔE/E        Δρ   ρ_lost/N    XC ρ  ρ_mix         accel   svMin    cfg       gap  
                            (1e+30)    (1e-06)   (1e-03)
----------------------------------------------------------------------------------------------------
  1      -7.108903646312   2.70e-01  -1.00e+00  5.12e-02  -1.02e-09      DM  PulD 0.45     Null    ----          6.69e-02  
  2      -7.111740379421   1.23e-02  -3.99e-04  2.60e-02  -1.38e-09      DM  PulD 0.45     Null    ----          6.83e-02  
  3      -7.113252674769   6.46e-03  -2.13e-04  1.34e-02  -1.68e-09      DM  PulD 0.45     Null    ----          6.90e-02  
  4      -7.114071951821   3.97e-03  -1.15e-04  7.18e-03  -1.68e-09      DM  PulD 0.45     Null    ----          6.93e-02  
----------------------------------------------------------------------------------------------------------
Energy    Breakdown  Total:  -7.11407195182  Kinetic:   4.36964927810  Potential: -11.48372122992
Potential Breakdown  Een  :  -1.14841263938  Eee    :   0.90563312895  Eex      :  -2.54558416749
Een       Breakdown  Eloc :  -2.88715582538  EenNL  :   1.73874318599
Virial               V/K  :  -2.62806474823  Eee/Exc:  -0.35576632685  
+-----------+-------------+----------+
| Occ/Degen |      ϵ (au) | Symmetry |
+-----------+-------------+----------+
|    2/2    | -0.11636597 | 2(0,0,0) |
|    6/6    |  0.32414685 | 3(0,0,0) |
+-----------+-------------+----------+
[gpw] NOT converged iters=4 charge=8.00000000000 Eelec=1.58128560017 Etot=-7.1140719518  (Ekin=4.36964927810 Een=-1.14841263938 Eee=0.90563312895 Exc=-2.54558416749 Enn=-8.40046478619 E_alphaZ=-0.29489276580)
)G";
}

//! The molecular layout (O2 HF triplet, polarized: m_site; virial; cfg flags; DIIS counts and svMin).
TEST(SCFTrace, MolecularConsoleIsByteIdentical)
{
    const std::string got=TraceBlock(CapturedTrace([]{
        Molecule o2;
        o2.Insert(new Atom(8, 0, Vector3D<double>(0, 0,  1.141)));
        o2.Insert(new Atom(8, 0, Vector3D<double>(0, 0, -1.141)));
        Calculation calc(MoleculeSetup::Resolve(o2, {.basis = "dzvp", .multiplicity = 3}),
                         {.NMaxIter = 80, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-8, .MinVirial = 1e2, .MinFD = 1e-6,
                          .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true});
        (void)calc.Result();
    }));
    EXPECT_EQ(got, std::string(kMolecularGolden));
}

//! The solid layout (Si Γ, 4 bounded iterations: ΔE/E, ρ_lost/N, XC ρ, the permanent gap, PulD mixing, the Null accelerator).
TEST(SCFTrace, SolidConsoleIsByteIdentical)
{
    const std::string got=TraceBlock(CapturedTrace([]{
        FCCUnitCell cell(10.26);
        cell.AddAtom(14, {0,0,0});
        cell.AddAtom(14, {0.25,0.25,0.25});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        auto mol = std::shared_ptr<const BasisSet::Real_BS>(
            BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                        BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        SCFParams par;
        par.NMaxIter=4; par.MinΔρPerElectron=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30;
        par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true;
        qchem::BoundedSCF bounded("the golden trace: 4 iterations print every solid column");
        SolidCalculation calc(SolidSetup::Resolve(lat, mol, {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}), par);
        (void)calc.Result();
    }));
    EXPECT_EQ(got, std::string(kSolidGolden));
}

// THE OBSERVER RECEIVES THE CONSOLE'S ROW (2026-10-31): the same columns, titled as the header titles them and in its order,
// and each cell's text where the console prints it -- built from the ONE column list, so a column cannot reach one and not the
// other.  On the molecular polarized run (the frame, the variable columns and m_site).
TEST(SCFTrace, TheObserverReceivesTheConsoleRow)
{
    std::vector<SCFIterator::IterationRow> rows;
    const std::string trace=CapturedTrace([&rows]{
        Molecule o2;
        o2.Insert(new Atom(8, 0, Vector3D<double>(0, 0,  1.141)));
        o2.Insert(new Atom(8, 0, Vector3D<double>(0, 0, -1.141)));
        Calculation calc(MoleculeSetup::Resolve(o2, {.basis = "dzvp", .multiplicity = 3}),
                         {.NMaxIter = 80, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-8, .MinVirial = 1e2, .MinFD = 1e-6,
                          .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true},
                         [&rows](const SCFIterator::IterationRow& r){ rows.push_back(r); });
        (void)calc.Result();
    });
    ASSERT_FALSE(rows.empty());
    std::istringstream is(trace); std::vector<std::string> lines; for (std::string l; std::getline(is,l); ) lines.push_back(l);
    const std::string header=HeaderLine(trace);
    // titles: in the header, in order
    size_t at=0;
    for (const auto& c : rows.front().cells)
    {
        const size_t p=header.find(c.title, at);
        ASSERT_NE(p, std::string::npos) << "cell '" << c.title << "' is not a header column (in order): " << header;
        at=p+c.title.size();
    }
    EXPECT_EQ(rows.front().cells.front().title, "Etotal");
    EXPECT_FALSE(std::isnan(rows.front().ValueOf("[F,D]")));
    EXPECT_TRUE (std::isnan(rows.front().ValueOf("ρ_lost/N"))) << "a molecular layout has no grid column";
    // every row's cell texts, in order, in that iteration's console line
    for (const auto& r : rows)
    {
        const auto line=std::find_if(lines.begin(), lines.end(), [&](const std::string& l)
            { std::istringstream ls(l); size_t n=0; return (ls >> n) && n==r.iteration && l.find("e")!=std::string::npos; });
        ASSERT_NE(line, lines.end()) << "no console row for iteration " << r.iteration;
        size_t pos=0;
        for (const auto& c : r.cells)
        {
            if (c.text.empty()) continue;                     // a blank cfg
            const size_t p=line->find(c.text, pos);
            ASSERT_NE(p, std::string::npos) << "iteration " << r.iteration << ": cell " << c.title << " '" << c.text << "' not in: " << *line;
            pos=p+c.text.size();
        }
        EXPECT_NEAR(r.ValueOf("Etotal"), r.eb.GetTotalEnergy(), 0.0);
    }
}
