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
#include <optional>   // std::nullopt: an OFF SCF gate
#include <sstream>
#include <ostream>
#include <iostream>
#include <string>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <type_traits>
import qchem.AtomCalculation;        // AtomCalculation, AtomType, BasisSetAccuracy, SpinGroup
import qchem.ElectronConfiguration;
import qchem.SCFCap;                 // BoundedSCF (the trace test runs a deliberately bounded SCF)
import qchem.SolidCalculation;       // SolidCalculation -- the NAMED periodic front door (Step 4)
import qchem.Structure;              // FCCUnitCell
import qchem.StructureData;          // GetMolecule("H2O") -- the Ladder hand-off gate
import qchem.Lattice_3D;             // Lattice_3D
import qchem.BasisSet;               // Real_BS
import qchem.BasisSet.Gaussian.Point.Factory;  // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.Types;                  // ivec3_t
import qchem.SCFIterator;            // SCFParams
import qchem.SCFIterator.Internal.Occupations;   // AnElectronMoved -- the cfg '*' comparison (unit tests may import Internal)
import qchem.MoleculeCalculation;            // the molecular facade (the golden molecular trace)
using namespace qchem;
using enum BasisSetAccuracy;

// V-STRUCTURE-CELL-PROBE: a solid iterator cannot be built without its cell, and a bare Structure is not one.
using SolidIterator=SCFIterator::SolidSCFIterator;
static_assert(!std::is_constructible_v<SolidIterator, const BasisSet::tBasisSet<dcmplx>*, const ElectronConfiguration*,
                                       Hamiltonian::tHamiltonian<dcmplx>*, SCFAccelerators::SCFAccelerator*>);
static_assert(!std::is_constructible_v<SolidIterator, const BasisSet::tBasisSet<dcmplx>*, const ElectronConfiguration*,
                                       Hamiltonian::tHamiltonian<dcmplx>*, SCFAccelerators::SCFAccelerator*, const Structure&>);

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
        AtomCalculation calc(2, 0, {.type = AtomType::Slater, .basisPool = Low},
            {.NMaxIter = 40, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = std::nullopt,
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

// A Hamiltonian with no real-space quadrature measures the S-metric size of ΔD and says so: no mesh exists for this column, and it is
// not the Δρ of a DFT run (user 2026-11-04: "HF should have no mesh!").
TEST(SCFTrace, HartreeFockShowsTheDensityMatrixNormNotDeltaRho)
{
    std::string trace = CapturedTrace([]{
        AtomCalculation calc(2, 0, {.type = AtomType::Slater, .basisPool = Low},   // the default model is HF
            {.NMaxIter = 40, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = std::nullopt,
             .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true});
        calc.Energy();
    });
    const std::string h = HeaderLine(trace);
    ASSERT_FALSE(h.empty()) << "no trace was printed at all:\n" << trace;
    EXPECT_TRUE(Has(h,"‖ΔD‖")) << h;
    EXPECT_FALSE(Has(h,"Δρ"))  << h;
    EXPECT_EQ(trace.find("residualMesh"), std::string::npos) << "a mesh built for the gate alone";
}

// Combo 2 of 4: MOLECULAR/ATOMIC x PP.  A pseudopotential is not Coulombic, so the virial theorem does
// NOT hold and the column must be ABSENT -- this is V1.27's landed half (IsVirialValid), pinned.  Before
// it, a molecular PP run printed a virial the tree's own docs call invalid for it.
TEST(SCFTrace, MolecularPP_HidesTheVirialColumn)
{
    std::string trace = CapturedTrace([]{
        AtomCalculation calc(14, 10, {.type = AtomType::Slater, .basisPool = Medium, .pseudopotential = true},
            {.NMaxIter = 60, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = std::nullopt,
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
        AtomCalculation calc(14, 10, {.type = AtomType::Slater, .basisPool = Medium, .pseudopotential = true},
            {.NMaxIter = 60, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-6, .MinVirial = std::nullopt,
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
        par.NMaxIter=2; par.MinΔρPerElectron=1e-3; par.MinΔE=1e-6; par.MinΔFD=std::nullopt;
        par.MinVirial=std::nullopt; par.MinFD=std::nullopt; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true;
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

// Combo 4 of 4, SOLID x NON-PP (an all-electron periodic run), has no test: no all-electron solid path exists to build one.  Its
// placeholder (an empty DISABLED_ body) was deleted 2026-11-09 under the DISABLED_ rule; its worry -- that SolidSCFIterator hard-codes
// the virial's absence -- is gone: SolidSCFIterator::AccumulateColumns calls the base layout, which adds the virial exactly when
// IsVirialValid().

// R2.7: a live [F,D] gate reading exactly 0 at convergence throws (the G3 shape: the Null accelerator's [F,D] = 0); a gate that is
// off may read 0, and a tiny but COMPUTED value is a verdict, not an absence.
TEST(SCFTrace, R2_7_ADeadGateThrows)
{
    using qchem::SCFIterator::RequireLiveGates;
    EXPECT_NO_THROW(RequireLiveGates(3e-7, 1e-5, 12));
    EXPECT_NO_THROW(RequireLiveGates(1e-300, 1e-5, 12));   // computed, just small
    EXPECT_NO_THROW(RequireLiveGates(0.0, std::nullopt, 12)); // [F,D] gate off
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
const char* const kMolecularGolden = R"G(  #               Etotal      [F,D]       ΔE/E     Δ[F,D]      ‖ΔD‖      2+V/K  ρ_mix         accel   svMin    cfg       m_site
                            (1e-06)               (1e-08)   (1e-06)    (1e+02)
----------------------------------------------------------------------------------------------------
  1    -150.860672827542   3.70e+00  -1.00e+00  -1.00e+10  2.15e-01   1.60e-01  Lin 0.50      ----    ----                 ----
  2    -149.253865566251   1.21e+00   1.08e-02  -2.49e+00  8.32e-02   5.66e-02  Lin 0.81      ----    ----                 ----
  3    -149.604239974138   2.96e-01  -2.34e-03  -9.11e-01  2.12e-02   5.55e-04  Lin 0.93      ----    ----                 ----
  4    -149.628364092134   7.86e-02  -1.61e-04  -2.17e-01  3.71e-03  -1.89e-03  Lin 0.97      DIIS:2  8.1e-02              ----
  5    -149.629387019536   2.30e-02  -6.84e-06  -5.56e-02  1.56e-03  -4.16e-03  Lin 0.99      DIIS:3  3.6e-03              ----
  6    -149.629544746148   4.15e-03  -1.05e-06  -1.88e-02  2.04e-04  -4.00e-03  Lin 1.00      DIIS:4  8.4e-05              ----
  7    -149.629546734869   1.37e-03  -1.33e-08  -2.77e-03  6.53e-05  -4.10e-03  Lin 1.00      DIIS:4  1.2e-05              ----
  8    -149.629546986082   2.67e-04  -1.68e-09  -1.11e-03  1.91e-05  -4.08e-03  Lin 1.00      DIIS:4  3.5e-07              ----
  9    -149.629547021205   8.69e-05  -2.35e-10  -1.80e-04  8.20e-06  -4.08e-03  Lin 1.00      DIIS:4  2.2e-08              ----
 10    -149.629547027597   1.02e-05  -4.27e-11  -7.67e-05  7.11e-07  -4.08e-03  Lin 1.00      DIIS:2  5.6e-09              ----
 11    -149.629547027645   1.95e-06  -3.18e-13  -8.25e-06  7.35e-08  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 12    -149.629547027645   3.08e-07  -3.42e-15  -1.64e-06  4.97e-08  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 13    -149.629547027645   2.15e-07  -2.66e-15  -9.26e-08  3.42e-08  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 14    -149.629547027645   1.55e-07   1.33e-15  -5.99e-08  2.38e-08  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 15    -149.629547027645   1.12e-07  -9.50e-16  -4.30e-08  1.67e-08  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 16    -149.629547027645   8.15e-08  -3.80e-16  -3.10e-08  1.19e-08  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 17    -149.629547027645   5.92e-08   1.90e-16  -2.23e-08  8.50e-09  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 18    -149.629547027645   4.31e-08   5.70e-16  -1.61e-08  6.12e-09  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 19    -149.629547027645   3.14e-08  -7.60e-16  -1.17e-08  4.43e-09  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
 20    -149.629547027645   2.29e-08   7.60e-16  -8.47e-09  3.22e-09  -4.08e-03  Lin 0.50      ----    1.0e+00              ----
[energy] final E_KS[D_out] = -149.6295470276  (the trace's first-order E_T: |E_T - E_KS| = 5.68e-14 Ha, second order)
----------------------------------------------------------------------------------------------------------
Energy    Breakdown  Total:    -149.629547028  Kinetic:     149.022225347  Potential:    -298.651772375
Potential Breakdown  Een  :    -410.569585387  Eee    :     100.155759296  Eex      :     -16.283520342
Virial               V/K  :    -2.00407537653  Eee/Exc:      -6.150743647  
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
[molecule] CONVERGED iters=20 charge=16.00000000000 Eelec=-177.67512108549 Etot=-149.6295470276  (Ekin=149.02222534747 Een=-410.56958538693 Eee=100.15575929635 Exc=-16.28352034239 Enn=28.04557405784)
)G";
const char* const kSolidGolden = R"G(  #               Etotal      [F,D]       ΔE/E        Δρ   ρ_lost/N    XC ρ    samp/τ  ρ_mix         accel   svMin    cfg       gap  
                                       (1e-06)   (1e-03)
----------------------------------------------------------------------------------------------------
  1      -7.116180722895   2.70e-01  -1.00e+00  5.12e-02   2.35e-10      DM      ----  PulD 0.45     ----    ----          6.69e-02  
  2      -7.115350767208   1.23e-02   1.17e-04  2.60e-02  -1.02e-09      DM      ----  PulD 0.45     ----    ----          6.83e-02  
  3      -7.115141435234   6.46e-03   2.94e-05  1.34e-02  -1.38e-09      DM      ----  PulD 0.45     ----    ----          6.90e-02  
  4      -7.115088757831   3.97e-03   7.40e-06  7.18e-03  -1.68e-09      DM      ----  PulD 0.45     ----    ----          6.93e-02  
[energy] final E_KS[D_out] = -7.115054455314  (the trace's first-order E_T: |E_T - E_KS| = 3.43e-05 Ha, second order)
----------------------------------------------------------------------------------------------------------
Energy    Breakdown  Total:  -7.11505445531  Kinetic:   4.36165986768  Potential: -11.47671432300
Potential Breakdown  Een  :  -1.13823286561  Eee    :   0.90114790101  Eex      :  -2.54427180641
Een       Breakdown  Eloc :  -2.86263085138  EenNL  :   1.72439798577
Virial               V/K  :  -2.63127219250  Eee/Exc:  -0.35418696176  
+-----------+-------------+----------+
| Occ/Degen |      ϵ (au) | Symmetry |
+-----------+-------------+----------+
|    2/2    | -0.11636597 | 2(0,0,0) |
|    6/6    |  0.32414685 | 3(0,0,0) |
+-----------+-------------+----------+
[gpw] NOT converged iters=4 charge=8.00000000000 Eelec=1.58030309668 Etot=-7.1150544553  (Ekin=4.36165986768 Een=-1.13823286561 Eee=0.90114790101 Exc=-2.54427180641 Enn=-8.40046478619 E_alphaZ=-0.29489276580)
[gpw raster pilot] SKIPPED: the SCF did not converge, so there is no D to audit -- the rasters stay the ball's own
)G";
}

//! The molecular layout (O2 HF triplet, polarized: m_site; virial; the cfg column; DIIS counts and svMin).  No cfg '*': one
//! no-symmetry block per spin always fills its lowest orbitals.  (Iterations 2 and 3 used to flag, an ↑ level crossing a ↓ level in an
//! energy-ordered label string -- OpenWork cfg row, 2026-11-03.)
TEST(SCFTrace, MolecularConsoleIsByteIdentical)
{
    const std::string got=TraceBlock(CapturedTrace([]{
        Molecule o2;
        o2.Insert(Site{8, Vector3D<double>(0, 0,  1.141)});
        o2.Insert(Site{8, Vector3D<double>(0, 0, -1.141)});
        MoleculeCalculation calc(MoleculeSetup::Resolve(o2, {.basis = "dzvp", .multiplicity = 3}),
                         {.NMaxIter = 80, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-8, .MinVirial = 1e2, .MinFD = 1e-6,
                          .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true});
        (void)calc.Result();
    }));
    EXPECT_EQ(got, std::string(kMolecularGolden));
}

//! The solid layout (Si Γ, 4 bounded iterations: ΔE/E, ρ_lost/N, XC ρ, samp/τ ("----": D mixing samples no ρ̃), the permanent gap,
//! PulD mixing, the Null accelerator).
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
        par.NMaxIter=4; par.MinΔρPerElectron=1e-3; par.MinΔE=1e-6; par.MinΔFD=std::nullopt;
        par.MinVirial=std::nullopt; par.MinFD=std::nullopt; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true;
        qchem::BoundedSCF bounded("the golden trace: 4 iterations print every solid column");
        SolidCalculation calc(SolidSetup::Resolve(lat, mol, {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}), par);
        (void)calc.Result();
    }));
    EXPECT_EQ(got, std::string(kSolidGolden));
}

// THE cfg '*' IS AN OCCUPATION CHANGE, NOT A REORDERING (OpenWork cfg row, 2026-11-03): CH₄ HF in C2v puts t2 on a₁ + b₁ + b₂, three
// exactly degenerate levels of DIFFERENT irreps, all full.  Their energy order reshuffles on noise every iteration; the occupations never
// move, so no row may flag.  (The energy-ordered label string flagged most rows here.)
TEST(SCFTrace, ADegenerateClosedShellAcrossIrrepsNeverFlagsCfg)
{
    const double a=1.18816;   // 1.089 Å
    Molecule ch4;
    ch4.Insert(Site{6, Vector3D<double>( 0, 0, 0)});
    ch4.Insert(Site{1, Vector3D<double>( a, a, a)});
    ch4.Insert(Site{1, Vector3D<double>(-a,-a, a)});
    ch4.Insert(Site{1, Vector3D<double>(-a, a,-a)});
    ch4.Insert(Site{1, Vector3D<double>( a,-a,-a)});
    std::vector<double> cfg;
    auto sink=[&cfg](const SCFIterator::IterationRow& r){ cfg.push_back(r.ValueOf("cfg")); };
    MoleculeCalculation calc(MoleculeSetup::Resolve(ch4, {.basis = "dzvp", .symmetry = true}), {.NMaxIter = 80, .Verbose = false}, sink);
    (void)calc.Result();
    ASSERT_GT(cfg.size(), 3u);
    for (size_t i=0;i<cfg.size();++i) EXPECT_EQ(cfg[i], 0.0) << "iteration " << i+1 << " flagged an occupation change";
}

// ...and a REAL occupation change still flags: NO (doublet, C2v) seeds from a cylindrically symmetric density, and its π* electron
// moves between the b₁ / b₂ blocks on iterations 2 and 3 (the aufbau's per-block counts change), so those rows flag.  Bounded: the
// change is the subject, not convergence.  (A metal does not serve: the solid fills every k-block with the same count.)
TEST(SCFTrace, AnElectronMovingBetweenIrrepsFlagsCfg)
{
    Molecule no;
    no.Insert(Site{7, Vector3D<double>(0, 0,  1.0865)});   // 1.15 Å
    no.Insert(Site{8, Vector3D<double>(0, 0, -1.0865)});
    qchem::BoundedSCF bounded("six iterations of NO: the π* electron's move between irreps is the subject");
    std::vector<double> cfg;
    auto sink=[&cfg](const SCFIterator::IterationRow& r){ cfg.push_back(r.ValueOf("cfg")); };
    MoleculeCalculation calc(MoleculeSetup::Resolve(no, {.basis = "dzvp", .multiplicity = 2, .symmetry = true}), {.NMaxIter = 6, .Verbose = false}, sink);
    (void)calc.Result();
    ASSERT_EQ(cfg.size(), 6u);
    EXPECT_EQ(cfg[0], 0.0) << "iteration 1 has no previous configuration to differ from";
    EXPECT_GT(std::count(cfg.begin(), cfg.end(), 1.0), 0) << "the π* electron moved between irreps and no row flagged it";
}

// ...but a SHARED level sitting at a half-integer is not a move (lane A, 2026-11-04): a pair of degenerate orbitals sharing one
// electron sits at ½ + ½, and noise of ~1e-13 tips each side of ½ from one iteration to the next.  Rounding each occupation (lround)
// flipped the pair between (1,0) and (0,1) and flagged on noise.  Driven on the comparison itself: no natural molecular run tried
// holds an exact ½ (smeared NO's π*↑ breaks the b₁ / b₂ pair self-consistently, kT 1e-3 to 5e-2).
TEST(SCFTrace, AHalfIntegerShareOnNoiseIsNotAMove)
{
    using SCFIterator::OccupiedConfiguration;
    const double eps=1e-13;
    const OccupiedConfiguration before{{0,{2.0,2.0,0.5+eps}}, {2,{2.0,0.5-eps}}};
    const OccupiedConfiguration after {{0,{2.0,2.0,0.5-eps}}, {2,{2.0,0.5+eps}}};
    EXPECT_FALSE(SCFIterator::AnElectronMoved(before, after)) << "a ½ + ½ share on 1e-13 noise flagged";
    // a smeared run's continuous drift below half an electron is not a move either
    const OccupiedConfiguration drift{{0,{2.0,1.99,0.7}}, {2,{2.0,0.31}}};
    EXPECT_FALSE(SCFIterator::AnElectronMoved(before, drift));
}

// ...and the same comparison still sees every real move: an electron between blocks, within a block, and into a block (or orbital)
// the other configuration does not list.
TEST(SCFTrace, AnElectronMovingIsAMoveWhereverItGoes)
{
    using SCFIterator::OccupiedConfiguration;
    const OccupiedConfiguration base{{0,{2.0,1.0}}, {2,{2.0}}};
    EXPECT_FALSE(SCFIterator::AnElectronMoved(base, base));
    EXPECT_TRUE (SCFIterator::AnElectronMoved(base, OccupiedConfiguration{{0,{2.0}},        {2,{2.0,1.0}}})) << "between blocks";
    EXPECT_TRUE (SCFIterator::AnElectronMoved(base, OccupiedConfiguration{{0,{2.0,0.0,1.0}}, {2,{2.0}}}))     << "within a block";
    EXPECT_TRUE (SCFIterator::AnElectronMoved(base, OccupiedConfiguration{{0,{2.0}},        {2,{2.0}}, {5,{1.0}}})) << "into a new block";
}

// THE ROW SINK'S EXCEPTION IS THE CALLER'S (2026-11-02): an observer that throws -- a GUI's cancel, a sink that fails to write --
// stops the run, and its exception leaves Converge (via the constructor's run, and Result) as itself: same type, same message.  No
// facade layer catches it, re-wraps it, or turns it into a failed Outcome.
namespace {
struct SinkStop : std::runtime_error { using std::runtime_error::runtime_error; };
}
TEST(SCFTrace, AnObserverExceptionLeavesConvergeUntouched)
{
    auto throwAt2=[](const SCFIterator::IterationRow& r){ if (r.iteration==2) throw SinkStop("the sink stopped at iteration 2"); };
    try
    {
        Molecule o2;
        o2.Insert(Site{8, Vector3D<double>(0, 0,  1.141)});
        o2.Insert(Site{8, Vector3D<double>(0, 0, -1.141)});
        MoleculeCalculation calc(MoleculeSetup::Resolve(o2, {.basis = "dzvp", .multiplicity = 3}), {.NMaxIter = 80, .Verbose = false}, throwAt2);
        (void)calc.Result();
        ADD_FAILURE() << "the molecular run did not stop";
    }
    catch (const SinkStop& e) { EXPECT_STREQ(e.what(), "the sink stopped at iteration 2"); }
    try
    {
        FCCUnitCell cell(10.26);
        cell.AddAtom(14, {0,0,0});
        cell.AddAtom(14, {0.25,0.25,0.25});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        auto mol = std::shared_ptr<const BasisSet::Real_BS>(
            BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                        BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        SCFParams par; par.NMaxIter=4; par.Verbose=false;
        SolidCalculation calc(SolidSetup::Resolve(lat, mol, {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}), par, throwAt2);
        (void)calc.Result();
        ADD_FAILURE() << "the solid run did not stop";
    }
    catch (const SinkStop& e) { EXPECT_STREQ(e.what(), "the sink stopped at iteration 2"); }
}

// THE MOLECULAR LADDER LADDERS (2026-11-02, the user's "Ladder is not laddering"): its tail hand-off was OFF (the json door's switchat 0),
// so a converging DIIS -- which never stalls -- kept the Ladder plain DIIS to the end (lane G: water and ethene LDA/dzvp, GDM never
// engaged).  The default now hands off at |ΔE/E| < 1e-6 (CLAUDE.md's rule): once a row reads |ΔE/E| below it, a GDM row follows.
TEST(SCFTrace, AMolecularLadderHandsOffToGDMBelowTheEnergyTail)
{
    std::vector<SCFIterator::IterationRow> rows;
    (void)CapturedTrace([&rows]{
        MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = qchem::Model::LDA}, {.type = "Ladder"}),
                         {.NMaxIter = 80, .Verbose = true},
                         [&rows](const SCFIterator::IterationRow& r){ rows.push_back(r); });
        (void)calc.Result();
    });
    auto accel=[](const SCFIterator::IterationRow& r)
    {
        for (const auto& c : r.cells) if (c.title=="accel") return c.text;
        return std::string();
    };
    size_t below=rows.size(), gdm=rows.size();
    for (size_t i=0; i<rows.size(); i++)
    {
        if (below==rows.size() && std::abs(rows[i].ValueOf("ΔE/E"))<1e-6) below=i;
        if (gdm==rows.size() && accel(rows[i]).rfind("GDM",0)==0) gdm=i;
    }
    ASSERT_LT(below, rows.size()) << "the run never reached |ΔE/E| < 1e-6";
    ASSERT_LT(gdm, rows.size()) << "no GDM row: the Ladder never handed off (" << rows.size() << " rows, |ΔE/E| < 1e-6 at row " << below+1 << ")";
    EXPECT_GT(gdm, below) << "the hand-off follows the energy tail, it does not precede it";
    for (size_t i=0; i<below; i++) EXPECT_EQ(accel(rows[i]).rfind("GDM",0), std::string::npos) << "row " << i+1 << " is GDM before the tail";
}

// ...and on a SYMMETRY-BLOCKED molecule whose electrons span several irreps (OpenWork, user 2026-11-03): CH4 C2v HF, a1 holding 3
// occupied orbitals, b1 and b2 one each, a2 none -- of ONE 10-electron reservoir.  GDM used to expect the reservoir's 5 orbitals in every
// block and declined everywhere (retreating to DIIS); then its trials refilled every block with all 10 electrons.  Now: it engages
// after the tail, never retreats, and lands on the DIIS fixed point.
TEST(SCFTrace, AMolecularLadderEngagesGDMOnASymmetryBlockedMolecule)
{
    const double a=1.18816;   // 1.089 Å
    Molecule ch4;
    ch4.Insert(Site{6, Vector3D<double>( 0, 0, 0)});
    ch4.Insert(Site{1, Vector3D<double>( a, a, a)});
    ch4.Insert(Site{1, Vector3D<double>(-a,-a, a)});
    ch4.Insert(Site{1, Vector3D<double>(-a, a,-a)});
    ch4.Insert(Site{1, Vector3D<double>( a,-a,-a)});
    auto run=[&ch4](const char* accel, std::vector<SCFIterator::IterationRow>& rows)
    {
        MoleculeCalculation calc(MoleculeSetup::Resolve(ch4, {.basis = "dzvp", .symmetry = true}, {.type = accel}),
                                 {.NMaxIter = 80, .Verbose = false}, [&rows](const SCFIterator::IterationRow& r){ rows.push_back(r); });
        auto r=calc.Result();
        EXPECT_TRUE(r) << accel << " did not converge";
        return r ? r->Energy() : 0.0;
    };
    std::vector<SCFIterator::IterationRow> ladder, diis;
    const double El=run("Ladder", ladder), Ed=run("DIIS", diis);
    auto accel=[](const SCFIterator::IterationRow& r)
    {
        for (const auto& c : r.cells) if (c.title=="accel") return c.text;
        return std::string();
    };
    size_t gdm=ladder.size();
    for (size_t i=0; i<ladder.size(); i++) if (gdm==ladder.size() && accel(ladder[i]).rfind("GDM",0)==0) gdm=i;
    ASSERT_LT(gdm, ladder.size()) << "the Ladder never handed off to GDM";
    for (size_t i=gdm; i<ladder.size(); i++) EXPECT_EQ(accel(ladder[i]).rfind("GDM",0), 0u) << "row " << i+1 << " retreated off GDM";
    for (size_t i=gdm+1; i<ladder.size(); i++)
        EXPECT_LE(ladder[i].eb.GetTotalEnergy(), ladder[i-1].eb.GetTotalEnergy()+1e-12) << "row " << i+1 << ": GDM must descend";
    EXPECT_NEAR(El, Ed, 1e-8) << "GDM and DIIS reach the same fixed point";
}

// THE OBSERVER RECEIVES THE CONSOLE'S ROW (2026-10-31): the same columns, titled as the header titles them and in its order,
// and each cell's text where the console prints it -- built from the ONE column list, so a column cannot reach one and not the
// other.  On the molecular polarized run (the frame, the variable columns and m_site).
TEST(SCFTrace, TheObserverReceivesTheConsoleRow)
{
    std::vector<SCFIterator::IterationRow> rows;
    const std::string trace=CapturedTrace([&rows]{
        Molecule o2;
        o2.Insert(Site{8, Vector3D<double>(0, 0,  1.141)});
        o2.Insert(Site{8, Vector3D<double>(0, 0, -1.141)});
        MoleculeCalculation calc(MoleculeSetup::Resolve(o2, {.basis = "dzvp", .multiplicity = 3}),
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
