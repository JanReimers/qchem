// File: IntegrationTests/ValenceBasisGen_UT.C
//
// The valence-basis generator (doc/GPWPlan.md sec 1): point it at an element + its GTH pseudopotential and
// get back a well-conditioned, GPW-usable valence Gaussian .bsd, validated by an atomic pseudo-atom SCF.
// These tests exercise F (from the F- closed shell -- the ionic state in NaF) and Na (neutral 3s^1), and
// PRINT the emitted .bsd so it can be captured into BasisSetData/.  The energies are "did-it-move" anchors
// (near the radial oracle, but oracle-matching is explicitly NOT the objective -- see doc/GPWPlan.md).
#include "gtest/gtest.h"
#include <map>
#include <cmath>
#include <iostream>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include "nlohmann/json.hpp"      // read report::GlobalReport() values (get<double>) across the module boundary

import qchem.ValenceBasisGen;
import qchem.AtomCalculation;            // the pseudo-atom vs CP2K ATOM gate
import qchem.SCFParams;
import qchem.Symmetry.Atom.Spherical;     // AtomicSymmetry::Getl
import qchem.Symmetry.Irrep;
import qchem.Structure;                  // Molecule, Atom
import qchem.BasisSet;                   // Real_BS
import qchem.BasisSet.Gaussian.Point.Factory;  // Gaussian::Factory, BasisSetData
import qchem.Types;
import qchem.Reporting;                  // report::GlobalReport -- inspect the recorded basis.usage

using namespace qchem;

TEST(ValenceBasisGen, Fluorine_q7)
{
    // F- (2s^2 2p^6, the ionic state in NaF): validate against 8 valence electrons over an s+p window.
    // Emit s over the full window; p over a DISJOINT, less-tight window (a p need not reach exp 40, and the
    // molecular reader would merge any shared exponent -- flagged bug).
    ValenceBasisRecipe r;
    r.element        = "F";
    r.Zion           = 7;
    r.electrons      = 8;                                   // F-
    r.shells         = { {0, EvenTemperedWindow(8, 0.12, 40.0)},
                         {1, EvenTemperedWindow(6, 0.14, 12.0)} };
    GeneratedBasis g = GenerateValenceBasis(r);
    std::cout << "[gen F] E=" << g.energy << " conv=" << g.converged << "\n" << g.block << std::endl;
    // ORACLE-PINNED (2026-08-06, after the KB radial-assembly fix): CP2K's numerically-exact ATOM code
    // gives -24.213540 for this F- pseudo-ion (deck ~/Code/cp2k-runs/fminus_q7.inp); we land 29 mHa above
    // it in this 8s+6p window (basis incompleteness).  The OLD bounds (-19 .. -22, "oracle ~ -20.93") were
    // set against the BROKEN 3-D-mesh KB route -- see A_PP.PerLKleinmanBylanderOracle.
    EXPECT_NEAR(g.energy, -24.2135, 0.10);   // vs the CP2K pseudo-ion oracle (window incompleteness)
    EXPECT_GT (g.energy, -24.3);             // not a variational collapse (below the complete-basis oracle)
}

TEST(ValenceBasisGen, Sodium_q1)
{
    // Na neutral (3s^1): diffuse s window (validated) + one diffuse p polarization shell (bonding in NaF).
    ValenceBasisRecipe r;
    r.element        = "Na";
    r.Zion           = 1;
    r.electrons      = 1;
    r.shells         = { {0, EvenTemperedWindow(5, 0.03, 2.0)},
                         {1, EvenTemperedWindow(2, 0.05, 0.3)} };   // p polarization
    GeneratedBasis g = GenerateValenceBasis(r);
    std::cout << "[gen Na] E=" << g.energy << " conv=" << g.converged << "\n" << g.block << std::endl;
    // ORACLE-PINNED (2026-08-06, after the KB radial-assembly fix): CP2K ATOM gives -0.184065 (deck
    // ~/Code/cp2k-runs/na_atom_q1.inp) and we land 100 uHa away -- Na q1 carries BOTH an l=0 and an l=1
    // projector, so it was doubly wrong before (l=1 deleted, l=0 leaking).  Old bounds assumed -0.1446.
    EXPECT_NEAR(g.energy, -0.184065, 5e-3);  // vs the CP2K pseudo-atom oracle
    EXPECT_GT (g.energy, -0.20);             // not a variational collapse
}

// V3.1/V3.2 (doc/CleanupCandidates.md): the fully-polarized ONE-electron pseudo-atom (Na q1 doublet,
// nUp=1 nDown=0) and the riding-along unoccupied p window.  The library Na entry was hand-constructed to
// dodge exactly this, so these are the regression anchors for the empty-minority-channel atom SCF.
TEST(ValenceBasisGen, SodiumSeedDensityUnpolarizedWithPolarizationShell)
{
    ValenceBasisRecipe r;
    r.element="Na"; r.Zion=1; r.electrons=1;
    r.shells={ {0, EvenTemperedWindow(5, 0.03, 2.0)}, {1, EvenTemperedWindow(2, 0.09, 0.3)} };
    GeneratedSeedDensity g = GenerateSeedDensity(r);
    EXPECT_NEAR(g.charge, 1.0, 0.01);
    EXPECT_EQ(g.moment, 0.0);
}

TEST(ValenceBasisGen, SodiumSeedDensitySpinResolved)             // V3.1: empty minority channel
{
    ValenceBasisRecipe r;
    r.element="Na"; r.Zion=1; r.electrons=1;
    r.shells={ {0, EvenTemperedWindow(5, 0.03, 2.0)} };
    r.spinResolved=true;
    GeneratedSeedDensity g = GenerateSeedDensity(r);
    EXPECT_NEAR(g.charge, 1.0, 0.01);
    EXPECT_NEAR(g.moment, 1.0, 0.05);   // one unpaired electron, all of it in the up channel
}

TEST(ValenceBasisGen, SodiumSeedDensitySpinResolvedWithPolarizationShell)   // V3.2 on the polarized path
{
    ValenceBasisRecipe r;
    r.element="Na"; r.Zion=1; r.electrons=1;
    r.shells={ {0, EvenTemperedWindow(5, 0.03, 2.0)}, {1, EvenTemperedWindow(2, 0.09, 0.3)} };
    r.spinResolved=true;
    GeneratedSeedDensity g = GenerateSeedDensity(r);
    EXPECT_NEAR(g.charge, 1.0, 0.01);
    EXPECT_NEAR(g.moment, 1.0, 0.05);
}

// SEED DENSITY generation (the offline library for IonicSAD): the SAME pseudo-atom SCF that makes the basis
// also emits a spherical rho(r) for the seed-density library.  THE POINT: an anion (F-) valence density is
// spatially DIFFUSE -- its <r> exceeds the neutral atom's -- which is exactly why a proper F- seed converges
// where the old neutral-density-scaled-x8/7 IonicSAD (too compact) did not (PlaneWaveDFTUT / doc/GPWPlan §0).
// This asserts that physics (charge conserved, F- more diffuse than neutral F) and PRINTS the F- library entry
// so it can be captured into atomic_valence_densities.json.
TEST(ValenceBasisGen, FluorineSeedDensityAnionIsDiffuse)
{
    auto window = []{ return std::vector<std::pair<int,std::vector<double>>>{
        {0, EvenTemperedWindow(8, 0.12, 40.0)}, {1, EvenTemperedWindow(6, 0.14, 12.0)} }; };
    ValenceBasisRecipe neutral; neutral.element="F"; neutral.Zion=7; neutral.electrons=7; neutral.shells=window();
    ValenceBasisRecipe anion;   anion.element  ="F"; anion.Zion  =7; anion.electrons  =8; anion.shells  =window();

    GeneratedSeedDensity n = GenerateSeedDensity(neutral);
    GeneratedSeedDensity a = GenerateSeedDensity(anion);
    std::cout << "[seed F ] neutral: charge="<<n.charge<<" <r>="<<n.meanR<<" conv="<<n.converged<<"\n"
              << "[seed F-] anion:   charge="<<a.charge<<" <r>="<<a.meanR<<" conv="<<a.converged<<std::endl;
    std::cout << "[F- seed entry] " << a.jsonEntry.substr(0, 180) << " ...rho[400]... }" << std::endl;

    EXPECT_NEAR(n.charge, 7.0, 0.1);          // neutral F: 7 valence e-
    EXPECT_NEAR(a.charge, 8.0, 0.1);          // F-: 8 valence e- (charge conserved by construction)
    EXPECT_GT(a.meanR, n.meanR);              // THE POINT: the anion density is more diffuse than the neutral
}

// Assemble the full valence_lowq.bsd (organised by TYPE, all elements in one file, per the BasisSetData
// convention).  Prints the file so it can be captured into BasisSetData/valence_lowq.bsd.  Grows one block
// per element (F, Na today; Si, Cs, I to follow).
TEST(ValenceBasisGen, AssembleValenceLowqFile)
{
    ValenceBasisRecipe f;
    f.element="F"; f.Zion=7; f.electrons=8;
    f.shells={ {0, EvenTemperedWindow(8, 0.12, 40.0)}, {1, EvenTemperedWindow(6, 0.14, 12.0)} };
    ValenceBasisRecipe na;
    na.element="Na"; na.Zion=1; na.electrons=1;
    na.shells={ {0, EvenTemperedWindow(5, 0.03, 2.0)}, {1, EvenTemperedWindow(2, 0.05, 0.3)} };
    ValenceBasisRecipe al;                          // Al metal (3s^2 3p^1), generated via CLIapps/valgen
    al.element="Al"; al.Zion=3; al.electrons=3;     // window tuned via the basis-usage heat map (0.75 mHa above the PP floor)
    al.shells={ {0, EvenTemperedWindow(6, 0.04, 4.0)}, {1, EvenTemperedWindow(5, 0.05, 2.5)} };

    std::vector<std::string> blocks = { GenerateValenceBasis(f).block, GenerateValenceBasis(na).block,
                                        GenerateValenceBasis(al).block };
    const std::string file = AssembleBasisFile(
        "Low-q GTH-pseudopotential valence basis, generated from atomic pseudo-atom SCFs (doc/GPWPlan sec 1)",
        blocks);
    std::cout << "===== BEGIN valence_lowq.bsd =====\n" << file << "===== END valence_lowq.bsd =====" << std::endl;
    EXPECT_NE(file.find(" F   0"),  std::string::npos);
    EXPECT_NE(file.find(" NA   0"), std::string::npos);
    EXPECT_NE(file.find(" AL   0"), std::string::npos);
}

// The basis-USAGE heat map (doc/GPWPlan1 §1): the pseudo-atom SCF records per-function occupation-weighted
// Mulliken populations P_i=(DS)_ii into the run report's basis.usage.  Their SUM over all functions is
// Tr(DS)=the valence electron count (the S-orthonormal orbitals make each occupied orbital contribute its
// occupation) -- the correctness anchor for the whole heat-map pipeline (Orbitals -> IrrepWF -> SCFIterator
// -> Reporting).  Al q3 (3s^2 3p^1) => 3 valence electrons.
TEST(ValenceBasisGen, BasisUsageSumsToElectronCount)
{
    namespace rpt = qchem::report;
    rpt::ClearGlobal();
    ValenceBasisRecipe r;
    r.element="Al"; r.Zion=3; r.electrons=3;
    r.shells={ {0, EvenTemperedWindow(6, 0.04, 4.0)}, {1, EvenTemperedWindow(5, 0.05, 2.5)} };
    GenerateValenceBasis(r);            // opens+closes an AtomCalculation run that records basis.usage

    double sum=0.0; size_t nrows=0; bool found=false;
    for (auto& [key, run] : rpt::GlobalReport().items())
        if (run.contains("basis") && run["basis"].contains("usage"))
        {
            found=true;
            for (auto& row : run["basis"]["usage"]) { sum += row["pop"].template get<double>(); nrows++; }
        }
    std::cout << "[basis.usage] found=" << found << " rows=" << nrows << " sum(pop)=" << sum << std::endl;
    EXPECT_TRUE(found);
    EXPECT_EQ(nrows, 11u);              // 6 s + 5 p functions
    EXPECT_NEAR(sum, 3.0, 1e-6);       // Sum_i (DS)_ii = Tr(DS) = N valence electrons
    rpt::ClearGlobal();
}

// End-to-end: the COMMITTED BasisSetData/valence_lowq.bsd parses through the molecular factory and yields the
// expected function counts (F: 8 s + 8 p = 32 Cartesian; Na: 5 s + 2 p = 11).  This closes the loop
// generator -> file -> loader, and guards the committed .bsd against drift from the generator above.
TEST(ValenceBasisGen, ValenceLowqFileLoads)
{
    using namespace qchem::BasisSet::Gaussian;
    auto nfun=[](int Z, BasisSetData d){
        Molecule m; m.Insert(new Atom(Z, 0.0, {0,0,0}));
        std::unique_ptr<qchem::BasisSet::Real_BS> bs(Factory(d, &m, Engine::MnD, Angular::Cartesian));
        return bs->GetNumFunctions();
    };
    std::cout << "[valence_lowq N] F="<<nfun(9,BasisSetData::VALENCE_LOWQ)
              << " Na="<<nfun(11,BasisSetData::VALENCE_LOWQ)
              << " Al="<<nfun(13,BasisSetData::VALENCE_LOWQ)
              << "   (calib: Si sipp="<<nfun(14,BasisSetData::SIPP)
              << " sipp_sr="<<nfun(14,BasisSetData::SIPP_SR)<<")" << std::endl;
    EXPECT_EQ(nfun(13, BasisSetData::VALENCE_LOWQ), 21u);   // Al (valgen): 6 s + 5 p = 6 + 15 Cartesian
    Molecule naf;
    naf.Insert(new Atom(9,  0.0, {0,0,0}));       // F
    naf.Insert(new Atom(11, 0.0, {3,0,0}));       // Na
    std::unique_ptr<qchem::BasisSet::Real_BS> bs(
        Factory(BasisSetData::VALENCE_LOWQ, &naf, Engine::MnD, Angular::Cartesian));
    std::cout << "[valence_lowq N] NaF total = " << bs->GetNumFunctions() << std::endl;
    EXPECT_EQ(bs->GetNumFunctions(), 26u + 11u);   // F(8s+6p=26) + Na(5s+2p=11), Cartesian (disjoint exponents)

}

// THE PSEUDO-ATOM AGAINST CP2K's ATOM CODE (DFT+U increment 3, 2026-09-21).  The atomic +U radial (and the
// UPF for the hp.x oracle) is our GTH pseudo-atom's own l orbital, so its spectrum must be the pseudo-atom's.
// Oracle: `cp2k.psmp` PROGRAM_NAME ATOM, Mn GTH-PADE-q7, CORE [Ar] 4s2 3d5, PADE LDA, 30 geometrical GTOs
// per l (scratch input `cp2katom/mn.inp`, 2026-09-21): E = -14.241357 Ha, eps(4s) = -0.194010 Ha,
// eps(3d) = -0.257385 Ha.  In a LARGE pool ours reproduces that to 3 / 0.3 / 0.1 mHa (the next test: the
// pseudopotential is right).  THIS test pins what the VA basis's OWN 7 s + 7 d exponents do to the FREE atom:
// the span was trimmed of the diffuse 0.18 d shell for the SOLID, so it cannot hold the free 3d -- E 46 mHa
// high, eps(4s) 34 mHa and eps(3d) 77 mHa too shallow (a too-compact 3d).  That is why the atomic +U radial
// is the complete-pool 3d PROJECTED onto the site's shells, never the pseudo-atom run in the shells themselves.
TEST(ValenceBasisGen, MnQ7PseudoAtomInTheVAExponentsShowsTheDiffuseTrim)
{
    AtomCalcOptions o;
    o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=7;
    o.exponentsByL={{0,{0.1,0.24928829,0.62144650,1.54919334,3.86195754,9.62740780,24.0}},
                    {2,{0.38369936,0.81791778,1.74352515,3.71660826,7.92255676,16.88822203,36.0}}};
    SCFParams p; p.MinVirial=1e30;
    AtomCalculation atom(25, 18, o, p);
    ASSERT_TRUE(atom.IsConverged());
    EXPECT_NEAR(atom.Energy(), -14.1949, 0.005) << "46 mHa above CP2K ATOM: the VA span's incompleteness for the free atom";
    std::map<int,double> eps;                                   // lowest occupied eigenvalue per l
    for (const Irrep& ir : atom.GetIrreps(Spin::None))
    {
        const auto* as=dynamic_cast<const Symmetry::Atom::AtomicSymmetry*>(ir.sym.get());
        ASSERT_TRUE(as);
        const auto* os=atom.Orbitals(ir);
        for (const auto* orb : os->Iterate())
            if (orb->IsOccupied()) { const int l=int(as->Getl()); if (!eps.count(l) || orb->GetEigenEnergy()<eps[l]) eps[l]=orb->GetEigenEnergy(); }
    }
    ASSERT_TRUE(eps.count(0) && eps.count(2));
    std::cout << "[Mn q7 pseudo-atom] E=" << atom.Energy() << " eps(4s)=" << eps[0] << " eps(3d)=" << eps[2]
              << "   (CP2K ATOM: -14.241357, -0.194010, -0.257385)" << std::endl;
    EXPECT_NEAR(eps[0], -0.1596, 0.005) << "4s: 34 mHa shallower than CP2K's -0.194";
    EXPECT_NEAR(eps[2], -0.1802, 0.005) << "3d: 77 mHa shallower than CP2K's -0.257 -- the trimmed 0.18 d shell";
}

// The same pseudo-atom in a LARGE even-tempered pool (16 s + 16 d, 0.05..200): separates basis
// incompleteness from the pseudopotential itself.  If eps(3d) stays 77 mHa above CP2K's here, the GTH d
// channel differs between the codes -- and that is a candidate for the 100 mHa MnO offset.
TEST(ValenceBasisGen, MnQ7PseudoAtomInALargePool)
{
    auto pool=[](int n, double emin, double emax){ std::vector<double> e; for (int i=0;i<n;i++) e.push_back(emin*std::pow(emax/emin, double(i)/(n-1))); return e; };
    AtomCalcOptions o;
    o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=7;
    o.exponentsByL={{0,pool(16,0.05,200.0)},{2,pool(16,0.08,200.0)}};
    SCFParams p; p.MinVirial=1e30;
    AtomCalculation atom(25, 18, o, p);
    ASSERT_TRUE(atom.IsConverged());
    std::map<int,double> eps;
    for (const Irrep& ir : atom.GetIrreps(Spin::None))
    {
        const auto* as=dynamic_cast<const Symmetry::Atom::AtomicSymmetry*>(ir.sym.get());
        for (const auto* orb : atom.Orbitals(ir)->Iterate())
            if (orb->IsOccupied()) { const int l=int(as->Getl()); if (!eps.count(l) || orb->GetEigenEnergy()<eps[l]) eps[l]=orb->GetEigenEnergy(); }
    }
    std::cout << "[Mn q7 pseudo-atom, 16+16 pool] E=" << atom.Energy() << " eps(4s)=" << eps[0] << " eps(3d)=" << eps[2]
              << "   (CP2K ATOM: -14.241357, -0.194010, -0.257385)" << std::endl;
    EXPECT_NEAR(atom.Energy(), -14.241357, 0.02);
    EXPECT_NEAR(eps[0], -0.194010, 0.01) << "4s";
    EXPECT_NEAR(eps[2], -0.257385, 0.01) << "3d";
}
