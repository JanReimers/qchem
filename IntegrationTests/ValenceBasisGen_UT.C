// File: IntegrationTests/ValenceBasisGen_UT.C
//
// The valence-basis generator (doc/GPWPlan.md sec 1): point it at an element + its GTH pseudopotential and
// get back a well-conditioned, GPW-usable valence Gaussian .bsd, validated by an atomic pseudo-atom SCF.
// These tests exercise F (from the F- closed shell -- the ionic state in NaF) and Na (neutral 3s^1), and
// PRINT the emitted .bsd so it can be captured into BasisSetData/.  The energies are "did-it-move" anchors
// (near the radial oracle, but oracle-matching is explicitly NOT the objective -- see doc/GPWPlan.md).
#include "gtest/gtest.h"
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include "nlohmann/json.hpp"      // read report::GlobalReport() values (get<double>) across the module boundary

import qchem.ValenceBasisGen;
import qchem.Structure;                  // Molecule, Atom
import qchem.BasisSet;                   // Real_BS
import qchem.BasisSet.Molecule.Factory;  // Molecule::Factory, BasisSetData
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
    EXPECT_LT(g.energy, -19.0);       // bound F- valence (oracle ~ -20.93); a sane, converged-enough basis
    EXPECT_GT(g.energy, -22.0);       // not a variational collapse
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
    EXPECT_LT(g.energy, -0.13);       // bound Na 3s^1 (oracle ~ -0.1446)
    EXPECT_GT(g.energy, -0.16);
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
    using namespace qchem::BasisSet::Molecule;
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
