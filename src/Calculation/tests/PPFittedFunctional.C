// File: Calculation/tests/PPFittedFunctional.C  PADE as a stated functional, and the PP / functional match policy (lane H, 2026-11-12).
//
// USER RULING 2026-11-12 (OpenWork PADE row, option a): OFFER PADE (libxc LDA_XC_TETER93, CP2K's XC_FUNCTIONAL PADE) as a functional a
// run states, and WARN AND RUN when a pseudopotential was fitted with another functional than the run's, naming both -- never an
// automatic switch.  LDA (Dirac + VWN5) stays the default.  Here: the names round-trip, the policy fires for a GTH-PADE PP under LDA and
// is silent under PADE (the helper, and a real pseudo-atom through the facade), and PADE is a different functional, not a relabelled LDA.
#include "gtest/gtest.h"
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

import qchem.AtomCalculation;
import qchem.MoleculeSetup;         // MoleculeSetup::PPFunctionalChecks (the pre-run query)
import qchem.SolidCalculation;      // SolidSetup::PPFunctionalChecks
import qchem.StructureData;         // GetMolecule
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData::SIPP_SR
import qchem.Types;
import qchem.Hamiltonian.Factory;   // Model, ModelToXC, Name, ModelFromName, Describe, GTHFittedName
import qchem.Policy;                // policy::PPFittedFunctional
import qchem.Reporting;
import qchem.SCFParams;
using namespace qchem;
using Hamiltonian::Model;

namespace
{
std::vector<report::json> PolicyRows(const report::json& j, const std::string& name)
{
    std::vector<report::json> out;
    std::function<void(const report::json&)> walk=[&](const report::json& x)
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
}
const report::json& FiledRun(const std::string& prefix)
{
    const auto& G=report::GlobalReport();
    const report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    if (!last) throw std::logic_error("no filed run "+prefix);
    return *last;
}
SCFParams Tight()
{
    return {.NMaxIter = 60, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-7, .MinVirial = std::nullopt, .MinFD = 1e-7,
            .StartingRelaxRo = 0.5, .MergeTol = 1e-7};
}
//! The Si pseudo-atom (GTH-PADE q4) under \a model; files its run as \a tag and returns its energy.
double SiPseudoAtom(Model model, const std::string& tag)
{
    report::Run run(tag);
    AtomCalculation calc(14, 10, {.type=AtomType::Slater, .basisPool=BasisSetAccuracy::Medium, .model=model, .pseudopotential=true}, Tight());
    EXPECT_TRUE(calc.IsConverged()) << tag;
    return calc.Energy();
}
} //anon

TEST(PPFittedFunctional, TheModelNamesRoundTripAndPadeIsLibxcTeter93)
{
    for (Model m : {Model::E1, Model::HF, Model::DE1, Model::DHF, Model::Xalpha, Model::LDA, Model::PADE})
    {
        const auto back=Hamiltonian::ModelFromName(Hamiltonian::Name(m));
        ASSERT_TRUE(back.has_value()) << Hamiltonian::Name(m);
        EXPECT_EQ(*back, m);
    }
    EXPECT_FALSE(Hamiltonian::ModelFromName("pade").has_value()) << "the spelling is exact";
    const auto pade=Hamiltonian::ModelToXC(Model::PADE);
    EXPECT_EQ(pade.kind, Hamiltonian::XC::LibXC);
    EXPECT_EQ(pade.libxcId, 20) << "libxc LDA_XC_TETER93";
    EXPECT_EQ(Hamiltonian::GTHFittedName(pade), "PADE");
    EXPECT_EQ(Hamiltonian::GTHFittedName(Hamiltonian::ModelToXC(Model::LDA)), "") << "no GTH table entry is fitted with Dirac + VWN5";
    EXPECT_THROW(Hamiltonian::ModelToXC(Model::HF), std::invalid_argument) << "a non-DFT model has no functional";
}

// The comparison is on LIBXC IDS, never on names across schemes (user 2026-11-12): the CP2K name only labels the fit.
TEST(PPFittedFunctional, TheFitAndTheRunCompareAsLibxcIDs)
{
    EXPECT_EQ(Hamiltonian::LibxcIDs(Hamiltonian::ModelToXC(Model::LDA)),  (std::vector<int>{1, 7}))  << "Dirac (LDA_X) + VWN5 (LDA_C_VWN)";
    EXPECT_EQ(Hamiltonian::LibxcIDs(Hamiltonian::ModelToXC(Model::PADE)), (std::vector<int>{20}))    << "LDA_XC_TETER93";
    EXPECT_EQ(Hamiltonian::LibxcIDsOfFitted("PADE"), (std::vector<int>{20}));
    EXPECT_EQ(Hamiltonian::LibxcIDsOfFitted("PBE"),  (std::vector<int>{101, 130}));
    EXPECT_EQ(Hamiltonian::LibxcIDsOfFitted("BP"),   (std::vector<int>{106, 132})) << "CP2K's BP is B88 + P86";
    EXPECT_THROW(Hamiltonian::LibxcIDsOfFitted("pade"), std::invalid_argument) << "an unmapped name is refused by name, never a silent match";
    EXPECT_THROW(Hamiltonian::PPFunctionalChecks({{"X","NOSUCHXC"}}, Hamiltonian::ModelToXC(Model::LDA)), std::invalid_argument);
    // Every fitted name the PP table carries is mapped (else a run with that PP would throw).
    for (const char* f : {"PADE","PBE","BLYP","BP","OLYP","PBESOL","HCTH120","HCTH407"})
        EXPECT_NO_THROW(Hamiltonian::LibxcIDsOfFitted(f)) << f;
    EXPECT_EQ(Hamiltonian::ModelForGTHFitted("PADE"), Model::PADE);
    EXPECT_FALSE(Hamiltonian::ModelForGTHFitted("PBE").has_value()) << "no model offered matches PBE until GGA lands";
}

TEST(PPFittedFunctional, ThePolicyNamesBothAndSaysWhatToDo)
{
    const auto lda=Hamiltonian::ModelToXC(Model::LDA);
    const auto cc=Hamiltonian::PPFunctionalChecks({{"Li","PADE"},{"O","PADE"}}, lda);
    ASSERT_EQ(cc.size(), 2u) << "every element is checked";
    for (const auto& c : cc)
    {
        EXPECT_FALSE(c.match) << c.element;
        EXPECT_EQ(c.fittedLabel, "PADE");
        EXPECT_EQ(c.fittedIDs, (std::vector<int>{20}));
        EXPECT_EQ(c.runIDs, (std::vector<int>{1, 7}));
        EXPECT_EQ(c.matchingModel, "PADE") << "a GTH-PADE fit is matched by stating model PADE";
    }
    {
        report::Run run("ppfitted helper differ");
        policy::PPFittedFunctional(cc);
    }
    const auto rows=PolicyRows(FiledRun("ppfitted helper differ@"), "PPFittedFunctional");
    ASSERT_EQ(rows.size(), 1u) << "ONE row naming every mismatched element";
    EXPECT_EQ(rows[0]["verdict"], "warn") << "warn and run, never a refusal or a switch";
    const std::string decided=rows[0]["decided"], reason=rows[0]["reason"];
    EXPECT_NE(decided.find("Dirac + VWN5"), std::string::npos) << decided;
    EXPECT_NE(decided.find("Li (fitted with PADE)"), std::string::npos) << decided;
    EXPECT_NE(reason.find("state model PADE"), std::string::npos) << "it says what to do: " << reason;
    EXPECT_EQ(rows[0]["detail"]["O"]["fitted"], "PADE");
    EXPECT_EQ(rows[0]["detail"]["O"]["fittedLibxcIDs"][0], 20);
    EXPECT_EQ(rows[0]["detail"]["O"]["matchingModel"], "PADE");
}

// Any mismatch, not only VWN5 vs PADE: a fit no offered model matches (a GTH-PBE PP under PADE) still warns, and its advice is the PPs.
TEST(PPFittedFunctional, AFitNoModelMatchesWarnsAndPointsAtThePPs)
{
    const auto cc=Hamiltonian::PPFunctionalChecks({{"Si","PBE"},{"O","PADE"}}, Hamiltonian::ModelToXC(Model::PADE));
    ASSERT_EQ(cc.size(), 2u);
    EXPECT_FALSE(cc[0].match);
    EXPECT_EQ(cc[0].matchingModel, "") << "no model offered matches PBE (until GGA lands)";
    EXPECT_TRUE(cc[1].match) << "O's PADE fit matches a PADE run";
    {
        report::Run run("ppfitted helper pbe");
        policy::PPFittedFunctional(cc);
    }
    const auto rows=PolicyRows(FiledRun("ppfitted helper pbe@"), "PPFittedFunctional");
    ASSERT_EQ(rows.size(), 1u);
    const std::string reason=rows[0]["reason"];
    EXPECT_NE(reason.find("use pseudopotentials fitted for PADE"), std::string::npos) << reason;
    EXPECT_FALSE(rows[0]["detail"].contains("O")) << "a matching element is not in the warning";
}

TEST(PPFittedFunctional, ThePolicyIsSilentWhenTheyMatch)
{
    const auto cc=Hamiltonian::PPFunctionalChecks({{"Li","PADE"},{"O","PADE"}}, Hamiltonian::ModelToXC(Model::PADE));
    for (const auto& c : cc) EXPECT_TRUE(c.match) << c.element;
    {
        report::Run run("ppfitted helper match");
        policy::PPFittedFunctional(cc);
    }
    EXPECT_TRUE(PolicyRows(FiledRun("ppfitted helper match@"), "PPFittedFunctional").empty());
}

// THE PRE-RUN QUERY (user 2026-11-12: "track this, warn and guide, mostly in the GUI"): a setup answers every element's check without
// running, on the molecular and the crystal path alike, so the GUI never computes it.
TEST(PPFittedFunctional, ASetupAnswersTheCheckBeforeTheRun)
{
    using StructureData::GetMolecule;
    const auto mLDA =MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA,  .pp="gth-lda-valence"});
    const auto mPADE=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::PADE, .pp="gth-lda-valence"});
    const auto mAE  =MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    const auto cc=mLDA->PPFunctionalChecks();
    ASSERT_EQ(cc.size(), 2u) << "H and O, each once";
    for (const auto& c : cc) {EXPECT_FALSE(c.match) << c.element; EXPECT_EQ(c.matchingModel, "PADE") << c.element;}
    for (const auto& c : mPADE->PPFunctionalChecks()) EXPECT_TRUE(c.match) << "GTH-PADE under PADE matches: " << c.element;
    EXPECT_TRUE(mAE->PPFunctionalChecks().empty()) << "an all-electron run has no PP to check";

    SolidCalcOptions o;
    const auto sLDA=SolidSetup::Resolve(Lattice_3D(UnitCell(*Materials::Get("Si_diamond").cell), ivec3_t(1,1,1)),
                                        GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o);
    const auto sc=sLDA->PPFunctionalChecks();
    ASSERT_EQ(sc.size(), 1u) << "Si, once";
    EXPECT_EQ(sc[0].element, "Si");
    EXPECT_FALSE(sc[0].match);
    EXPECT_EQ(sc[0].matchingModel, "PADE");
}

// Through the facade: a GTH-PADE pseudo-atom files the row under LDA and not under PADE, and the two functionals give different
// energies (PADE is not LDA relabelled; on Li q3 the two differ by 0.6 mHa, all of it E_xc).
TEST(PPFittedFunctional, APseudoAtomWarnsUnderLDAAndIsSilentUnderPADE)
{
    const double eLDA =SiPseudoAtom(Model::LDA,  "ppfitted si lda");
    const double ePADE=SiPseudoAtom(Model::PADE, "ppfitted si pade");
    const auto warned=PolicyRows(FiledRun("ppfitted si lda@"), "PPFittedFunctional");
    ASSERT_EQ(warned.size(), 1u);
    EXPECT_EQ(warned[0]["detail"]["Si"]["fitted"], "PADE");
    EXPECT_TRUE(PolicyRows(FiledRun("ppfitted si pade@"), "PPFittedFunctional").empty()) << "a GTH-PADE PP under PADE matches";
    std::printf("  [Si q4 pseudo-atom] LDA %.8f  PADE %.8f  diff %.3e Ha\n", eLDA, ePADE, ePADE-eLDA);
    EXPECT_GT(std::fabs(ePADE-eLDA), 1e-5) << "PADE and Dirac + VWN5 are different functionals";
    EXPECT_LT(std::fabs(ePADE-eLDA), 5e-3) << "...but both are the LDA: mHa apart, not Ha";
}
