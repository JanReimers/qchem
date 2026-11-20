// File: Calculation/tests/GGARoute.C  GGA stage 1 (user ruling 2026-11-20): the GGA models are named, map to their libxc ids and to the
// GTH PPs fitted with them, and EVERY XC route refuses them by policy (GGARoute) until it samples the density gradient -- never an LDA in
// their place.
#include <gtest/gtest.h>
#include <string>
#include <tuple>
#include <vector>
import qchem.MoleculeCalculation;
import qchem.AtomCalculation;
import qchem.SolidCalculation;
import qchem.Hamiltonian.Factory;   // Model, ModelToXC, LibxcIDs, ModelForGTHFitted, ConsumesDensityGradient
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData::SIPP_SR
import qchem.Structure;
import qchem.Types;

using namespace qchem;

TEST(GGARoute, TheGGAModelsAreTheirLibxcPairsAndMatchTheirGTHFits)
{
    using H=Hamiltonian::Model;
    for (const auto& [m, name, ids] : std::vector<std::tuple<H,std::string,std::vector<int>>>{
             {H::PBE, "PBE", {101,130}}, {H::PBEsol, "PBESOL", {116,133}}, {H::BLYP, "BLYP", {106,131}}})
    {
        EXPECT_TRUE(Hamiltonian::IsDFT(m));
        EXPECT_EQ(Hamiltonian::ModelFromName(Hamiltonian::Name(m)), m) << "a deck spells it back";
        const auto xc=Hamiltonian::ModelToXC(m);
        EXPECT_EQ(Hamiltonian::LibxcIDs(xc), ids);
        EXPECT_EQ(Hamiltonian::GTHFittedName(xc), name);
        EXPECT_EQ(Hamiltonian::ModelForGTHFitted(name), m) << "a GTH-" << name << " PP now names the model that matches it";
        EXPECT_TRUE(Hamiltonian::ConsumesDensityGradient(xc));
    }
    for (H m : {H::LDA, H::PADE, H::Xalpha}) EXPECT_FALSE(Hamiltonian::ConsumesDensityGradient(Hamiltonian::ModelToXC(m)));
}

TEST(GGARoute, EveryRouteRefusesAGGAUntilItSamplesTheGradient)
{
    Molecule he; he.Insert(Site{2, rvec3_t(0,0,0)});
    const auto mol=MoleculeSetup::Admit(he, {.model=Hamiltonian::Model::PBE, .molecule={.basis="dzvp"}});
    ASSERT_FALSE(mol) << "a molecule";
    EXPECT_EQ(mol.Error().policy, "GGARoute");

    RunOptions atom; atom.model=Hamiltonian::Model::BLYP;
    const auto a=AtomSetup::Admit(Atom(8), atom);
    ASSERT_FALSE(a) << "an atom";
    EXPECT_EQ(a.Error().policy, "GGARoute");
    RunOptions over; over.atom.xc=Hamiltonian::ModelToXC(Hamiltonian::Model::PBEsol);
    const auto ao=AtomSetup::Admit(Atom(8), over);
    ASSERT_FALSE(ao) << "an atom's xc override";
    EXPECT_EQ(ao.Error().policy, "GGARoute");

    const auto si=Materials::Get("Si_diamond");
    RunOptions o; o.Nelec=8; o.solid.densityEcut=20.0; o.model=Hamiltonian::Model::PBE;
    const auto s=SolidSetup::Admit(Lattice_3D(*si.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o);
    ASSERT_FALSE(s) << "a periodic cell";
    EXPECT_EQ(s.Error().policy, "GGARoute");
}
