// File: src/Structure/tests/StructureDataUT.C  qchem.StructureData: the pre-defined molecules and cells
// (D-STRUCTDATA).  Structure only -- valence counts are the Calculation library's (src/Calculation/tests/Materials.C).
#include "gtest/gtest.h"
#include <stdexcept>
#include <string>
#include <vector>

import qchem.StructureData;
import qchem.Types;

using namespace qchem;
namespace SD=qchem::StructureData;

TEST(StructureData, WaterIsTheExperimentalGeometryInBohr)
{
    const Molecule w=SD::GetMolecule("H2O");
    ASSERT_EQ(w.GetNumAtoms(), 3u);
    EXPECT_EQ(SD::MoleculeNames().size(), 4u);
    EXPECT_EQ(SD::KindOf("H2O"), SD::Kind::Molecule);
}

TEST(StructureData, CellsLoadWithTheirSpinDecoration)
{
    const UnitCell mno=SD::GetCell("MnO_AFM2");
    EXPECT_EQ(mno.GetNumAtoms(), 4u);
    std::vector<int> Z; std::vector<bool> flip;
    mno.ForEachSite([&](int z, const rvec3_t&, bool f){ Z.push_back(z); flip.push_back(f); });
    EXPECT_EQ(Z,    (std::vector<int>{25,25,8,8}));
    EXPECT_EQ(flip, (std::vector<bool>{false,true,false,false}));
    for (const std::string& n : SD::CellNames()) EXPECT_EQ(SD::KindOf(n), SD::Kind::Cell) << n << " is in the cell pick-list, so it must be a cell";
    EXPECT_FALSE(SD::CellNames().empty());
    // The lattice-constant override scales the cell: a=2x => 8x volume.
    EXPECT_NEAR(SD::GetCell("Si_diamond", 2*10.26).GetCellVolume(), 8*SD::GetCell("Si_diamond").GetCellVolume(), 1e-9);
}

TEST(StructureData, BoxesAreTheBankedNamedBoxes)
{
    const UnitCell box=SD::AtomInBox("Mn", 16.0), file=SD::GetCell("Mn_box16");
    EXPECT_NEAR(box.GetCellVolume(), file.GetCellVolume(), 1e-9);
    EXPECT_EQ(SD::DimerInBox("Na", 16.0, 5.8, true).GetNumAtoms(), 2u);
}

TEST(StructureData, AKindMismatchOrAMissThrows)
{
    EXPECT_EQ(SD::KindOf("Si_diamond"), SD::Kind::Cell);
    EXPECT_THROW(SD::GetCell("H2O"),            std::runtime_error) << "a molecule is not a cell";
    EXPECT_THROW(SD::GetMolecule("Si_diamond"), std::runtime_error) << "a cell is not a molecule";
    EXPECT_THROW(SD::GetMolecule("C60"),        std::runtime_error);
    EXPECT_THROW(SD::KindOf("_doc"),            std::runtime_error) << "the documentation key is not an entry";
}
