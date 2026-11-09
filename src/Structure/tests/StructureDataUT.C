// File: src/Structure/tests/StructureDataUT.C  qchem.StructureData: the pre-defined molecules and cells
// (D-STRUCTDATA).  Structure only -- valence counts are the MoleculeCalculation library's (src/Calculation/tests/Materials.C).
#include "gtest/gtest.h"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <unistd.h>
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

// ------------------------------------------------------------------------------------- THE LAYERS (OpenWork §2, ruled 2026-11-09)
// The shipped file stays curated; imports go to the user's library (here a temp file through QCHEM_USER_MATERIALS), a workspace
// adds its own on top, and one name is one structure across all of them.
namespace
{
struct UserLibrary
{
    std::string path;
    UserLibrary()
    {
        path=(std::filesystem::temp_directory_path() / ("qchem-ut-materials-"+std::to_string(::getpid())+"-"+std::to_string(++n)+".json")).string();
        std::filesystem::remove(path);
        ::setenv("QCHEM_USER_MATERIALS", path.c_str(), 1);
    }
    ~UserLibrary() { std::filesystem::remove(path); ::setenv("QCHEM_USER_MATERIALS", "/nonexistent/qchem-test-user-materials.json", 1); }
    static inline int n=0;
};
const char* kRocksalt=R"({"comment": "an import", "lattice": {"vectors_bohr": [[0,4,4],[4,0,4],[4,4,0]]},
                          "atoms": [{"el":"Mg","frac":[0,0,0]},{"el":"O","frac":[0.5,0.5,0.5]}], "provenance": {"source":"x.cif"}})";
std::string Thrown(const std::function<void()>& f)
{
    try { f(); } catch (const std::exception& e) { return e.what(); }
    return "";
}
}

TEST(StructureData, ASavedMaterialIsAUserLayerCellBesideTheShippedOnes)
{
    UserLibrary lib;
    const size_t shipped=SD::CellNames().size();
    SD::SaveUserMaterial("MgO_imported", kRocksalt);
    ASSERT_EQ(SD::CellNames().size(), shipped+1);
    EXPECT_EQ(SD::CellNames().back(), "MgO_imported") << "the user's layer follows the shipped one";
    const SD::CellEntry e=SD::GetCellEntry("MgO_imported");
    EXPECT_EQ(e.layer, SD::Layer::User);
    EXPECT_EQ(e.file, lib.path);
    EXPECT_EQ(SD::KindOf("MgO_imported"), SD::Kind::Cell);
    EXPECT_EQ(SD::GetCell("MgO_imported").GetNumAtoms(), 2u);
    EXPECT_NE(e.entryJson.find("provenance"), std::string::npos) << "the entry is saved whole";
    // a second save of the name is refused unless it says replace; a replace keeps the entry's place
    EXPECT_NE(Thrown([]{ SD::SaveUserMaterial("MgO_imported", kRocksalt); }).find("replace it, or choose another name"), std::string::npos);
    SD::SaveUserMaterial("Second", kRocksalt);
    SD::SaveUserMaterial("MgO_imported", kRocksalt, /*replace*/true);
    EXPECT_EQ(SD::CellNames()[shipped], "MgO_imported");
}

TEST(StructureData, AUserEntryMayNotReuseAShippedNameAndTheShippedFileIsUntouched)
{
    UserLibrary lib;
    EXPECT_NE(Thrown([]{ SD::SaveUserMaterial("NaF_rocksalt", kRocksalt, true); }).find("already in the shipped materials"), std::string::npos);
    EXPECT_NE(Thrown([]{ SD::SaveUserMaterial("H2O", kRocksalt); }).find("shipped molecule name"), std::string::npos);
    EXPECT_NE(Thrown([]{ SD::SaveUserMaterial("_x", kRocksalt); }).find("must be non-empty, not start with '_'"), std::string::npos);
    EXPECT_NE(Thrown([]{ SD::SaveUserMaterial("bad", R"({"lattice":{"vectors_bohr":[[1,0,0],[0,1,0],[0,0,1]]},"atoms":[{"el":"Qq","frac":[0,0,0]}]})"); })
              .find("unknown element symbol 'Qq'"), std::string::npos) << "the parser refuses before the file is written";
    EXPECT_FALSE(std::filesystem::exists(lib.path)) << "nothing refused reached the file";
    // a hand-edited library that reuses a shipped name is refused the moment it is read, naming the file
    { std::ofstream f(lib.path); f << R"({"Si_diamond": )" << kRocksalt << "}"; }
    const std::string why=Thrown([]{ (void)SD::CellNames(); });
    EXPECT_NE(why.find(lib.path), std::string::npos) << why;
    EXPECT_NE(why.find("'Si_diamond' is already a name in the shipped materials"), std::string::npos) << why;
}

TEST(StructureData, AWorkspaceLibraryLayersOnTopAndMayNotCollideWithTheUsers)
{
    UserLibrary lib;
    SD::SaveUserMaterial("Mine", kRocksalt);
    const std::string ws=lib.path+".workspace.json";
    { std::ofstream f(ws); f << R"({"Theirs": )" << kRocksalt << R"(, "Mine": )" << kRocksalt << "}"; }
    EXPECT_NE(Thrown([&]{ SD::AddWorkspaceLibrary(ws); }).find("'Mine' is already a name in your materials library"), std::string::npos);
    { std::ofstream f(ws); f << R"({"Theirs": )" << kRocksalt << "}"; }
    SD::AddWorkspaceLibrary(ws);
    EXPECT_EQ(SD::GetCellEntry("Theirs").layer, SD::Layer::Workspace);
    EXPECT_EQ(SD::CellNames().back(), "Theirs");
    EXPECT_NE(Thrown([]{ SD::SaveUserMaterial("Theirs", kRocksalt); }).find("already in a workspace library"), std::string::npos);
    std::filesystem::remove(ws);   // a missing file is an empty layer, so the rest of this process is unaffected
}
