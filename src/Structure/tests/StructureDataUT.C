// File: src/Structure/tests/StructureDataUT.C  qchem.StructureData: the pre-defined molecules and cells
// (D-STRUCTDATA).  Structure only -- valence counts are the MoleculeCalculation library's (src/Calculation/tests/Materials.C).
#include "gtest/gtest.h"
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
// The shipped file stays curated; imports go to the user's library and a workspace adds its own on top -- each passed in as an
// explicit PATH (the library reads no environment, controller ruling 2026-11-09) -- and one name is one structure across them.
namespace
{
struct TempLibrary
{
    std::string path;
    SD::LibraryFile file;
    explicit TempLibrary(SD::Layer l=SD::Layer::User)
        : path((std::filesystem::temp_directory_path() / ("qchem-ut-materials-"+std::to_string(::getpid())+"-"+std::to_string(++n)+".json")).string())
        , file{l, path}
    { std::filesystem::remove(path); }
    ~TempLibrary() { std::filesystem::remove(path); }
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

TEST(StructureData, ASavedMaterialIsAUserLayerEntryBesideTheShippedOnes)
{
    TempLibrary mine;
    const size_t shipped=SD::CellEntries().size();
    EXPECT_EQ(shipped, SD::CellNames().size());
    EXPECT_EQ(SD::CellEntries({mine.file}).size(), shipped) << "a missing library is an empty layer";
    SD::SaveMaterial(mine.file, "MgO_imported", kRocksalt);
    const auto all=SD::CellEntries({mine.file});
    ASSERT_EQ(all.size(), shipped+1);
    EXPECT_EQ(all.back().name, "MgO_imported") << "the user's layer follows the shipped one";
    const SD::CellEntry e=SD::GetCellEntry("MgO_imported", {mine.file});
    EXPECT_EQ(e.layer, SD::Layer::User);
    EXPECT_EQ(e.file, mine.path);
    EXPECT_EQ(SD::CellFromEntry(e.entryJson, e.name).GetNumAtoms(), 2u);
    EXPECT_NE(e.entryJson.find("provenance"), std::string::npos) << "the entry is saved whole";
    EXPECT_EQ(SD::CellNames().size(), shipped) << "the shipped list (what a deck may name) is unchanged";
    EXPECT_NE(Thrown([]{ (void)SD::KindOf("MgO_imported"); }).find("stated INLINE"), std::string::npos) << "a deck never names it";
    // a second save of the name is refused unless it says replace; a replace keeps the entry's place
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial(mine.file, "MgO_imported", kRocksalt); }).find("replace it, or choose another name"), std::string::npos);
    SD::SaveMaterial(mine.file, "Second", kRocksalt);
    SD::SaveMaterial(mine.file, "MgO_imported", kRocksalt, /*replace*/true);
    EXPECT_EQ(SD::CellEntries({mine.file})[shipped].name, "MgO_imported");
}

TEST(StructureData, AUserEntryMayNotReuseAShippedNameAndTheShippedFileIsNeverWritten)
{
    TempLibrary mine;
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial(mine.file, "NaF_rocksalt", kRocksalt, true); }).find("already in the shipped materials"), std::string::npos);
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial(mine.file, "H2O", kRocksalt); }).find("shipped molecule name"), std::string::npos);
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial(mine.file, "_x", kRocksalt); }).find("must be non-empty, not start with '_'"), std::string::npos);
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial({SD::Layer::Shipped, mine.path}, "x", kRocksalt); }).find("CURATED"), std::string::npos);
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial(mine.file, "bad", R"({"lattice":{"vectors_bohr":[[1,0,0],[0,1,0],[0,0,1]]},"atoms":[{"el":"Qq","frac":[0,0,0]}]})"); })
              .find("unknown element symbol 'Qq'"), std::string::npos) << "the parser refuses before the file is written";
    EXPECT_FALSE(std::filesystem::exists(mine.path)) << "nothing refused reached the file";
    // a hand-edited library that reuses a shipped name is refused the moment it is read, naming the file
    { std::ofstream f(mine.path); f << R"({"Si_diamond": )" << kRocksalt << "}"; }
    const std::string why=Thrown([&]{ (void)SD::CellEntries({mine.file}); });
    EXPECT_NE(why.find(mine.path), std::string::npos) << why;
    EXPECT_NE(why.find("'Si_diamond' is already a name in the shipped materials"), std::string::npos) << why;
}

TEST(StructureData, AWorkspaceLibraryLayersOnTopAndMayNotCollideWithTheUsers)
{
    TempLibrary mine, ws(SD::Layer::Workspace);
    SD::SaveMaterial(mine.file, "Mine", kRocksalt);
    { std::ofstream f(ws.path); f << R"({"Theirs": )" << kRocksalt << R"(, "Mine": )" << kRocksalt << "}"; }
    EXPECT_NE(Thrown([&]{ (void)SD::CellEntries({mine.file, ws.file}); }).find("'Mine' is already a name in your materials library"), std::string::npos);
    { std::ofstream f(ws.path); f << R"({"Theirs": )" << kRocksalt << "}"; }
    EXPECT_EQ(SD::GetCellEntry("Theirs", {mine.file, ws.file}).layer, SD::Layer::Workspace);
    EXPECT_EQ(SD::CellEntries({mine.file, ws.file}).back().name, "Theirs");
    EXPECT_NE(Thrown([&]{ SD::SaveMaterial(mine.file, "Theirs", kRocksalt, false, {ws.file}); }).find("already in a workspace library"), std::string::npos);
}
