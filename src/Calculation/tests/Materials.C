// File: src/Calculation/tests/Materials.C  The pre-defined materials/molecules data (doc/OpenWork.md row MD).
//
// The data file is what every SCF-running test builds its cell from, so what is pinned here is the CONTRACT
// of the file: every entry loads, the electron counts derive to the numbers the anchors were banked at, the
// one re-based cell (MnO AFM-II) reproduces the nine literals the campaign wrote by hand, and the decoration
// reaches the atoms.  A lattice constant in the file that drifts moves an anchor -- these tests do not pin
// the constants themselves (the anchored SCF tests do), only that the file says what the tests expect.
#include <filesystem>
#include "gtest/gtest.h"
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

import qchem.Materials;
import qchem.StructureData;   // SaveMaterial, GetCellEntry: the user layer
import qchem.CrystalFiles;   // an imported entry IS an entry: it runs through Inline
import qchem.Pseudopotential.Set;
import qchem.Types;

using namespace qchem;
namespace M=qchem::Materials;

TEST(Materials, EveryEntryLoadsAndDerivesItsElectronCount)
{
    struct Expect { const char* name; size_t atoms; int Nelec; };
    for (Expect e : { Expect{"Si_diamond",2,8}, {"Al_fcc",1,3}, {"Na_fcc",1,1}, {"NaF_rocksalt",2,8},
                      {"CsI_cscl",2,8}, {"MnO_AFM2",4,26}, {"NiO_AFM2",4,32}, {"Si_box16",1,4}, {"Na_box16",1,1},
                      {"Mn_box16",1,7}, {"O2_box16",2,12}, {"Mn2_box7",2,14}, {"Na2_box16",2,2},
                      {"LiMn2O4_spinel",14,78}, {"MnO2_lambda_spinel",12,76},
                      {"Si_diamond_2x1x1",4,16}, {"Si_diamond_2x2x1",8,32}, {"Si_diamond_2x2x2",16,64} })
    {
        M::Material m=M::Get(e.name);
        EXPECT_EQ(m.name, e.name);
        EXPECT_EQ(m.cell->GetNumAtoms(), e.atoms) << e.name;
        EXPECT_EQ(m.Nelec(), e.Nelec) << e.name;
        EXPECT_FALSE(m.species.empty()) << e.name;
    }
    // EVERY entry in the file -- not only the ones with a row above -- loads, names itself, has atoms and a species vocabulary, and derives a positive electron
    // count.  (No pinned entry COUNT: adding a material must not break a test; a broken entry must.)
    const std::vector<std::string> names=M::Names();
    EXPECT_FALSE(names.empty());
    std::set<std::string> seen;
    for (const std::string& n : names)
    {
        EXPECT_TRUE(seen.insert(n).second) << n << " appears twice in the pick-list";
        M::Material m=M::Get(n);
        EXPECT_EQ(m.name, n);
        EXPECT_GT(m.cell->GetNumAtoms(), 0u) << n;
        EXPECT_FALSE(m.species.empty()) << n;
        EXPECT_GT(m.Nelec(), 0) << n;
        // The vocabulary comes FROM the named set (D-STRUCTDATA step 2): every species is the set's own variant for that element.
        EXPECT_FALSE(m.pp.empty()) << n;
        const auto set=qchem::Pseudopotential::PseudopotentialSet::Get(m.pp);
        for (const auto& [el,q] : m.species) EXPECT_EQ(q, set.Valence(el)) << n << ' ' << el;
    }
}

// The re-based cell: CubicF at a=8.40 under T=[[0,1,1],[1,0,1],[1,1,0]] is the AFM-II rhombohedral cell
// GPW_SCF_UT.C's RunMnO spelled as (a, a/2, a/2; a/2, a, a/2; a/2, a/2, a) -- and the +1/-1 site spins
// land on the two Mn as the spin-flip bits the seed reads, with the O unflipped.
TEST(Materials, MnO_AFM2_IsTheHandWrittenCellWithItsDecoration)
{
    M::Material m=M::Get("MnO_AFM2");
    const double a=8.40;
    const Matrix3D<double> hand(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) EXPECT_NEAR(m.cell->GetCellMatrix()(i,j), hand(i,j), 1e-14);
    std::vector<int> Z; std::vector<bool> flip;
    m.cell->ForEachSite([&](int z, const rvec3_t&, bool f){ Z.push_back(z); flip.push_back(f); });
    EXPECT_EQ(Z,    (std::vector<int>{25,25,8,8}));
    EXPECT_EQ(flip, (std::vector<bool>{false,true,false,false})) << "-m on the second Mn, O undecorated";
    EXPECT_EQ(m.species, (std::vector<std::pair<std::string,int>>{{"Mn",7},{"O",6}}));
}

// NiO AFM-II is the SAME construction as MnO one lattice constant apart, and its geometry is not ours to
// choose: a = 7.88 bohr is the cell QE's hp.x NiO benchmark uses, and the matched-UPF hp.x run that banked
// U(Ni 3d) = 5.27 eV (doc/OpenWork.md step 5) ran it -- an estimator compared against an oracle must stand
// on the oracle's cell.  The QE deck (IntegrationTests/QE/NiOg.scf.1.in) writes the same cell as
// CELL_PARAMETERS (1,1/2,1/2; 1/2,1,1/2; 1/2,1/2,1) x alat with ATOMIC_POSITIONS {alat} Ni1 (0,0,0),
// Ni2 (1/2,1/2,0), O (1/2,0,0), O (1,1/2,0); in fractional coordinates of that cell those are Ni 0 and 1/2,
// O 3/4 and 1/4 -- the same four sites this entry lists, with the two O in the other order.
TEST(Materials, NiO_AFM2_IsTheQEBenchmarkCell)
{
    M::Material m=M::Get("NiO_AFM2");
    const double a=7.88;                                    // QE celldm(1), the hp.x NiO benchmark cell
    const Matrix3D<double> hand(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) EXPECT_NEAR(m.cell->GetCellMatrix()(i,j), hand(i,j), 1e-14);
    std::vector<int> Z; std::vector<bool> flip;
    m.cell->ForEachSite([&](int z, const rvec3_t&, bool f){ Z.push_back(z); flip.push_back(f); });
    EXPECT_EQ(Z,    (std::vector<int>{28,28,8,8}));
    EXPECT_EQ(flip, (std::vector<bool>{false,true,false,false})) << "-m on the second Ni, O undecorated";
    EXPECT_EQ(m.species, (std::vector<std::pair<std::string,int>>{{"Ni",10},{"O",6}}));
}

// Si diamond IS the cell every Si anchor was banked on: FCC a=10.26, atoms at 0 and 1/4 -- in Cartesian.
TEST(Materials, SiDiamondAtomsAreWhereTheAnchorsExpect)
{
    M::Material m=M::Get("Si_diamond");
    std::vector<rvec3_t> R;
    m.cell->ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
    ASSERT_EQ(R.size(), 2u);
    EXPECT_NEAR(R[0].x, 0.0, 1e-14);
    EXPECT_NEAR(R[1].x, 10.26/4, 1e-12) << "frac (1/4,1/4,1/4) of the FCC primitive cell = (a/4,a/4,a/4): the diamond bond along [111]";
    EXPECT_NEAR(R[1].y, R[1].x, 1e-14); EXPECT_NEAR(R[1].z, R[1].x, 1e-14);
    // The lattice-constant override -- a ladder point -- scales the cell and the atoms with it.
    M::Material big=M::Get("Si_diamond", 2*10.26);
    EXPECT_NEAR(big.cell->GetCellVolume(), 8*m.cell->GetCellVolume(), 1e-9);
}

TEST(Materials, BoxesAndMoleculesAndMisses)
{
    M::Material mn=M::AtomInBox("Mn", 7, 16.0);
    EXPECT_EQ(mn.Nelec(), 7);
    EXPECT_NEAR(mn.cell->GetCellVolume(), 16.0*16.0*16.0, 1e-9);
    M::Material na2=M::DimerInBox("Na", 1, 16.0, 5.8, /*afm*/true);
    EXPECT_EQ(na2.Nelec(), 2);
    std::vector<bool> flip; na2.cell->ForEachSite([&](int, const rvec3_t&, bool f){ flip.push_back(f); });
    EXPECT_EQ(flip, (std::vector<bool>{false,true}));
    // ...and the file's Na2_box16 is that dimer at the banked size, atom for atom.
    M::Material file=M::Get("Na2_box16");
    std::vector<rvec3_t> Rf, Rd;
    file.cell->ForEachSite([&](int, const rvec3_t& r, bool){ Rf.push_back(r); });
    na2 .cell->ForEachSite([&](int, const rvec3_t& r, bool){ Rd.push_back(r); });
    for (size_t i=0;i<2;i++) EXPECT_NEAR(Rf[i].x, Rd[i].x, 1e-12);

    EXPECT_THROW(M::Get("Unobtainium"), std::runtime_error);
    EXPECT_THROW(M::Get("_doc"),        std::runtime_error) << "the documentation key is not an entry";
    EXPECT_THROW(M::Get("H2O"),         std::runtime_error) << "a molecule is not a material";
}

// AN IMPORT IS AN ENTRY (OpenWork §2, ruled 2026-11-09): what a CIF translates to, provenance and all, is a material the deck path
// (Inline) takes as it is -- rocksalt NaF from a two-site F m -3 m asymmetric unit lands on the shipped NaF_rocksalt cell.
TEST(Materials, AnImportedCrystalIsAnInlineMaterial)
{
    const std::string cif="data_NaF\n_cell_length_a 4.619717051183\n_cell_length_b 4.619717051183\n_cell_length_c 4.619717051183\n"
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
        "loop_\n_symmetry_equiv_pos_as_xyz\nx,y,z\nx,1/2+y,1/2+z\n1/2+x,y,1/2+z\n1/2+x,1/2+y,z\n-x,-y,-z\n"
        "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nNa 0 0 0\nF 0.5 0.5 0.5\n";
    const auto imp=CrystalFiles::Read(cif, CrystalFiles::Format::CIF, "NaF.cif", "2026-11-09");
    const Materials::Material m=Materials::Inline(imp.entry, "gth-lda-valence");
    const Materials::Material shipped=Materials::Get("NaF_rocksalt");
    EXPECT_EQ(m.cell->GetNumAtoms(), 2u);
    EXPECT_EQ(m.Nelec(), shipped.Nelec());
    EXPECT_NEAR(m.cell->GetCellVolume(), shipped.cell->GetCellVolume(), 1e-9) << "a = 8.73 bohr, the shipped cell";
}

// A SAVED MATERIAL RUNS INLINE, NEVER BY NAME (controller ruling 2026-11-09): the deck is the one source of a run's numbers, so a
// material of the user's own enters it as the entry itself -- which Inline takes, provenance and all -- and Get refuses its name.
TEST(Materials, AUserMaterialRunsInlineNeverByName)
{
    const StructureData::LibraryFile mine{StructureData::Layer::User,
                                          (std::filesystem::temp_directory_path() / "qchem-ut-calc-materials.json").string()};
    std::filesystem::remove(mine.path);
    StructureData::SaveMaterial(mine, "NaF_mine", R"({"lattice":{"vectors_bohr":[[0,4.365,4.365],[4.365,0,4.365],[4.365,4.365,0]]},"atoms":[{"el":"Na","frac":[0,0,0]},{"el":"F","frac":[0.5,0.5,0.5]}],"provenance":{"source":"NaF.cif"}})");
    try { (void)Materials::Get("NaF_mine"); ADD_FAILURE() << "a user material ran by name"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("a user material must be inlined in the deck"), std::string::npos) << e.what(); }
    const Materials::Material m=Materials::Inline(StructureData::GetCellEntry("NaF_mine", {mine}).entryJson, "gth-lda-valence");
    EXPECT_EQ(m.Nelec(), Materials::Get("NaF_rocksalt").Nelec());
    std::filesystem::remove(mine.path);
}
