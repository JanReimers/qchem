// File: src/Calculation/tests/Materials.C  The pre-defined materials/molecules data (doc/OpenWork.md row MD).
//
// The data file is what every SCF-running test builds its cell from, so what is pinned here is the CONTRACT
// of the file: every entry loads, the electron counts derive to the numbers the anchors were banked at, the
// one re-based cell (MnO AFM-II) reproduces the nine literals the campaign wrote by hand, and the decoration
// reaches the atoms.  A lattice constant in the file that drifts moves an anchor -- these tests do not pin
// the constants themselves (the anchored SCF tests do), only that the file says what the tests expect.
#include "gtest/gtest.h"
#include <stdexcept>
#include <string>
#include <vector>

import qchem.Materials;
import qchem.Types;

using namespace qchem;
namespace M=qchem::Materials;

TEST(Materials, EveryEntryLoadsAndDerivesItsElectronCount)
{
    struct Expect { const char* name; size_t atoms; int Nelec; };
    for (Expect e : { Expect{"Si_diamond",2,8}, {"Al_fcc",1,3}, {"Na_fcc",1,1}, {"NaF_rocksalt",2,8},
                      {"CsI_cscl",2,8}, {"MnO_AFM2",4,26}, {"Si_box16",1,4}, {"Na_box16",1,1},
                      {"Mn_box16",1,7}, {"O2_box16",2,12}, {"Mn2_box7",2,14}, {"Na2_box16",2,2} })
    {
        M::Material m=M::Get(e.name);
        EXPECT_EQ(m.name, e.name);
        EXPECT_EQ(m.cell->GetNumAtoms(), e.atoms) << e.name;
        EXPECT_EQ(m.Nelec(), e.Nelec) << e.name;
        EXPECT_FALSE(m.species.empty()) << e.name;
    }
    const std::vector<std::string> names=M::Names();
    EXPECT_EQ(names.size(), 12u) << "the file has exactly the entries this test knows; add a row here when you add one";
    EXPECT_EQ(names.front(), "Si_diamond") << "file order is the pick-list order";
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

    Molecule w=M::GetMolecule("H2O");
    EXPECT_EQ(w.GetNumAtoms(), 3u);
    EXPECT_EQ(M::MoleculeNames().size(), 4u);

    EXPECT_THROW(M::Get("Unobtainium"), std::runtime_error);
    EXPECT_THROW(M::Get("_doc"),        std::runtime_error) << "the documentation key is not an entry";
    EXPECT_THROW(M::GetMolecule("C60"), std::runtime_error);
}
