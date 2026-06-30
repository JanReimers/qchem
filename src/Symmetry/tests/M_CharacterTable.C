// File: UnitTests/M_CharacterTable.C  Abelian point-group character tables (stage 3a-i).
#include <gtest/gtest.h>
#include <string>
#include <vector>
import qchem.Symmetry.CharacterTable;
using namespace qchem;

using namespace qchem::Symmetry;

// Every abelian table is square (h irreps x h ops), its first irrep is totally symmetric,
// and its rows are orthogonal with norm^2 = h (the great orthogonality theorem, 1-D irreps).
TEST(CharacterTable, all_abelian_orthonormal)
{
    for (std::string s : {"C1","Ci","Cs","C2","C2h","C2v","D2","D2h"})
    {
        auto t = AbelianCharacterTable(s);
        size_t h = t.order();
        EXPECT_EQ(t.nIrreps(), h);                       // abelian: #irreps == group order
        for (const auto& row : t.chi) EXPECT_EQ(row.size(), h);
        for (int c : t.chi[0]) EXPECT_EQ(c, 1);          // first irrep totally symmetric

        for (size_t r1=0;r1<t.nIrreps();++r1)
            for (size_t r2=0;r2<t.nIrreps();++r2)
            {
                int dot=0; for (size_t k=0;k<h;k++) dot += t.chi[r1][k]*t.chi[r2][k];
                EXPECT_EQ(dot, (r1==r2) ? (int)h : 0);
            }
    }
}

TEST(CharacterTable, c2v_labels)
{
    auto t = AbelianCharacterTable("C2v");
    EXPECT_EQ(t.symbol, "C2v");
    EXPECT_EQ(t.order(), 4u);
    EXPECT_EQ(t.irreps, (std::vector<std::string>{"A1","A2","B1","B2"}));
    EXPECT_EQ(t.opTags, (std::vector<std::string>{"E","C2z","sxz","syz"}));
}

TEST(CharacterTable, d2h_is_8x8)
{
    auto t = AbelianCharacterTable("D2h");
    EXPECT_EQ(t.order(), 8u);
    EXPECT_EQ(t.nIrreps(), 8u);
    EXPECT_EQ(t.irreps[0], "Ag");
    EXPECT_EQ(t.irreps[4], "Au");
}
