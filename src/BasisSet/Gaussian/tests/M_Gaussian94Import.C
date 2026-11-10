// File: M_Gaussian94Import.C  A basis from the Basis Set Exchange, translated into the BasisSetData layout, is what the
// Gaussian94Reader builds; what the library cannot represent is refused by name (OpenWork §2 BSE row, lane G 2026-11-10).
//
// The two inputs are the Exchange's own "gaussian94" output, verbatim (cc-pVDZ H and O; 6-31G O), fetched 2026-11-10 from
// https://www.basissetexchange.org/api/basis/<name>/format/gaussian94/ -- the D exponents, the mixed-case and "!" header, and
// the Pople SP shells whose S and P coefficients differ are exactly what the translator exists for.
#include "gtest/gtest.h"
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94Import;
import qchem.BasisSet.Gaussian.Point.Factory;   // FactoryFromFile: the translated file builds
import qchem.BasisSet;
import qchem.Structure;
import qchem.Types;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;

namespace
{
const char* kCCpVDZ_HO = R"(!----------------------------------------------------------------------
! Basis Set Exchange
! Version 0.12
! https://www.basissetexchange.org
!----------------------------------------------------------------------
!   Basis set: cc-pVDZ
! Description: cc-pVDZ
!        Role: orbital
!     Version: 1  (Data from ccRepo/Grant Hill)
!----------------------------------------------------------------------


H     0
S    4   1.00
      1.301000D+01           1.968500D-02
      1.962000D+00           1.379770D-01
      4.446000D-01           4.781480D-01
      1.220000D-01           5.012400D-01
S    1   1.00
      1.220000D-01           1.000000D+00
P    1   1.00
      7.270000D-01           1.0000000
****
O     0
S    9   1.00
      1.172000D+04           7.100000D-04
      1.759000D+03           5.470000D-03
      4.008000D+02           2.783700D-02
      1.137000D+02           1.048000D-01
      3.703000D+01           2.830620D-01
      1.327000D+01           4.487190D-01
      5.025000D+00           2.709520D-01
      1.013000D+00           1.545800D-02
      3.023000D-01          -2.585000D-03
S    9   1.00
      1.172000D+04          -1.600000D-04
      1.759000D+03          -1.263000D-03
      4.008000D+02          -6.267000D-03
      1.137000D+02          -2.571600D-02
      3.703000D+01          -7.092400D-02
      1.327000D+01          -1.654110D-01
      5.025000D+00          -1.169550D-01
      1.013000D+00           5.573680D-01
      3.023000D-01           5.727590D-01
S    1   1.00
      3.023000D-01           1.000000D+00
P    4   1.00
      1.770000D+01           4.301800D-02
      3.854000D+00           2.289130D-01
      1.046000D+00           5.087280D-01
      2.753000D-01           4.605310D-01
P    1   1.00
      2.753000D-01           1.000000D+00
D    1   1.00
      1.185000D+00           1.0000000
****

)";

const char* k631G_O = R"(!----------------------------------------------------------------------
! Basis Set Exchange
!   Basis set: 6-31G
!----------------------------------------------------------------------


O     0
S    6   1.00
      0.5484671660D+04       0.1831074430D-02
      0.8252349460D+03       0.1395017220D-01
      0.1880469580D+03       0.6844507810D-01
      0.5296450000D+02       0.2327143360D+00
      0.1689757040D+02       0.4701928980D+00
      0.5799635340D+01       0.3585208530D+00
SP   3   1.00
      0.1553961625D+02      -0.1107775495D+00       0.7087426823D-01
      0.3599933586D+01      -0.1480262627D+00       0.3397528391D+00
      0.1013761750D+01       0.1130767015D+01       0.7271585773D+00
SP   1   1.00
      0.2700058226D+00       0.1000000000D+01       0.1000000000D+01
****
)";

std::filesystem::path Written(const std::string& name, const std::string& text)
{
    const auto p=std::filesystem::temp_directory_path() / ("qchem_UT_Gaussian94Import_"+name+".bsd");
    std::ofstream(p) << text;
    return p;
}

size_t Functions(const std::filesystem::path& file, const Structure& st)
{
    std::unique_ptr<qchem::BasisSet::tBasisSet<double>> bs(FactoryFromFile(file.string(), &st, Engine::MnD, Angular::Cartesian));
    return bs->GetNumFunctions();
}

Molecule Water()
{
    Molecule w;
    w.Insert(Site{8, Vector3D(0.0, 0.0, 0.0)});
    w.Insert(Site{1, Vector3D(0.0, 1.431, 1.107)});
    w.Insert(Site{1, Vector3D(0.0, -1.431, 1.107)});
    return w;
}
}

TEST(Gaussian94Import, AnExchangeCCpVDZBuildsWaterWithEveryShellOfTheFile)
{
    const auto r=ImportGaussian94(kCCpVDZ_HO, "cc-pVDZ", {"source: Basis Set Exchange", "version: 1"});
    ASSERT_TRUE(r) << r.Error();
    EXPECT_EQ(r->elements, (std::vector<std::string>{"H", "O"}));
    EXPECT_EQ(r->maxL, 2);
    EXPECT_NE(r->bsd.find("BASIS=\"cc-pVDZ\""), std::string::npos);
    EXPECT_NE(r->bsd.find("version: 1"), std::string::npos) << "the provenance rides in the preamble";
    const std::string blocks=r->bsd.substr(r->bsd.find("BASIS="));
    EXPECT_EQ(blocks.find("D+"), std::string::npos) << "a Fortran D exponent survived";
    EXPECT_EQ(blocks.find("D-"), std::string::npos) << "a Fortran D exponent survived";
    EXPECT_NE(blocks.find("1.301000E+01"), std::string::npos) << "the digits are the file's own, D -> E only";
    // Cartesian: O 3s 2p 1d = 3 + 6 + 6, each H 2s 1p = 2 + 3
    EXPECT_EQ(Functions(Written("ccpvdz", r->bsd), Water()), 25u);
}

TEST(Gaussian94Import, APopleSPShellIsItsSAndItsPContractionsOnTheSameExponents)
{
    const auto r=ImportGaussian94(k631G_O, "6-31G");
    ASSERT_TRUE(r) << r.Error();
    ASSERT_EQ(r->notes.size(), 2u);
    EXPECT_NE(r->notes[0].find("2 shared-exponent shell(s) (SP)"), std::string::npos) << r->notes[0];
    Molecule o; o.Insert(Site{8, Vector3D(0.0, 0.0, 0.0)});
    // 1s(6) + SP(3) + SP(1): s 3, p 2 -> 3 + 6 Cartesian functions, read without the reader's shared-radial refusal
    EXPECT_EQ(Functions(Written("631g", r->bsd), o), 9u);
}

TEST(Gaussian94Import, WhatTheLibraryCannotRepresentIsRefusedByName)
{
    auto why=[](const std::string& text){ auto r=ImportGaussian94(text, "t"); return r ? std::string("ACCEPTED") : r.Error(); };
    const std::string h="H     0\nS    1   1.00\n      1.0D+00   1.0D+00\n****\n";
    EXPECT_NE(why(h+"FE     0\nFE-ECP     3     10\n").find("ECP"), std::string::npos);
    EXPECT_NE(why("H     0\nS    1   2.00\n      1.0   1.0\n****\n").find("scale factor"), std::string::npos);
    EXPECT_NE(why("H     0\nK    1   1.00\n      1.0   1.0\n****\n").find("angular label 'K'"), std::string::npos);
    EXPECT_NE(why(h+h).find("appears twice"), std::string::npos);
    EXPECT_NE(why("H     0\nS    1   1.00\n      1.0   1.0\n").find("terminator"), std::string::npos);
    EXPECT_NE(why("Xx     0\n").find("not an element"), std::string::npos);
    EXPECT_NE(why("H     0\nS    2   1.00\n      1.0   1.0\n****\n").find("needs 2 numbers"), std::string::npos) << "a short shell is named";
    EXPECT_NE(why("! only comments\n").find("no element block"), std::string::npos);
    EXPECT_EQ(why(h), "ACCEPTED");
}
