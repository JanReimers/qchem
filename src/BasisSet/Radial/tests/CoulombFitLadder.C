// File: src/BasisSet/Radial/tests/CoulombFitLadder.C  The radial basis's Coulomb-fit ladder (the atom pilot's Dunlap axis).
#include "gtest/gtest.h"
#include <memory>
#include <stdexcept>
#include <vector>
#include <cmath>

import qchem.BasisSet.Radial.Factory;
import qchem.BasisSet;
using namespace qchem;
using namespace qchem::BasisSet::Radial;

namespace
{
std::unique_ptr<BasisSet::Real_BS> HighPool(Type t, size_t Z) {return std::unique_ptr<BasisSet::Real_BS>(Factory(BasisSetAccuracy::High, t, Z));}
}

// Rung 0 is the basis's own Coulomb fit (CreateCDFitBasisSet, the s block x 2), so a run on rung 0 is the run before the ladder existed.
// Rungs >= 1 add functions (Slater: each block's exponents x 2 at its own power 2l; Gaussian: every block's merged into one s set), then
// insert 2^{r-1} - 1 points between neighbours: NESTED, so the Dunlap error can only fall.  Gaussian's sizes say so exactly,
// n_r = n_1 + (n_1 - 1)(2^{r-1} - 1); Slater's grow block by block.
TEST(CoulombFitLadder, RungZeroIsTheBasissOwnFitAndTheRungsAboveAreNested)
{
    for (Type t : {Type::Slater, Type::Gaussian})
        for (size_t Z : {18ul, 54ul, 92ul})
        {
            const auto bs=HighPool(t, Z);
            const auto* ladder=dynamic_cast<const CoulombFitLadder*>(bs.get());
            ASSERT_NE(ladder, nullptr) << "a Slater / Gaussian NR basis carries the ladder";
            const std::unique_ptr<BasisSet::rFIT_CD_ABS> own(bs->CreateCDFitBasisSet(nullptr));
            const std::unique_ptr<BasisSet::rFIT_CD_ABS> r0(ladder->MakeCoulombFitRung(0));
            EXPECT_EQ(r0->GetNumFunctions(), own->GetNumFunctions()) << "Z " << Z;
            const std::unique_ptr<BasisSet::rFIT_CD_ABS> r1(ladder->MakeCoulombFitRung(1));
            const size_t n1=r1->GetNumFunctions();
            EXPECT_GE(n1, r0->GetNumFunctions()) << "the merged set holds at least the s block's span, Z " << Z;
            size_t below=n1;
            for (int r=2; r<CoulombFitLadder::kCoulombFitRungs; r++)
            {
                const std::unique_ptr<BasisSet::rFIT_CD_ABS> fr(ladder->MakeCoulombFitRung(r));
                const size_t m=(size_t(1)<<(r-1))-1;
                if (t==Type::Gaussian) EXPECT_EQ(fr->GetNumFunctions(), n1+(n1-1)*m) << "rung " << r << ", Z " << Z;
                EXPECT_GT(fr->GetNumFunctions(), below);
                below=fr->GetNumFunctions();
            }
        }
}

// A Slater rung >= 1 fits each l-shell with its own shape, r^{2l} e^{-zeta r}: U's rung 1 adds the p, d and f blocks' functions to the s
// block's, and every one is a normalized positive charge distribution (a power-m function with a wrong norm or power shows here first).
TEST(CoulombFitLadder, ASlaterRungAboveZeroCarriesEachShellsOwnPower)
{
    const auto bs=HighPool(Type::Slater, 92);
    const auto* ladder=dynamic_cast<const CoulombFitLadder*>(bs.get());
    ASSERT_NE(ladder, nullptr);
    const std::unique_ptr<BasisSet::rFIT_CD_ABS> r0(ladder->MakeCoulombFitRung(0)), r1(ladder->MakeCoulombFitRung(1));
    EXPECT_GT(r1->GetNumFunctions(), r0->GetNumFunctions()) << "the p, d and f blocks' own functions are added";
    const auto* nonOrtho=dynamic_cast<const BasisSet::FIT_CD_NonOrtho*>(r1.get());
    ASSERT_NE(nonOrtho, nullptr);
    const auto q=nonOrtho->Charge();
    for (size_t i=0;i<q.size();i++) { EXPECT_GT(q[i], 0.0); EXPECT_TRUE(std::isfinite(q[i])); }
}

TEST(CoulombFitLadder, ARungOffTheLadderThrows)
{
    const auto bs=HighPool(Type::Slater, 54);
    const auto* ladder=dynamic_cast<const CoulombFitLadder*>(bs.get());
    ASSERT_NE(ladder, nullptr);
    EXPECT_THROW(delete ladder->MakeCoulombFitRung(-1), std::out_of_range);
    EXPECT_THROW(delete ladder->MakeCoulombFitRung(CoulombFitLadder::kCoulombFitRungs), std::out_of_range);
}

// A basis with no Coulomb fit (B-spline: no DFT block) carries no ladder, so the atom pilot can tell without a concrete cast.
TEST(CoulombFitLadder, ABSplineBasisCarriesNone)
{
    const auto bs=HighPool(Type::BSpline6, 18);
    EXPECT_EQ(dynamic_cast<const CoulombFitLadder*>(bs.get()), nullptr);
}
