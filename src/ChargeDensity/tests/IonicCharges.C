// File: src/ChargeDensity/tests/IonicCharges.C  The IonicSAD electronegativity heuristic.
//
// IonicFormalCharges assigns charge-conserving integer formal charges from Pauling electronegativity + the
// per-atom valence count: electrons flow least- -> most-electronegative, capped by valence/octet.  These are
// the seed charges (Na+ + F-) that pre-bake the ~1 e- transfer the neutral SAD leaves for the SCF.
#include "gtest/gtest.h"
#include <vector>
#include <utility>
#include <numeric>

import qchem.ChargeDensity.Seed;   // IonicFormalCharges
using namespace qchem;

using namespace qchem::ChargeDensity;
using P = std::pair<int,int>;      // {Z, Nval}

static int Sum(const std::vector<int>& v){ return std::accumulate(v.begin(),v.end(),0); }

TEST(IonicFormalCharges, NaF)
{
    // Na (Z=11, Zion=1) + F (Z=9, Zion=7): F is far more electronegative -> Na+, F-.
    auto q = IonicFormalCharges({{11,1},{9,7}});
    ASSERT_EQ(q.size(), 2u);
    EXPECT_EQ(q[0], +1);          // Na+
    EXPECT_EQ(q[1], -1);          // F-
    EXPECT_EQ(Sum(q), 0);         // neutral cell (charge conserved)
}

TEST(IonicFormalCharges, CsI)
{
    // Cs (Z=55, Zion=1) + I (Z=53, Zion=7): I more electronegative -> Cs+, I-.
    auto q = IonicFormalCharges({{55,1},{53,7}});
    EXPECT_EQ(q[0], +1);          // Cs+
    EXPECT_EQ(q[1], -1);          // I-
    EXPECT_EQ(Sum(q), 0);
}

TEST(IonicFormalCharges, Na2O_multiplicity)
{
    // 2 Na (Zion=1) + 1 O (Zion=6): O fills its octet by taking 2 -> 2 Na+ + O2-.  Counts must balance.
    auto q = IonicFormalCharges({{11,1},{11,1},{8,6}});
    EXPECT_EQ(q[0], +1);
    EXPECT_EQ(q[1], +1);
    EXPECT_EQ(q[2], -2);          // O2- (accepted 2, one from each Na)
    EXPECT_EQ(Sum(q), 0);
}

TEST(IonicFormalCharges, SingleSpeciesNeutral)
{
    // Si-Si (same element): no electronegativity gap -> no transfer -> neutral, falls back to neutral SAD.
    auto q = IonicFormalCharges({{14,4},{14,4}});
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
}

TEST(IonicFormalCharges, NearCovalentNoTransfer)
{
    // C-Si: EN gap (2.55 vs 1.90 = 0.65) is modest; still ionic-ish but small.  Mainly a guard that a tiny
    // gap below the threshold (here use near-equal P and S, gap 0.39 < 0.5) yields NO transfer.
    auto q = IonicFormalCharges({{15,5},{16,6}});   // P (2.19) vs S (2.58): gap 0.39 < 0.5
    EXPECT_EQ(q[0], 0);
    EXPECT_EQ(q[1], 0);
}
