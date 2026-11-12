// File: BasisSet/Gaussian/tests/M_BoysTable.C
//
// The Boys F_j(T) lookup table (AuxillaryFJ::theLookUp[121][23]) serves jmax <= 16.  Its bound is checked ONCE, at
// construction, and an over-range jmax THROWS (CLAUDE.md "Where to check, and with what"; OpenWork row "Fixed-size
// tables indexed by a runtime L are guarded by ASSERT only", 2026-11-04).  Before, an assert (dead in Release, and
// looser than the table: it let jmax reach 20) was the only guard, so four h shells would have read past the table.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>
#include <vector>
import qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM;
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Index3;
import qchem.Parallel;   // OuterFor + ScopedWorkerThreads (the serial-vs-threaded bit gate)
import qchem.Types;
using namespace qchem;
using BasisSet::Gaussian::Evaluators::Internal::Boys::AuxillaryFJ;
using BasisSet::Gaussian::Evaluators::Internal::MnD::RNLM;
using BasisSet::Gaussian::Evaluators::Internal::MnD::Index3;

TEST(M_BoysTable, RefusesJBeyondTheTable)
{
    const int jMax=AuxillaryFJ::thejMax;
    EXPECT_NO_THROW(AuxillaryFJ{jMax});
    EXPECT_THROW   (AuxillaryFJ{jMax+1}, std::out_of_range);
    EXPECT_THROW   (AuxillaryFJ{-1}    , std::out_of_range);
    // The evaluator's entry: an RNLM with Ltotal beyond the table is refused before the static Rj table is indexed.
    EXPECT_THROW   (RNLM(jMax+1, 1.0, rvec3_t(0.1,0.2,0.3)), std::out_of_range);
}

TEST(M_BoysTable, ValuesAtTheBound)
{
    // F_j(0) = 1/(2j+1) exactly; at large T, F_j(T) -> (2j-1)!!/2^(j+1) sqrt(π/T^(2j+1)).
    const int jMax=AuxillaryFJ::thejMax;
    AuxillaryFJ fj(jMax);
    rvec_t F(jMax+1,0.0);
    fj.GetFjAt(0.0,F);
    for (int j=0;j<=jMax;j++) EXPECT_NEAR(F[j],1.0/(2*j+1),1e-12) << "j=" << j;
    // At T = 5 (the Taylor-table branch) against the defining integral F_j(T) = ∫₀¹ u^(2j) exp(-T u²) du, Simpson.
    const double T=5.0;
    fj.GetFjAt(T,F);
    const int n=4000;
    for (int j=0;j<=jMax;j++)
    {
        double s=0.0;
        for (int i=0;i<=n;i++)
        {
            double u=double(i)/n, w=(i==0||i==n) ? 1.0 : (i%2 ? 4.0 : 2.0);
            s+=w*std::pow(u,2*j)*std::exp(-T*u*u);
        }
        s/=3.0*n;
        EXPECT_NEAR(F[j]/s,1.0,1e-9) << "j=" << j;
    }
}

// RNLM::Build re-fills ONE scratch RNLM in place (the cold ERI4 group kernel rebuilds it per primitive quartet instead of caching
// 3.9 M of them, OpenWork §4).  THE CLAIM: a scratch rebuilt through a sequence of L that grows AND shrinks is, element for element and
// bit for bit, the RNLM a fresh constructor gives -- no stale storage or reach survives a rebuild.
TEST(M_BoysTable, ARebuiltScratchRNLMIsAFreshOne)
{
    struct Arg {int L; double α; rvec3_t dR;};
    const Arg args[]={{8,0.7,{0.3,-1.1,0.4}}, {2,13.0,{0.0,0.0,0.0}}, {4,0.05,{2.0,1.5,-3.0}}, {0,1.0,{0.1,0.1,0.1}}, {6,2.5,{-0.4,0.9,0.0}}};
    RNLM scratch;
    for (const Arg& a : args)
    {
        scratch.Build(a.L, a.α, a.dR);
        const RNLM fresh(a.L, a.α, a.dR);
        EXPECT_TRUE (scratch.CheckLMax(a.L));
        EXPECT_FALSE(scratch.CheckLMax(a.L+1)) << "a rebuild must not keep a larger earlier reach";
        size_t differ=0;
        for (int n=0;n<=a.L;n++) for (int l=0;l<=a.L-n;l++) for (int m=0;m<=a.L-n-l;m++)
            if (scratch(Index3{n,l,m})!=fresh(Index3{n,l,m})) ++differ;
        EXPECT_EQ(differ, 0u) << "L=" << a.L;
    }
}

// RNLMs BUILT CONCURRENTLY are the serial ones, bit for bit.  The build recurses through a scratch of Triangle3D layers; while that
// scratch was function-STATIC (shared by every thread) two concurrent builds overwrote each other's layers.  No production path builds an
// RNLM inside a parallel region today (every MnD integral is built in a serial phase: the cold Fock pass, a refresh; OpenWork §3), so this
// gate drives the builder directly through qchem::OuterFor at 4 workers against the same list built serially.
TEST(M_BoysTable, RNLMsBuiltConcurrentlyAreTheSerialOnes)
{
    struct Arg {int L; double α; rvec3_t dR;};
    std::vector<Arg> args;
    for (int i=0;i<4000;i++)   // every L in [0,8], α across six decades, T = α|dR|² across the table's three Boys branches
        args.push_back({i%9, std::pow(10.0, -2.0+5.0*((i*37)%101)/100.0), rvec3_t(0.3+0.01*(i%17), -0.2*(i%5), 0.05*(i%23))});
    auto build=[&](std::vector<std::vector<double>>& out)
    {
        out.assign(args.size(), {});
        qchem::OuterFor(args.size(), [&](size_t i)
        {
            const Arg& a=args[i];
            RNLM r; r.Build(a.L, a.α, a.dR);   // half through Build (the ERI4 scratch), half through the constructor (the cache's entries)
            const RNLM c(a.L, a.α, a.dR);
            for (int n=0;n<=a.L;n++) for (int l=0;l<=a.L-n;l++) for (int m=0;m<=a.L-n-l;m++)
            {
                out[i].push_back(r(Index3{n,l,m}));
                out[i].push_back(c(Index3{n,l,m}));
            }
        });
    };
    std::vector<std::vector<double>> serial, threaded;
    { qchem::ScopedWorkerThreads one(1);  build(serial);   }
    { qchem::ScopedWorkerThreads four(4); build(threaded); }
    size_t differ=0;
    for (size_t i=0;i<args.size();i++) if (serial[i]!=threaded[i]) ++differ;
    EXPECT_EQ(differ, 0u) << differ << " of " << args.size() << " RNLMs built at 4 workers differ from the serial ones";
}
