// File: BasisSet/Gaussian/tests/M_BoysTable.C
//
// The Boys F_j(T) lookup table (AuxillaryFJ::theLookUp[121][23]) serves jmax <= 16.  Its bound is checked ONCE, at
// construction, and an over-range jmax THROWS (CLAUDE.md "Where to check, and with what"; OpenWork row "Fixed-size
// tables indexed by a runtime L are guarded by ASSERT only", 2026-11-04).  Before, an assert (dead in Release, and
// looser than the table: it let jmax reach 20) was the only guard, so four h shells would have read past the table.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>
import qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM;
import qchem.Types;
using namespace qchem;
using BasisSet::Gaussian::Evaluators::Internal::Boys::AuxillaryFJ;
using BasisSet::Gaussian::Evaluators::Internal::MnD::RNLM;

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
