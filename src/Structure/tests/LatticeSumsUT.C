// File: src/Structure/tests/LatticeSumsUT.C  Structure's two lattice-sum faces (the quotient range audit's bounds, lane H 2026-11-27).
//
// The oracle is POISSON SUMMATION, an identity independent of either loop:
//   (1/Omega) Sum_{all G} e^{-G^2/(4 alpha)} = (alpha/pi)^{3/2} Sum_{all R} e^{-alpha R^2},
// so the reciprocal sum (G != 0) must equal the image sum (R != 0, SumAroundSite on one atom) plus the two G = 0 / R = 0 terms.  A skewed,
// non-commensurate cell, so a shell ordering that only works for a cube cannot pass.
#include "gtest/gtest.h"
#include <cmath>
#include <optional>

import qchem.UnitCell;
import qchem.Types;

using namespace qchem;

TEST(LatticeSums, TheReciprocalSumAndTheImageSumAgreeByPoissonSummation)
{
    UnitCell cell(5.3, 6.1, 4.7, 77.0, 84.0, 101.0);
    cell.AddAtom(8, rvec3_t(0.13, 0.27, 0.41));
    const double Omega=cell.GetCellVolume();
    for (double alpha : {0.05, 0.3, 1.5})
    {
        const std::optional<double> recip=cell.SumOverReciprocalLattice([alpha](const rvec3_t& G){return std::exp(-(G*G)/(4.0*alpha));});
        ASSERT_TRUE(recip) << "a periodic structure has a reciprocal lattice";
        const double images=cell.SumAroundSite(0, [alpha](int, double d){return std::exp(-alpha*d*d);});
        const double poisson=std::pow(alpha/M_PI, 1.5)*(1.0+images) - 1.0/Omega;
        EXPECT_NEAR(*recip, poisson, 1e-12*std::max(1.0, std::fabs(poisson))) << "alpha " << alpha;
    }
}
