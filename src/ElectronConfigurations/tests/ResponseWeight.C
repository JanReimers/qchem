// File: ElectronConfigurations/tests/ResponseWeight.C  The occupancy rules' RESPONSE WEIGHTS
// (doc/LinearResponsePlan.md §3 row E1, stage R0).
//
// \f$(f_n-f_m)/(\varepsilon_n-\varepsilon_m)\f$ is the level-pair weight of every independent-particle response.
// What these pin: (1) the integer rule's zero between equally-occupied levels and its REFUSAL of an unresolved
// degenerate pair (a broken invariant: the caller gates on the response gap first); (2) the Fermi rule's
// degenerate branch is the derivative f'(ε) the quotient tends to, so the weight is CONTINUOUS through the
// branch switch; (3) the factory builds the rule the policy uses from the same configuration value.
#include "gtest/gtest.h"
#include <cmath>
#include <memory>
#include <stdexcept>
import qchem.ElectronConfiguration.OccupationPolicy;

using namespace qchem;

namespace {
double Fermi(double e, double mu, double kT) {return 1.0/(1.0+std::exp((e-mu)/kT));}
}

TEST(ResponseWeight, IntegerZeroBetweenEqualOccupations)
{
    IntegerOccupancy r;
    EXPECT_EQ(r.ResponseWeight(-0.5, 1.0, -0.3, 1.0), 0.0);   // two full levels
    EXPECT_EQ(r.ResponseWeight( 0.2, 0.0,  0.4, 0.0), 0.0);   // two empty levels
    EXPECT_EQ(r.ResponseWeight(-0.1, 1.0, -0.1, 1.0), 0.0);   // the n = m diagonal of a full level
    EXPECT_EQ(r.ResponseWeight( 0.1, 0.5,  0.1, 0.5), 0.0);   // inside one partially-filled degenerate level
}

TEST(ResponseWeight, IntegerQuotientIsNegativeAcrossAGap)
{
    IntegerOccupancy r;
    // occupied at -0.2, empty at +0.3: (1-0)/(-0.5) -- the static response is negative-definite.
    EXPECT_DOUBLE_EQ(r.ResponseWeight(-0.2, 1.0, 0.3, 0.0), -2.0);
    EXPECT_DOUBLE_EQ(r.ResponseWeight( 0.3, 0.0,-0.2, 1.0), -2.0);   // symmetric in the pair
    EXPECT_TRUE(r.RequiresResolvedGap());
}

TEST(ResponseWeight, IntegerThrowsOnAnUnresolvedDegeneratePair)
{
    IntegerOccupancy r;
    EXPECT_THROW(r.ResponseWeight(0.1, 1.0, 0.1, 0.0), std::logic_error);
}

TEST(ResponseWeight, FermiIsTheQuotientAwayFromDegeneracy)
{
    const double kT=0.01, mu=0.0;
    FermiOccupancy r(kT);
    const double en=-0.013, em=0.021;
    EXPECT_DOUBLE_EQ(r.ResponseWeight(en,Fermi(en,mu,kT),em,Fermi(em,mu,kT)),
                     (Fermi(en,mu,kT)-Fermi(em,mu,kT))/(en-em));
    EXPECT_FALSE(r.RequiresResolvedGap());
}

TEST(ResponseWeight, FermiDegenerateBranchIsTheDerivativeAndContinuous)
{
    const double kT=0.01, mu=0.0;
    FermiOccupancy r(kT);
    for (double e : {-0.03, -0.005, 0.0, 0.008, 0.04})
    {
        const double f=Fermi(e,mu,kT);
        const double exact=-f*(1.0-f)/kT;                         // f'(ε)
        EXPECT_NEAR(r.ResponseWeight(e,f,e,f), exact, 1e-14) << "e=" << e;
        // Just outside the 1e-6 kT branch switch the quotient must already agree to ~1e-6 relative.
        const double d=3e-6*kT, fd=Fermi(e+d,mu,kT);
        EXPECT_NEAR(r.ResponseWeight(e,f,e+d,fd), exact, 1e-5*std::fabs(exact)+1e-12) << "e=" << e;
        EXPECT_LE(std::fabs(r.ResponseWeight(e,f,e,f)), 0.25/kT+1e-12);   // bounded by 1/4kT
    }
}

TEST(ResponseWeight, FactoryBuildsTheRuleThePolicyUses)
{
    auto integer=MakeOccupancyRule(OccupationConfig{});
    auto fermi  =MakeOccupancyRule(OccupationConfig{.kT=0.02});
    EXPECT_TRUE (integer->RequiresResolvedGap());
    EXPECT_EQ   (integer->kT(), 0.0);
    EXPECT_FALSE(fermi->RequiresResolvedGap());
    EXPECT_EQ   (fermi->kT(), 0.02);
}
