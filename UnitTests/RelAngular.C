// File: UnitTests/RelAngular.C  Tests for relativistic Wigner 3j symbols.
#include "gtest/gtest.h"
#include "wignerSymbols/wignerSymbols-cpp.h"
#include <cmath>
#include <iostream>
import qchem.BasisSet.Atom.Evaluators.Internal.RelWigner3j;

using std::cout;
using std::endl;

class RelWigner3jTests : public ::testing::Test {};

// Check parity symbol against direct WignerSymbols call for all κ, k.
TEST_F(RelWigner3jTests, ParitySymbol)
{
    const int LMax=4, KMax=LMax+1;
    for (int κa=-(KMax); κa<=KMax; κa++)
    {
        if (κa==0) continue;
        double ja = κa>0 ? κa-0.5 : -κa-0.5;
        for (int κb=-(KMax); κb<=KMax; κb++)
        {
            if (κb==0) continue;
            double jb = κb>0 ? κb-0.5 : -κb-0.5;
            for (int k=0; k<=2*LMax; k++)
            {
                double expected = WignerSymbols::wigner3j(ja, k, jb, 0.5, 0.0, -0.5);
                double got      = RelWigner3j::w3j(κa, κb, k);
                EXPECT_NEAR(got, expected, 1e-14)
                    << "κa=" << κa << " κb=" << κb << " k=" << k;
            }
        }
    }
}

// Check m-dependent symbol against direct WignerSymbols call for all κ, k, mj.
TEST_F(RelWigner3jTests, MSymbol)
{
    const int LMax=4, KMax=LMax+1;
    for (int κa=-(KMax); κa<=KMax; κa++)
    {
        if (κa==0) continue;
        double ja = κa>0 ? κa-0.5 : -κa-0.5;
        for (int κb=-(KMax); κb<=KMax; κb++)
        {
            if (κb==0) continue;
            double jb = κb>0 ? κb-0.5 : -κb-0.5;
            for (int k=0; k<=2*LMax; k++)
            {
                for (double mja=-ja; mja<=ja; mja+=1.0)
                for (double mjb=-jb; mjb<=jb; mjb+=1.0)
                {
                    double expected = WignerSymbols::wigner3j(ja, k, jb, -mja, mja-mjb, mjb);
                    double got      = RelWigner3j::w3j(κa, κb, k, mja, mjb);
                    EXPECT_NEAR(got, expected, 1e-14)
                        << "κa=" << κa << " κb=" << κb << " k=" << k
                        << " mja=" << mja << " mjb=" << mjb;
                }
            }
        }
    }
}

// Completeness sum rule: Σ_k (2k+1) * (ja k jb / ½ 0 -½)^2 = 1  [for fixed ja==jb]
// Follows from completeness of CG coefficients; k runs up to ja+jb = 2*LMax+1.
TEST_F(RelWigner3jTests, SumRule)
{
    const int LMax=4, KMax=LMax+1, KkMax=2*LMax+1;
    for (int κ=-(KMax); κ<=KMax; κ++)
    {
        if (κ==0) continue;
        double j = κ>0 ? κ-0.5 : -κ-0.5;
        double sum=0.0;
        for (int k=0; k<=KkMax; k++)
        {
            double w = RelWigner3j::w3j(κ, κ, k);
            sum += (2*k+1)*w*w;
        }
        EXPECT_NEAR(sum, 1.0, 1e-13) << "κ=" << κ << " j=" << j;
    }
}
