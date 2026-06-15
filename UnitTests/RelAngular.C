// File: UnitTests/RelAngular.C  Tests for relativistic Wigner 3j symbols.
#include "gtest/gtest.h"
#include "wignerSymbols/wignerSymbols-cpp.h"
#include <blaze/Math.h>
#include <cmath>
#include <iostream>
import qchem.BasisSet.Atom.Evaluators.Internal.RelWigner3j;
import qchem.BasisSet.Atom.Evaluators.Internal.RelAngularIntegrals;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;

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

class RelAngularIntegralsTests : public ::testing::Test {};

// For s1/2 (κ=-1, l=0): ml=0 always, all CG=1.
// Direct : RelDirect  == NR_Direct (0,0,0,0) for ALL (mja,mjc) pairs.
// Exchange: RelExchange includes the spin δ, so it is non-zero only when mja==mjc
//           (same-spin), where it equals NR_Exchange(0,0,0,0); zero otherwise.
TEST_F(RelAngularIntegralsTests, S12ExactlyNR)
{
    rvec11_t nr_coulomb  = AngularIntegrals::Direct  (0, 0, 0, 0);
    rvec11_t nr_exchange = AngularIntegrals::Exchange(0, 0, 0, 0);
    for (double mja : {-0.5, 0.5})
    for (double mjc : {-0.5, 0.5})
    {
        rvec11_t rc = RelAngularIntegrals::Direct  (-1,-1,mja,mjc);
        rvec11_t re = RelAngularIntegrals::Exchange(-1,-1,mja,mjc);
        rvec11_t re_expected = (mja==mjc) ? nr_exchange : rvec11_t(0.0);
        for (size_t k=0; k<rc.size(); k++)
        {
            EXPECT_NEAR(rc[k], nr_coulomb[k],    1e-13) << "Direct   mja=" << mja << " mjc=" << mjc << " k=" << k;
            EXPECT_NEAR(re[k], re_expected[k],   1e-13) << "Exchange mja=" << mja << " mjc=" << mjc << " k=" << k;
        }
    }
}

// Hermitian symmetry: Direct (κa,κc,mja,mjc) == Direct (κc,κa,mjc,mja)
TEST_F(RelAngularIntegralsTests, DirectSymmetry)
{
    const int LMax=3, KMax=LMax+1;
    for (int κa=-(KMax); κa<=KMax; κa++) { if (κa==0) continue;
    for (int κc=-(KMax); κc<=KMax; κc++) { if (κc==0) continue;
        double ja=κa>0?κa-0.5:-κa-0.5, jc=κc>0?κc-0.5:-κc-0.5;
        for (double mja=-ja; mja<=ja; mja+=1.0)
        for (double mjc=-jc; mjc<=jc; mjc+=1.0)
        {
            rvec11_t ac=RelAngularIntegrals::Direct (κa,κc,mja,mjc);
            rvec11_t ca=RelAngularIntegrals::Direct (κc,κa,mjc,mja);
            for (size_t k=0; k<ac.size(); k++)
                EXPECT_NEAR(ac[k], ca[k], 1e-13)
                    << "κa=" << κa << " κc=" << κc << " mja=" << mja << " mjc=" << mjc << " k=" << k;
        }
    }}
}

// Full-shell sum rule via CG completeness:
// Summing RelDirect  over ALL (κa,κb) pairs in the l-shell (including cross terms)
// and over all mja,mjb gives exactly 4 × NR sum, because the independent CG sums
// for a and c each collapse to 1 (completeness), and the two independent ms sums
// each contribute a factor of 2.
// Σ_{κa,κb∈{κm,κp}} Σ_{mja,mjc} RelDirect (κa,κb,mja,mjc) = 4 × Σ_{ma,mc} NR_Direct (l,l,ma,mc)
TEST_F(RelAngularIntegralsTests, CombinedShellsEqualNR)
{
    struct LShell { int l, κminus, κplus; }; // κminus=l (j=l-½), κplus=-(l+1) (j=l+½)
    for (auto [l,κm,κp] : std::initializer_list<LShell>{{1,1,-2},{2,2,-3},{3,3,-4}})
    {
        rvec11_t Ak_rel(0.0);
        for (int κa : {κm, κp})
        for (int κb : {κm, κp})  // all cross-κ pairs
        {
            double ja=κa>0?κa-0.5:-κa-0.5;
            double jb=κb>0?κb-0.5:-κb-0.5;
            for (double mja=-ja; mja<=ja; mja+=1.0)
            for (double mjc=-jb; mjc<=jb; mjc+=1.0)
                Ak_rel += RelAngularIntegrals::Direct (κa,κb,mja,mjc);
        }
        rvec11_t Ak_nr(0.0);
        for (int ma=-l; ma<=l; ma++)
        for (int mc=-l; mc<=l; mc++)
            Ak_nr += AngularIntegrals::Direct (l,l,ma,mc);

        for (size_t k=0; k<Ak_rel.size(); k++)
            EXPECT_NEAR(Ak_rel[k], 4.0*Ak_nr[k], 1e-8) << "l=" << l << " k=" << k;
    }
}
