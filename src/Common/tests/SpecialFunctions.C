// File: src/Common/tests/SpecialFunctions.C  qchem.SpecialFunctions::SphericalBessel against a hard-coded reference table.
//
// The table is the project's own oracle for j_l (user, 2026-11-24: "UT our own SphericalBessel against a hard coded table of known
// values"); libc++ ships no std::sph_bessel, so no test may lean on the standard library's.
//
// SOURCE of the table: the power series j_l(x) = x^l * sum_k (-x^2/2)^k / (k! (2l+2k+1)!!), summed in Python's decimal module at
// 500 significant digits, at the EXACT binary value of each double x (Decimal(float)), until a term falls below 1e-60 of the sum
// with k > x.  The series converges for every x; 500 digits absorb its cancellation (terms reach ~1e40 at x=100).  It shares nothing
// with the downward recurrence under test.  Generator, verbatim:
//     from decimal import Decimal, getcontext;  getcontext().prec=500
//     def jl(l,xf):
//         x=Decimal(xf); df=Decimal(1)
//         for m in range(1,2*l+2,2): df*=m
//         term=x**l/df; s=term; k=0; y=-(x*x)/2
//         while True:
//             k+=1; term=term*y/(k*(2*l+2*k+1)); s+=term
//             if abs(term)<abs(s)*Decimal(10)**-60 and k>x: break
//         return s
//     print("%.17e" % float(jl(l,x)))
// Cross-checked against scipy.special.spherical_jn (pyscf-env): worst relative difference 5.7e-15, over x <= 50.
//
// The x are chosen to exercise each regime: the series regime (1e-3, 0.1), l > x and l < x, the zeros of j_0 (pi, 2pi) and of j_1
// (4.4934...), where a normalisation on j_0 alone fails, and large x (33.3, 50, 100), where the recurrence must start well above x.
#include "gtest/gtest.h"
#include <cmath>
#include <vector>
import qchem.SpecialFunctions;
import qchem.Types;

namespace
{
struct Row { double x; double j[11]; };   // j_0(x) .. j_10(x)
const std::vector<Row> theTable=
{
    {0.001, {9.99999833333341637e-01, 3.33333300000001188e-04, 6.66666619047620395e-08, 9.52380899470900688e-12, 1.05820101010101121e-15, 9.62000925000925732e-20, 7.40000715334049159e-24, 4.93333812157341867e-28, 2.90196360990994783e-32, 1.52734927220159963e-36, 7.27309178744673128e-41}},
    {0.1, {9.98334166468281548e-01, 3.33000119025575747e-02, 6.66190608445568766e-04, 9.51851972086556938e-06, 1.05772015020987336e-07, 9.61631023291644802e-10, 7.39754109358770945e-12, 4.93188747573197526e-14, 2.90120010253019116e-16, 1.52698569349482037e-18, 7.27151099671367484e-21}},
    {0.5, {9.58851077208406011e-01, 1.62537030636066560e-01, 1.63711066079934124e-02, 1.17403544386755722e-03, 6.53896061523897090e-05, 2.97746687545744574e-06, 1.14665107674094210e-07, 3.82592406900384716e-09, 1.12614396021212883e-10, 2.96539571739077667e-12, 7.06412396366187782e-14}},
    {1.0, {8.41470984807896505e-01, 3.01168678939756795e-01, 6.20350520113738604e-02, 9.00658111711251548e-03, 1.01101580841375270e-03, 9.25611586112581576e-05, 7.15693631008708582e-06, 4.79013419873948862e-07, 2.82649880221472958e-08, 1.49137650255514560e-09, 7.11655264004731351e-11}},
    {2.5, {2.39388857641582603e-01, 4.16212989275406509e-01, 2.60066729488905246e-01, 1.03920469702403942e-01, 3.09105856778257998e-02, 7.35763873776893590e-03, 1.46302476835752067e-03, 2.50090057690171717e-04, 3.75155777835096512e-05, 5.01587123769388205e-06, 6.05043622963854017e-07}},
    {3.141592653589793, {3.89817183251937549e-17, 3.18309886183790691e-01, 3.03963550927013326e-01, 1.65463130314201651e-01, 6.47163003184774710e-02, 1.99354133832935765e-02, 5.08573049721547897e-03, 1.10948446119762370e-03, 2.11667591282020651e-04, 3.59056159285872631e-05, 5.48554659014209661e-06}},
    {4.493409457909064, {-2.17233628211221663e-01, -7.21830072942819814e-18, 2.17233628211221663e-01, 2.41724719554389150e-01, 1.59334110102148735e-01, 7.74108954729203769e-02, 3.01700644402632424e-02, 9.87490418485911965e-03, 2.79456167725370590e-03, 6.97813249046818902e-04, 1.56082339519064678e-04}},
    {5.0, {-1.91784854932627702e-01, -9.50894080791707952e-02, 1.34731210085125230e-01, 2.29820618164295998e-01, 1.87017655344889205e-01, 1.06811161456504544e-01, 4.79668998594207988e-02, 1.79027781779895305e-02, 5.74143467454779101e-03, 1.61809971547296105e-03, 4.07344244249460414e-04}},
    {6.283185307179586, {-3.89817183251937549e-17, -1.59154943091895346e-01, -7.59908877317533038e-02, 9.86833160296463341e-02, 1.85932450559478490e-01, 1.67645301259927659e-01, 1.07564911638643318e-01, 5.49080348170017177e-02, 2.35183657101724203e-02, 8.72405583868514703e-03, 2.86268989002016333e-03}},
    {10.0, {-5.44021110889369788e-02, 7.84669417987515494e-02, 7.79421936285624450e-02, -3.94958449844703269e-02, -1.05589285117691667e-01, -5.55345116214521831e-02, 4.45013223340942760e-02, 1.13386230655774731e-01, 1.25578023649567827e-01, 1.00096409548490570e-01, 6.46051544925642651e-02}},
    {20.0, {4.56472625363813847e-02, -1.81217399638505285e-02, -4.83655235309589648e-02, 6.03035908111078985e-03, 5.04761492093477385e-02, 1.66839080630956926e-02, -4.12999997746451089e-02, -4.35289079166150117e-02, 8.65331883718384841e-03, 5.08842289282212812e-02, 3.96866986446263725e-02}},
    {33.3, {2.85684305941976795e-02, 1.01125030394986441e-02, -2.76573942843329375e-02, -1.42652649440531392e-02, 2.46586899417391842e-02, 2.09297757391177855e-02, -1.77449502080966742e-02, -2.78572337783146856e-02, 5.19664670435131780e-03, 3.05101765402958393e-02, 1.22115621324541158e-02}},
    {50.0, {-5.24749707407857601e-03, -1.94042705113238387e-02, 4.08324084339914579e-03, 1.98125945956637503e-02, -1.30947760000622035e-03, -2.00483005636648698e-02, -3.10114852400005143e-03, 1.92420019474248589e-02, 8.87374910822750877e-03, -1.62249272506275050e-02, -1.50392214634659602e-02}},
    {100.0, {-5.06365641109758798e-03, -8.67382528698781505e-03, 4.80344165248795373e-03, 8.91399736961221226e-03, -4.17946183661509815e-03, -9.29014893490757131e-03, 3.15754545377526551e-03, 9.70062984389835625e-03, -1.70245097719051220e-03, -9.99004651002074334e-03, -1.95657859713429014e-04}},
};

// |delta| <= 1e-13 |j| + 1e-16.  The absolute floor is the conditioning limit: |j_l| <= 1, and near a zero of j_l its double value
// is fixed only to ~1e-16 absolute by the rounding of x itself.
::testing::AssertionResult Near(double got, double want, int l, double x, int lmax)
{
    const double tol=1e-13*std::fabs(want)+1e-16;
    if (std::fabs(got-want)<=tol) return ::testing::AssertionSuccess();
    return ::testing::AssertionFailure() << "j_" << l << "(" << x << ") with lmax=" << lmax << ": got " << got << ", want " << want
                                         << ", |delta|=" << std::fabs(got-want) << " > " << tol;
}
}

// Every lmax from 0 to 10, not just the largest: the recurrence's starting index depends on lmax, and a small lmax at large x was
// the case that went wrong (j_2(50) off by 4e-10 with lmax=3).
TEST(SpecialFunctions, SphericalBesselMatchesTheReferenceTable)
{
    for (const Row& r : theTable)
        for (int lmax=0; lmax<=10; ++lmax)
        {
            const rvec_t j=qchem::SpecialFunctions::SphericalBessel(lmax, r.x);
            ASSERT_EQ(j.size(), size_t(lmax+1));
            for (int l=0; l<=lmax; ++l) EXPECT_TRUE(Near(j[l], r.j[l], l, r.x, lmax));
        }
}

// Tiny x with a large lmax: the trial values grow like prod (2l+1)/x and overflow a double without rescaling.
TEST(SpecialFunctions, SphericalBesselStaysFiniteAtTinyX)
{
    const rvec_t j=qchem::SpecialFunctions::SphericalBessel(10, 1e-11);
    for (size_t l=0; l<j.size(); ++l) EXPECT_TRUE(std::isfinite(j[l])) << "l=" << l;
    EXPECT_TRUE(Near(j[0], 1.00000000000000000e+00, 0, 1e-11, 10));
    EXPECT_TRUE(Near(j[1], 3.33333333333333313e-12, 1, 1e-11, 10));
    EXPECT_NEAR(j[2]/6.66666666666666542e-24, 1.0, 1e-13);   // below the absolute floor: check it relatively
}

// j(0) = delta_{l0}, the documented special case.
TEST(SpecialFunctions, SphericalBesselAtZeroIsTheDelta)
{
    const rvec_t j=qchem::SpecialFunctions::SphericalBessel(4, 0.0);
    EXPECT_EQ(j[0], 1.0);
    for (size_t l=1; l<j.size(); ++l) EXPECT_EQ(j[l], 0.0);
}
