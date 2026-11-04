// File: src/Common/tests/NUFFT.C  The type-2 NUFFT (qchem.NUFFT) against the direct series it replaces.
#include "gtest/gtest.h"
#include <complex>
#include <cmath>
#include <random>
#include <vector>
#include <utility>
#include <iostream>
import qchem.NUFFT;
import qchem.Types;
using namespace qchem;

namespace
{
using Terms=std::vector<std::pair<ivec3_t,dcmplx>>;
// A deterministic random coefficient set on the box [lo, hi] (every k present) and its l1 norm.
Terms RandomTerms(const ivec3_t& lo, const ivec3_t& hi, double& l1)
{
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> u(-1.0,1.0);
    Terms t; l1=0.0;
    for (int x=lo.x; x<=hi.x; x++) for (int y=lo.y; y<=hi.y; y++) for (int z=lo.z; z<=hi.z; z++)
    {
        const dcmplx c(u(gen),u(gen));
        t.push_back({ivec3_t(x,y,z),c}); l1+=std::abs(c);
    }
    return t;
}
dcmplx Direct(const Terms& t, const rvec3_t& s)
{
    dcmplx f(0.0);
    for (const auto& [k,c] : t) f+=c*std::exp(dcmplx(0.0, 2.0*M_PI*(k.x*s.x+k.y*s.y+k.z*s.z)));
    return f;
}
// The max error over a point set, relative to ||c||_1 (the norm the width rule is stated against).
double MaxRelError(double tol, const ivec3_t& lo, const ivec3_t& hi)
{
    double l1=0.0;
    const Terms t=RandomTerms(lo,hi,l1);
    const NUFFT::Type2 plan(lo,hi,tol);
    const cvec_t g=plan.Grid(t);
    std::mt19937 gen(777);
    std::uniform_real_distribution<double> u(-1.5,2.5);   // points outside [0,1): the series is periodic
    double err=0.0;
    for (int p=0; p<200; p++)
    {
        const rvec3_t s(u(gen),u(gen),u(gen));
        err=std::max(err, std::abs(plan.Interpolate(g,s)-Direct(t,s)));
    }
    return err/l1;
}
}

// The width rule delivers its tolerance at every width, on an ASYMMETRIC box (lo != -hi, odd and even extents),
// against ||c||_1 -- the box the grid sizing once got wrong (sized from the extent, not max|k|: 12x short at w 8).
// The table is printed so a record can cite the delivered error per width.
TEST(NUFFT, Type2DeliversItsToleranceAgainstTheDirectSeries)
{
    const ivec3_t lo(-5,-3,-7), hi(4,6,7);
    for (double tol : {1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13})
    {
        const double e=MaxRelError(tol,lo,hi);
        std::cout << "[NUFFT width] tol " << tol << "  w " << NUFFT::Type2::WidthFor(tol) << "  delivered " << e << std::endl;
        EXPECT_LT(e, tol) << "tol " << tol;
    }
}
// A one-term box: the plan reproduces a single plane wave (the deconvolution of the m = 0 image is exact).
TEST(NUFFT, Type2ReproducesOnePlaneWave)
{
    const ivec3_t k(3,-2,5);
    const NUFFT::Type2 plan(k,k,1e-13);
    const cvec_t g=plan.Grid(Terms{{k,dcmplx(0.25,-0.5)}});
    for (const rvec3_t s : {rvec3_t(0.1,0.2,0.3), rvec3_t(-0.77,1.31,0.05), rvec3_t(0.999,0.0,0.5)})
        EXPECT_NEAR(std::abs(plan.Interpolate(g,s)-Direct(Terms{{k,dcmplx(0.25,-0.5)}},s)), 0.0, 1e-13);
}
