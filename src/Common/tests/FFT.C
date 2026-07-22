// File: src/Common/tests/FFT.C  Tests for the radix-2 Cooley-Tukey FFT (qchem.FFT).
#include "gtest/gtest.h"
#include <complex>
#include <cmath>
import qchem.FFT;
import qchem.Types;
import qchem.Blaze;   // scalar /= on cvec_t (the blaze operators live in Blaze's global fragment)
import qchem.Math;
using namespace qchem;

using qchem::FFT::FFT1D;
using qchem::FFT::FFT3D;
using qchem::FFT::IsPow2;
using qchem::FFT::NextPow2;

namespace
{
// O(N^2) reference 1-D DFT (unnormalised), same sign convention as FFT1D.
cvec_t DirectDFT1(const cvec_t& x, int sign)
{
    size_t n=x.size();
    cvec_t X(n);
    for (size_t k=0; k<n; ++k)
    {
        dcmplx s(0,0);
        for (size_t j=0; j<n; ++j)
        {
            double a=sign*2.0*Pi*double(j*k)/double(n);
            s += x[j]*dcmplx(std::cos(a),std::sin(a));
        }
        X[k]=s;
    }
    return X;
}

cvec_t Ramp(size_t n)   // a deterministic, non-symmetric test signal
{
    cvec_t x(n);
    for (size_t i=0; i<n; ++i) x[i]=dcmplx(0.3*double(i)-1.0, 0.7-0.2*double(i));
    return x;
}

double MaxAbsDiff(const cvec_t& a, const cvec_t& b)
{
    double m=0;
    for (size_t i=0; i<a.size(); ++i) m=std::max(m, std::abs(a[i]-b[i]));
    return m;
}
} //anon

TEST(FFT, Pow2Helpers)
{
    EXPECT_TRUE (IsPow2(1)); EXPECT_TRUE(IsPow2(2)); EXPECT_TRUE(IsPow2(16));
    EXPECT_FALSE(IsPow2(0)); EXPECT_FALSE(IsPow2(3)); EXPECT_FALSE(IsPow2(13));
    EXPECT_EQ(NextPow2(1), 1u); EXPECT_EQ(NextPow2(13), 16u); EXPECT_EQ(NextPow2(16), 16u);
    EXPECT_EQ(NextPow2(17), 32u);
}

TEST(FFT, Forward1DMatchesDirectDFT)
{
    for (size_t n : {2u,4u,8u,16u,32u})
    {
        cvec_t x=Ramp(n), X=x;
        FFT1D(X,-1);
        EXPECT_LT(MaxAbsDiff(X, DirectDFT1(x,-1)), 1e-10) << "n="<<n;
    }
}

TEST(FFT, RoundTrip1D)
{
    for (size_t n : {2u,8u,64u})
    {
        cvec_t x=Ramp(n), y=x;
        FFT1D(y,-1);          // forward
        FFT1D(y,+1);          // inverse (unnormalised)
        y/=double(n);         // normalise
        EXPECT_LT(MaxAbsDiff(y,x), 1e-10) << "n="<<n;
    }
}

TEST(FFT, RoundTrip3D)
{
    ivec3_t n(4,8,2);
    size_t N=size_t(n.x)*n.y*n.z;
    cvec_t x=Ramp(N);
    cvec_t y=FFT3D(x, n, -1);
    y=FFT3D(y, n, +1);
    y/=double(N);
    EXPECT_LT(MaxAbsDiff(y,x), 1e-10);
}

// 3-D forward FFT == the separable product of 1-D DFTs (checked element-wise against a triple direct sum).
TEST(FFT, Forward3DMatchesDirect)
{
    ivec3_t n(4,4,2);
    size_t n0=n.x, n1=n.y, n2=n.z, N=n0*n1*n2;
    cvec_t x=Ramp(N);
    cvec_t X=FFT3D(x, n, -1);

    for (size_t k0=0;k0<n0;++k0) for (size_t k1=0;k1<n1;++k1) for (size_t k2=0;k2<n2;++k2)
    {
        dcmplx s(0,0);
        for (size_t j0=0;j0<n0;++j0) for (size_t j1=0;j1<n1;++j1) for (size_t j2=0;j2<n2;++j2)
        {
            double a=-2.0*Pi*(double(j0*k0)/n0 + double(j1*k1)/n1 + double(j2*k2)/n2);
            s += x[(j0*n1+j1)*n2+j2]*dcmplx(std::cos(a),std::sin(a));
        }
        EXPECT_LT(std::abs(X[(k0*n1+k1)*n2+k2]-s), 1e-10);
    }
}

// ---------------- Mixed-radix (PocketFFT) path: arbitrary N (doc/GPWPlan.md mixed-radix increment) -----

TEST(FFT, Next5Smooth)
{
    using qchem::FFT::Next5Smooth;
    EXPECT_EQ(Next5Smooth(1),   1u);
    EXPECT_EQ(Next5Smooth(16), 16u);   // pow2 stays itself
    EXPECT_EQ(Next5Smooth(17), 18u);   // 2*3^2
    EXPECT_EQ(Next5Smooth(33), 36u);   // 2^2*3^2  (the CP2K NaF raster)
    EXPECT_EQ(Next5Smooth(37), 40u);   // 2^3*5
    EXPECT_EQ(Next5Smooth(41), 45u);   // 3^2*5
    EXPECT_EQ(Next5Smooth(129),135u);  // 3^3*5  (vs pow2 padding's 256 -- the raster-volume lever)
    EXPECT_EQ(Next5Smooth(73), 75u);   // 3*5^2
}

// The PocketFFT path (any non-pow2 axis) against the same brute-force reference the radix-2 path is held
// to -- the two engines are cross-validated through one oracle, and the convention (sign, unnormalised)
// is pinned identical.
TEST(FFT, Forward3DMixedRadixMatchesDirect)
{
    ivec3_t n(6,10,9);                 // 2*3, 2*5, 3^2 -- all three axes non-pow2
    size_t n0=n.x, n1=n.y, n2=n.z, N=n0*n1*n2;
    cvec_t x=Ramp(N);
    cvec_t X=FFT3D(x, n, -1);
    for (size_t k0=0;k0<n0;++k0) for (size_t k1=0;k1<n1;++k1) for (size_t k2=0;k2<n2;++k2)
    {
        dcmplx s(0,0);
        for (size_t j0=0;j0<n0;++j0) for (size_t j1=0;j1<n1;++j1) for (size_t j2=0;j2<n2;++j2)
        {
            double a=-2.0*Pi*(double(j0*k0)/n0 + double(j1*k1)/n1 + double(j2*k2)/n2);
            s += x[(j0*n1+j1)*n2+j2]*dcmplx(std::cos(a),std::sin(a));
        }
        EXPECT_LT(std::abs(X[(k0*n1+k1)*n2+k2]-s), 1e-9);
    }
}

TEST(FFT, RoundTrip3DMixedRadix)
{
    for (ivec3_t n : {ivec3_t(36,36,36), ivec3_t(12,45,8), ivec3_t(6,10,15)})
    {
        size_t N=size_t(n.x)*n.y*n.z;
        cvec_t x=Ramp(N);
        cvec_t y=FFT3D(x, n, -1);
        y=FFT3D(y, n, +1);
        y/=double(N);
        EXPECT_LT(MaxAbsDiff(y,x), 1e-9) << "n=("<<n.x<<","<<n.y<<","<<n.z<<")";
    }
}

// A MIXED shape with pow2 axes present (only one axis non-pow2) still routes the WHOLE call through
// PocketFFT -- pin that the per-axis convention agrees with the radix-2 engine by comparing a pow2-only
// transform (radix-2 path) against the same data zero-extended... simplest honest check: a (8,12,8)
// shape against the direct sum (covers pocketfft's pow2 kernels alongside radix-3).
TEST(FFT, Forward3DMixedPow2AxesMatchesDirect)
{
    ivec3_t n(8,12,8);
    size_t n0=n.x, n1=n.y, n2=n.z, N=n0*n1*n2;
    cvec_t x=Ramp(N);
    cvec_t X=FFT3D(x, n, -1);
    for (size_t k0=0;k0<n0;++k0) for (size_t k1=0;k1<n1;++k1) for (size_t k2=0;k2<n2;++k2)
    {
        dcmplx s(0,0);
        for (size_t j0=0;j0<n0;++j0) for (size_t j1=0;j1<n1;++j1) for (size_t j2=0;j2<n2;++j2)
        {
            double a=-2.0*Pi*(double(j0*k0)/n0 + double(j1*k1)/n1 + double(j2*k2)/n2);
            s += x[(j0*n1+j1)*n2+j2]*dcmplx(std::cos(a),std::sin(a));
        }
        EXPECT_LT(std::abs(X[(k0*n1+k1)*n2+k2]-s), 1e-9);
    }
}
