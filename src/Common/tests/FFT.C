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
