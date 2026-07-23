// File: Math/FFT.C  Complex FFT for the plane-wave transforms: radix-2 Cooley-Tukey (power-of-two N,
// the historical engine, kept BIT-IDENTICAL) + PocketFFT (submodules/pocketfft, header-only BSD-3 --
// the engine NumPy/SciPy ship) for ARBITRARY N.
//
// WHY two engines (doc/GPWPlan.md, mixed-radix increment 2026-07-22): the pow2 raster padding costs GPW
// dearly -- at NaF's Ecut=160 the alias-free raster is 128^3 where CP2K's mixed-radix engine runs 36^3
// (45x the points, multiplying through every N^3-scaled cost: pair boxes, stream caches, sweeps).  A
// finer menu of admissible N (5-smooth numbers 2^a 3^b 5^c) needs an arbitrary-N FFT; PocketFFT also
// positions the future FFT-rate-limited workloads (PW + ultrasoft PPs).  The pow2 path DISPATCHES to the
// original radix-2 so every existing grid result stays bit-identical until the AutoGrid policy flip
// (a separate, anchor-re-pinning increment) deliberately changes N.
//
// Conventions: the transforms are UNNORMALISED.  FFT1D/FFT3D with
//   sign=-1  compute  X_k = sum_j x_j e^{-2*pi*i*j.k/N}   (forward),
//   sign=+1  compute  x_j = sum_k X_k e^{+2*pi*i*j.k/N}   (inverse, caller divides by prod(N)).
// PocketFFT's c2c uses the same sign/normalisation convention (FORWARD=-1, no scaling at fct=1).
module;
#include <complex>   // std::complex free operators (* + -) for dcmplx -- not leaked by qchem.Types
#include <cmath>
#include <cassert>
#include "pocketfft_hdronly.h"   // submodules/pocketfft (header-only; include path on qcMath)
export module qchem.FFT;
import qchem.Types;   // cvec_t, ivec3_t
import qchem.CMath;    // Pi

export namespace qchem::FFT
{

//! True if \a n is a power of two (and non-zero).
inline bool IsPow2(size_t n) {return n!=0 && (n & (n-1))==0;}

//! Smallest power of two >= \a n (>=1).
inline size_t NextPow2(size_t n)
{
    size_t p=1;
    while (p<n) p<<=1;
    return p;
}

//! Smallest 5-SMOOTH number (\f$2^a3^b5^c\f$) >= \a n (>=1) -- the mixed-radix FFT grid menu (PocketFFT
//! handles any N, but 5-smooth sizes keep the transform in its fast radix kernels; CP2K restricts its
//! grids the same way).  The AutoGrid policy flip replaces \c NextPow2 with this (36 instead of 64,
//! 45 instead of 128 per axis -- the raster-volume lever).
inline size_t Next5Smooth(size_t n)
{
    if (n<=1) return 1;
    size_t best=NextPow2(n);                       // a valid 5-smooth bound to beat
    for (size_t p5=1; p5<best; p5*=5)
        for (size_t p35=p5; p35<best; p35*=3)
        {
            size_t s=p35;
            while (s<n) s<<=1;                     // s = p35 * 2^a >= n
            if (s<best) best=s;
        }
    return best;
}

//! In-place radix-2 Cooley-Tukey FFT of \a a (size MUST be a power of two).  \a sign = -1 forward,
//! +1 inverse; unnormalised (see the file header).  The HISTORICAL engine -- kept verbatim so every
//! pow2-grid result in the suite stays bit-identical; arbitrary-N callers go through FFT3D's dispatch.
inline void FFT1D(cvec_t& a, int sign)
{
    const size_t n=a.size();
    assert(IsPow2(n));
    // Decimation-in-time bit-reversal permutation.
    for (size_t i=1, j=0; i<n; ++i)
    {
        size_t bit=n>>1;
        for (; j & bit; bit>>=1) j^=bit;
        j^=bit;
        if (i<j) {dcmplx t=a[i]; a[i]=a[j]; a[j]=t;}
    }
    // Butterflies.
    for (size_t len=2; len<=n; len<<=1)
    {
        const double  ang=sign*2.0*Pi/double(len);
        const dcmplx  wlen(std::cos(ang), std::sin(ang));
        for (size_t i=0; i<n; i+=len)
        {
            dcmplx w(1.0,0.0);
            for (size_t k=0; k<len/2; ++k)
            {
                dcmplx u=a[i+k], v=a[i+k+len/2]*w;
                a[i+k]        =u+v;
                a[i+k+len/2]  =u-v;
                w*=wlen;
            }
        }
    }
}

//! 3-D FFT of a grid stored row-major as a flat vector (index = (i0*n1 + i1)*n2 + i2).  \a sign /
//! normalisation as FFT1D.  DISPATCH: all-pow2 \a n runs the historical per-axis radix-2 (bit-identical
//! to every committed grid anchor); any other \a n runs PocketFFT's mixed-radix c2c (any N, all three
//! axes in one call).
inline cvec_t FFT3D(cvec_t data, const ivec3_t& n, int sign)
{
    assert(data.size()==size_t(n.x)*n.y*n.z);
    const size_t n0=n.x, n1=n.y, n2=n.z;
    if (!(IsPow2(n0) && IsPow2(n1) && IsPow2(n2)))
    {
        // PocketFFT mixed-radix path.  Row-major C layout; strides in BYTES; unnormalised (fct=1),
        // matching the radix-2 convention exactly (verified by the pow2 cross-check in the unit tests).
        const pocketfft::shape_t  shape{n0,n1,n2};
        const ptrdiff_t cs=ptrdiff_t(sizeof(dcmplx));
        const pocketfft::stride_t stride{ptrdiff_t(n1*n2)*cs, ptrdiff_t(n2)*cs, cs};
        const pocketfft::shape_t  axes{0,1,2};
        pocketfft::c2c(shape, stride, stride, axes, sign==-1 /*==pocketfft::FORWARD*/,
                       data.data(), data.data(), 1.0);
        return data;
    }
    // Axis 2 (innermost, contiguous stride 1).
    for (size_t i0=0; i0<n0; ++i0)
        for (size_t i1=0; i1<n1; ++i1)
        {
            cvec_t line(n2);
            const size_t base=(i0*n1+i1)*n2;
            for (size_t i2=0; i2<n2; ++i2) line[i2]=data[base+i2];
            FFT1D(line, sign);
            for (size_t i2=0; i2<n2; ++i2) data[base+i2]=line[i2];
        }
    // Axis 1 (stride n2).
    for (size_t i0=0; i0<n0; ++i0)
        for (size_t i2=0; i2<n2; ++i2)
        {
            cvec_t line(n1);
            const size_t base=i0*n1*n2+i2;
            for (size_t i1=0; i1<n1; ++i1) line[i1]=data[base+i1*n2];
            FFT1D(line, sign);
            for (size_t i1=0; i1<n1; ++i1) data[base+i1*n2]=line[i1];
        }
    // Axis 0 (outermost, stride n1*n2).
    for (size_t i1=0; i1<n1; ++i1)
        for (size_t i2=0; i2<n2; ++i2)
        {
            cvec_t line(n0);
            const size_t base=i1*n2+i2, stride=n1*n2;
            for (size_t i0=0; i0<n0; ++i0) line[i0]=data[base+i0*stride];
            FFT1D(line, sign);
            for (size_t i0=0; i0<n0; ++i0) data[base+i0*stride]=line[i0];
        }
    return data;
}

} //namespace
