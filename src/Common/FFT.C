// File: Common/FFT.C  Self-contained complex FFT (radix-2 Cooley-Tukey) for the plane-wave transforms.
//
// Lives in Common for now (a later pass reorganises the Common math: Math / SpecialFunctions / Blaze /
// FFT).  Replaces the O(Ngrid^2) direct real-space<->G-space DFT sums in the plane-wave basis with an
// O(Ngrid log Ngrid) transform.  Radix-2 only (lengths must be powers of two -- plane-wave FFT grids are
// padded up to the next power of two, which still resolves the difference set without aliasing).
//
// Conventions: the transforms are UNNORMALISED.  FFT1D/FFT3D with
//   sign=-1  compute  X_k = sum_j x_j e^{-2*pi*i*j.k/N}   (forward),
//   sign=+1  compute  x_j = sum_k X_k e^{+2*pi*i*j.k/N}   (inverse, caller divides by prod(N)).
module;
#include <complex>   // std::complex free operators (* + -) for dcmplx -- not leaked by qchem.Types
#include <cmath>
#include <cassert>
export module qchem.FFT;
import qchem.Types;   // cvec_t, ivec3_t
import qchem.Math;    // Pi

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

//! In-place radix-2 Cooley-Tukey FFT of \a a (size MUST be a power of two).  \a sign = -1 forward,
//! +1 inverse; unnormalised (see the file header).
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

//! 3-D FFT of a grid stored row-major as a flat vector (index = (i0*n1 + i1)*n2 + i2), applying FFT1D
//! along each axis.  Every \a n component must be a power of two.  \a sign / normalisation as FFT1D.
inline cvec_t FFT3D(cvec_t data, const ivec3_t& n, int sign)
{
    assert(data.size()==size_t(n.x)*n.y*n.z);
    assert(IsPow2(n.x) && IsPow2(n.y) && IsPow2(n.z));
    const size_t n0=n.x, n1=n.y, n2=n.z;

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
