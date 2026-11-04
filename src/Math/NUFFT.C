// File: Math/NUFFT.C  Type-2 non-uniform FFT on an integer index box: f(s) = Sum_k c_k e^{2 pi i k.s}.
//
// WHY (doc/OpenWork.md §4, the rho-tilde sampler, 2026-11-04): a density held as Fourier coefficients on a
// G-ball is sampled at every point of an atom-centred (Becke) mesh EVERY Kerker/Pulay iteration.  The direct
// sum is O(N_G * N_r) and grows with the ball though the mesh does not (MnO geometry run: 1415 + 689 s of a
// 4100 s SCF).  This is the SAME band-limited series evaluated to a STATED tolerance in
// O(N_G + n^3 log n + N_r w^3): no lattice sum is cut (doc/Pins.md pin 1) -- the only error is the
// spreading kernel's Fourier alias tail, which the width w controls.
//
// The algorithm (Dutt-Rokhlin / Greengard-Lee, in FINUFFT's "exponential of semicircle" form):
//   1. DECONVOLVE  b_k = c_k / psi^(2 pi k/n)              (psi^ = the kernel's Fourier transform)
//   2. FFT         g_j = Sum_k b_k e^{+2 pi i k.j/n}       (an oversampled grid, n >= sigma * box)
//   3. INTERPOLATE f(s) = Sum_j g_j psi(n s - j)            (the w^3 grid nodes nearest n s)
// Poisson summation makes step 3 reproduce e^{2 pi i k.s} psi^(2 pi k/n) exactly up to the alias images
// psi^(2 pi (k/n + m)), m != 0; step 1 divides the m = 0 image out.  The kernel
//   psi(x) = exp(beta (sqrt(1 - (2x/w)^2) - 1)),  |x| <= w/2 (grid units),  beta = 2.30 w,  sigma = 2,
// leaves a relative error ~10^{1-w} against ||c||_1 (Barnett, Magland & af Klinteberg, SIAM J. Sci. Comput.
// 41 (2019) C479; the width rule w = ceil(log10(1/tol)) + 1 is theirs).  A kernel of compact support is a
// property of the kernel, not a cut of a field: the error is its alias tail, bounded and reported.
//
// Pure math, no threading: Grid() runs once per coefficient set, Interpolate() is const and per point, so a
// caller threads the point loop itself (and each point is one thread's gather -- bit-identical at any count).
module;
#include <complex>
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
export module qchem.NUFFT;
import qchem.Types;   // cvec_t, ivec3_t, rvec3_t, dcmplx
import qchem.FFT;     // FFT3D, Next5Smooth
import qchem.CMath;   // Pi

namespace qchem::NUFFT
{
//! Component \a a of a fixed 3-vector (Vector3D has no indexing).
template <class T> T At(const Vector3D<T>& v, int a) {return a==0 ? v.x : (a==1 ? v.y : v.z);}
}

export namespace qchem::NUFFT
{

//! A type-2 plan on the index box [lo, hi] (inclusive, per axis) at a requested relative tolerance.
class Type2
{
public:
    static constexpr double sigma = 2.0;   //!< the oversampling ratio the width rule is calibrated for
    //! The kernel width for a relative tolerance \a tol at sigma = 2 (FINUFFT's rule, capped at 16 -- double
    //! precision floors near 1e-14 whatever the width).
    static int WidthFor(double tol)
    {
        assert(tol>0.0);
        return std::clamp(int(std::ceil(std::log10(1.0/tol)))+1, 2, 16);
    }

    Type2(const ivec3_t& lo, const ivec3_t& hi, double tol)
        : itsLo(lo), itsHi(hi), itsW(WidthFor(tol)), itsBeta(2.30*itsW)
    {
        for (int a=0; a<3; a++)
        {
            // The alias error is governed by max|k|/n, not by the box's extent, so size from the SYMMETRIC hull
            // (an asymmetric box sized by its extent runs below sigma on that axis: measured 12x short at w 8).
            assert(At(hi,a)>=At(lo,a));
            const int nk=2*std::max(std::abs(At(lo,a)),std::abs(At(hi,a)))+1;
            n[a]=int(FFT::Next5Smooth(size_t(std::max(int(std::ceil(sigma*nk)), 2*itsW))));
        }
        itsN=ivec3_t(n[0],n[1],n[2]);
        for (int a=0; a<3; a++) itsDeconv[a]=MakeDeconvolution(a);
    }

    int Width() const {return itsW;}
    const ivec3_t& GridSize() const {return itsN;}

    //! Steps 1+2: deconvolve and inverse-FFT the coefficients onto the oversampled grid.  \a terms are
    //! (k, c_k) with every k inside the box; a k absent from \a terms is a zero coefficient.
    template <class Terms> cvec_t Grid(const Terms& terms) const
    {
        const size_t n0=itsN.x, n1=itsN.y, n2=itsN.z;
        cvec_t g(n0*n1*n2, dcmplx(0.0));
        for (const auto& [k, c] : terms)
        {
            assert(k.x>=itsLo.x && k.x<=itsHi.x && k.y>=itsLo.y && k.y<=itsHi.y && k.z>=itsLo.z && k.z<=itsHi.z);
            const size_t j0=Wrap(k.x,itsN.x), j1=Wrap(k.y,itsN.y), j2=Wrap(k.z,itsN.z);
            const double d=itsDeconv[0][k.x-itsLo.x]*itsDeconv[1][k.y-itsLo.y]*itsDeconv[2][k.z-itsLo.z];
            g[(j0*n1+j1)*n2+j2] += d*dcmplx(c);
        }
        return FFT::FFT3D(std::move(g), itsN, +1);
    }

    //! Step 3: f(s) at one point, \a s in FRACTIONAL coordinates (any real value; the series is 1-periodic).
    dcmplx Interpolate(const cvec_t& g, const rvec3_t& s) const
    {
        assert(g.size()==size_t(itsN.x)*itsN.y*itsN.z);
        constexpr int wMax=16;
        int    idx[3][wMax];
        double wt [3][wMax];
        for (int a=0; a<3; a++)
        {
            const double sa=At(s,a), y=n[a]*(sa-std::floor(sa));        // in [0, n)
            const int j0=int(std::ceil(y-0.5*itsW));
            for (int t=0; t<itsW; t++)
            {
                const int j=j0+t;
                idx[a][t]=int(Wrap(j,n[a]));
                wt [a][t]=Kernel(y-j);
            }
        }
        const size_t n1=itsN.y, n2=itsN.z;
        dcmplx f(0.0);
        for (int t0=0; t0<itsW; t0++)
        {
            dcmplx f1(0.0);
            for (int t1=0; t1<itsW; t1++)
            {
                const size_t base=(size_t(idx[0][t0])*n1+size_t(idx[1][t1]))*n2;
                dcmplx f2(0.0);
                for (int t2=0; t2<itsW; t2++) f2 += wt[2][t2]*g[base+idx[2][t2]];
                f1 += wt[1][t1]*f2;
            }
            f += wt[0][t0]*f1;
        }
        return f;
    }

private:
    static size_t Wrap(int j, int n) {const int m=j%n; return size_t(m<0 ? m+n : m);}
    //! psi(x), x in grid units: the ES kernel on |x| <= w/2, zero outside.
    double Kernel(double x) const
    {
        const double z=2.0*x/itsW, q=1.0-z*z;
        return q>0.0 ? std::exp(itsBeta*(std::sqrt(q)-1.0)) : 0.0;
    }
    //! 1/psi^(2 pi k/n_a) for every k in the box along axis \a a.  psi^(xi) = 2 Int_0^{w/2} psi(x) cos(xi x) dx
    //! by Gauss-Legendre: psi is smooth on the open interval and ~e^{-beta} at its ends, so a rule of a few
    //! times w nodes integrates it to machine precision.
    std::vector<double> MakeDeconvolution(int a) const
    {
        const int q=4*itsW+20;
        std::vector<double> x, wq;
        GaussLegendre(q, x, wq);                                    // on [-1, 1]
        const double half=0.5*itsW;
        const int lo=At(itsLo,a), hi=At(itsHi,a);
        std::vector<double> r(hi-lo+1);
        for (int k=lo; k<=hi; k++)
        {
            const double xi=2.0*Pi*k/n[a];
            double s=0.0;
            for (int i=0; i<q; i++)
            {
                const double xx=half*0.5*(x[i]+1.0);                // [0, w/2]
                s += wq[i]*Kernel(xx)*std::cos(xi*xx);
            }
            const double hat=2.0*half*0.5*s;                        // 2 * (w/4) * Sum
            assert(hat>0.0);
            r[k-lo]=1.0/hat;
        }
        return r;
    }
    //! Gauss-Legendre nodes/weights on [-1, 1] (Newton on the three-term recurrence).
    static void GaussLegendre(int q, std::vector<double>& x, std::vector<double>& w)
    {
        x.assign(q,0.0); w.assign(q,0.0);
        for (int i=0; i<q; i++)
        {
            double z=std::cos(Pi*(i+0.75)/(q+0.5)), dp=0.0;
            for (int it=0; it<100; it++)
            {
                double p0=1.0, p1=z;
                for (int l=2; l<=q; l++) {const double p2=((2*l-1)*z*p1-(l-1)*p0)/l; p0=p1; p1=p2;}
                dp=q*(z*p1-p0)/(z*z-1.0);                 // P_q' from P_q, P_{q-1}
                const double dz=p1/dp;
                z-=dz;
                if (std::abs(dz)<1e-16) break;
            }
            x[i]=z; w[i]=2.0/((1.0-z*z)*dp*dp);
        }
    }

    ivec3_t itsLo, itsHi, itsN;
    int     n[3];        //!< itsN as an array (the per-axis loops)
    int     itsW;
    double  itsBeta;
    std::vector<double> itsDeconv[3];
};

} //namespace
