// File: Common/Math.C
module;
#include <cmath>
#include <vector>
export module qchem.Math;
export import qchem.IntPow;
export import qchem.Constants;
export import qchem.Factorials;

export
{
    using std::sqrt;
    using std::fabs;
    using std::abs;
    using std::floor;
    using std::ceil;
    using std::pow;
    using std::exp;
    using std::log;
    using std::log10;
    using std::sin;
    using std::cos;
    using std::acos;
    using std::isfinite;
    using std::max;
    using std::min;
    using std::lround;

}

// Special functions, namespaced (qchem::Math, matching the module) so they never collide with the
// C math library -- e.g. spherical j1 vs the library's cylindrical ::j1.
export namespace qchem::Math
{
    //! Spherical Bessel functions \f$j_0(x)\dots j_{lmax}(x)\f$ via the stable downward (Miller)
    //! recurrence, normalised to \f$j_0=\sin x/x\f$.  Robust for all x (upward recurrence is unstable
    //! for \f$x<l\f$).  \f$j_l(0)=\delta_{l0}\f$.
    inline std::vector<double> SphericalBessel(int lmax, double x)
    {
        std::vector<double> j(lmax+1, 0.0);
        if (std::fabs(x)<1e-12) { j[0]=1.0; return j; }
        int top=lmax+15+static_cast<int>(x);
        std::vector<double> t(top+2, 0.0);
        t[top]=1e-30;
        for (int l=top; l>=1; --l) t[l-1]=(2*l+1)/x*t[l]-t[l+1];
        double scale=(std::sin(x)/x)/t[0];
        for (int l=0; l<=lmax; ++l) j[l]=t[l]*scale;
        return j;
    }

    //! Spherical Bessel \f$j_1(x)\f$ in closed form (convenience scalar; avoids a vector allocation).
    inline double SphericalBessel1(double x)
    {
        if (std::fabs(x)<1e-6) return x/3.0;
        return std::sin(x)/(x*x) - std::cos(x)/x;
    }

    //! Derivatives \f$j_0'(x)\dots j_{lmax}'(x)\f$ from a precomputed j array (\a x>0):
    //! \f$j_l'=j_{l-1}-\tfrac{l+1}{x}j_l\f$, with \f$j_{-1}(x)=\cos x/x\f$.
    inline std::vector<double> SphericalBesselPrime(int lmax, double x, const std::vector<double>& j)
    {
        std::vector<double> jp(lmax+1, 0.0);
        for (int l=0; l<=lmax; ++l)
        {
            double jprev=(l==0) ? std::cos(x)/x : j[l-1];
            jp[l]=jprev - (l+1)/x*j[l];
        }
        return jp;
    }

    //! Legendre polynomials \f$P_0(x)\dots P_{lmax}(x)\f$ via the (stable) upward recurrence.
    inline std::vector<double> LegendreP(int lmax, double x)
    {
        std::vector<double> P(lmax+1, 0.0);
        P[0]=1.0; if (lmax>=1) P[1]=x;
        for (int l=1; l<lmax; ++l) P[l+1]=((2*l+1)*x*P[l]-l*P[l-1])/(l+1);
        return P;
    }
}

