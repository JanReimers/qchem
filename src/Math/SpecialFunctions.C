// File: Common/SpecialFunctions.C  Special functions (spherical Bessel, Legendre) returning rvec_t.
//
// These live in their own module rather than qchem.Math because they return rvec_t: qchem.Math is
// below qchem.Types in the module DAG (qchem.Vector3D re-imports qchem.Math, qchem.Types imports
// qchem.Vector3D), so Math cannot depend on rvec_t without a cycle.  qchem.SpecialFunctions sits
// above Types and so returns rvec_t directly -- callers avoid std::vector<->rvec_t conversions.
module;
#include <cmath>

export module qchem.SpecialFunctions;
import qchem.Types;   // rvec_t

export namespace qchem::SpecialFunctions
{
    //! Spherical Bessel functions \f$j_0(x)\dots j_{lmax}(x)\f$ via the stable downward (Miller)
    //! recurrence, normalised to \f$j_0=\sin x/x\f$.  Robust for all x.  \f$j_l(0)=\delta_{l0}\f$.
    inline rvec_t SphericalBessel(int lmax, double x)
    {
        rvec_t j(lmax+1, 0.0);
        if (std::fabs(x)<1e-12) { j[0]=1.0; return j; }
        int top=lmax+15+static_cast<int>(x);
        rvec_t t(top+2, 0.0);
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
    inline rvec_t SphericalBesselPrime(int lmax, double x, const rvec_t& j)
    {
        rvec_t jp(lmax+1, 0.0);
        for (int l=0; l<=lmax; ++l)
        {
            double jprev=(l==0) ? std::cos(x)/x : j[l-1];
            jp[l]=jprev - (l+1)/x*j[l];
        }
        return jp;
    }

    //! Legendre polynomials \f$P_0(x)\dots P_{lmax}(x)\f$ via the (stable) upward recurrence.
    inline rvec_t LegendreP(int lmax, double x)
    {
        rvec_t P(lmax+1, 0.0);
        P[0]=1.0; if (lmax>=1) P[1]=x;
        for (int l=1; l<lmax; ++l) P[l+1]=((2*l+1)*x*P[l]-l*P[l-1])/(l+1);
        return P;
    }
}
