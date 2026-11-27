// File: Common/SpecialFunctions.C  Special functions (spherical Bessel, Legendre) returning rvec_t.
//
// These live in their own module rather than qchem.Math because they return rvec_t: qchem.Math is
// below qchem.Types in the module DAG (qchem.Vector3D re-imports qchem.Math, qchem.Types imports
// qchem.Vector3D), so Math cannot depend on rvec_t without a cycle.  qchem.SpecialFunctions sits
// above Types and so returns rvec_t directly -- callers avoid std::vector<->rvec_t conversions.
module;

export module qchem.SpecialFunctions;
import std.compat;   // std as a module; std.compat because the code uses ::size_t, ::abs
import qchem.Types;   // rvec_t

export namespace qchem::SpecialFunctions
{
    //! Spherical Bessel functions \f$j_0(x)\dots j_{lmax}(x)\f$ via the stable downward (Miller) recurrence.  \f$j_l(0)=\delta_{l0}\f$.
    //! Three details, each earned by a wrong value against the reference table (src/Common/tests/SpecialFunctions.C, 2026-11-24):
    //!   - NORMALISE on whichever of \f$j_0=\sin x/x\f$, \f$j_1=\sin x/x^2-\cos x/x\f$ is larger in magnitude: they have no common
    //!     zero.  Normalising on \f$j_0\f$ alone gave every \f$j_{l\ge1}\f$ 40-90% wrong at \f$x=n\pi\f$.  \f$j_0\f$ itself is then
    //!     set from its closed form (the recurrence gets it by cancellation near its zeros).
    //!   - START at \f$\max(lmax,x)+20+4\sqrt{x}\f$: the recurrence must begin well above \f$x\f$, where it is stable.  The old
    //!     \f$lmax+15+x\f$ left \f$j_2(50)\f$ wrong by 4e-10 (relative) when lmax=3.
    //!   - RESCALE as the trial values grow (\f$\prod (2l+1)/x\f$ overflows a double for tiny \f$x\f$ and large lmax).
    inline rvec_t SphericalBessel(int lmax, double x)
    {
        rvec_t j(lmax+1, 0.0);
        if (std::fabs(x)<1e-12) { j[0]=1.0; return j; }
        const int top=std::max(lmax,static_cast<int>(std::fabs(x)))+20+static_cast<int>(4.0*std::sqrt(std::fabs(x)));
        rvec_t t(top+2, 0.0);
        t[top]=1e-30;
        for (int l=top; l>=1; --l)
        {
            t[l-1]=(2*l+1)/x*t[l]-t[l+1];
            if (std::fabs(t[l-1])>1e250)
                for (int k=l-1; k<=top; ++k) t[k]*=1e-250;
        }
        const double j0=std::sin(x)/x, j1=std::sin(x)/(x*x)-std::cos(x)/x;
        const double scale = std::fabs(j0)>=std::fabs(j1) ? j0/t[0] : j1/t[1];
        for (int l=0; l<=lmax; ++l) j[l]=t[l]*scale;
        j[0]=j0;
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
