// File: Polarization.C  Structure describing just the polarization portion of a basis function.
module;
#include <iosfwd>
export module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Types;
import qchem.Math.Angular;                                     // Monomial (the shared (n,l,m) index skeleton)
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Index3;  // Cartesian->Hermite index seam

export namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{

//! \brief A Cartesian monomial exponent triple \f$(n,l,m)\f$ that ALSO evaluates itself in real space and
//! maps to the generic MnD Hermite index.  The pure index skeleton -- the \a n,l,m members, axis access,
//! ordering (\c operator<) and equality -- lives in the shared \c qchem::Math::Monomial base; Polarization
//! adds only the evaluator-coupled behaviour (arithmetic that must return \c Polarization, real-space
//! evaluation, and the \c Index3 conversion).
class Polarization : public qchem::Math::Monomial
{
public:
    Polarization(                    ) : Monomial{0,0,0} {}
    Polarization(int _n,int _l,int _m) : Monomial{_n,_l,_m} {}
    // Implicit lift from a bare Monomial: a Polarization adds no data, only evaluator behaviour, so a
    // Monomial (e.g. a qcMath CartTerm's exponents) can stand in wherever a Polarization is expected -- the
    // spherical evaluator feeds SphericalShell()'s Monomial terms straight into the Cartesian ERI kernels.
    Polarization(const qchem::Math::Monomial& mono) : Monomial(mono) {}

    // Arithmetic MUST return Polarization (not Monomial) so the rnlm(pa+pb)->Index3 call sites are unchanged.
    Polarization  operator+ (const Polarization& p) const {return Polarization(n+p.n, l+p.l, m+p.m);}
    Polarization  operator- (const Polarization& p) const {return Polarization(n-p.n, l-p.l, m-p.m);}
    bool          operator> (const Polarization& p) const {return n>p.n || l>p.l || m>p.m;}

    //! \f$x^n y^l z^m\f$.  A VALUE LEAF of the Enzyme-registered basis gradient (GaussianRF::Gradient), so it is inline (the AD unit
    //! differentiates its body) and branch-free: \c intpow returns 0 at x = 0 through a branch, and the derivative of a branch is the
    //! branch's, so AD through it would read \f$\partial x/\partial x = 0\f$ at every mesh point on an axis through the centre.
    inline double operator()(const rvec3_t& r) const   // explicit inline: a module's in-class body is NOT inline, so not importable
    {
        return Power(r.x,n)*Power(r.y,l)*Power(r.z,m);
    }
    //! \f$x^k\f$, k >= 0, as a plain product (see operator()).
    static inline double Power(double x, int k) {double p=1.0; for (int i=0;i<k;i++) p*=x; return p;}

    // A Cartesian polarization IS a Hermite (N,L,M) index to the generic MnD core (RNLM etc.).  Implicit so
    // existing rnlm(pa+pb) call sites are unchanged.
    operator Evaluators::Internal::MnD::Index3() const {return {n,l,m};}

    int   GetTotalL  () const;
    int   GetMaximumL() const;
    int   GetSign    () const; //(-1)^(n+l+m)

    friend std::ostream& operator<<(std::ostream&, const Polarization&);
};

inline int Polarization::GetTotalL() const
{
    return n+l+m;
}

inline int   Polarization::GetSign    () const
{
    return GetTotalL()%2==0 ? 1 : -1;
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD

