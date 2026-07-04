// File: Math/Angular.C  Foundational angular / polynomial basis-math shared by qcSymmetry and the basis
// evaluators.  Lives in qcMath (a PUBLIC dependency of both) so neither has to depend on the other -- in
// particular the low-level Cartesian/spherical integral evaluators must NOT depend on group theory
// (qcSymmetry), and qcSymmetry's rep builders must not depend on the evaluators.  Angular math (Cartesian
// monomials today; Wigner/Clebsch-Gordan, solid harmonics later) is genuinely foundational and sits below
// both -- hence the general name "Angular" rather than "Monomial".
module;
#include <vector>
#include <cassert>
#include <ostream>
export module qchem.Math.Angular;

export namespace qchem::Math
{

//! \brief The exponents \f$(n_x,n_y,n_z)\f$ of a Cartesian monomial \f$x^{n_x}y^{n_y}z^{n_z}\f$.
//!
//! A plain aggregate (so \c Monomial{2,0,0} works) carrying only the pure index operations the two consumers
//! share: axis access (the rep builders' \c e[i] / \c te[j]+=1), equality, and a strict-weak order for use as
//! a \c std::map key.  The order is LEXICOGRAPHIC in (n,l,m); for valid exponents \f$\in[0,64)\f$ this is
//! identical to the old \c std::array order (qcSymmetry rep builders) AND to \c Polarization's former
//! \f$n\cdot64^2+l\cdot64+m\f$ radix order (the PG_Cart_MnD Hermite index cache), so both maps' iteration
//! order is preserved. \c Polarization derives from this, adding the evaluator-coupled parts (arithmetic,
//! real-space evaluation, the Hermite-index conversion).
struct Monomial
{
    int n = 0, l = 0, m = 0;   //!< the three Cartesian exponents \f$(n_x,n_y,n_z)\f$

    int  operator[](int i) const {return (&n)[i];}   //!< axis access, 0->n(x) 1->l(y) 2->m(z)
    int& operator[](int i)       {return (&n)[i];}

    bool operator==(const Monomial& p) const {return n==p.n && l==p.l && m==p.m;}
    bool operator!=(const Monomial& p) const {return !(*this==p);}

    //! Lexicographic ordering (n, then l, then m).  See the class note: byte-compatible with the two
    //! former orderings for valid (non-negative, <64) exponents.
    bool operator<(const Monomial& p) const
    {
        if (n != p.n) return n < p.n;
        if (l != p.l) return l < p.l;
        return m < p.m;
    }
};

//! Compact unambiguous stream form \c xNyLzM (used e.g. in geometry-aware integral-cache keys).
inline std::ostream& operator<<(std::ostream& os, const Monomial& p)
{
    return os << 'x' << p.n << 'y' << p.l << 'z' << p.m;
}

//! \brief One term of a function's Cartesian-monomial expansion: a monomial and its coefficient.  An
//! aggregate, so \c CartTerm{ {n,l,m}, c } (and the c2s literals in the qcSymmetry rep tests) work.
struct CartTerm
{
    Monomial p;   //!< the Cartesian monomial \f$x^{p.n} y^{p.l} z^{p.m}\f$
    double   c;   //!< its coefficient
};

//! \brief A "Cartesian-to-spherical" map: a list of harmonics, each expressed as its Cartesian-monomial
//! expansion.  \c c2s[m] is the m-th harmonic of a shell, ordered by the basis's own m-convention.
using HarmonicC2S = std::vector<std::vector<CartTerm>>;

//! \brief The \f$2l+1\f$ real regular solid harmonics \f$R_{l,m}\f$ (homogeneous degree-\f$l\f$ harmonic
//! polynomials, \f$\nabla^2 R=0\f$) for shell \a l, ordered \f$m=-l..+l\f$, each as its Cartesian-monomial
//! expansion (\c l=0..3 = s,p,d,f; \c l>=4 asserts).  RAW coefficients: only the RELATIVE weights within
//! each \f$m\f$ matter (the overall per-harmonic scale is fixed elsewhere by unit self-overlap).  Shared by
//! the spherical-Gaussian basis evaluator (which folds in normalisation) and the qcSymmetry spherical
//! operation-rep builder (which projects the Cartesian rep through this map).
inline HarmonicC2S SphericalShell(int l)
{
    using M = Monomial;   // M{a,b,c} == monomial x^a y^b z^c
    switch (l)
    {
    case 0: return {
        { {M{0,0,0}, 1.0} },                                              // s
    };
    case 1: return {
        { {M{0,1,0}, 1.0} },                                             // m=-1  y
        { {M{0,0,1}, 1.0} },                                             // m= 0  z
        { {M{1,0,0}, 1.0} },                                             // m=+1  x
    };
    case 2: return {
        { {M{1,1,0}, 1.0} },                                            // m=-2  xy
        { {M{0,1,1}, 1.0} },                                            // m=-1  yz
        { {M{0,0,2}, 1.0}, {M{2,0,0},-0.5}, {M{0,2,0},-0.5} },          // m= 0  z^2 - (x^2+y^2)/2
        { {M{1,0,1}, 1.0} },                                            // m=+1  xz
        { {M{2,0,0}, 1.0}, {M{0,2,0},-1.0} },                           // m=+2  x^2 - y^2
    };
    case 3: return {
        { {M{2,1,0}, 3.0}, {M{0,3,0},-1.0} },                          // m=-3  3x^2 y - y^3
        { {M{1,1,1}, 1.0} },                                           // m=-2  xyz
        { {M{0,1,2}, 4.0}, {M{2,1,0},-1.0}, {M{0,3,0},-1.0} },         // m=-1  y(4z^2 - x^2 - y^2)
        { {M{0,0,3}, 2.0}, {M{2,0,1},-3.0}, {M{0,2,1},-3.0} },         // m= 0  z(2z^2 - 3x^2 - 3y^2)
        { {M{1,0,2}, 4.0}, {M{3,0,0},-1.0}, {M{1,2,0},-1.0} },         // m=+1  x(4z^2 - x^2 - y^2)
        { {M{2,0,1}, 1.0}, {M{0,2,1},-1.0} },                          // m=+2  z(x^2 - y^2)
        { {M{3,0,0}, 1.0}, {M{1,2,0},-3.0} },                          // m=+3  x(x^2 - 3y^2)
    };
    default:
        assert(false && "qchem::Math::SphericalShell: only l=0..3 (s,p,d,f) supported");
        return {};
    }
}

} //namespace qchem::Math
