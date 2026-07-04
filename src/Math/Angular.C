// File: Math/Angular.C  Foundational angular / polynomial basis-math shared by qcSymmetry and the basis
// evaluators.  Lives in qcMath (a PUBLIC dependency of both) so neither has to depend on the other -- in
// particular the low-level Cartesian/spherical integral evaluators must NOT depend on group theory
// (qcSymmetry), and qcSymmetry's rep builders must not depend on the evaluators.  Angular math (Cartesian
// monomials today; Wigner/Clebsch-Gordan, solid harmonics later) is genuinely foundational and sits below
// both -- hence the general name "Angular" rather than "Monomial".
module;
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

} //namespace qchem::Math
