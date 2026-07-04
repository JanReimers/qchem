// File: Polarization.C  Structure describing just the polarization portion of a basis function.
module;
#include <iosfwd>
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Types;
import qchem.IntPow;
import qchem.Math.Angular;                                     // Monomial (the shared (n,l,m) index skeleton)
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD.Index3;  // Cartesian->Hermite index seam

export namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
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

    // Arithmetic MUST return Polarization (not Monomial) so the rnlm(pa+pb)->Index3 call sites are unchanged.
    Polarization  operator+ (const Polarization& p) const {return Polarization(n+p.n, l+p.l, m+p.m);}
    Polarization  operator- (const Polarization& p) const {return Polarization(n-p.n, l-p.l, m-p.m);}
    bool          operator> (const Polarization& p) const {return n>p.n || l>p.l || m>p.m;}

    double operator   ()(const rvec3_t& r) const
    {
        return intpow(r.x,n)*intpow(r.y,l)*intpow(r.z,m);
    }
    rvec3_t  Gradient     (const rvec3_t& r) const;

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

} //namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD

