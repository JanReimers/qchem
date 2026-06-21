// File: Polarization.C  Structure describing just the polarization portion of a basis function.
module;
#include <iosfwd>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.Types;
import qchem.IntPow;
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD.Index3;  // Cartesian->Hermite index seam

export namespace BasisSet::Molecule::PolarizedGaussian
{

class Polarization
{
public:
    Polarization(                    ) : n(  0),l(  0),m(  0) {};
    Polarization(int _n,int _l,int _m) : n( _n),l( _l),m( _m) {};
    Polarization(const Polarization& p  ) : n(p.n),l(p.l),m(p.m) {};

    Polarization& operator= (const Polarization& p)
    {
        n=p.n;
        l=p.l;
        m=p.m;
        return *this;
    }
    Polarization  operator+ (const Polarization& p) const
    {
        return Polarization(n+p.n,l+p.l,m+p.m);
    }
    Polarization  operator- (const Polarization& p) const
    {
        return Polarization(n-p.n,l-p.l,m-p.m);
    }
    bool          operator==(const Polarization& p) const
    {
        return n==p.n&&l==p.l&&m==p.m;
    }
    bool          operator!=(const Polarization& p) const
    {
        return !(*this==p);
    }
    bool          operator>(const Polarization& p) const
    {
        return n>p.n|| l>p.l || m>p.m;
    }
    bool          operator<(const Polarization& p) const
    {
        int i1=  n*LMax*LMax+  l*LMax+  m;
        int i2=p.n*LMax*LMax+p.l*LMax+p.m;
        return i1<i2;
    }

    double operator   ()(const rvec3_t& r) const
    {
        return intpow(r.x,n)*intpow(r.y,l)*intpow(r.z,m);
    }
    rvec3_t  Gradient     (const rvec3_t& r) const;

    // A Cartesian polarization IS a Hermite (N,L,M) index to the generic MnD core (RNLM etc.); the
    // monomial above is the Cartesian-only part.  Implicit so existing rnlm(pa+pb) call sites are unchanged.
    operator Evaluators::Internal::MnD::Index3() const {return {n,l,m};}

    int   GetTotalL  () const;
    int   GetMaximumL() const;
    int   GetSign    () const; //(-1)^(n+l+m)

    friend std::ostream& operator<<(std::ostream&, const Polarization&);

    int n,l,m;
    const static int LMax=64;
};

inline int Polarization::GetTotalL() const
{
    return n+l+m;
}

inline int   Polarization::GetSign    () const
{
    return GetTotalL()%2==0 ? 1 : -1;
}

} //namespace BasisSet::Molecule::PolarizedGaussian

