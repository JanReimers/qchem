// File: Polarization.C  Structure describing just the polarization portion of a basis function.
module;
#include <iosfwd>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.Types;
import Common.IntPow;
import oml.Vector3D;

export namespace PolarizedGaussian
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

    double operator   ()(const RVec3& r) const
    {
        return intpow(r.x,n)*intpow(r.y,l)*intpow(r.z,m);
    }
    RVec3  Gradient     (const RVec3& r) const;

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

} //namespace PolarizedGaussian

