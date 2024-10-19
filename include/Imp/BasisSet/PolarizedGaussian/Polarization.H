// File: Polarization.H  Structure describing just the polarization portion of a basis function.
#ifndef _Polarization_H_
#define _Polarization_H_

#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
#include <iosfwd>

typedef Vector3D<double> RVec3;
namespace PolarizedGaussian
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

    double operator   ()(const RVec3& r) const
    {
        return intpow(r.x,n)*intpow(r.y,l)*intpow(r.z,m);
    }
    RVec3  Gradient     (const RVec3& r) const;

    int   GetTotalL  () const;
    int   GetMaximumL() const;
    int   GetSign    () const; //(-1)^(n+l+m)

    friend std::ostream& operator<<(std::ostream&, const Polarization&);
    friend std::istream& operator>>(std::istream&,       Polarization&);

    int n,l,m;
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

#endif // _Polarization_H_
