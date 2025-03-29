// File: Polarization.C  Structure describing just the polarization portion of a basis function.



#include "Imp/BasisSet/Molecule/PolarizedGaussian/Polarization.H"
#include "oml/imp/stream.h"
#include "oml/imp/binio.h"

namespace PolarizedGaussian
{

int Polarization::GetMaximumL() const
{
    int t=l>m ? l : m;
    return t>n ? t : n;
}



inline double DiffIntPow(double x,int n)
{
    return n==0 ? 0 : ((double)n)*intpow(x,n-1);
}

RVec3 Polarization::Gradient(const RVec3& r) const
{
    double x=intpow(r.x,n);
    double y=intpow(r.y,l);
    double z=intpow(r.z,m);
    return RVec3(DiffIntPow(r.x,n)*y*z, DiffIntPow(r.y,l)*x*z, DiffIntPow(r.z,m)*x*y );
}

std::ostream& operator<<(std::ostream& os, const Polarization& p)
{
    if ( StreamableObject::Binary())
    {
        BinaryWrite(p.n,os);
        BinaryWrite(p.l,os);
        BinaryWrite(p.m,os);
    }
    if ( StreamableObject::Ascii()) os << p.n << " " << p.l << " " << p.m << " ";
    if ( StreamableObject::Pretty())
    {
        static char SPDF[]="SPDFGHIJ";
        int n=1;
        os << SPDF[p.GetTotalL()];
        for (int i=0; i<p.n; i++,n++) os << "x";
        for (int i=0; i<p.l; i++,n++) os << "y";
        for (int i=0; i<p.m; i++,n++) os << "z";
//    for (;n<=5;n++) os << " ";
        os << " ";
    }
    return os;
}

std::istream& operator>>(std::istream& is,       Polarization& p)
{
    if ( StreamableObject::Binary())
    {
        BinaryRead(p.n,is);
        BinaryRead(p.l,is);
        BinaryRead(p.m,is);
    }
    if ( StreamableObject::Ascii())
        is >> p.n >> p.l >> p.m;

    return is;
}

} //namespace PolarizedGaussian
