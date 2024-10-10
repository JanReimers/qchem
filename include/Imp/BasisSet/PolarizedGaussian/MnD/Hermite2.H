// File: Hermite2.H  Class for managing 2 function Hermite coefficients
#ifndef _Hermite2_H_
#define _Hermite2_H_

#include "../Polarization.H"
#include "oml/vector.h"
#include "oml/vector3d.h"
#include <iosfwd>

#define LMAX 3

//
//  Encapsulate calculation and storage of Hermite expansion coefficients
//    _   _   _
//   nn  ll  mm
//  d   e   f
//   N   L   M
//
//  The op() function returns the product of all three.  This is more efficient
//  because this function can do fast checks to see if any of the three are
//  zero before doing a more expensive lookup.
//

typedef Vector3D<double> RVec3;
class Polarization;
namespace PolarizedGaussian
{

class Hermite2
{
public:
    Hermite2();
    Hermite2(double AlphaP, const RVec3& PA, const RVec3& PB, int LA, int LB);

    double operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const;

    double Getd(int N,int na,int nb) const
    {
        return GetAny(d,N,na,nb);
    }
    double Gete(int L,int la,int lb) const
    {
        return GetAny(e,L,la,lb);
    }
    double Getf(int M,int ma,int mb) const
    {
        return GetAny(f,M,ma,mb);
    }

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual Hermite2* Clone(        ) const;

private:
    double GetAny(const Vector<double>& def, int N,int na,int nb) const;

    static int      theIndex[2*LMAX+1][LMAX+1][LMAX+1];
    static int      theArraySizes[LMAX+1][LMAX+1];
    int             itsLA, itsLB;
    Vector<double>  d,e,f;
};

} //namespace PolarizedGaussian

#endif //_Hermite2_H_
