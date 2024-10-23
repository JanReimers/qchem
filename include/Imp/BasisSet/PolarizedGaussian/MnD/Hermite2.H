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

    
    double Getd1(int N,int na,int nb) const
    {
        return GetAny(d1,N,na,nb);
    }
    double Gete1(int L,int la,int lb) const
    {
        return GetAny(e1,L,la,lb);
    }
    double Getf1(int M,int ma,int mb) const
    {
        return GetAny(f1,M,ma,mb);
    }
    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual Hermite2* Clone(        ) const;

private:
    double GetAny(const std::vector<double>& def, int N,int na,int nb) const;
    size_t GetIndex(size_t N, size_t na, size_t nb) const 
    {
        assert(N>=0);
        assert(N<=LA+LB);
        assert(na>=0);
       // assert(na<=LA);
        assert(nb>=0);
        //assert(nb<=LB);
        assert(N<=na+nb);
        size_t index=N*(LA+1)*(LB+1)+na*(LB+1)+nb;
        if (d1.size()>0 && index>d1.size())
            std::cout << N << " " << na << " " << nb << " " << index << " " << d1.size() << std::endl;
        assert(d1.size()==0 || index<=d1.size());
        return index;
    }
    size_t GetSize() const {return (LA+LB)*(LA+1)*(LB+1)+LA*(LB+1)+LB+1;}


    int             LA, LB;
    std::vector<double> d1,e1,f1;
};

} //namespace PolarizedGaussian

#endif //_Hermite2_H_
