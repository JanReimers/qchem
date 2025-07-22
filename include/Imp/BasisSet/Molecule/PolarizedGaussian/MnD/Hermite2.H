// File: Hermite2.H  Class for managing 2 function Hermite coefficients
#ifndef _Hermite2_H_
#define _Hermite2_H_

#include <vector>
#include <iosfwd>
#include <cassert>

#include "../Polarization.H"
import oml;
import qchem.Types;
//#define USE_CACHE

#ifdef USE_CACHE
#include <map>
#endif

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

class Polarization;
namespace PolarizedGaussian
{

class Hermite2
{
public:
    //Hermite2();
    Hermite2(double AlphaP, const RVec3& PA, const RVec3& PB, int LA, int LB);
    ~Hermite2();
    double operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual Hermite2* Clone(        ) const;

private:
    Hermite2(const Hermite2&) {};
    Hermite2& operator=(const Hermite2&) {return *this;}
    RVec3 Get(int N,int na,int nb) const 
    {
        size_t index=GetIndex(N,na,nb);
        return RVec3(d[index],e[index],f[index]);
    } 
    void Assign(int N,int na,int nb,const RVec3& a);
    //
    //  This has very high hits and total time in profiling.
    //
    inline size_t GetIndex(int N, int na, int nb) const 
    {
        assert(N>=0);
        assert(N<=LA+LB);
        assert(na>=0);
        assert(na<=LA);
        assert(nb>=0);
        assert(nb<=LB);
        assert(N<=na+nb);
        return N*LAB+na*LB1+nb;
    }
    size_t GetSize() const {return (LA+LB)*(LA+1)*(LB+1)+LA*(LB+1)+LB+1;}


    int LA, LB;
    int LAB,LB1;
    std::vector<double> d,e,f;
    #ifdef USE_CACHE
    std::map<Polarization,size_t> indexCache;
    #endif
};

} //namespace PolarizedGaussian

#endif //_Hermite2_H_
