// File: PolarizedGaussian/MnD/Hermite2.C  Class for managing 2 function Hermite coefficients
module;

#include <vector>
#include <iosfwd>
#include <cassert>

export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite2;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
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

export
namespace PolarizedGaussian
{

class Hermite2
{
public:
    //Hermite2();
    Hermite2(double AlphaP, const rvec3_t& PA, const rvec3_t& PB, int LA, int LB);
    ~Hermite2();
    double operator()(const Polarization& P,const Polarization& Pa,const Polarization& Pb) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual Hermite2* Clone(        ) const;

private:
    Hermite2(const Hermite2&) {};
    Hermite2& operator=(const Hermite2&) {return *this;}
    rvec3_t Get(int N,int na,int nb) const 
    {
        size_t index=GetIndex(N,na,nb);
        return rvec3_t(d[index],e[index],f[index]);
    } 
    void Assign(int N,int na,int nb,const rvec3_t& a);
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

