// File: BSpline/BFGrouper.H  Group BSpline basis functions by support start positions.
#ifndef _BSpline_BFGrouper_H_
#define _BSpline_BFGrouper_H_

#include <bspline/Core.h>
#include <vector>
#include <map>

class GLCache;
// 
// Keep a list of unique exponents for Group Slater or Gaussian basis functions.
// For each unique exponent also store am index and the maximum l angular momentum used
// for that exponent.  The idea is to share exponents between different l-irreps 
// and work out the radial integral tables up to LMax for each exponent combination. 
// THis class should be working together with the charge distribution caching mechanism.
//  
namespace BSpline
{

template <class T, size_t K> struct cmpSplines {
    bool operator()(const bspline::Spline<T,K>& a, const bspline::Spline<T,K>& b) const
    {
        double af=a.getSupport().front(), bf=b.getSupport().front();
        if (af!=bf) return af<bf;
        return a.getSupport().back()<b.getSupport().back();
    }
};

template <size_t K> class IrrepIEClient;

template <size_t K> class BFGrouper
{
protected:
    typedef bspline::Spline<double, K> spline_t;
    void Append(IrrepIEClient<K>*);
    size_t LMax(size_t ia, size_t ib, size_t ic, size_t id) const;
    const GLCache* GetGL(size_t l) const;
    //! Linear array of unique splines.
    std::vector<spline_t> unique_spv; 
private:
    std::map<spline_t,size_t,cmpSplines<double,K>> unique_sp; //Unique splines.
    //! For each unique spline, store the maximum l value.
    std::vector<size_t> maxls;
    std::map<size_t,const GLCache*> itsGLs;
};

}


#endif // _BFGrouper_H_
