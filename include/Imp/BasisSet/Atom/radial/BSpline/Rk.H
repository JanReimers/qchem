// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 
#ifndef _BSpline_RkEngine_H_
#define _BSpline_RkEngine_H_

#include "Imp/Integrals/Cache4.H"
#include <bspline/Core.h>
#include "oml/vector.h"

namespace BSpline
{
//
// Store Slater repulsions integral for 0 <= k <= 2*LMax.
//
template <size_t K> class RkEngine : public virtual Cacheable
{
    typedef bspline::Spline<double,K> sp_t;
public:
    RkEngine(const sp_t& a,const sp_t& b,const sp_t& c,const sp_t& d,size_t LMax);
    double Coulomb_R0() const; //R_0(la,la,lc,lc);
    Vector<double> Coulomb_Rk(size_t la,size_t lc) const; 
    Vector<double> ExchangeRk(size_t la,size_t lb) const; 
    //! R_k(la,lb,la,lb) with |Ala-Alb| <= k <= Ala+Alb
    // Vector<double> ExchangeRk() const; 
    
private:
    size_t LMax;
    Vector<double> Rabcd_k;
};

} //namespace

#endif // _BSpline_RkEngine_H_
