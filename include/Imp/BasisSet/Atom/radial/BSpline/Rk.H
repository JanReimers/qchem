// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 
#ifndef _BSpline_RkEngine_H_
#define _BSpline_RkEngine_H_

#include "../../../Cache4.H"
#include <bspline/Core.h>
import oml;

class GLCache;
namespace BSpline
{
//
//  Calculate and store all off diagonal and diagonal cell slater integrals.
//
template <size_t K> class RkCache
{
    typedef bspline::Spline<double,K> sp_t;
    typedef std::vector<double> dv_t;
public:
    RkCache(const std::vector<sp_t>& splines,const GLCache& gl,size_t lmax);
    const dv_t& find_plus (size_t ia,size_t ib) const {return find(ia,ib,itsMomentsPlus);}
    const dv_t& find_minus(size_t ia,size_t ib) const {return find(ia,ib,itsMomentsMinus);}
private:
    typedef std::pair<size_t,size_t> id2_t; //convention id_1 < id_2
    typedef std::map<id2_t,dv_t> moment_t;
    static const dv_t& find(size_t ia,size_t ib,const moment_t&);


    moment_t itsMomentsPlus;  //<r^(2+k) B1*B2>
    moment_t itsMomentsMinus;  //<r^(1-k) B1*B2>
};
//
// Store Slater repulsions integral for 0 <= k <= 2*LMax.
//
template <size_t K> class RkEngine : public virtual Cacheable
{
    typedef bspline::Spline<double,K> sp_t;
public:
    RkEngine(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t LMax, const GLCache& gl, const RkCache<K>&);
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
