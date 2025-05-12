// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 

#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"
#include "oml/vector.h"

namespace BSpline
{
//
//  Ranges:  
//    0 <= k <= 2LMax  in steps of 2
//    3 <= Lab_p=la+lb+3+k <= 4LMax+3
//    1 <= Lcd_m=lc+ld+1-k <= 2LMax+1
//    3 <= Lab_m=la+lb+1-k <= 4LMax+3
//    1 <= Lcd_p=lc+ld+3+k <= 2LMax+1
//
//  Build up the derivative look up tables.
//
RkEngine::RkEngine(double _eab, double _ecd, size_t _LMax)
 : eab(_eab), ecd(_ecd), LMax(_LMax)
 {
        
 }
 
 
double RkEngine::Coulomb_R0(size_t la,size_t lc) const
{
    assert(la<=LMax);
    assert(lc<=LMax);
   
    return 0.0;
}

Vector<double> RkEngine::Coulomb_Rk(size_t la,size_t lc) const
{
    assert(la>=0);
    assert(lc>=0);
    assert(la<=LMax);
    assert(lc<=LMax);
    Vector<double> ret(la+lc+1,0.0);
    
    return ret;
}

Vector<double> RkEngine::ExchangeRk(size_t la,size_t lb) const
{
    assert(la>=0);
    assert(lb>=0);
    assert(la<=LMax);
    assert(lb<=LMax);
    size_t kmin=std::abs((int)la-(int)lb);
    size_t kmax=la+lb;
    size_t N=(kmax-kmin)/2+1;
    Vector<double> ret(N,0.0);
    
    return ret;
}

Vector<double> RkEngine::ExchangeRk(size_t Ala,size_t Alb, size_t la,size_t lb) const
{
    assert (Ala>=la);
    assert (Alb>=lb);
    assert(la>=0);
    assert(lb>=0);
    assert(la<=LMax);
    assert(lb<=LMax);
    size_t kmin=std::abs((int)Ala-(int)Alb);
    size_t kmax=Ala+Alb;
    size_t N=(kmax-kmin)/2+1;
    Vector<double> ret(N,0.0);
    
    return ret;
}

} //namespace