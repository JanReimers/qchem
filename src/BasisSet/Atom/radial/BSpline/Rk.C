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
template <size_t K> RkEngine<K>::RkEngine(const sp_t& a,const sp_t& b,const sp_t& c,const sp_t& d, size_t _LMax)
 : LMax(_LMax)
 {
    assert(c.getSupport().hasSameGrid(d));
    assert(a.getSupport().hasSameGrid(b));
    for (size_t k=0;k<2*LMax;k++)
    {
        // GLCache glcd1(c.getSupport().getGrid(),K+1+k);
        // GLCache glcd2(c.getSupport().getGrid(),K+1);
        // std::function< double (double)> wcd1 = [k](double r2){return pow(r2,k);};
        // std::function< double (double)> wcd2 = [k](double r2){return pow(r2,-(k+1));};
        // // std::function< double (double)> Yk1 = [glcd1,c,d](double r1)
        // {
        //     // return 
        // }

        
    }
 }
 
 

} //namespace