// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 

#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"
#include "Imp/BasisSet/Atom/radial/BSpline/GLQuadrature.H"
#include "oml/vector.h"

using std::cout;
using std::endl;

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
 : LMax(_LMax), Rabcd_k(VecLimits(0,2*LMax))
 {
    auto& sa=a.getSupport();
    auto& sb=b.getSupport();
    auto& sc=c.getSupport();
    auto& sd=d.getSupport();
    
    assert(sc.hasSameGrid(sd));
    assert(sa.hasSameGrid(sb));
    bspline::Grid grid=sa.getGrid();
    size_t rmin=grid.front(),rmax=grid.back();
    for (size_t k=0;k<=2*LMax;k++)
    {
        GLCache glcd1(sc.getGrid(),K+1+k);
        GLCache glcd2(sc.getGrid(),K+1);
        std::function< double (double)> wcd1 = [k](double r2)
        {
            return pow(r2,k);
        };
        std::function< double (double)> wcd2 = [k](double r2)
        {
            assert(r2>0);
            assert(std::isfinite(pow(1.0/r2,k+1)));
            return pow(1.0/r2,k+1);
        };
        std::function< double (double)> Yk1 = [glcd1,wcd1,c,d,rmin](double r1)
        {
            return glcd1.Integrate(wcd1,c,d,rmin,r1);
        };
        std::function< double (double)> Yk2 = [glcd2,wcd2,c,d,rmax](double r1)
        {
            return glcd2.Integrate(wcd2,c,d,r1,rmax);
        };

        std::function< double (double)> wab = [k,grid,Yk1,Yk2] (double r1)
        {
            assert(r1>0);
            return pow(1.0/r1,k+1)*Yk1(r1)+pow(r1,k)*Yk2(r1);
        };

        Rabcd_k(k)=glcd1.Integrate(wab,a,b);

        
    }
 }

 template <size_t K> Vector<double> RkEngine<K>::Coulomb_Rk() const
 {
    return Rabcd_k;
 }
 
 template class RkEngine<6>;

} //namespace


