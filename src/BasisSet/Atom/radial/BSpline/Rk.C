// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 

#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"
#include "Imp/BasisSet/Atom/radial/BSpline/GLQuadrature.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector.h"

using std::cout;
using std::endl;

namespace BSpline
{
//
//  Calculate and store 2 electron radial repulsion (Slater) integrals for all valules of k.
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
        GLCache glcd1(sc.getGrid(),K+1+k+2);
        GLCache glcd2(sc.getGrid(),K+2);
        std::function< double (double)> wcd1 = [k](double r2)
        {
            return intpow(r2,k+2);
        };
        std::function< double (double)> wcd2 = [k](double r2)
        {
            assert(r2>0);
            return intpow(r2,1-k);
        };
        std::function< double (double)> Yk1 = [glcd1,wcd1,c,d,rmin](double r1)
        {
            return glcd1.Integrate(wcd1,c,d,rmin,r1);
        };
        std::function< double (double)> Yk2 = [glcd2,wcd2,c,d,rmax](double r1)
        {
            return glcd2.Integrate(wcd2,c,d,r1,rmax);
        };

        std::function< double (double)> wab = [k,Yk1,Yk2] (double r1)
        {
            assert(r1>0);
            return intpow(r1,1-k)*Yk1(r1)+intpow(r1,k+1)*Yk2(r1);
        };

        Rabcd_k(k)=glcd1.Integrate(wab,a,b);

        
    }
 }

 template <size_t K> Vector<double> RkEngine<K>::Coulomb_Rk(size_t la,size_t lc) const
 {
    assert(la>=0);
    assert(lc>=0);
    assert(la<=LMax);
    assert(lc<=LMax);
    Vector<double> Rk(std::min(la,lc)+1,0.0);
    size_t i=1;
    for (size_t k=0;k<=2*std::min(la,lc);k+=2)
    {
        Rk(i++)=Rabcd_k(k); 
    }
    return Rk;
 }
 template <size_t K> Vector<double> RkEngine<K>::ExchangeRk(size_t la,size_t lb) const
 {
    assert(la>=0);
    assert(lb>=0);
    assert(la<=LMax);
    assert(lb<=LMax);
    int kmin=std::abs((int)la-(int)lb);
    int kmax=la+lb;
    int N=(kmax-kmin)/2+1;
    Vector<double> Rk(N);
    int i=1;
    for (int k=kmin;k<=kmax;k+=2)
    {
        assert((k+la+lb)%2==0);
        Rk(i++)=Rabcd_k(k); 
    }
    return Rk;
 }
 
 template class RkEngine<6>;

} //namespace


