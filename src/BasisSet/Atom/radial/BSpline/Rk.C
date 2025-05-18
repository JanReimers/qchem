// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 

#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"
#include "Imp/BasisSet/Atom/radial/BSpline/GLQuadrature.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector.h"

using std::cout;
using std::endl;

namespace BSpline
{
template <size_t K> RkCache<K>::RkCache(const std::vector<sp_t>& splines,const GLCache& gl, size_t lmax)
{
    for (size_t ia=0;ia<splines.size();ia++)
        for (size_t ib=ia;ib<splines.size();ib++)
        {
            std::vector<double> mp,mm;
            for (size_t k=0;k<=2*lmax;k++)
            {
                std::function< double (double)> wp = [k](double r) {return intpow(r,k+2);};
                std::function< double (double)> wm = [k](double r) {return intpow(r,1-k);};
                mp.push_back(gl.Integrate(wp,splines[ia],splines[ib]));
                mm.push_back(gl.Integrate(wm,splines[ia],splines[ib]));

            }
            itsMomentsPlus [std::make_pair(ia,ib)]=mp;
            itsMomentsMinus[std::make_pair(ia,ib)]=mm;
        }

}

template <size_t K> const typename RkCache<K>::dv_t& RkCache<K>::find(size_t ia,size_t ib,const moment_t& mm)
{
    if (ia>ib) std::swap(ia,ib);
    auto i=mm.find(std::make_pair(ia,ib));
    assert(i!=mm.end());
    return i->second;
}

template class RkCache<6>;
//
//  Calculate and store 2 electron radial repulsion (Slater) integrals for all valules of k.
//
template <size_t K> RkEngine<K>::RkEngine(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t _LMax, const GLCache& gl, const RkCache<K>& rkcache)
 : LMax(_LMax), Rabcd_k(VecLimits(0,2*LMax))
 {
    sp_t a=splines[ia];
    sp_t b=splines[ib];
    sp_t c=splines[ic];
    sp_t d=splines[id];
    auto& sa=a.getSupport();
    auto& sb=b.getSupport();
    auto& sc=c.getSupport();
    auto& sd=d.getSupport();

    auto sab=sa.calcIntersection(sb);
    auto scd=sc.calcIntersection(sd);
    assert(sab.containsIntervals());
    assert(scd.containsIntervals());
    //
    //  THis is not the right condition
    // auto sabcd=sab.calcIntersection(scd);
    // assert(sabcd.containsIntervals());

    assert(sc.hasSameGrid(sd));
    assert(sa.hasSameGrid(sb));
    bspline::Grid grid=sa.getGrid();
    // double rmin=std::max(sab.front(),scd.front()),rmax=std::min(sab.back(),scd.back());
    for (size_t k=0;k<=2*LMax;k++)
    {
        std::function< double (double)> wcd1 = [k](double r2)
        {
            return intpow(r2,k+2);
        };
        std::function< double (double)> wcd2 = [k](double r2)
        {
            assert(r2>=0);
            return intpow(r2,1-k);
        };
        std::function< double (double)> Yk1 = [gl,wcd1,c,d,scd](double r1)
        {
            if (r1<scd.front()) return 0.0; 
            double rmax=std::min(r1,scd.back());
            return gl.Integrate(wcd1,c,d,scd.front(),rmax);
        };
        std::function< double (double)> Yk2 = [gl,wcd2,c,d,scd](double r1)
        {
            if (r1>scd.back()) return 0.0;
            double rmin=std::max(r1,scd.front());
            return gl.Integrate(wcd2,c,d,rmin,scd.back());
        };

        std::function< double (double)> wab = [k,Yk1,Yk2] (double r1)
        {
            assert(r1>=0);
            return intpow(r1,1-k)*Yk1(r1)+intpow(r1,k+2)*Yk2(r1);
        };

        if (sab.calcIntersection(scd).containsIntervals())
            Rabcd_k(k)=gl.Integrate(wab,a,b,sab.front(),sab.back()); //Diagonal
        else if (sab.back()<=scd.front())
        {
            std::vector<double> mp=rkcache.find_plus(ia,ib);
            std::vector<double> mm=rkcache.find_minus(ic,id);
            Rabcd_k(k)=mp[k]*mm[k];
        }
        else
        {
            assert(sab.front()>=scd.back());
            std::vector<double> mm=rkcache.find_minus(ia,ib);
            std::vector<double> mp=rkcache.find_plus(ic,id);
            Rabcd_k(k)=mp[k]*mm[k];
        }
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


