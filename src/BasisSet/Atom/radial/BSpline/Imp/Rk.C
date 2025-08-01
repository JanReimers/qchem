// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 
module;
#include <iostream>
#include <vector>
#include <functional>
#include <cassert>
#include <bspline/Core.h>
module qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.Basisset.Atom.radial.BSpline.GLQuadrature;

import Common.IntPow;
import oml;

using std::cout;
using std::endl;

namespace BSpline
{
template <size_t K> RkCache<K>::RkCache(const std::vector<sp_t>& splines,const GLCache& gl, size_t lmax)
{
    for (size_t ia=0;ia<splines.size();ia++)
        for (size_t ib=ia;ib<splines.size();ib++)
        {
            // TODO skip for zero overlap
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
    assert(i!=mm.end()); //If not found return zero, there was no support overlap.
    return i->second;
}

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

    bspline::Support<double> sab=sa.calcIntersection(sb);
    bspline::Support<double> scd=sc.calcIntersection(sd);
    assert(sab.containsIntervals());
    assert(scd.containsIntervals());
    auto Sabcd=sab.calcIntersection(scd);
    
    assert(sc.hasSameGrid(sd));
    assert(sa.hasSameGrid(sb));
    bspline::Grid grid=sa.getGrid();
    // cout << "grid = ";
    // for (auto r:grid) cout << r << ", ";
    // cout << endl;
    // // double rmin=std::max(sab.front(),scd.front()),rmax=std::min(sab.back(),scd.back());
    for (size_t k=0;k<=2*LMax;k++)
    {
        std::function< double (double)> wp = [k](double r2)
        {
            return intpow(r2,k+2);
        };
        std::function< double (double)> wm = [k](double r2)
        {
            assert(r2>=0);
            return intpow(r2,1-k);
        };
        
        if (Sabcd.containsIntervals())
        {
            cout.precision(6);
            // cout << sab.getStartIndex() << " " << sab.getEndIndex() << " " << sab.numberOfIntervals() << " " << scd.getStartIndex() << " " << scd.getEndIndex() << endl;
            double RkOff=0.0, RkDiag=0.0;
            // #pragma omp parallel for collapse(1) 
            for (size_t iab=sab.getStartIndex();iab<sab.getEndIndex()-1;iab++)
            {
                double rab=grid[iab],rab1=grid[iab+1];
                double Iab_p=gl.IntegrateIndex(wp,a,b,iab);
                double Iab_m=gl.IntegrateIndex(wm,a,b,iab);
                double Icd_p=gl.IntegrateIndex(wp,c,d,scd.getStartIndex(),iab);
                double Icd_m=gl.IntegrateIndex(wm,c,d,iab+1,scd.getEndIndex()-1);
                std::function< double (double)> Yk1_diag = [&gl,&wp,&c,&d,rab](double r1)
                {
                    assert(rab<=r1);
                    const GLQuadrature& gl1=gl.find(rab,r1);
                    std::function< double (double)> f=[&wp,&c,&d](double r) {return wp(r)*c(r)*d(r);};
                    return gl1.Integrate(f);
                };
                std::function< double (double)> Yk2_diag = [&gl,&wm,&c,&d,rab1](double r1)
                {
                    assert(r1<=rab1);
                    const GLQuadrature& gl1=gl.find(r1,rab1);
                    std::function< double (double)> f=[&wm,&c,&d](double r) {return wm(r)*c(r)*d(r);};
                    return gl1.Integrate(f);
                };

                std::function< double (double)> wab_diag = [k,&Yk1_diag,&Yk2_diag] (double r1)
                {
                    assert(r1>=0);
                    return intpow(r1,1-k)*Yk1_diag(r1)+intpow(r1,k+2)*Yk2_diag(r1);
                };
                double Idiag=gl.IntegrateIndex(wab_diag,a,b,iab);
                // # pragma omp critical
                RkOff+=Iab_m*Icd_p + Iab_p*Icd_m;
                RkDiag+=Idiag;
            }
            Rabcd_k(k)=RkDiag + RkOff;
        }
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

 template <size_t K> double RkEngine<K>::Coulomb_R0() const
 {
    return Rabcd_k(0);
 }

 template <size_t K> Vector<double> RkEngine<K>::Coulomb_Rk(size_t la,size_t lc) const
 {
    assert(la>=0);
    assert(lc>=0);
    assert(la<=LMax);
    assert(lc<=LMax);
    Vector<double> Rk(VecLimits(1,la+lc+1),0.0);
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
 
#define INSTANCEk(k) template class RkCache<k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class RkEngine<k>;
#include "../Instance.hpp"

} //namespace


