// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 
module;
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <cassert>
#include <bspline/Core.h>
module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;

using std::cout;
using std::endl;

namespace qchem::BSpline
{

template <size_t K> RkCache<K>::RkCache(const std::vector<sp_t>& splines,const GLCache1D& gl1, size_t lmax, const func_t& wp, const func_t& wm)
: itsMomentsPlus (splines.size())
, itsMomentsMinus(splines.size())
{
    for (size_t ia=0;ia<splines.size();ia++)
    {
        for (size_t ib=ia;ib<splines.size();ib++)
        {
            std::vector<double> mp,mm;
            for (size_t k=0;k<=2*lmax;k++)
            {
                mp.push_back(gl1.Integrate(wp,splines[ia],splines[ib],k)); //gl1 knows how to skip zero overlap.
                mm.push_back(gl1.Integrate(wm,splines[ia],splines[ib],k));
            }
            itsMomentsPlus (ia,ib)=mp;
            itsMomentsMinus(ia,ib)=mm;
        }
    }
}

template <size_t K>  size_t RkCache<K>::RAMsize() const
{
    size_t ndoubles=0;
    for (size_t ia=0;ia<itsMomentsPlus.rows();ia++)
        for (size_t ib=0;ib<itsMomentsPlus.columns();ib++)
        {
            ndoubles+=itsMomentsPlus (ia,ib).size();
            ndoubles+=itsMomentsMinus(ia,ib).size();
        }
    return ndoubles;
}

//
//  Calculate and store 2 electron radial repulsion (Slater) integrals for all valules of k.
//
template <size_t K> RkEngine<K>::RkEngine(const std::vector<sp_t>& splines, size_t ia, size_t ib, size_t ic, size_t id, size_t _LMax
    , const GLCache1D& gl1,const GLCache2D& gl2, const RkCache<K>& rkcache, const func_t& wp, const func_t& wm)
 : itsLMax(_LMax), Rabcd_k(2*itsLMax+1,0.0)
 {
    sp_t   a=splines[ia]   ,  b=splines[ib]   ,  c=splines[ic]   ,  d=splines[id];
    auto &sa=a.getSupport(),&sb=b.getSupport(),&sc=c.getSupport(),&sd=d.getSupport();

    bspline::Support<double> sab=sa.calcIntersection(sb), scd=sc.calcIntersection(sd);
    if (!sab.containsIntervals())
    {
        std::cerr << "BSpline::RkEngine<K>::RkEngine no ab overlap!" << std::endl;
        return;
    }
    if (!scd.containsIntervals())
    {
        std::cerr << "BSpline::RkEngine<K>::RkEngine no cd overlap!" << std::endl;
        return;
    }
    assert(sab.containsIntervals());
    assert(scd.containsIntervals());
    auto Sabcd=sab.calcIntersection(scd);
    
    assert(sc.hasSameGrid(sd));
    assert(sa.hasSameGrid(sb));
    bspline::Grid grid=sa.getGrid();
    for (size_t k=0;k<=2*itsLMax;k++)
    {
        // THis is hard/hot loop for the whole process.  abcd all overlap and no factoring is possible.
        if (Sabcd.containsIntervals())
        {
            std::function< double (double)> wpab=[k,&wp,&a,&b](double r1) {return wp(r1,k)*a(r1)*b(r1);};
            std::function< double (double)> wmab=[k,&wm,&a,&b](double r1) {return wm(r1,k)*a(r1)*b(r1);};
            std::function< double (double)> wpcd=[k,&wp,&c,&d](double r2) {return wp(r2,k)*c(r2)*d(r2);};
            std::function< double (double)> wmcd=[k,&wm,&c,&d](double r2) {return wm(r2,k)*c(r2)*d(r2);};
            for (size_t iab=sab.getStartIndex();iab<sab.getEndIndex()-1;iab++)
            {
                double rab=grid[iab],rab1=grid[iab+1];
                double Iab_p=gl1[iab].Integrate(wpab);
                double Iab_m=gl1[iab].Integrate(wmab);
                double Icd_p=gl1.IntegrateIndex(wpcd,scd.getStartIndex(),iab);
                double Icd_m=gl1.IntegrateIndex(wmcd,iab+1,scd.getEndIndex()-1);
                // These need to be in the loop because they capture iab
                std::function< double (double,size_t)> Yk1_diag = [&gl2,&wpcd,rab,iab](double r1, size_t i1)
                {
                    // double ret=gl2.find_grid_gl(iab,i1).Integrate(wpcd);
                    // auto gl=GLQuadrature(rab,r1,K+1);
                    // // cout << gl.Integrate(wpcd) << " " << ret << " " << gl.Integrate(wpcd)-ret << endl;
                    // assert(gl.Integrate(wpcd)==ret);
                    // return ret;
                    // return gl.Integrate(wpcd);
                    return gl2.find_grid_gl(iab,i1).Integrate(wpcd);
                };
                std::function< double (double,size_t)> Yk2_diag = [&gl2,&wmcd,rab1,iab](double r1, size_t i1)
                {
                    // double ret=gl2.find_gl_grid(i1,iab+1).Integrate(wmcd);
                    // auto gl=GLQuadrature(r1,rab1,K+1);
                    // assert(gl.Integrate(wmcd)==ret);
                    // return ret;
                    // return gl.Integrate(wmcd);
                    return gl2.find_gl_grid(i1,iab+1).Integrate(wmcd);
                };
                std::function< double (double,size_t)> wab_diag1 = [k,&a,&b,&wp,&wm,&Yk1_diag,&Yk2_diag](double r1, size_t i1)
                {
                    return (wm(r1,k)*Yk1_diag(r1,i1)+wp(r1,k)*Yk2_diag(r1,i1))*a(r1)*b(r1);
                };
                double Idiag=gl2.GL1d()[iab].Integrate(wab_diag1);
                //  std::cout  << "iab=" << std::setw(2) << iab << " Idiag=" << std::setw(12) << Idiag 
                // << " Iab_m=" << std::setw(12) << Iab_m << " Icd_p=" << std::setw(12) << Icd_p 
                // << " Iab_p=" << std::setw(12) << Iab_p << " Icd_m=" << std::setw(12) << Icd_m << std::endl;
                Rabcd_k[k]+=Idiag + Iab_m*Icd_p + Iab_p*Icd_m;
            }
            
        }
        // (ab) has no overlap (cd) and 2D integral factors.  Use look up tables in rkcache to
        // integral factores.
        else if (sab.back()<=scd.front())
        {
            std::vector<double> mp=rkcache.plus (ia,ib);
            std::vector<double> mm=rkcache.minus(ic,id);
            Rabcd_k[k]=mp[k]*mm[k];
        }
        else
        {
            assert(sab.front()>=scd.back());
            std::vector<double> mm=rkcache.minus(ia,ib);
            std::vector<double> mp=rkcache.plus (ic,id);
            Rabcd_k[k]=mp[k]*mm[k];
        }
    }
 }

 template <size_t K> double RkEngine<K>::DirectR0  () const
 {
    return Rabcd_k[0];
 }

 template <size_t K> double RkEngine<K>::DirectRk  (size_t la,size_t lc, const rvec11_t& Ak) const
 {
    assert(la>=0);
    assert(lc>=0);
    assert(la<=itsLMax);
    assert(lc<=itsLMax);
    double Rk(0.0);
    for (size_t k=0;k<=2*std::min(la,lc);k+=2)
    {
        Rk+=Rabcd_k[k]*Ak[k]; 
    }
    return Rk;
 }
 template <size_t K> double RkEngine<K>::ExchangeRk(size_t la,size_t lb, const rvec11_t& Ak) const
 {
    assert(la>=0);
    assert(lb>=0);
    assert(la<=itsLMax);
    assert(lb<=itsLMax);
    int kmin=std::abs((int)la-(int)lb);
    int kmax=la+lb;
    double Rk(0.0);
    for (int k=kmin;k<=kmax;k+=2)
    {
        assert((k+la+lb)%2==0);
        Rk+=Rabcd_k[k]*Ak[k]; 
    }
    return Rk;
 }

template <size_t K> size_t RkEngine<K>::RAMsize() const
{
    return Rabcd_k.size();
}

 
#define INSTANCEk(k) template class RkCache<k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class RkEngine<k>;
#include "../Instance.hpp"

} //namespace


