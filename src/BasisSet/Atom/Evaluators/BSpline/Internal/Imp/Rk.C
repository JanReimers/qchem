// File: BSpline::RkEngine.H  4 electron Charge distribution of BSpline orbitals. 
module;
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <algorithm>
#include <cassert>
#include <bspline/Core.h>
module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.Blaze;   // rmat_t + blazem::column() for the pre-sampled B-spline value tables

using std::cout;
using std::endl;

namespace qchem::BSpline
{

template <size_t K> RkCache<K>::RkCache(const std::vector<sp_t>& splines,const GLCache1D& gl1, size_t lmax, const func_t& wp, const func_t& wm, const bspline::Grid<double>& grid)
: itsMomentsPlus (splines.size())
, itsMomentsMinus(splines.size())
{
    for (size_t ia=0;ia<splines.size();ia++)
    {
        if (!(splines[ia].getSupport().getGrid()==grid)) continue; //off-grid spline: not covered by this gl1
        for (size_t ib=ia;ib<splines.size();ib++)
        {
            if (!(splines[ib].getSupport().getGrid()==grid)) continue;
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

    //  The support topology (which of the three cases applies) is independent of k, so decide it ONCE and
    //  put the k loop innermost.  In the hard "all four overlap" case the B-spline VALUES at the fixed GL
    //  nodes are also k-independent (only the weight powers wp/wm carry k), so we sample each spline at its
    //  nodes ONCE and reuse the samples across all k via GLQuadrature::Integrate(w,k,u,v).  Spline
    //  evaluation (binary search for the knot span + polynomial) is the dominant cost, so this removes a
    //  ~(2*LMax+1)x redundancy; the sampled Integrate reproduces the old Integrate(f=w*u*v) term-for-term,
    //  so results stay bit-identical.
    // Sample a spline at every node of a quadrature into a (contiguous, column-major) matrix column.
    auto sampleCol=[](const sp_t& s, const GLQuadrature& g, auto col)
    {
        for (size_t j=0;j<g.size();j++) col[j]=s(g.node(j));
    };
    if (Sabcd.containsIntervals())
    {
        // THis is hard/hot loop for the whole process.  abcd all overlap and no factoring is possible.
        const size_t ab0=sab.getStartIndex(), ab1=sab.getEndIndex()-1;   // outer (ab) interval range
        const size_t cd0=scd.getStartIndex(), cd1=scd.getEndIndex()-1;   // (cd) support interval range
        // Sample matrices: rows index the quadrature node, columns index the iterated interval/node, so a
        // per-index sample is one contiguous column (rmat_t is column-major).  Sampled once, reused across
        // all k -- spline evaluation (the binary knot-span search) is the dominant cost.
        const size_t n1=gl1[ab0].size(), nO=gl2.GL1d()[ab0].size(), n2=gl2.find_grid_gl(ab0,0).size();
        rmat_t A1(n1,ab1-ab0), B1(n1,ab1-ab0), AO(nO,ab1-ab0), BO(nO,ab1-ab0);  // a,b on gl1 and gl2 outer
        for (size_t iab=ab0;iab<ab1;iab++)
        {
            const size_t t=iab-ab0;
            sampleCol(a,gl1[iab],       blazem::column(A1,t)); sampleCol(b,gl1[iab],       blazem::column(B1,t));
            sampleCol(a,gl2.GL1d()[iab],blazem::column(AO,t)); sampleCol(b,gl2.GL1d()[iab],blazem::column(BO,t));
        }
        rmat_t C1(n1,cd1-cd0), D1(n1,cd1-cd0);                                    // c,d on gl1 over (cd)
        for (size_t i=cd0;i<cd1;i++) { sampleCol(c,gl1[i],blazem::column(C1,i-cd0)); sampleCol(d,gl1[i],blazem::column(D1,i-cd0)); }

        rmat_t CGG(n2,nO), DGG(n2,nO), CGg(n2,nO), DGg(n2,nO);  // c,d on the 2D inner sub-grids (per iab)
        for (size_t iab=ab0;iab<ab1;iab++)
        {
            const size_t t=iab-ab0;
            const GLQuadrature& glab =gl1[iab];
            const GLQuadrature& glout=gl2.GL1d()[iab];
            auto a1=blazem::column(A1,t); auto b1=blazem::column(B1,t); auto aO=blazem::column(AO,t); auto bO=blazem::column(BO,t);
            for (size_t i1=0;i1<glout.size();i1++)
            {
                sampleCol(c,gl2.find_grid_gl(iab,i1),  blazem::column(CGG,i1)); sampleCol(d,gl2.find_grid_gl(iab,i1),  blazem::column(DGG,i1));
                sampleCol(c,gl2.find_gl_grid(i1,iab+1),blazem::column(CGg,i1)); sampleCol(d,gl2.find_gl_grid(i1,iab+1),blazem::column(DGg,i1));
            }
            for (size_t k=0;k<=2*itsLMax;k++)
            {
                double Iab_p=glab.Integrate(wp,k,a1,b1);
                double Iab_m=glab.Integrate(wm,k,a1,b1);
                // Icd_p sums gl1 intervals [cd0,iab); Icd_m sums (iab,cd1).  The old IntegrateIndex added
                // each interval's OWN Integrate() result, so we accumulate per interval (+= its Integrate)
                // to round identically -- not one flat sum over all nodes.  Intervals outside (cd)'s support
                // contribute exactly 0 (c or d vanishes), so skipping them is the same as adding 0.0 and
                // keeps C1/D1 indexing inside [cd0,cd1).
                double Icd_p=0.0;
                for (size_t i=cd0;i<std::min(iab,cd1);i++) Icd_p+=gl1[i].Integrate(wp,k,blazem::column(C1,i-cd0),blazem::column(D1,i-cd0));
                double Icd_m=0.0;
                for (size_t i=std::max(iab+1,cd0);i<cd1;i++) Icd_m+=gl1[i].Integrate(wm,k,blazem::column(C1,i-cd0),blazem::column(D1,i-cd0));
                double Idiag=glout.Integrate([&](double r1, size_t i1)
                {
                    double Yk1=gl2.find_grid_gl(iab,i1)  .Integrate(wp,k,blazem::column(CGG,i1),blazem::column(DGG,i1));
                    double Yk2=gl2.find_gl_grid(i1,iab+1).Integrate(wm,k,blazem::column(CGg,i1),blazem::column(DGg,i1));
                    return (wm(r1,k)*Yk1+wp(r1,k)*Yk2)*aO[i1]*bO[i1];
                });
                Rabcd_k[k]+=Idiag + Iab_m*Icd_p + Iab_p*Icd_m;
            }
        }
    }
    // (ab) has no overlap with (cd): the 2D integral factors -- use the precomputed rkcache moment tables.
    else if (sab.back()<=scd.front())
        for (size_t k=0;k<=2*itsLMax;k++) Rabcd_k[k]=rkcache.plus (ia,ib)[k]*rkcache.minus(ic,id)[k];
    else
    {
        assert(sab.front()>=scd.back());
        for (size_t k=0;k<=2*itsLMax;k++) Rabcd_k[k]=rkcache.minus(ia,ib)[k]*rkcache.plus (ic,id)[k];
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


