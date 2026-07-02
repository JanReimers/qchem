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
    //  nodes ONCE and reuse the samples across all k.  Spline evaluation (binary search for the knot span +
    //  polynomial) is the dominant cost, so this removes a ~(2*LMax+1)x redundancy.  A hand-rolled
    //  sum_i weight(i)*wp(node(i),k)*u(i)*v(i) with cached u=spline(node) matches the old Integrate() calls
    //  term-for-term (same order and associativity), so results stay bit-identical.
    auto sample=[](const sp_t& s, const GLQuadrature& g)
    {
        rvec_t v(g.size());
        for (size_t j=0;j<g.size();j++) v[j]=s(g.node(j));
        return v;
    };
    if (Sabcd.containsIntervals())
    {
        // THis is hard/hot loop for the whole process.  abcd all overlap and no factoring is possible.
        const size_t ab0=sab.getStartIndex(), ab1=sab.getEndIndex()-1;   // outer (ab) interval range
        const size_t cd0=scd.getStartIndex(), cd1=scd.getEndIndex()-1;   // (cd) support interval range
        // Engine-level samples on the 1D grids (independent of iab): a,b on gl1 and on the gl2 outer grid;
        // c,d on gl1 over (cd)'s support.
        std::vector<rvec_t> A1(ab1-ab0), B1(ab1-ab0), AO(ab1-ab0), BO(ab1-ab0);
        for (size_t iab=ab0;iab<ab1;iab++)
        {
            const size_t t=iab-ab0;
            A1[t]=sample(a,gl1[iab]);          B1[t]=sample(b,gl1[iab]);
            AO[t]=sample(a,gl2.GL1d()[iab]);   BO[t]=sample(b,gl2.GL1d()[iab]);
        }
        std::vector<rvec_t> C1(cd1-cd0), D1(cd1-cd0);
        for (size_t i=cd0;i<cd1;i++) { C1[i-cd0]=sample(c,gl1[i]); D1[i-cd0]=sample(d,gl1[i]); }

        for (size_t iab=ab0;iab<ab1;iab++)
        {
            const size_t t=iab-ab0;
            const GLQuadrature& glab =gl1[iab];
            const GLQuadrature& glout=gl2.GL1d()[iab];
            const rvec_t& a1=A1[t]; const rvec_t& b1=B1[t]; const rvec_t& aO=AO[t]; const rvec_t& bO=BO[t];
            // Per-iab samples of c,d on the 2D inner sub-grids (these depend on iab and the outer node i1).
            std::vector<rvec_t> CGG(glout.size()), DGG(glout.size()), CGg(glout.size()), DGg(glout.size());
            for (size_t i1=0;i1<glout.size();i1++)
            {
                const GLQuadrature& gg=gl2.find_grid_gl(iab,i1);
                CGG[i1]=sample(c,gg); DGG[i1]=sample(d,gg);
                const GLQuadrature& gr=gl2.find_gl_grid(i1,iab+1);
                CGg[i1]=sample(c,gr); DGg[i1]=sample(d,gr);
            }
            for (size_t k=0;k<=2*itsLMax;k++)
            {
                // Every term below reproduces the old GLQuadrature::Integrate(f) sum exactly: weight(j)
                // times the SAME integrand f (e.g. wp(x,k)*a*b, evaluated left-to-right) -- the extra
                // parentheses keep weight*(wp*a*b) rather than (weight*wp)*(a*b), which would round
                // differently and, through the nonlinear SCF, drift the converged energy.
                double Iab_p=0.0, Iab_m=0.0;
                for (size_t j=0;j<glab.size();j++)
                {
                    const double x=glab.node(j), w=glab.weight(j);
                    Iab_p+=w*(wp(x,k)*a1[j]*b1[j]);
                    Iab_m+=w*(wm(x,k)*a1[j]*b1[j]);
                }
                // Icd_p sums gl1 intervals [cd0,iab); Icd_m sums (iab,cd1).  The old IntegrateIndex added
                // each interval's OWN Integrate() result, so we must accumulate per interval (an inner sum
                // per interval, then += it) to round identically -- not one flat sum over all nodes.
                // Intervals outside (cd)'s support contribute exactly 0 (c or d vanishes), so skipping them
                // is the same as adding 0.0 and keeps C1/D1 indexing inside [cd0,cd1).
                double Icd_p=0.0;
                for (size_t i=cd0;i<std::min(iab,cd1);i++)
                {
                    const GLQuadrature& g=gl1[i]; const rvec_t& cu=C1[i-cd0]; const rvec_t& du=D1[i-cd0];
                    double s=0.0;
                    for (size_t j=0;j<g.size();j++) s+=g.weight(j)*(wp(g.node(j),k)*cu[j]*du[j]);
                    Icd_p+=s;
                }
                double Icd_m=0.0;
                for (size_t i=std::max(iab+1,cd0);i<cd1;i++)
                {
                    const GLQuadrature& g=gl1[i]; const rvec_t& cu=C1[i-cd0]; const rvec_t& du=D1[i-cd0];
                    double s=0.0;
                    for (size_t j=0;j<g.size();j++) s+=g.weight(j)*(wm(g.node(j),k)*cu[j]*du[j]);
                    Icd_m+=s;
                }
                double Idiag=0.0;
                for (size_t i1=0;i1<glout.size();i1++)
                {
                    double Yk1=0.0;
                    { const GLQuadrature& gg=gl2.find_grid_gl(iab,i1); const rvec_t& cu=CGG[i1]; const rvec_t& du=DGG[i1];
                      for (size_t j=0;j<gg.size();j++) Yk1+=gg.weight(j)*(wp(gg.node(j),k)*cu[j]*du[j]); }
                    double Yk2=0.0;
                    { const GLQuadrature& gr=gl2.find_gl_grid(i1,iab+1); const rvec_t& cu=CGg[i1]; const rvec_t& du=DGg[i1];
                      for (size_t j=0;j<gr.size();j++) Yk2+=gr.weight(j)*(wm(gr.node(j),k)*cu[j]*du[j]); }
                    const double r1=glout.node(i1);
                    Idiag+=glout.weight(i1)*((wm(r1,k)*Yk1+wp(r1,k)*Yk2)*aO[i1]*bO[i1]);
                }
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


