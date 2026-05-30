// File: src/BasisSet1/Atom/Evaluators/BSpline/Imp/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <functional>
#include <sstream>

module qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.BasisSet.Atom.Evaluators.IBS;
import Common.Constants;
// import Common.IntPow;

namespace BasisSet::Atom::Evaluators::BSpline
{

using namespace bspline::operators; 
using namespace bspline::integration; 

template<size_t K> using spline_t = bspline::Spline<double, K>;

    // Alternate Overlap without operators.
    // std::function< double (double)> x2 = [](double r) {return r*r;};
    // return gl.Integrate(x2,a,b)*FourPi;

    // Alternate Grad2
    // std::function< double (double)> x1 = [](double r) {return r;};    
    // static const auto T = -X<2>{} * Dx<2>{};
    // assert(la==lb);
    // auto dbdx=transformSpline(bspline::operators::Dx<1>{},b);
    // double Iadb=gl.Integrate(x1,a,dbdx);
    // return (BilinearForm{T}(a,b) - 2*Iadb)*FourPi;

    // Alternate Inv_r1
    // std::function< double (double)> x1 = [](double r) {return r;};
    // return gl.Integrate(x1,a,b)*FourPi;

    // Alternate Inv_r2
    // std::function< double (double)> x0 = [](double r){return 1.0;};
    // return gl.Integrate(x0,a,b)*FourPi;    

//---------------------------------------------------------------------------
//
//  Start member functions.
//

template <size_t K> Evaluator<K>::Evaluator(size_t Ngrid, double rmin, double rmax,const Irrep_QNs::sym_t& ylm) 
: Internal::EvaluatorCommon<K>(Ngrid,rmin,rmax,ylm)
{
    // using l=IBS_Evaluator::l;
    for (size_t n=0;n<=3-l;n++) splines.pop_back(); //For s orbital the last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    // itsGrid=splines[0].getSupport().getGrid();
    ns=norms();
    assert(size()==splines.size());
    assert(size()==ns.size());
};


template <size_t K> std::string Evaluator<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << ">";
    return os.str();
}



std::ostream& operator<<(std::ostream& os, const bspline::Grid<double>& grid)
{
    os << "{";
    for (auto g:grid) os << g << ",";

    return os << "}";
}


template <size_t K> Cache4*    Evaluator<K>::MakeCache4() const
{
    return new BSpline_Cache4<K>(itsGrid);
}

template <size_t K> rvec_t Evaluator<K>::norms() const
{
    size_t N=splines.size();
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=Norm(i); 
    return ret;
}


template <size_t K> rvec_t Evaluator<K>::operator() (const rvec3_t& r) const
{
    rvec_t ret(size());
    double mr=norm(r);
    size_t i=0;
    for (auto s:splines) 
    {
        ret[i]=ns[i]*s(mr);
        ++i;
    }
    return ret;
}

template <size_t K> rvec3vec_t Evaluator<K>::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0) 
    {
        
        ret=rvec3_t(0,0,0);
        return ret; //Cusp at the origin so grad is undefined.
    }
    assert(mr>0);
    ret=r/mr;
    size_t i=0;
    for (auto s:splines) 
    {
        auto dsdx=transformSpline(bspline::operators::Dx<1>{},s);
        ret[i]*=ns[i]*dsdx(mr);
        ++i;
    }
    return ret;
}

template <size_t K> BSpline_Cache4<K>::BSpline_Cache4(const bspline::Grid<double>& grid) 
    : wp([](double r2,size_t k) {return intpow(r2,k+2);})
    , wm([](double r2,size_t k) {return intpow(r2,1-k);})
    , itsMaxl(0)
    , itsGL1D(grid,K+3)
    , itsGL2D(grid,2*K+3,K+3)
    , itsRkCache(0) 
    {
    };

template <size_t K> void BSpline_Cache4<K>::Register(Cache4_Client * eval)
{
    assert(eval);
    Evaluator<K>* geval=dynamic_cast<Evaluator<K>*>(eval);
    geval->Register(&grouper);
    if (geval->Getl()>itsMaxl) itsMaxl=geval->Getl();
    //
    //  At this point we need sweep through all Cacheable* (Rks) in Cache4::cache_t
    //  and check if geval is supported (geval.l <= Rk.LMax).
    //  All unsupport Rks will be removed.  These will then automatically be recreated next time
    //  loop_4 is called.
    //
    Cache4::Register(eval);

    delete itsRkCache;
    itsRkCache=new ::BSpline::RkCache<K>(grouper.unique_spv,itsGL1D, itsMaxl,wp,wm);
}

template <size_t K> Rk*  BSpline_Cache4<K>::Create (size_t ia,size_t ic,size_t ib,size_t id) const
{
    assert(itsRkCache);
    size_t lmax=grouper.LMax(ia,ib,ic,id);
    return new ::BSpline::RkEngine(grouper.unique_spv,ia,ib,ic,id,lmax,itsGL1D,itsGL2D,*itsRkCache,wp,wm);
}

template <size_t K>  size_t BSpline_Cache4<K>::RAMsize() const
{
    size_t ndoubles=Cache4::RAMsize();
    ndoubles+=itsGL1D.RAMsize();
    ndoubles+=itsGL2D.RAMsize();
    ndoubles+=itsRkCache->RAMsize();
    return ndoubles;
}


#define INSTANCEk(k) template class Evaluator<k>;
#include "../Internal/Instance.hpp"
#define INSTANCEk(k) template class BSpline_Cache4<k>;
#include "../Internal/Instance.hpp"

} //namespace