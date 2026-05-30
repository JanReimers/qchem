// File: BasisSet1/Atom/Evaluators/BSpline/Imp/IBS_r_Evaluator.C
module;
#include <bspline/Core.h>
#include <InvPosition.H> // 1/x^n operator are not provided in bspline package, so a roll our own.
#include <cmath>
#include <cassert>
#include <iostream>
#include <functional>
#include <sstream>

module qchem.BasisSet.Atom.Evaluators.BSpline.IBS_r;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import Common.Constants;

namespace BasisSet::Atom::Evaluators::BSpline
{

using namespace bspline::operators; 
using namespace bspline::integration; 

template<size_t K> using spline_t = bspline::Spline<double, K>;


template <size_t K> BSpline_r_IBS_Evaluator<K>::BSpline_r_IBS_Evaluator(size_t Ngrid, double rmin, double rmax,const Irrep_QNs::sym_t& ylm) 
: Internal::EvaluatorCommon<K>(Ngrid,rmin,rmax,ylm)
{
    splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    splines.pop_back(); //Last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    itsGL1D.reset(new GLCache1D(itsGrid,K+1));
    ns=norms();
    assert(size()==splines.size());
};

//  template <size_t K> std::vector<double> BSpline_r_IBS_Evaluator<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
// {
//     assert(Ngrid>1);
//     std::vector<double> knots;
//     size_t numberOfZeros = 1;

//     if (K + 1 > l)  numberOfZeros = K + 1 - l;

//     for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

//     // logarithmic step
//     const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid-1));
//     for (size_t i = 0; i < Ngrid-1; i++) //Skip 0.0 and rmax
//         knots.push_back(rmin * pow(step, i));
    
//     // cout << Ngrid << " " << l << " " << numberOfZeros << " ";
//     if (numberOfZeros>Ngrid-numberOfZeros) numberOfZeros=Ngrid-numberOfZeros;
//     if (numberOfZeros<1) numberOfZeros=1;
//     // cout << numberOfZeros << endl;

//     for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
//     // std::cout << knots << std::endl;
//     return knots;
// }

template <size_t K> std::string BSpline_r_IBS_Evaluator<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << "> 1/r ";
    return os.str();
}

template <size_t K> Cache4*    BSpline_r_IBS_Evaluator<K>::MakeCache4() const
{
    return new BSpline_r_Cache4<K>(itsGrid);
}

template <size_t K> rvec_t BSpline_r_IBS_Evaluator<K>::norms() const
{
    size_t N=splines.size();
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=Norm(i); 
    return ret;
}

template <size_t K> rvec_t BSpline_r_IBS_Evaluator<K>::operator() (const rvec3_t& r) const
{
    rvec_t ret(size());
    double mr=norm(r);
    size_t i=0;
    for (auto s:splines) 
    {
        ret[i]=ns[i]*s(mr)/mr;
        ++i;
    }
    return ret;
}

template <size_t K> rvec3vec_t BSpline_r_IBS_Evaluator<K>::Gradient(const rvec3_t& r) const
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

template <size_t K> BSpline_r_Cache4<K>::BSpline_r_Cache4(const bspline::Grid<double>& grid) 
: wp([](double r2,size_t k) {return intpow(r2,k);})
, wm([](double r2,size_t k) {return intpow(r2,-1-k);})
, itsMaxl(0)
, itsGL1D(grid,K+1)
, itsGL2D(grid,2*K+1,K+3)
, itsRkCache(0) 
{};

template <size_t K>  void BSpline_r_Cache4<K>::Register(Cache4_Client * eval)
{
    assert(eval);
    BSpline_r_IBS_Evaluator<K>* geval=dynamic_cast<BSpline_r_IBS_Evaluator<K>*>(eval);
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
template <size_t K>  Rk*  BSpline_r_Cache4<K>::Create (size_t ia,size_t ic,size_t ib,size_t id) const
{
        assert(itsRkCache);
    // std::cout << "ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
    size_t lmax=grouper.LMax(ia,ib,ic,id);
    return new ::BSpline::RkEngine(grouper.unique_spv,ia,ib,ic,id,lmax,itsGL1D,itsGL2D,*itsRkCache,wp,wm);
}
template <size_t K>  size_t BSpline_r_Cache4<K>::RAMsize() const
{
    size_t ndoubles=Cache4::RAMsize();
    ndoubles+=itsGL1D.RAMsize();
    ndoubles+=itsGL2D.RAMsize();
    ndoubles+=itsRkCache->RAMsize();
    return ndoubles;
}

#define INSTANCEk(k) template class BSpline_r_IBS_Evaluator<k>;
#include "../Internal/Instance.hpp"
#define INSTANCEk(k) template class BSpline_r_Cache4<k>;
#include "../Internal/Instance.hpp"

} //namespace