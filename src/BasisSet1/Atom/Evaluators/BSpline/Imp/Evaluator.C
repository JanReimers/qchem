// File: src/BasisSet1/Atom/Evaluators/BSpline/Imp/Evaluator.C
module;
#include <bspline/Core.h>
#include <cassert>
#include <iostream>
#include <functional>
#include <sstream>

module qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.Symmetry.Spherical;
import qchem.Math;

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

template <size_t K> Evaluator<K>::Evaluator(size_t Ngrid, double rmin, double rmax,const sym_t& ylm)
: Internal::EvaluatorCommon<K>(Ngrid,rmin,rmax,ylm)
, NR_Angular(ylm)
{
    for (size_t n=0;n<=3-Getl();n++) splines.pop_back(); //For s orbital the last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    ns=norms();
    assert(ns.size()==splines.size());
    assert(ns.size()==ns.size());
};


template <size_t K> std::string Evaluator<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << ">";
    return os.str();
}




template <size_t K> Cache4*    Evaluator<K>::MakeCache4() const
{
    return new Internal::Cache4<K>
    (
        itsGrid,
        [](double r2,size_t k) {return intpow(r2,k+2);},
        [](double r2,size_t k) {return intpow(r2,1-k);},
        3
    );
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



#define INSTANCEk(k) template class Evaluator<k>;
#include "../Internal/Instance.hpp"

} //namespace