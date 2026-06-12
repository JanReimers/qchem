// File: BasisSet1/Atom/Evaluators/BSpline/Imp/IBS_r_Evaluator.C
module;
#include <bspline/Core.h>
#include <cassert>
#include <iostream>
#include <functional>
#include <sstream>

module qchem.BasisSet.Atom.Evaluators.BSpline.IBS_r;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;
import qchem.Symmetry.Spherical;
import qchem.Math;

namespace BasisSet::Atom::Evaluators::BSpline
{

using namespace bspline::operators; 
using namespace bspline::integration; 

template<size_t K> using spline_t = bspline::Spline<double, K>;


template <size_t K> Evaluator_r<K>::Evaluator_r(size_t Ngrid, double rmin, double rmax,const sym_t& ylm)
: Evaluators::Evaluator(Symmetry::Getl(ylm))
, Internal::EvaluatorCommon<K>(Ngrid,rmin,rmax,ylm)
, NR_Angular(Symmetry::Getmls(ylm))
{
    splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    splines.pop_back(); //Last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    itsGL1D.reset(new GLCache1D(itsGrid,K+1));
    ns=norms();
    assert(ns.size()==splines.size());
};

//  template <size_t K> std::vector<double> Evaluator_r<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
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

template <size_t K> std::string Evaluator_r<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << "> 1/r ";
    return os.str();
}

template <size_t K> Cache4* Evaluator_r<K>::MakeCache4() const
{
    return new Internal::Cache4<K>
    (
        itsGrid,
        [](double r2,size_t k) {return intpow(r2,k);},
        [](double r2,size_t k) {return intpow(r2,-1-k);},
        1
    );
}

template <size_t K> rvec_t Evaluator_r<K>::norms() const
{
    size_t N=splines.size();
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=Norm(i); 
    return ret;
}

template <size_t K> rvec_t Evaluator_r<K>::operator() (const rvec3_t& r) const
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

template <size_t K> rvec3vec_t Evaluator_r<K>::Gradient(const rvec3_t& r) const
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

#define INSTANCEk(k) template class Evaluator_r<k>;
#include "../Internal/Instance.hpp"

} //namespace