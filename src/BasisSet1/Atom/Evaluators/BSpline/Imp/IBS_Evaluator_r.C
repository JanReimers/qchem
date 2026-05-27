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


template <size_t K> BSpline_r_IBS_Evaluator<K>::BSpline_r_IBS_Evaluator(size_t Ngrid, double _rmin, double _rmax,const Irrep_QNs::sym_t& ylm) 
: IBS_Evaluator(ylm), rmin(_rmin), rmax(_rmax) , itsGrid({0,1})
{
    knots=MakeLogKnots(Ngrid,rmin,rmax);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    splines.pop_back(); //Last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    itsGrid=splines[0].getSupport().getGrid();
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    itsGL1D.reset(new GLCache1D(itsGrid,K+1));
    ns=norms();
    assert(size()==splines.size());
};

 template <size_t K> std::vector<double> BSpline_r_IBS_Evaluator<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
{
    assert(Ngrid>1);
    std::vector<double> knots;
    size_t numberOfZeros = 1;

    if (K + 1 > l)  numberOfZeros = K + 1 - l;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

    // logarithmic step
    const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid-1));
    for (size_t i = 0; i < Ngrid-1; i++) //Skip 0.0 and rmax
        knots.push_back(rmin * pow(step, i));
    
    // cout << Ngrid << " " << l << " " << numberOfZeros << " ";
    if (numberOfZeros>Ngrid-numberOfZeros) numberOfZeros=Ngrid-numberOfZeros;
    if (numberOfZeros<1) numberOfZeros=1;
    // cout << numberOfZeros << endl;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
    // std::cout << knots << std::endl;
    return knots;
}

template <size_t K> void BSpline_r_IBS_Evaluator<K>::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<SplineGrouper<K>*>(_grouper);
    assert(grouper);
    for (auto s:splines) es_indices.push_back(grouper->Insert(s,l));
}


template <size_t K> std::string BSpline_r_IBS_Evaluator<K>::Name () const
{
    std::ostringstream os;
    os << "BSpline<" << K << "> 1/r ";
    return os.str();
}

template <size_t K> std::string BSpline_r_IBS_Evaluator<K>::RadialID () const
{
    assert(splines.size()>0);
    std::ostringstream os;
    bspline::Grid<double> grid=splines[0].getSupport().getGrid();
    os << Name() << " grid: N=" << grid.size() << " {";
    assert(grid.size()>2);
    os << grid[0] << "," << grid[1] << "," << grid[2] << " ... " << grid[grid.size()-1];
    os << "}";
    return os.str();
}
template <size_t K> std::string BSpline_r_IBS_Evaluator<K>::RadialType() const
{
    std::ostringstream os;
    os << "BS<" << K << "> 1/r grid=" << splines[0].getSupport().getGrid();;
    return os.str();
}
template <size_t K> Cache4*    BSpline_r_IBS_Evaluator<K>::MakeCache4() const
{
    return new BSpline_r_Cache4<K>(itsGrid);
}

template <size_t K> rvec_t BSpline_r_IBS_Evaluator<K>::norms() const
{
    size_t N=splines.size();
    // std::cout << "BSpline_IBS_Evaluator<K>::norms() N=" << N << std::endl;
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(BilinearForm{IdentityOperator{}}(splines[i],splines[i])*FourPi); 
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

template <size_t K> std::ostream&  BSpline_r_IBS_Evaluator<K>::Write(std::ostream& os) const
{
    return os << " with " << size() << " basis functions, {" << rmin << " ... " << rmax << "}" << std::endl;
}


#define INSTANCEk(k) template class BSpline_r_IBS_Evaluator<k>;
#include "../Internal/Instance.hpp"

} //namespace