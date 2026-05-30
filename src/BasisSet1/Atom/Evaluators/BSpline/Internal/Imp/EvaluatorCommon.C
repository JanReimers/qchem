// File: BasisSet1/Atom/Evaluators/BSpline/Internal/Imp/IBS_EvaluatorCommon.C
module;
#include <vector>
#include <bspline/Core.h>
#include <cmath>
#include <cassert>

module qchem.BasisSet.Atom.Evaluators.BSpline.Internal.Common;
import qchem.BasisSet.Atom.Evaluators.BSpline.Internal.SplineGrouper;

namespace BasisSet::Atom::Evaluators::BSpline::Internal
{

template <size_t K> EvaluatorCommon<K>::EvaluatorCommon(size_t Ngrid, double _rmin, double _rmax,const Irrep_QNs::sym_t& ylm) 
: IBS_Evaluator(ylm)
, rmin(_rmin), rmax(_rmax) , itsGrid({0,1})
{
    knots=MakeLogKnots(Ngrid,rmin,rmax);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    itsGrid=splines[0].getSupport().getGrid();
    // // splines.erase(splines.begin()); //First spline has B(0)=1.0 with violates B(0)=0 boundary condition for 1/r prefactor.
    // for (size_t n=0;n<=3-l;n++) splines.pop_back(); //For s orbital the last spline has B(R)=1.0 with violates B(R)=0 boundary condition for 1/r prefactor.
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    // ns=norms();
    // std::cout << "Evaluator<K>::Evaluator size=" << size() << std::endl;
    // assert(size()==splines.size());
    // assert(size()==ns.size());
};

template <size_t K> std::vector<double> EvaluatorCommon<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
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
    
    // std::cout << Ngrid << " " << l << " " << numberOfZeros << " ";
    // if (numberOfZeros>Ngrid-numberOfZeros) numberOfZeros=Ngrid-numberOfZeros;
    // if (numberOfZeros<1) numberOfZeros=1;
    //     std::cout << numberOfZeros << std::endl;

     for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(rmax);
    // std::cout << knots << std::endl;
    return knots;
}


template <size_t K> void EvaluatorCommon<K>::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<SplineGrouper<K>*>(_grouper);
    assert(grouper);
    for (auto s:splines) es_indices.push_back(grouper->Insert(s,l));
}

template <size_t K> std::string EvaluatorCommon<K>::RadialID () const
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

template <size_t K> std::string EvaluatorCommon<K>::RadialType() const
{
    assert(grid.size()>2);
    std::ostringstream os;
    os << Name() << " grid=<" << itsGrid[0] << "," << itsGrid[1] << " ... " << itsGrid[itsGrid.size()-1] << "}";
    return os.str();
}

template <size_t K> std::ostream&  EvaluatorCommon<K>::Write(std::ostream& os) const
{
    return os << " N= " << size() << " basis functions, {" << rmin << " ... " << rmax << "}" << std::endl;
}

#define INSTANCEk(k) template class EvaluatorCommon<k>;
#include "../../Internal/Instance.hpp"


} //namespace