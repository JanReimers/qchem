// File: BSpline/IEC.C Common IE client code for all atom BSpline basis sets and IEs.
module;
#include <vector>
#include <cassert>
#include <cmath>
#include <bspline/Core.h>

module qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.Basisset.Atom.radial.BSpline.GLQuadrature;
import qchem.stl_io;

namespace BSpline
{
template <size_t K> IrrepIEClient<K>::IrrepIEClient( size_t Ngrid, double _rmin, double _rmax, size_t _l)
: AtomIrrepIEClient(Ngrid)
,rmin(_rmin), rmax(_rmax), itsGL(0)
{
    l=_l;
    std::vector<double> knots=MakeLogKnots(Ngrid,rmin,rmax);
    // std::cout << "Knots=" << knots << std::endl;
    splines=bspline::generateBSplines<K>(knots);
    auto grid=splines[0].getSupport().getGrid();
    // std::cout << "Grid = " << grid.size() << "    ";
    // for (auto r:grid) std::cout << r << ",";
    // std::cout << std::endl;
    ns.SetLimits(splines.size());
    itsGL=new GLCache(grid,K+2*l);
    assert(itsGL);
}

template <size_t K> IrrepIEClient<K>::IrrepIEClient( size_t Ngrid, double _rmin, double _rmax, size_t _l, const std::vector<int>& _ml)
: IrrepIEClient(Ngrid,_rmin,_rmax,_l)
{
    ml=_ml;
}


template <size_t K> IrrepIEClient<K>::~IrrepIEClient()
{
    delete itsGL;
}

template <size_t K> std::vector<double> IrrepIEClient<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
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

#define INSTANCEk(k) template class IrrepIEClient<k>;
#include "../Instance.hpp"

} //namespace


