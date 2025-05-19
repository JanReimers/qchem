// File: BSpline/IEC.C Common IE client code for all atom BSpline basis sets and IEs.

#include "Imp/BasisSet/Atom/radial/BSpline/IEC.H"
#include "Imp/Containers/stl_io.h"

namespace BSpline
{
template <size_t K> IrrepIEClient<K>::IrrepIEClient( size_t Ngrid, double _rmin, double _rmax, size_t _l, int _m)
: rmin(_rmin), rmax(_rmax), l(_l), m(_m), itsGL(0)
{
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
    for (size_t i = 0; i < Ngrid; i++) 
        knots.push_back(rmin * pow(step, i));
    return knots;
}

template class IrrepIEClient<6>;

} //namespace