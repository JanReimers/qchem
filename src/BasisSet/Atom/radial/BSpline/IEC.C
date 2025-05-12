// File: BSpline/IEC.C Common IE client code for all atom BSpline basis sets and IEs.

#include "Imp/BasisSet/Atom/radial/BSpline/IEC.H"

namespace BSpline
{
template <size_t K> IrrepIEClient<K>::IrrepIEClient( size_t Ngrid, double _rmin, double _rmax, size_t _l, int _m)
: rmin(_rmin), rmax(_rmax), l(_l), m(_m)
{
    std::vector<double> knots=MakeLogKnots(Ngrid,rmin,rmax);
    splines=bspline::generateBSplines<K>(knots);
    ns.SetLimits(splines.size());
}

template <size_t K> IrrepIEClient<K>::~IrrepIEClient()
{
    std::cout << "BSpline::IrrepIEClient<K>::~IrrepIEClient()" << std::endl;
}

template <size_t K> std::vector<double> IrrepIEClient<K>::MakeLogKnots(size_t Ngrid, double rmin, double rmax)
{
    std::vector<double> knots;
    size_t numberOfZeros = 1;

    if (K + 1 > l)  numberOfZeros = K + 1 - l;

    for (size_t i = 0; i < numberOfZeros; i++) knots.push_back(0.0);

    // logarithmic step
    const double step =  pow(rmax / rmin, 1 / static_cast<double>(Ngrid));
    for (size_t i = 0; i <= Ngrid; i++) 
        knots.push_back(rmin * pow(step, i));
    return knots;
}

template class IrrepIEClient<6>;

} //namespace