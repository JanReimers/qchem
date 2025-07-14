#ifndef SLATERINTEGRALS_H
#define SLATERINTEGRALS_H

#include "Common/IntPower.H"
#include <cmath>
#include <cassert>

import Common.Constants;
import Common.Factorials;

namespace Slater
{
//##############################################################################
//      /
//  4Pi |  r^l exp(-a*r) r^2 dr
//     /
inline double Integral(double a, int l)
{
    assert(l>=-2);
    return FourPi* qchem::Fact[l+2]/pow(a,l+3);
}

//--------------------------------------------------------------------
//
//  1/sqrt(<a|a>)
//
inline double Norm(double a, int n)
{
    return 1.0/std::sqrt(Slater::Integral(2*a,2*n-2));
}

} //namespace

#endif // SLATERINTEGRALS_H
