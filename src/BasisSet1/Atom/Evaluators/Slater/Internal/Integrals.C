// File: BasisSet1/Atom/Evaluators/Slater/Internal/Integrals.C   Slater radial integral functions.
module;
#include <cassert>
export module qchem.BasisSet.Atom.Evaluators.Slater.Internal.Integrals; 
import qchem.Math;

export namespace Slater
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
    return 1.0/sqrt(Slater::Integral(2*a,2*n-2));
}

} //namespace

