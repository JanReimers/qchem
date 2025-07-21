// File: BasisSet/Atom/radial/Gaussian/Integrals.C     Gaussian radial integral functions.
module;
#include <cassert>
#include <cmath>
export module qchem.BasisSet.Atom.radial.GaussianIntegrals;
import Common.Constants;
import Common.Factorials;
//###################################################################
//
//  Hand coded Gaussian integrals.
//
namespace Gaussian
{
double Oddl (double a, int n)
{
    return n==0 ? 1.0/(2*a) : Oddl(a,n-1)*n/a;
}

double Evenl(double a, int n)
{
    return n==0 ? sqrt(Pi/a)/2.0 : Evenl(a,n-1)*(2*n-1)/(2.0*a);
}

//##############################################################################
//      /
//  4Pi |  r^l exp(-a*r^2) r^2 dr
//     /
export double Integral(double a, int l)
{
    assert(l>=-2);
    return FourPi*( l%2 ? Oddl(a,(l+1)/2) : Evenl(a,(l+2)/2 ));
}
//--------------------------------------------------------------------
//
//  1/sqrt(<a|a>)
//
export double Norm(double a, int l) {return 1.0/sqrt(Integral(2*a,2*l));}

} //namespace
