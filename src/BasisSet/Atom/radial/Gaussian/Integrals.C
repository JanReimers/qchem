// File: Gaussian::Integrals.C  General gaussian integral.


#include <iostream>
#include <cassert>
#include <cmath>
#include "radial/Gaussian/Integrals.H"
import qchem.BasisSet.Atom.AngularIntegrals;;


import Common.Constants;
import Common.Factorials;

namespace Gaussian
{

double Oddl (double exp, int n);
double Evenl(double exp, int n);

//##############################################################################
//      /
//  4Pi |  r^l exp(-a*r^2) r^2 dr
//     /
double Integral(double a, int l)
{
    assert(l>=-2);
    return FourPi*( l%2 ? Oddl(a,(l+1)/2) : Evenl(a,(l+2)/2 ));
}

double Oddl (double a, int n)
{
    return n==0 ? 1.0/(2*a) : Oddl(a,n-1)*n/a;
}

double Evenl(double a, int n)
{
    return n==0 ? sqrt(Pi/a)/2.0 : Evenl(a,n-1)*(2*n-1)/(2.0*a);
}

double Norm(double a, int l)
{
    return 1.0/sqrt(Integral(2*a,2*l));
}


};