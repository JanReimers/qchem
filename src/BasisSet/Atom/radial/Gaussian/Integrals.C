// File: guass.cpp  General gaussian integral.


#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Factorials.H"

#include <iostream>
#include <cassert>

double Oddl (double exp, int n);
double Evenl(double exp, int n);

//##############################################################################
//      /
//  4Pi |  r^l exp(-a*r^2) r^2 dr
//     /
double GaussianIntegral(double a, int l)
{
    return 4*Pi*( l%2 ? Oddl(a,(l+1)/2) : Evenl(a,(l+2)/2 ));
}

double Oddl (double a, int n)
{
    return n==0 ? 1.0/(2*a) : Oddl(a,n-1)*n/a;
}

double Evenl(double a, int n)
{
    return n==0 ? sqrt(Pi/a)/2.0 : Evenl(a,n-1)*(2*n-1)/(2.0*a);
}

double GaussianNorm(double a, int l)
{
    return 1.0/sqrt(GaussianIntegral(2*a,2*l));
}


