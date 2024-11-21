// File: guass.cpp  General gaussian integral.


#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Integrals/Factorials.H"

#include <iostream>
#include <cassert>

double Oddl (double exp, int n);
double Evenl(double exp, int n);


//##############################################################################
//      /
//  4Pi |  r^l exp(-er^2) r^2 dr
//     /
double GaussianIntegral(double exp, int l)
{
//    return 4*Pi*( l%2 ? Oddl(exp,(l+1)/2) : Evenl(exp,l/2+1) );
    return 4*Pi*( l%2 ? Oddl(exp,(l+1)/2) : Evenl(exp,(l+2)/2 ));
}

double Oddl (double exp, int n)
{
    return n==0 ? 1.0/(2*exp) : Oddl(exp,n-1)*n/exp;
}

double Evenl(double exp, int n)
{
    return n==0 ? sqrt(Pi/exp)/2.0 : Evenl(exp,n-1)*(2*n-1)/(2.0*exp);
//    return n==0 ? sqrt(Pi/exp)/2.0 : Evenl(exp,n-1)*(2*n+1)/(2.0*exp);
}

double GaussianNorm(double e, int l)
{
    return 1.0/sqrt(GaussianIntegral(2*e,2*l));
}

inline void Swap(double& a, double& b, int& la, int& lb)
{
    double t=a;
    a=b;
    b=t;
    int    it=la;
    la=lb;
    lb=it;
}



