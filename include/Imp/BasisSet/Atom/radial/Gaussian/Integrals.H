// File: _Gaussian::Integrals.H   Gaussian integral functions.
#ifndef  _GaussianIntegrals_H_
#define  _GaussianIntegrals_H_

namespace Gaussian
{
//###################################################################
//
//  Hand coded Gaussian integrals.
//
// inline double GaussianIntegralm2(double exp)
// {
//     return 2.0*Pi*sqrt(Pi/exp);
// }

// inline double GaussianIntegralm1(double exp)
// {
//     return 2.0*Pi/exp;
// }

// inline double GaussianIntegral0(double exp)
// {
//     double t=Pi/exp;
//     return t*sqrt(t);
// }
// inline double GaussianIntegral2(double exp)
// {
//     double den=1.0/exp;
//     double t=den*Pi;
//     return 3.0/2.0*Pi*Square(den)*sqrt(t);
// }

// inline double GaussianIntegral4(double exp)
// {
//     double den=1.0/exp;
//     double t=den*Pi;
//     return 15.0/4.0*Pi*Cube(den)*sqrt(t);
// }

// inline double GaussianRepulsionIntegral(double e1, double e2)
// {
//     return 2.0*Pi52/(e1*e2)/sqrt(e1+e2);
// }

//--------------------------------------------------------------------
//      /
//  4Pi |  r^l exp(-e*r^2) r^2 dr
//     /
double Integral(double exp, int l);

//--------------------------------------------------------------------
//
//  1/sqrt(<a|a>)
//
double Norm(double e, int l);

} //namespace
#endif //_GaussianIntegrals_H_
