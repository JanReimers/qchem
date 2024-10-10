// File: _GaussianIntegrals_H_   Gaussian integral functions.
#ifndef  _GaussianIntegrals_H_
#define  _GaussianIntegrals_H_



#include <cmath>
#include "Misc/DFTDefines.H"

inline double TwoLPlusOne(int L)
{
    return 1.0;
//    return 1.0/(2*L+1);
}

//###################################################################
//
//  Hand coded Gaussian integrals.
//
inline double GaussianIntegralm2(double exp)
{
    return 2.0*Pi*sqrt(Pi/exp);
}

inline double GaussianIntegralm1(double exp)
{
    return 2.0*Pi/exp;
}

inline double GaussianIntegral0(double exp)
{
    double t=Pi/exp;
    return t*sqrt(t);
}
inline double GaussianIntegral2(double exp)
{
    double den=1.0/exp;
    double t=den*Pi;
    return 3.0/2.0*Pi*Square(den)*sqrt(t);
}

inline double GaussianIntegral4(double exp)
{
    double den=1.0/exp;
    double t=den*Pi;
    return 15.0/4.0*Pi*Cube(den)*sqrt(t);
}

inline double GaussianRepulsionIntegral(double e1, double e2)
{
    return 2.0*Pi52/(e1*e2)/sqrt(e1+e2);
}

//--------------------------------------------------------------------
//      /
//  4Pi |  r^l exp(-e*r^2) r^2 dr
//     /
double GaussianIntegral(double exp, int l);

//--------------------------------------------------------------------
//
//  1/sqrt(<a|a>)
//
double GaussianNorm(double e, int l);
//--------------------------------------------------------------------
//      /
//  4Pi |  r1^l1 exp(-e1*r1^2) * 1/r12 * r2^l2 exp(-e2*r2^2) * r^2 dr
//     /
double GaussianRepulsionIntegral(double e1, double e2, int l1, int l2);
double GaussianRepulsionIntegral(double ea, double eb, double ec, double ed, int la, int lb, int lc, int ld);
double GaussianExchangeIntegral (double ea, double eb, double ec, double ed, int la, int lb, int lc, int ld);

#endif //_GaussianIntegrals_H_
