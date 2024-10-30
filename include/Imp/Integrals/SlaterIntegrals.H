#ifndef SLATERINTEGRALS_H
#define SLATERINTEGRALS_H

#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Misc/IntPower.H"
#include "Imp/Integrals/Factorials.H"

#include <cmath>
#include <cassert>



//class GaussianRadialIntegrals
//{
//public:
//    GaussianRadialIntegrals(double e_ab, double e_cd);
//
//    double operator()   (int l,int la, int lb, int lc, int ld) const;
//    double DoExchangeSum(      int la, int lb, int lc, int ld) const; //Sum over l
//   
//private:
//    double Dab(int m, int n) const;
//    double Dcd(int m, int n) const;
//    double D(double _a, int m, int n) const;
//    double I(double _a, int m, int n) const;
//
//
//
//    double a,c; //Compact notation for exponents a+b and c+d.
//};


//--------------------------------------------------------------------
//      /
//  4Pi |  r^n exp(-a*r) dr
//     /
inline double SlaterIntegral(double a, int n)
{
    assert(n>=0);
    return 4*Pi* qchem::Fact[n]/pow(a,n+1);
}

//--------------------------------------------------------------------
//
//  1/sqrt(<a|a>)
//
inline double SlaterNorm(double a, int n)
{
    return 1.0/std::sqrt(SlaterIntegral(2*a,2*n));
}

class SlaterRadialIntegrals
{
public:
    SlaterRadialIntegrals(double _eab, double _ecd);
                   
//    double operator()   (int k,int la, int lb, int lc, int ld) const;
//    double operator()   (int lab, int lcd) const;
    double Coulomb      (int lab, int lcd) const;
    double Coulomb      (int la, int lb, int lc, int ld) const;
    double DoExchangeSum(int la, int lb, int lc, int ld) const; //Sum over k
    double Coulomb      (int la, int lb, int lc, int ld, int ma, int mb, int mc, int md) const; //Sum over k, Ylm version
    double DoExchangeSum(int la, int lb, int lc, int ld, int ma, int mb, int mc, int md) const; //Sum over k, Ylm version
    
private:
    double R(int k,int la, int lb, int lc, int ld) const;
    double D(double _a, int m, int n) const;
    double fk(double _a, int m, int n) const;
    
    double eab,ecd;
};


#endif // SLATERINTEGRALS_H
