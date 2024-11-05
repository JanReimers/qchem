#ifndef SLATERINTEGRALS_H
#define SLATERINTEGRALS_H

//#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Misc/IntPower.H"
#include "Imp/Integrals/Factorials.H"

#include <cmath>
#include <cassert>

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

    double R0(int lab, int lcd) const {return R(0,lab,lab,lcd,lcd);}
    double R0(int la, int lb, int lc, int ld) const {return R(0,la,lb,lc,ld);}

        
private:
    double R(int k,int la, int lb, int lc, int ld) const;
    double D(double _a, int m, int n) const;
    double fk(double _a, int m, int n) const;
    
    double eab,ecd;
    static double FourPi2;
};

template <class T> class Vector;

Vector<double> Rk(int la,int lc,double eab,double ecd);


#endif // SLATERINTEGRALS_H
