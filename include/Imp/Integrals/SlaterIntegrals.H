#ifndef SLATERINTEGRALS_H
#define SLATERINTEGRALS_H

#include "Misc/IntPower.H"
#include <cmath>



class SlaterIntegrals
{
public:
    SlaterIntegrals(double e_ab, double e_cd);

    double operator()   (int l,int la, int lb, int lc, int ld) const;
    double DoExchangeSum(      int la, int lb, int lc, int ld) const; //Sum over l
    static double DFact[]; //Double factorials 1,3,3*5,3*5*7 etc. lookup table.
    static double Fact[]; //Factorials 1,2,2*3,2*3*4 etc. lookup table.
    static double Twon[];  //2^n lookup table.
private:
    double Dab(int m, int n) const;
    double Dcd(int m, int n) const;
    double D(double _a, int m, int n) const;
    double I(double _a, int m, int n) const;



    double a,c; //Compact notation for exponents a+b and c+d.
};

//--------------------------------------------------------------------
//      /
//  4Pi |  r^n exp(-a*r) r^2 dr
//     /
inline double SlaterIntegral(double a, int n)
{
//    return 4*Pi*( l%2 ? Oddl(exp,(l+1)/2) : Evenl(exp,l/2+1) );
    return 4*Pi* SlaterIntegrals::Fact[n+2]/intpow(a,n+1);
}

//--------------------------------------------------------------------
//
//  1/sqrt(<a|a>)
//
inline double SlaterNorm(double a, int n)
{
    return 1.0/std::sqrt(SlaterIntegral(2*a,2*n));
}



#endif // SLATERINTEGRALS_H
