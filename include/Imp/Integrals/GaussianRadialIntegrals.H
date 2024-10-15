#ifndef _GaussianRadialIntegrals_H_
#define _GaussianRadialIntegrals_H_

#include "Imp/Integrals/Wigner3j.H"
#include "Imp/Integrals/Factorials.H"
#include "Misc/IntPower.H"
#include <cmath>



class GaussianRadialIntegrals
{
public:
    GaussianRadialIntegrals(double e_ab, double e_cd);

    double operator()   (int k,int la, int lb, int lc, int ld) const;
    double DoExchangeSum(      int la, int lb, int lc, int ld) const; //Sum over l
private:
    double Dab(int m, int n) const;
    double Dcd(int m, int n) const;
    double D(double _a, int m, int n) const;
    double I(double _a, int m, int n) const;

    double a,c; //Compact notation for exponents a+b and c+d.
};

/*
//--------------------------------------------------------------------
//      /
//  4Pi |  r^n exp(-a*r) r^2 dr
//     /
inline double SlaterIntegral(double a, int n)
{
    return 4*Pi* qchem::Fact[n+2]/intpow(a,n+1);
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
    SlaterRadialIntegrals(double _ea, double _eb, double _ec, double _ed
                     , int _na, int _nb , int _nc, int _nd
                     , int _la, int _lb , int _lc, int _ld)
    : ea(_ea), eb(_eb), ec(_ec), ed(_ed)
    , na(_na), nb(_nb), nc(_nc), nd(_nd)
    , la(_la), lb(_lb), lc(_lc), ld(_ld)
    , nac(na+nc), nbd(nb+nd), L(la+lb+lc+ld)
    , eac(ea+ec), ebd(eb+ed), eabcd(ea+eb+ec+ed)
    {};
    
    double Coulomb() const;
    double Exchange() const;
private:
    double R(int k) const;
    double Habcd(int k) const;
    double Hbadc(int k) const;
    
    double ea,eb,ec,ed;
    int na,nb,nc,nd;
    int la,lb,lc,ld;
    int nac, nbd, L;
    double eac,ebd,eabcd;
    Wigner3j itsW3j;
};
*/

#endif // _GaussianRadialIntegrals_H_
