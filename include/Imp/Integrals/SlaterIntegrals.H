#ifndef SLATERINTEGRALS_H
#define SLATERINTEGRALS_H

#include "Imp/Misc/IntPower.H"
#include "Imp/Integrals/Factorials.H"

#include <cmath>

//--------------------------------------------------------------------
//      /
//  4Pi |  r^n exp(-a*r) dr
//     /
inline double SlaterIntegral(double a, int n)
{
    if (n==-1) 
    {
        assert(a>0.0);
        return -4*Pi*log(a);
    }
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


#endif // SLATERINTEGRALS_H
