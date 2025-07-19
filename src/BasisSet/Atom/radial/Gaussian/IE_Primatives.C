// File: Gaussian/IE_Primatives.C get all calculation of primative integrals in one place.

#include <cassert>
#include "radial/Gaussian/IE_Primatives.H"
#include "radial/Gaussian/Integrals.H"
#include "radial/Gaussian/Rk.H"

import Common.Constants;

namespace Gaussian
{

double  IE_Primatives::Overlap(double ea , double eb,size_t l_total) const
{
    return Gaussian::Integral(ea+eb,l_total); //Already has 4*Pi
}
    
double IE_Primatives::Grad2(double ea , double eb,size_t l, size_t lb) const
{
    assert(l==lb);
    double t=ea+eb;
    size_t l1=l+1;
    return  l1*l1     * Gaussian::Integral(t,2*l-2)
            -2*l1 * t * Gaussian::Integral(t,2*l  )
            +4*ea*eb  * Gaussian::Integral(t,2*l+2);
        
}

double IE_Primatives::Inv_r1(double ea , double eb,size_t l_total) const
{
    return Gaussian::Integral(ea+eb,l_total-1); //Already has 4*Pi
}
double IE_Primatives::Inv_r2(double ea , double eb,size_t l_total) const
{
    return Gaussian::Integral(ea+eb,l_total-2); //Already has 4*Pi
}
double IE_Primatives::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    Gaussian::RkEngine cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

double IE_Primatives::Charge(double ea, size_t l) const
{
    return Gaussian::Integral(ea,l);
}

} //namespace
