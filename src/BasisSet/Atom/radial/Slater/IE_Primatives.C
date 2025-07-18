// File: Slater/IE_Primatives.C get all calculation of primative integrals in one place.

#include <cmath>
#include "radial/Slater/IE_Primatives.H"
#include "radial/Slater/Integrals.H"
#include "radial/Slater/Rk.H"

import Common.Constants;

namespace Slater
{
    double IE_Primatives::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    Slater::RkEngine cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

double  IE_Primatives::Overlap(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

double IE_Primatives::Grad2(double ea , double eb,size_t la, size_t lb) const
{
    assert(la==lb);
    double ab=ea+eb;
    int l=la; //Safer to do formulas with int.
    // int ll=l*(l+1);
    double Term1=(l+1)*(l+1)*Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
    double Term2=-(l+1)*ab* Slater::Integral(ab,2*l-1);
    double Term3=ea*eb*Slater::Integral(ab,2*l);
    return Term1+Term2+Term3;
}

double IE_Primatives::Inv_r1(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}
double IE_Primatives::Inv_r2(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total-2); //Already has 4*Pi
}

double IE_Primatives::Charge(double ea, size_t l) const
{
    return ::Slater::Integral(ea,l);
}


} //namespace
