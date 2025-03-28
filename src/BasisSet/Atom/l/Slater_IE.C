// File: Atom/l/Slater_IE.H Integral Engine for Slater atom basis functions.


#include "Imp/BasisSet/Atom/l/Slater_IE.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"
#include "Imp/BasisSet/Atom/radial/Slater/Rk.H"


namespace Slater
{

double IE_Primatives::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    Slater::RkEngine cd(eab,ec,std::max(la,lc));
    return 4*4*Pi*Pi*cd.Coulomb_R0(la,lc);
}
// double Fit_IE::Repulsion(double eab, double ec,size_t la,size_t lc) const
// {    
//     Slater::RkEngine cd(eab,ec,std::max(la,lc));
//     return 4*4*Pi*Pi*cd.Coulomb_R0(la,lc);
// }

double  IE_Primatives::Overlap(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

double IE_Primatives::Grad2(double ea , double eb,size_t la, size_t lb) const
{
    assert(la==lb);
    double ab=ea+eb;
    int l=la; //Safer to do formulas with int.
    int ll=l*(l+1);
    double Term1=((l+1)*(l+1)+ll)*Slater::Integral(ab,2*l-2); //SlaterIntegral already has 4*Pi
    double Term2=-(l+1)*ab* Slater::Integral(ab,2*l-1);
    double Term3=ea*eb*Slater::Integral(ab,2*l);
    return Term1+Term2+Term3;
}

double IE_Primatives::Nuclear(double ea , double eb,size_t l_total) const
{
    return Slater::Integral(ea+eb,l_total-1); //Already has 4*Pi
}

double IE_Primatives::Charge(double ea, size_t l) const
{
    return ::Slater::Integral(ea,l);
}

} //namespace

double Atoml::Slater::Fit_IE::Charge(double ea, size_t l) const
{
    return ::Slater::Integral(ea,l);
}
