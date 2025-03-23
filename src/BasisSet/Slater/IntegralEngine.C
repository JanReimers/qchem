// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/Integrals/SlaterCD.H"
#include "Imp/Integrals/SlaterIntegrals.H"


namespace Slater
{

double IE_Primatives::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    SlaterCD cd(eab,ec,std::max(la,lc));
    return 4*4*Pi*Pi*cd.Coulomb_R0(la,lc);
}
// double Fit_IE::Repulsion(double eab, double ec,size_t la,size_t lc) const
// {    
//     SlaterCD cd(eab,ec,std::max(la,lc));
//     return 4*4*Pi*Pi*cd.Coulomb_R0(la,lc);
// }

double  IE_Primatives::Overlap(double ea , double eb,size_t l_total) const
{
    return SlaterIntegral(ea+eb,l_total); //Already has 4*Pi and r^2 from dr.
}

double IE_Primatives::Grad2(double ea , double eb,size_t la, size_t lb) const
{
    assert(la==lb);
    double ab=ea+eb;
    int l=la; //Safer to do formulas with int.
    int ll=l*(l+1);
    double Term1=((l+1)*(l+1)+ll)*SlaterIntegral(ab,2*l-2); //SlaterIntegral already has 4*Pi
    double Term2=-(l+1)*ab* SlaterIntegral(ab,2*l-1);
    double Term3=ea*eb*SlaterIntegral(ab,2*l);
    return 0.5*(Term1+Term2+Term3);
}

double IE_Primatives::Nuclear(double ea , double eb,size_t l_total) const
{
    return SlaterIntegral(ea+eb,l_total-1); //Already has 4*Pi
}

double Fit_IE::Charge(double ea, size_t l) const
{
    return SlaterIntegral(ea,l);
}


} //namespace
