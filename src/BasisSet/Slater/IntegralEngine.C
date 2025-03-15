// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/Integrals/SlaterCD.H"
#include "Imp/Integrals/SlaterIntegrals.H"

namespace Slater
{

double IE_Common::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    SlaterCD cd(eab,ec,std::max(la,lc));
    return 4*4*pi*pi*cd.Coulomb_R0(la,lc);
}

double  IE_Common::Overlap(double ea , double eb,size_t l_total) const
{
    return SlaterIntegral(ea+eb,l_total+2); //Already has 4*Pi
}

double IE_Common::Kinetic(double ea , double eb,size_t l, size_t lb) const
{
    assert(l==lb);
    double ab=ea+eb;
    int na=l+1,nb=l+1;
    size_t ll=l*(l+1);
    int n=na+nb;
    double Term1=0.5*(na*nb+ll)*SlaterIntegral(ab,n-2); //SlaterIntegral already has 4*Pi
    double Term2=-0.5*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
    double Term3=0.5*ea*eb*SlaterIntegral(ab,n);
    //cout << "Slater::IntegralEngine::Kinetic Terms 1,2,3=" << Term1 << " " << Term2 << " " << Term3 << endl;

    return Term1+Term2+Term3;
}

double IE_Common::Nuclear(double ea , double eb,size_t l_total) const
{
    return SlaterIntegral(ea+eb,l_total+1); //Already has 4*Pi
}

double Fit_IE::Charge(double ea, size_t l) const
{
    return SlaterIntegral(ea,l+2);
}


} //namespace
