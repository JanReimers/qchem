// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater/IEClient.H" 
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"

namespace Slater
{

double IntegralEngine::FourPi2=4*4*pi*pi;

double IntegralEngine::Overlap(double ea, double eb,size_t l) const
{
    return SlaterIntegral(ea+eb,l+2);
}

double IntegralEngine::Kinetic(double ea, double eb,size_t l) const
{
    double ab=ea+eb;
    int na=l+1,nb=l+1;
    int n=na+nb;
    double Term1=0.5*(na*nb+l*(l+1))*SlaterIntegral(ab,n-2);
    double Term2=-0.5*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
    double Term3=0.5*ea*eb*SlaterIntegral(ab,n);
    return Term1+Term2+Term3;
}

double IntegralEngine::Nuclear(double ea, double eb,size_t l) const
{
    return SlaterIntegral(ea+eb,2*l+1);
}

double IntegralEngine::Charge (double ea,           size_t l) const
{
    return SlaterIntegral(ea,l+2);
}

double IntegralEngine::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    SlaterCD cd(eab,ec,std::max(la,lc));
    return FourPi2*cd.Coulomb_R0(la,lc);
}

IntegralEngine::RVec IntegralEngine::Coulomb_AngularIntegrals(size_t la, size_t lc, int, int) const
{
    return AngularIntegrals::Coulomb(la,lc);
}

IntegralEngine::RVec IntegralEngine::ExchangeAngularIntegrals(size_t la, size_t lb, int, int) const
{
    return AngularIntegrals::Exchange(la,lb);
}



} //namespace
