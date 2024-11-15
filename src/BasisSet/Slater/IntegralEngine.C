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


void IntegralEngine::Make4C(ERI4& J, ERI4& K,const ::IEClient* iec) const
{
    const IEClient* sg=dynamic_cast<const IEClient*>(iec);

    for (size_t ia:sg->es.indices())
    {
        sg->loop_1(ia); //Start a cache for SlaterCD*
        for (size_t ic:sg->es.indices())
        {
            sg->loop_2(ic);
            int la=sg->Ls(ia), lc=sg->Ls(ic);
            //No angular integrals required for the Coulomb part.
            for (size_t ib:sg->indices(la)) //Only loop over indices with l=la
            {
                if (ib<ia) continue;
                sg->loop_3(ib);
                for (size_t id:sg->indices(lc)) //Only loop over indices with l=lc
                {
                    if (id<ic) continue;
                    const SlaterCD* cd1=sg->loop_4(id);
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    J(ia,ib,ic,id)=FourPi2*cd1->Coulomb_R0(la,lc)*norm;
                 }
            }
        }
    }

                
    for (size_t ia:sg->es.indices())
    {
        sg->loop_1(ia);
        for (size_t ib:sg->es.indices(ia))
        {
            int la=sg->Ls(ia), lb=sg->Ls(ib);
            RVec Akab=AngularIntegrals::Exchange(la,lb);
            for (size_t ic:sg->indices(la))
            {
                if (ic<ia) continue;
                sg->loop_2(ic);
                sg->loop_3(ib);
                for (size_t id:sg->indices(lb))
                {
                    if (id<ic) continue;
                    const SlaterCD* cd1=sg->loop_4(id);
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    K(ia,ib,ic,id)=FourPi2*Akab*cd1->ExchangeRk(la,lb)*norm;
                }
            }
        }
    }
    
}


} //namespace
