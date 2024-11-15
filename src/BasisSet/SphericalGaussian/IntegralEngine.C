// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/IEClient.H" 
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/GaussianRadialIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include <Cluster.H>
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Imp/Containers/ERI4.H"

using std::cout;
using std::endl;

namespace SphericalGaussian
{

double IntegralEngine::FourPi2=4*4*Pi*Pi;

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
AnalyticIE<double>* IntegralEngine::Clone() const
{
    return new IntegralEngine(*this);
}

const IrrepIEClient* IntegralEngine::dcast(iec_t* iea)
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient*>(iea);
    assert(a);
    return a;
}


//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
double IntegralEngine::Overlap(double ea, double eb,size_t l) const
{
    return GaussianIntegral(ea+eb,l);
}

double IntegralEngine::Kinetic(double ea, double eb,size_t l) const
{
    double t=ea+eb;
    size_t l1=l+1;
    return 0.5*(
               (l1*l1 + l*l1) * GaussianIntegral(t,2*l-2)
               -2*l1 * t      * GaussianIntegral(t,2*l  )
               +4*ea*eb       * GaussianIntegral(t,2*l+2)
           );
}

double IntegralEngine::Nuclear(double ea, double eb,size_t l) const
{
    return GaussianIntegral(ea+eb,2*l-1);
}

double IntegralEngine::Charge (double ea,           size_t l) const
{
    return GaussianIntegral(ea,l);
}



//----------------------------------------------------------------------------------------
//
//  Repulsion type integrals
//


double IntegralEngine::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    GaussianRadialIntegrals R(eab,ec);
    return R.Coulomb(la,la,lc,0); 
}

//

void IntegralEngine::Make4C(ERI4& J, ERI4& K,const ::IEClient* iec) const
{
    const IEClient* sg=dynamic_cast<const IEClient*>(iec);

    for (index_t ia:sg->es.indices())
    {
        sg->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (index_t ic:sg->es.indices())
        {
            sg->loop_2(ic);
            int la=sg->Ls(ia), lc=sg->Ls(ic);
            for (index_t ib:sg->indices(la))
            {
                if (ib<ia) continue;
                sg->loop_3(ib);
                for (index_t id:sg->indices(lc))
                {
                    if (id<ic) continue;
                    const SphericalGaussianCD* cd1=sg->loop_4(id);
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    J(ia,ib,ic,id)=FourPi2*cd1->Coulomb_R0(la,lc)*norm;
                }
            }
        }
    }
    
    for (index_t ia:sg->es.indices())
    {
        sg->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (index_t ib:sg->es.indices(ia))
        {
            int la=sg->Ls(ia), lb=sg->Ls(ib);
            RVec Ak=AngularIntegrals::Exchange(la,lb);

            for (index_t ic:sg->indices(la))
            {
                if (ic<ia) continue;
                sg->loop_2(ic);
                sg->loop_3(ib);
                
                for (index_t id:sg->indices(lb))
                {
                    if (id<ic) continue;
                    const SphericalGaussianCD* cd1=sg->loop_4(id);
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    K(ia,ib,ic,id)=FourPi2*Ak*cd1->ExchangeRk(la,lb)*norm;
                }
            }
        }
    }
    
}

////
//
////----------------------------------------------------------------------------------------
////
////  Special integrals
////

void IntegralEngine::Report(std::ostream& os) const
{
    os << "Spherical Gaussian integral engine cache:" << std::endl;
    os << "    No cache." << std::endl;
}


} //namespace
