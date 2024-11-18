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

double IntegralEngine::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    GaussianRadialIntegrals R(eab,ec);
    return R.Coulomb(la,la,lc,0); 
}

IntegralEngine::RVec IntegralEngine::Coulomb_AngularIntegrals(size_t la, size_t lc, int, int) const
{
    return AngularIntegrals::Coulomb(la,lc);
}

IntegralEngine::RVec IntegralEngine::ExchangeAngularIntegrals(size_t la, size_t lb, int, int) const
{
    return AngularIntegrals::Exchange(la,lb);
}

const Cacheable* IntegralEngine::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}

Vector<double>  IntegralEngine::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double>  IntegralEngine::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
    return cd->ExchangeRk(la,lc);
}

} //namespace
