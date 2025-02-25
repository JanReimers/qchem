// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater/IEClient.H" 
#include "Imp/Integrals/SlaterCD.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"

using std::cout;
using std::endl;

namespace Slater
{

double IntegralEngine::Overlap(double ea, double eb,size_t l) const
{
    return SlaterIntegral(ea+eb,l+2); //Already has 4*Pi
}

double IntegralEngine::Kinetic(double ea, double eb,size_t l) const
{
    double ab=ea+eb;
    int na=l+1,nb=l+1;
    int n=na+nb;
    double Term1=0.5*(na*nb+l*(l+1))*SlaterIntegral(ab,n-2); //SlaterIntegral already has 4*Pi
    double Term2=-0.5*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
    double Term3=0.5*ea*eb*SlaterIntegral(ab,n);
    //cout << "Slater::IntegralEngine::Kinetic Terms 1,2,3=" << Term1 << " " << Term2 << " " << Term3 << endl;

    return Term1+Term2+Term3;
}
double IntegralEngine::Kinetic(double ea, double eb,size_t la, size_t lb) const
{
    double ab=ea+eb;
    int na=la+1,nb=lb+1;
    size_t ll=(la*(la+1)+lb*(lb+1))/2;
    int n=na+nb;
    double Term1=0.5*(na*nb+ll)*SlaterIntegral(ab,n-2); //SlaterIntegral already has 4*Pi
    double Term2=-0.5*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
    double Term3=0.5*ea*eb*SlaterIntegral(ab,n);
    //cout << "Slater::IntegralEngine::Kinetic Terms 1,2,3=" << Term1 << " " << Term2 << " " << Term3 << endl;

    return Term1+Term2+Term3;
}

double IntegralEngine::Nuclear(double ea, double eb,size_t l) const
{
    return SlaterIntegral(ea+eb,2*l+1); //Already has 4*Pi
}

double IntegralEngine::Charge (double ea,           size_t l) const
{
    return SlaterIntegral(ea,l+2);
}

double IntegralEngine::Repulsion(double eab, double ec,size_t la,size_t lc) const
{    
    SlaterCD cd(eab,ec,std::max(la,lc));
    return 4*4*pi*pi*cd.Coulomb_R0(la,lc);
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
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}


Vector<double> IntegralEngine::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double> IntegralEngine::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->ExchangeRk(la,lc);
}

double IntegralEngine1::Integral(qchem::IType type,double ea, double eb,size_t l) const
{
    switch(type)
    {
        case qchem::Overlap1: return SlaterIntegral(ea+eb,2*l+2); //Already has 4*Pi
        case qchem::Kinetic1:
        {
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
        case qchem::Nuclear1:  return SlaterIntegral(ea+eb,2*l+1); //Already has 4*Pi
        default: assert(false);
    }
    return 0.0;
}

} //namespace
