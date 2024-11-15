// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/IEClient.H" 
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"
#include <iomanip>

namespace Slater_m
{

double IntegralEngine::FourPi2=4*4*pi*pi;

using std::cout;
using std::endl;

void IntegralEngine::Make4C(ERI4& J, ERI4& K,const ::IEClient* iec) const
{
    const IEClient* sg=dynamic_cast<const IEClient*>(iec);

    for (size_t ia:sg->es.indices())
    {
        sg->loop_1(ia); //Start a cache for SlaterCD*
        for (size_t ic:sg->es.indices(ia))
        {
            sg->loop_2(ic);
            int la=sg->Ls(ia), lc=sg->Ls(ic);
            int ma=sg->Ms(ia), mc=sg->Ms(ic);
            RVec Akac=AngularIntegrals::Coulomb(la,lc,ma,mc);
            //cout << std::setprecision(6) << "Akac=" << Akac << endl;
            for (size_t ib:sg->indices(la))
            {
                if (ib<ia) continue;
                sg->loop_3(ib);
                for (size_t id:sg->indices(lc))
                {
                    if (id<ic) continue;
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    const SlaterCD* cd1=sg->loop_4(id);
                    //cout << "cd.Rk=" << cd.Rk(la,lc) << endl;
                    J(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lc+1)*Akac*cd1->Coulomb_Rk(la,lc)*norm;
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
            int ma=sg->Ms(ia), mb=sg->Ms(ib);
            RVec Akab=AngularIntegrals::Exchange(la,lb,ma,mb);
            //cout << std::setprecision(6) << "Akab=" << Akab << endl;
            for (size_t ic:sg->indices(la))
            {
                if (ic<ia) continue;
                sg->loop_2(ic);
                sg->loop_3(ib);

                for (size_t id:sg->indices(lb))
                {
                    if (id<ic) continue;
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    const SlaterCD* cd1=sg->loop_4(id);
                    //cout << "cd.ExchangeRk=" << cd.ExchangeRk(la,lb) << endl;
                    K(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lb+1)*Akab*cd1->ExchangeRk(la,lb)*norm;                        
                }
           }
        }
    }
         
    
}



} //namespace
