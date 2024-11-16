// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SphericalGaussian_m/IntegralEngine.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/GaussianRadialIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/Containers/ERI4.H"

using std::cout;
using std::endl;

namespace SphericalGaussian_m
{
    
double IntegralEngine::FourPi2=4*4*Pi*Pi;

IntegralEngine::RVec IntegralEngine::Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const
{
    return AngularIntegrals::Coulomb(la,lc,ma,mc);
}

IntegralEngine::RVec IntegralEngine::ExchangeAngularIntegrals(size_t la, size_t lb, int ma, int mb) const
{
    return AngularIntegrals::Exchange(la,lb,ma,mb);
}

void IntegralEngine::Make4C(ERI4& J, ERI4& K,const ::IEClient* iec) const
{
    const IEClient* sg=dynamic_cast<const IEClient*>(iec);

    for (index_t ia:sg->es.indices())
    {
        sg->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (index_t ic:sg->es.indices(ia))
        {
            sg->loop_2(ic);
            int la=sg->Ls(ia), lc=sg->Ls(ic);
            int ma=sg->Ms(ia), mc=sg->Ms(ic);
            RVec Akac=AngularIntegrals::Coulomb(la,lc,ma,mc);
            for (index_t ib:sg->indices(la))
            {
                if (ib<ia) continue;
                sg->loop_3(ib);
                for (index_t id:sg->indices(lc))
                {
                    if (id<ic) continue;
                    if (ia==ic && id<ib) continue;
                    if (J(ia,ib,ic,id)!=0.0)
                    {
                        cout << "overwritting Jold(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << J(ia,ib,ic,id) << endl;    
                        assert(false);
                    }
                    const SphericalGaussianCD* cd1=sg->loop_4(id);
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    J(ia,ib,ic,id)=FourPi2*Akac*cd1->Coulomb_Rk(la,lc)*norm;
//                    cout << "J(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                    cout << std::setprecision(8) << J(ia,ib,ic,id) << endl; 
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
            int ma=sg->Ms(ia), mb=sg->Ms(ib);
            RVec Akab=AngularIntegrals::Exchange(la,lb,ma,mb);

            for (index_t ic:sg->indices(la))
            {
                if (ic<ia) continue;
                sg->loop_2(ic);
                sg->loop_3(ib);
                
                for (index_t id:sg->indices(lb))
                {
                    if (id<ic) continue;
                    if (ia==ic && id<ib) continue;
                    if (K(ia,ib,ic,id)!=0.0)
                    {
                        cout << "overwritting Kold(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << std::setprecision(8) << K(ia,ib,ic,id) << endl;    
                        assert(false);
                    }
                    const SphericalGaussianCD* cd1=sg->loop_4(id);
                    double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                    K(ia,ib,ic,id)=FourPi2*Akab*cd1->ExchangeRk(la,lb)*norm;
//                    cout << "Kold(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                    cout << std::setprecision(8) << K(ia,ib,ic,id) << endl; 
                }
            }
        }
    }
    
}

//#define SymmetryCheck

//
//#ifdef SymmetryCheck
//    double tol=1e-12;
//    typedef std::tuple<int,int,int,int> i4_t;
//    std::map<double,std::vector<i4_t> > Jsym;
//    for (size_t ia:a->indices())
//        for (size_t ib:a->indices())
//            for (size_t ic:c->indices())
//                for (size_t id:c->indices())
//                {
//                    double key= J(ia,ib,ic,id);
//                    auto il=Jsym.lower_bound(key-tol);
//                    auto iu=Jsym.upper_bound(key+tol);
//                    if (il==Jsym.end() || il==iu)
//                    {
//                        std::vector<i4_t> indices;
//                        indices.push_back(std::make_tuple(ia,ib,ic,id));
//                        Jsym[key]=indices;
//                    }
//                    else
//                        il->second.push_back(std::make_tuple(ia,ib,ic,id));
//                }
//
//    if (a==c) 
//        cout << "J Irreps are equal" << endl;
//    else
//        cout << "J Irreps are NOT equal" << endl;
//
//    for (auto k:Jsym)
//    {
//        cout << std::setprecision(14) << k.first;
//        for (auto i:k.second)
//            cout << " (" << std::get<0>(i) << "," << std::get<1>(i) << "," <<std::get<2>(i) << "," <<std::get<3>(i) << ") ";
//        cout << std::endl;
//    }
//#endif
//    
//    
//    return J;
//}

//
//    #ifdef SymmetryCheck    
//    double tol=1e-12;
//    typedef std::tuple<int,int,int,int> i4_t;
//    std::map<double,std::vector<i4_t> > Ksym;
//    for (size_t ia:a->indices())
//        for (size_t ib:b->indices())
//            for (size_t ic:a->indices())
//                for (size_t id:b->indices())
//                {
//                    double key= K(ia,ic,ib,id);
//                    auto il=Ksym.lower_bound(key-tol);
//                    auto iu=Ksym.upper_bound(key+tol);
//                    if (il==Ksym.end() || il==iu)
//                    {
//                        std::vector<i4_t> indices;
//                        indices.push_back(std::make_tuple(ia,ic,ib,id));
//                        Ksym[key]=indices;
//                    }
//                    else
//                        il->second.push_back(std::make_tuple(ia,ic,ib,id));
//                }
//
//    if (a==b) 
//        cout << "Irreps are equal" << endl;
//    else
//        cout << "Irreps are NOT equal" << endl;
//
//    for (auto k:Ksym)
//    {
//        cout << std::setprecision(14) << k.first;
//        for (auto i:k.second)
//            cout << " (" << std::get<0>(i) << "," << std::get<1>(i) << "," <<std::get<2>(i) << "," <<std::get<3>(i) << ") ";
//        cout << std::endl;
//    }
//#endif
//
//    return K;
//}
//



} //namespace
