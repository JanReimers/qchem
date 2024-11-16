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
