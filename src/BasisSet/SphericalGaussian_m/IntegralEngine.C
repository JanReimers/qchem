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
                    J(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lc+1)*Akac*cd1->Coulomb_Rk(la,lc)*norm;
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
                    K(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lb+1)*Akab*cd1->ExchangeRk(la,lb)*norm;
//                    cout << "Kold(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                    cout << std::setprecision(8) << K(ia,ib,ic,id) << endl; 
                }
            }
        }
    }
    
}

//#define SymmetryCheck

ERIJ IntegralEngine::MakeDirect(const IrrepIEClient* a, const IrrepIEClient* c,const IEClient* iec) const
{
    size_t Na=a->size(), Nc=c->size();
    ERIJ J(Na,Nc);
    for (size_t ia:a->indices())
    {
        iec->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            iec->loop_2(ic);
            int la=a->Ls(ia), lc=c->Ls(ic);
            int ma=a->Ms(ia), mc=c->Ms(ic);
            RVec Akac=AngularIntegrals::Coulomb(la,lc,ma,mc);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                iec->loop_3(ib);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    const SphericalGaussianCD* cd=iec->loop_4(id);
                    assert(la==a->Ls(ib));
                    assert(lc==c->Ls(id));
                    if (J(ia,ib,ic,id)!=0.0)
                    {
                        cout << "overwritting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << J(ia,ib,ic,id) << endl;    
                        assert(false);
                    }
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    J(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lc+1)*Akac*cd->Coulomb_Rk(la,lc)*norm;
                }
            }
        }
    }

#ifdef SymmetryCheck
    double tol=1e-12;
    typedef std::tuple<int,int,int,int> i4_t;
    std::map<double,std::vector<i4_t> > Jsym;
    for (size_t ia:a->indices())
        for (size_t ib:a->indices())
            for (size_t ic:c->indices())
                for (size_t id:c->indices())
                {
                    double key= J(ia,ib,ic,id);
                    auto il=Jsym.lower_bound(key-tol);
                    auto iu=Jsym.upper_bound(key+tol);
                    if (il==Jsym.end() || il==iu)
                    {
                        std::vector<i4_t> indices;
                        indices.push_back(std::make_tuple(ia,ib,ic,id));
                        Jsym[key]=indices;
                    }
                    else
                        il->second.push_back(std::make_tuple(ia,ib,ic,id));
                }

    if (a==c) 
        cout << "J Irreps are equal" << endl;
    else
        cout << "J Irreps are NOT equal" << endl;

    for (auto k:Jsym)
    {
        cout << std::setprecision(14) << k.first;
        for (auto i:k.second)
            cout << " (" << std::get<0>(i) << "," << std::get<1>(i) << "," <<std::get<2>(i) << "," <<std::get<3>(i) << ") ";
        cout << std::endl;
    }
#endif
    
    
    return J;
}

ERIK IntegralEngine::MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b,const IEClient* iec) const
{
    size_t Na=a->size(), Nb=b->size();
    ERIK K(Na,Nb);
    for (size_t ia:a->indices())
    {
        iec->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (size_t ib:b->indices())
        {
            int la=a->Ls(ia), lb=b->Ls(ib);
            int ma=a->Ms(ia), mb=b->Ms(ib);
            RVec Akab=AngularIntegrals::Exchange(la,lb,ma,mb);

            for (size_t ic:a->indices())
            {
                if (ic<ia) continue;
                iec->loop_2(ic);
                iec->loop_3(ib);
                
                for (size_t id:b->indices())
                {
//                    if (id<ic) continue;
                    if (ia==ic && id<ib) continue;
                    //if (id<ib) continue;
                    const SphericalGaussianCD* cd=iec->loop_4(id);
                    assert(la==a->Ls(ic));
                    assert(lb==b->Ls(id));
                    if (K(ia,ic,ib,id)!=0.0)
                    {
                        cout << "overwritting Knew(" << ia << " " << ic << " " << ib << " " << id << ")="; 
                        cout << K(ia,ic,ib,id) << endl;    
                        assert(false);
                    }
                    double norm=a->ns(ia)*b->ns(ib)*a->ns(ic)*b->ns(id);
                    K(ia,ic,ib,id)=FourPi2*(2*la+1)*(2*lb+1)*Akab*cd->ExchangeRk(la,lb)*norm; 
                    if (ia==ic) K(ia,ic,id,ib)=K(ia,ic,ib,id); //ERIK container does support this symmetry yet.
//                    cout << "Knew(" << ia << " " << ic << " " << ib << " " << id << ")="; 
//                    cout << std::setprecision(8) << K(ia,ic,ib,id) << endl;    

                }
            }
        }
    }

    #ifdef SymmetryCheck    
    double tol=1e-12;
    typedef std::tuple<int,int,int,int> i4_t;
    std::map<double,std::vector<i4_t> > Ksym;
    for (size_t ia:a->indices())
        for (size_t ib:b->indices())
            for (size_t ic:a->indices())
                for (size_t id:b->indices())
                {
                    double key= K(ia,ic,ib,id);
                    auto il=Ksym.lower_bound(key-tol);
                    auto iu=Ksym.upper_bound(key+tol);
                    if (il==Ksym.end() || il==iu)
                    {
                        std::vector<i4_t> indices;
                        indices.push_back(std::make_tuple(ia,ic,ib,id));
                        Ksym[key]=indices;
                    }
                    else
                        il->second.push_back(std::make_tuple(ia,ic,ib,id));
                }

    if (a==b) 
        cout << "Irreps are equal" << endl;
    else
        cout << "Irreps are NOT equal" << endl;

    for (auto k:Ksym)
    {
        cout << std::setprecision(14) << k.first;
        for (auto i:k.second)
            cout << " (" << std::get<0>(i) << "," << std::get<1>(i) << "," <<std::get<2>(i) << "," <<std::get<3>(i) << ") ";
        cout << std::endl;
    }
#endif

    return K;
}

void IntegralEngine::MakeDirect(erij_t& Jac, const ::IEClient* iec) const
{
    Jac.clear();
    const IEClient& sg=*dynamic_cast<const IEClient*>(iec);
    size_t NIrrep=sg.GetNumIrreps();
    for (size_t ia=1;ia<=NIrrep;ia++)
        for (size_t ic=1;ic<=NIrrep;ic++) //TODO run from ia n
        {
            const IrrepIEClient* a=sg[ia];
            const IrrepIEClient* c=sg[ic];
            Jac[ia][ic]=MakeDirect(a,c,&sg);
        }

}

void IntegralEngine::MakeExchange(erik_t& Kab, const ::IEClient* iec) const
{
    Kab.clear();
    const IEClient& sg=*dynamic_cast<const IEClient*>(iec);
    size_t NIrrep=sg.GetNumIrreps();
    for (size_t ia=1;ia<=NIrrep;ia++)
        for (size_t ib=1;ib<=NIrrep;ib++) //TODO run from ib 
        {
            const IrrepIEClient* a=sg[ia];
            const IrrepIEClient* b=sg[ib];
            Kab[ia][ib]=MakeExchange(a,b,&sg);
        }
    
}



} //namespace
