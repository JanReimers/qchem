// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SphericalGaussian_m/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian_m/IEClient.H" 
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

namespace SphericalGaussian_m
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
    return GaussianIntegral(ea+eb,2*l);
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


IntegralEngine::SMat IntegralEngine::MakeOverlap(iec_t* ieab, const bf_tuple& c) const
{    
    auto ab=dcast(ieab);;
    size_t N=ab->size();
    int Nc,Lc,Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=GaussianIntegral(ab->es(i)+ab->es(j)+ec,ab->Ls(i)+ab->Ls(j)+Lc)*ab->ns(i)*ab->ns(j)*nc;
    return s;
}

IntegralEngine::ERI3 IntegralEngine::MakeOverlap3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);;
   
    ERI3 s3;
    for (auto i:c->es.indices()) s3.push_back(MakeOverlap(ieab,(*c)(i)));
    return s3;
}

//----------------------------------------------------------------------------------------
//
//  Repulsion type integrals
//
IntegralEngine::SMat IntegralEngine::MakeRepulsion(iec_t* iea ) const
{
    auto a=dcast(iea);;
    assert(a);
    size_t N=a->size();
    SMat r(N,N);
    for (auto i:r.rows())
        for (auto j:r.cols(i))
            r(i,j)=GaussianRepulsionIntegral(a->es(i),a->es(j),a->Ls(i),a->Ls(j))*a->ns(i)*a->ns(j);

    return r;
}

IntegralEngine::Mat IntegralEngine::MakeRepulsion(iec_t* iea,iec_t* ieb) const
{
    auto a=dcast(iea);;
    auto b=dcast(ieb);;
    size_t Na=a->es.size(), Nb=b->es.size();
    Mat s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=GaussianRepulsionIntegral(a->es(i),b->es(j),a->Ls(i),b->Ls(j))*a->ns(i)*b->ns(j);

    return s;
}

//
IntegralEngine::SMat IntegralEngine::MakeRepulsion(iec_t* ieab,const bf_tuple& c) const
{    
    auto ab=dcast(ieab);;
    size_t N=ab->size();
    int Nc,Lc,Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
        {
            GaussianRadialIntegrals R(ab->es(i)+ab->es(j),ec);
            s(i,j)=R.Coulomb(ab->Ls(i),ab->Ls(j),Lc,0)*ab->ns(i)*ab->ns(j)*nc;
        }
    return s;
}


IntegralEngine::ERI3 IntegralEngine::MakeRepulsion3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);;
    
    ERI3 s3;
    for (auto i:c->es.indices()) s3.push_back(MakeRepulsion(ieab,(*c)(i)));
    return s3;
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

////
//
////----------------------------------------------------------------------------------------
////
////  Special integrals
////


//IntegralEngine::SMat IntegralEngine::MakeKinetic(iec_t* iea) const
//{
//    auto a=dcast(iea);;
//    size_t N=a->size();
//    SMatrix<double> Hk(N);
//    for (auto i:Hk.rows())
//        for (auto j:Hk.cols(i))
//        {
//            assert(a->Ls(i)==a->Ls(j));
//            double t=a->es(i)+a->es(j);
//            int L=a->Ls(i),L1=L+1;
//            Hk(i,j)=0.5*a->ns(i)*a->ns(j)*
//                   (
//                       (L1*L1 + L*L1) * GaussianIntegral(t,2*L-2)
//                       -2*L1 * t      * GaussianIntegral(t,2*L  )
//                       +4*a->es(i)*a->es(j) * GaussianIntegral(t,2*L+2)
//                   );
//        }
//
//    return Hk;
//}
////

//IntegralEngine::SMat IntegralEngine::MakeNuclear(iec_t* iea,const Cluster& cl) const
//{
//    auto a=dcast(iea);;
//    size_t N=a->size(),L=a->Ls(1);
//    SMatrix<double> Hn(N);
//    double Z=-cl.GetNuclearCharge();
//    for (auto i:Hn.rows())
//        for (auto j:Hn.cols(i))
//            Hn(i,j)= Z*GaussianIntegral(a->es(i)+a->es(j),2*L-1)*a->ns(i)*a->ns(j);
//
//    return Hn;
//}

//IntegralEngine::RVec IntegralEngine::MakeNormalization(iec_t* iea) const
//{
//
//    auto a=dcast(iea);;
//    RVec n(a->size());
//    for (auto i:a->es.indices())  n(i)=GaussianNorm(a->es(i),a->Ls(i));
//    return n;
//}
double IntegralEngine::Charge (double ea,           size_t l) const
{
    return GaussianIntegral(ea,l);
}

//IntegralEngine::RVec IntegralEngine::MakeCharge(iec_t* iea) const
//{
//    auto a=dcast(iea);;
//    RVec c(a->size());
//    for (auto i:a->es.indices())  c(i)=GaussianIntegral(a->es(i),a->Ls(i))*a->ns(i);
//    return c;
//}

void IntegralEngine::Report(std::ostream& os) const
{
    os << "Spherical Gaussian integral engine cache:" << std::endl;
    os << "    No cache." << std::endl;
}


} //namespace
