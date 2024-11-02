// File: SlaterIE.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/IEClient.H" 
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/Factorials.H"
#include <Cluster.H>
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Imp/Containers/ERI4.H"

namespace Slater_m
{

double IntegralEngine::FourPi2=4*4*Pi*Pi;

IntegralEngine::IntegralEngine()
{
}

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

using std::cout;
using std::endl;
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
// <a|b>
IntegralEngine::SMat IntegralEngine::MakeOverlap(iec_t* iea ) const
{
    auto  a=dcast(iea);
    size_t N=a->size();
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=SlaterIntegral(a->es(i)+a->es(j),a->Ns(i)+a->Ns(j))*a->ns(i)*a->ns(j);

    return s;
}

//
//
// <ab|c>
IntegralEngine::SMat IntegralEngine::MakeOverlap(iec_t* ieab, const bf_tuple& c) const
{    
    auto ab=dcast(ieab);
    size_t N=ab->size();
    int Nc,Lc,Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=SlaterIntegral(ab->es(i)+ab->es(j)+ec,ab->Ns(i)+ab->Ns(j)+Lc)*ab->ns(i)*ab->ns(j)*nc;
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
        {
            SlaterRadialIntegrals R(a->es(i),a->es(j));
            r(i,j)=FourPi2*R.R0(a->Ls(i),a->Ls(j))*a->ns(i)*a->ns(j);
        }

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
        {
            SlaterRadialIntegrals R(a->es(i),b->es(j));
            s(i,j)=FourPi2*R.R0(a->Ls(i),b->Ls(j))*a->ns(i)*b->ns(j);
        }

    return s;
}

//
IntegralEngine::SMat IntegralEngine::MakeRepulsion(iec_t* ieab,const bf_tuple& c) const
{    
    auto ab=dcast(ieab);;
    size_t N=ab->size();
    int Nc,Lc, Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
        {
            SlaterRadialIntegrals R(ab->es(i)+ab->es(j),ec);
            s(i,j)=FourPi2*R.R0(ab->Ls(i),ab->Ls(j),Lc,0)*ab->ns(i)*ab->ns(j)*nc;
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

void IntegralEngine::Make4C(ERI4& J,const ::IEClient* iec) const
{
    Make4C(&J,NULL,iec);
}
void IntegralEngine::Make4C(ERI4& J, ERI4& K,const ::IEClient* iec) const
{
    Make4C(&J,&K,iec);
}

bool Jmatch(const IEClient& iec, index_t ia, index_t ib, index_t ic, index_t id)
{
    bool l=iec.Ls(ia)==iec.Ls(ib) && iec.Ls(ic)==iec.Ls(id); 
    //bool m=iec.Ms(ia)==iec.Ms(ib) && iec.Ms(ic)==iec.Ms(id); 
    return l;// && m;
}

bool Kmatch(const IEClient& iec, index_t ia, index_t ib, index_t ic, index_t id)
{
    bool lac=iec.Ls(ia)==iec.Ls(ic) && iec.Ls(ib)==iec.Ls(id); 
    bool mac=iec.Ms(ia)==iec.Ms(ic) && iec.Ms(ib)==iec.Ms(id); 
    bool lad=iec.Ls(ia)==iec.Ls(id) && iec.Ls(ib)==iec.Ls(ic); 
    bool mad=iec.Ms(ia)==iec.Ms(id) && iec.Ms(ib)==iec.Ms(ic); 
    return (lac && mac) || (lad && mad);
}


void IntegralEngine::Make4C(ERI4* J, ERI4* K,const ::IEClient* iec) const
{
    const IEClient* sg=dynamic_cast<const IEClient*>(iec);
    size_t N=sg->size();
    J->SetSize(N,0.0);
    if (K)
        K->SetSize(N,0.0);
    
    for (index_t ia:sg->es.indices())
        for (index_t ib:sg->es.indices())
            for (index_t ic:sg->es.indices())
                for (index_t id:sg->es.indices())
                {
                    bool doJ = Jmatch(*sg,ia,ib,ic,id) ;
                    bool doK = K && Kmatch(*sg,ia,ib,ic,id) ;
                    if (doJ || doK)
                    {
                        double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                        SlaterRadialIntegrals S(sg->es(ia)+sg->es(ib),sg->es(ic)+sg->es(id));
                        if (doJ)
                        {
                            (*J)(ia,ib,ic,id)=S.Coulomb(sg->Ls(ia),sg->Ls(ib),sg->Ls(ic),sg->Ls(id),sg->Ms(ia),sg->Ms(ib),sg->Ms(ic),sg->Ms(id))*norm;
                        //                           std::cout << "L=(" << sg->Ls(ia) << "," << sg->Ls(ib) << "," << sg->Ls(ic) << "," << sg->Ls(id) 
                        //                            << ") abcd=(" << ia << "," << ib << "," << ic << "," << id << ")  J/norm=" << J(ia,ib,ic,id)/norm << std::endl;

                        }
                        if (doK)
                        {
//                           std::cout << "L=(" << sg->Ls(ia) << "," << sg->Ls(ib) << "," << sg->Ls(ic) << "," << sg->Ls(id) 
//                            << ") m=(" << sg->Ms(ia) << "," << sg->Ms(ib) << "," << sg->Ms(ic) << "," << sg->Ms(id) 
//                            << ") abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                            (*K)(ia,ib,ic,id)=S.DoExchangeSum(sg->Ls(ia),sg->Ls(ib),sg->Ls(ic),sg->Ls(id),sg->Ms(ia),sg->Ms(ib),sg->Ms(ic),sg->Ms(id))*norm;
//                           std::cout << "L=(" << sg->Ls(ia) << "," << sg->Ls(ib) << "," << sg->Ls(ic) << "," << sg->Ls(id) 
//                            << ") abcd=(" << ia << "," << ib << "," << ic << "," << id << ")  K/norm=" << K(ia,ib,ic,id)/norm << std::endl;
                            
                        }
//                        else
//                            if (K) (*K)(ia,ib,ic,id)=0.0;
                     }
                }
    
}

////
//
////----------------------------------------------------------------------------------------
////
////  Special integrals
////
IntegralEngine::SMat IntegralEngine::MakeKinetic(iec_t* iea) const
{
    auto a=dcast(iea);;
    size_t N=a->size();
    SMatrix<double> Hk(N);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols(i))
        {
            assert(a->Ls(i)==a->Ls(j));
            double ea=a->es(i), eb=a->es(j);
            double ab=ea+eb;
            int la=a->Ls(i),lb=a->Ls(j);
            int na=la+1,nb=lb+1;
            int n=a->Ns(i)+a->Ns(j);
            int l=a->Ls(i);
            assert(la==lb);
            double Term1=0.5*a->ns(i)*a->ns(j)*(na*nb+l*(l+1))*SlaterIntegral(ab,n-2);
            double Term2=-0.5*a->ns(i)*a->ns(j)*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
            double Term3=0.5*a->ns(i)*a->ns(j)*ea*eb*SlaterIntegral(ab,n);
            Hk(i,j)=Term1+Term2+Term3;
            
        }

    return Hk;
}
//
IntegralEngine::SMat IntegralEngine::MakeNuclear(iec_t* iea,const Cluster& cl) const
{
    auto a=dcast(iea);;
    size_t N=a->size();
    SMatrix<double> Hn(N);
    double Z=-cl.GetNuclearCharge();
    for (auto i:Hn.rows())
        for (auto j:Hn.cols(i))
            Hn(i,j)= Z*SlaterIntegral(a->es(i)+a->es(j),a->Ns(i)+a->Ns(j)-1)*a->ns(i)*a->ns(j);

    return Hn;
}

IntegralEngine::RVec IntegralEngine::MakeNormalization(iec_t* iea) const
{

    auto a=dcast(iea);;
    RVec n(a->size());
    for (auto i:a->es.indices())  n(i)=SlaterNorm(a->es(i),a->Ls(i));
    return n;
}

IntegralEngine::RVec IntegralEngine::MakeCharge(iec_t* iea) const
{
    auto a=dcast(iea);;
    RVec c(a->size());
    for (auto i:a->es.indices())  c(i)=SlaterIntegral(a->es(i),a->Ns(i)+1)*a->ns(i);
    return c;
}

void IntegralEngine::Report(std::ostream& os) const
{
    os << "Spherical Gaussian integral engine cache:" << std::endl;
    os << "    No cache." << std::endl;
}


} //namespace
