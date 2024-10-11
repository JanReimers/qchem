// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/IEClient.H" 
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include <Cluster.H>
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Imp/Containers/ERI4.H"

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
IntegralEngine::SMat IntegralEngine::MakeOverlap(iec_t* iea ) const
{
    auto  a=dcast(iea);
    size_t N=a->size();
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=GaussianIntegral(a->es(i)+a->es(j),2*a->Ls(i))*a->ns(i)*a->ns(j);

    return s;
}


IntegralEngine::SMat IntegralEngine::MakeOverlap(iec_t* ieab, const bf_tuple& c) const
{    
    auto ab=dcast(ieab);;
    size_t N=ab->size();
    int Lc;
    double ec,nc;
    std::tie(Lc,ec,nc)=c;
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
    int Lc;
    double ec,nc;
    std::tie(Lc,ec,nc)=c;
    SMat s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
        {
            SlaterIntegrals R(ab->es(i)+ab->es(j),ec);
            s(i,j)=FourPi2*R(0,ab->Ls(i),ab->Ls(j),Lc,0)*ab->ns(i)*ab->ns(j)*nc;
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
    size_t N=sg->size();
    J.SetSize(N,-1.0);
    K.SetSize(N,-1.0);
    std::cout << N << " " << J.itsData.size() <<" " << K.itsData.size() << std::endl;

    for (index_t ia:sg->es.indices())
        for (index_t ib:sg->es.indices(ia))
            for (index_t ic:sg->es.indices())
                for (index_t id:sg->es.indices(ic))
                {
                    bool doJ = sg->Ls(ia)==sg->Ls(ib) && sg->Ls(ic)==sg->Ls(id) && J(ia,ib,ic,id)==-1.0;
                    bool doK = sg->Ls(ia)==sg->Ls(ic) && sg->Ls(ib)==sg->Ls(id) && K(ia,ib,ic,id)==-1.0;
                    if (doJ || doK)
                    {
                        double norm=sg->ns(ia)*sg->ns(ib)*sg->ns(ic)*sg->ns(id);
                        SlaterIntegrals R(sg->es(ia)+sg->es(ib),sg->es(ic)+sg->es(id));
                        if (doJ)
                            J(ia,ib,ic,id)=FourPi2*R(0,sg->Ls(ia),sg->Ls(ib),sg->Ls(ic),sg->Ls(id))*norm;
                        if (doK)
                            K(ia,ib,ic,id)=FourPi2*R.DoExchangeSum(sg->Ls(ia),sg->Ls(ib),sg->Ls(ic),sg->Ls(id))*norm;
                        else
                            K(ia,ib,ic,id)=0.0;
//                        std::cout << "L=(" << sg->Ls(ia) << "," << sg->Ls(ib) << "," << sg->Ls(ic) << "," << sg->Ls(id) 
//                        << ") abcd=(" << ia << "," << ib << "," << ic << "," << id << ")  J=" << J(ia,ib,ic,id) << std::endl;
                                
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
            double t=a->es(i)+a->es(j);
            int L=a->Ls(i),L1=L+1;
            Hk(i,j)=0.5*a->ns(i)*a->ns(j)*
                   (
                       (L1*L1 + L*L1) * GaussianIntegral(t,2*L-2)
                       -2*L1 * t      * GaussianIntegral(t,2*L  )
                       +4*a->es(i)*a->es(j) * GaussianIntegral(t,2*L+2)
                   );
        }

    return Hk;
}
//
IntegralEngine::SMat IntegralEngine::MakeNuclear(iec_t* iea,const Cluster& cl) const
{
    auto a=dcast(iea);;
    size_t N=a->size(),L=a->Ls(1);
    SMatrix<double> Hn(N);
    double Z=-cl.GetNuclearCharge();
    for (auto i:Hn.rows())
        for (auto j:Hn.cols(i))
            Hn(i,j)= Z*GaussianIntegral(a->es(i)+a->es(j),2*L-1)*a->ns(i)*a->ns(j);

    return Hn;
}

IntegralEngine::RVec IntegralEngine::MakeNormalization(iec_t* iea) const
{

    auto a=dcast(iea);;
    RVec n(a->size());
    for (auto i:a->es.indices())  n(i)=GaussianNorm(a->es(i),a->Ls(i));
    return n;
}

IntegralEngine::RVec IntegralEngine::MakeCharge(iec_t* iea) const
{
    auto a=dcast(iea);;
    RVec c(a->size());
    for (auto i:a->es.indices())  c(i)=GaussianIntegral(a->es(i),a->Ls(i))*a->ns(i);
    return c;
}

std::ostream& IntegralEngine::Write(std::ostream& os) const
{
    return os ;
}
std::istream& IntegralEngine::Read (std::istream& is)
{
    return is ;
}



} //namespace
