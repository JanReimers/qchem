// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include <Cluster.H>
#include "oml/smatrix.h"

const AtomIEClient* AtomIE::dcast(iec_t* iea)
{
    const AtomIEClient* a=dynamic_cast<const AtomIEClient*>(iea);
    assert(a);
    return a;
}


AtomIE::SMat AtomIE::MakeOverlap(iec_t* iea ) const
{
    auto  a=dcast(iea);
    size_t N=a->size();
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=Overlap(a->es(i),a->es(j),2*a->Ls(i))*a->ns(i)*a->ns(j);

    return s;
}

AtomIE::SMat AtomIE::MakeKinetic(iec_t* iea) const
{
    auto a=dcast(iea);;
    size_t N=a->size();
    SMatrix<double> Hk(N);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols(i))
        {
            Hk(i,j)=Kinetic(a->es(i),a->es(j),a->Ls(i))*a->ns(i)*a->ns(j);
            assert(a->Ls(i)==a->Ls(j));
        }

    return Hk;
}

AtomIE::SMat AtomIE::MakeNuclear(iec_t* iea,const Cluster& cl) const
{
    auto a=dcast(iea);;
    size_t N=a->size(),L=a->Ls(1);
    SMatrix<double> Hn(N);
    double Z=-cl.GetNuclearCharge();
    for (auto i:Hn.rows())
        for (auto j:Hn.cols(i))
            Hn(i,j)= Z*Nuclear(a->es(i),a->es(j),L)*a->ns(i)*a->ns(j);

    return Hn;
}

AtomIE::RVec AtomIE::MakeCharge(iec_t* iea) const
{
    auto a=dcast(iea);;
    RVec c(a->size());
    for (auto i:a->es.indices())  c(i)=Charge(a->es(i),a->Ls(i))*a->ns(i);
    return c;
}

AtomIE::SMat AtomIE::MakeRepulsion(iec_t* iea ) const
{
    auto a=dcast(iea);;
    assert(a);
    size_t N=a->size();
    SMat r(N,N);
    for (auto i:r.rows())
        for (auto j:r.cols(i))
            //r(i,j)=GaussianRepulsionIntegral(a->es(i),a->es(j),a->Ls(i),a->Ls(j))*a->ns(i)*a->ns(j);
            r(i,j)=Repulsion(a->es(i),a->es(j),a->Ls(i),a->Ls(j))*a->ns(i)*a->ns(j);

    return r;
}

AtomIE::Mat AtomIE::MakeRepulsion(iec_t* iea,iec_t* ieb) const
{
    auto a=dcast(iea);;
    auto b=dcast(ieb);;
    size_t Na=a->es.size(), Nb=b->es.size();
    Mat s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=Repulsion(a->es(i),b->es(j),a->Ls(i),b->Ls(j))*a->ns(i)*a->ns(j);
//            s(i,j)=GaussianRepulsionIntegral(a->es(i),b->es(j),a->Ls(i),b->Ls(j))*a->ns(i)*b->ns(j);

    return s;
}


AtomIE::SMat AtomIE::MakeOverlap(iec_t* ieab, const bf_tuple& c) const
{    
    auto ab=dcast(ieab);;
    size_t N=ab->size();
    int Nc,Lc,Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
        {
            assert(ab->Ls(i)==ab->Ls(j));
            assert(Lc==0); //Non-polarized fit basis
            //assert(ab->Ls(i)==Lc); //TODO what going on here?
            s(i,j)=Overlap(ab->es(i)+ab->es(j),ec,ab->Ls(i)+ab->Ls(j)+Lc)*ab->ns(i)*ab->ns(j)*nc;            
        }
    return s;
}


AtomIE::ERI3 AtomIE::MakeOverlap3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);;
   
    ERI3 s3;
    for (auto i:c->es.indices()) s3.push_back(MakeOverlap(ieab,(*c)(i)));
    return s3;
}


AtomIE::SMat AtomIE::MakeRepulsion(iec_t* ieab,const bf_tuple& c) const
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
            assert(ab->Ls(i)==ab->Ls(j));
            s(i,j)=Repulsion(ab->es(i)+ab->es(j),ec,ab->Ls(i),Lc)*ab->ns(i)*ab->ns(j)*nc;            
        }
    return s;
}


AtomIE::ERI3 AtomIE::MakeRepulsion3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);;
    
    ERI3 s3;
    for (auto i:c->es.indices()) s3.push_back(MakeRepulsion(ieab,(*c)(i)));
    return s3;
}


