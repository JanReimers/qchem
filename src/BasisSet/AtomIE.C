// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include "Imp/Containers/ERI4.H"
#include <Cluster.H>

using std::cout;
using std::endl;

void AtomIE::Append(const IrrepIEClient* ciec)
{
    AnalyticIE<double>::Append(ciec);
    IrrepIEClient* iec=const_cast<IrrepIEClient*>(ciec);
    AtomIrrepIEClient* aiec=dynamic_cast<AtomIrrepIEClient*>(iec);
    BFGrouper::Append(aiec);
}


const AtomIrrepIEClient* AtomIE::dcast(iec_t* iea)
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(iea);
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
            s(i,j)=Overlap(a->es(i),a->es(j),2*a->l)*a->ns(i)*a->ns(j);

    return s;
}

AtomIE::SMat AtomIE::MakeKinetic(iec_t* iea) const
{
    auto a=dcast(iea);;
    size_t N=a->size();
    SMatrix<double> Hk(N);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols(i))
            Hk(i,j)=Kinetic(a->es(i),a->es(j),a->l)*a->ns(i)*a->ns(j);

    return Hk;
}

AtomIE::Mat AtomIE::MakeKinetic(iec_t* iea, iec_t* ieb) const
{
    auto a=dcast(iea);;
    auto b=dcast(ieb);;
    size_t Na=a->size();
    size_t Nb=b->size();
    Matrix<double> Hk(Na,Nb);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols())
            Hk(i,j)=Kinetic(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*b->ns(j);

    return Hk;
}

AtomIE::SMat AtomIE::MakeNuclear(iec_t* iea,const Cluster& cl) const
{
    auto a=dcast(iea);;
    size_t N=a->size(),L=a->l;
    SMatrix<double> Hn(N);
    double Z=-cl.GetNuclearCharge();
    for (auto i:Hn.rows())
        for (auto j:Hn.cols(i))
            Hn(i,j)= Z*Nuclear(a->es(i),a->es(j),L)*a->ns(i)*a->ns(j);

    return Hn;
}

AtomIE::SMat AtomIE::MakeRestMass(iec_t* iea) const
{
    assert(false);
    return SMat();
}

AtomIE::RVec AtomIE::MakeCharge(iec_t* iea) const
{
    auto a=dcast(iea);;
    RVec c(a->size());
    for (auto i:a->es.indices())  c(i)=Charge(a->es(i),a->l)*a->ns(i);
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
            r(i,j)=Repulsion(a->es(i),a->es(j),a->l,a->l)*a->ns(i)*a->ns(j);

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
            s(i,j)=Repulsion(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*a->ns(j);

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
            s(i,j)=Overlap(ab->es(i)+ab->es(j),ec,ab->l+ab->l+Lc)*ab->ns(i)*ab->ns(j)*nc;            

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
            s(i,j)=Repulsion(ab->es(i)+ab->es(j),ec,ab->l,Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}


AtomIE::ERI3 AtomIE::MakeRepulsion3C(iec_t* ieab,iec_t* iec) const
{
    auto c=dcast(iec);;
    
    ERI3 s3;
    for (auto i:c->es.indices()) s3.push_back(MakeRepulsion(ieab,(*c)(i)));
    return s3;
}

ERI4 AtomIE::MakeDirect  (const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient* >(_a);
    const AtomIrrepIEClient* c=dynamic_cast<const AtomIrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            RVec Akac=Coulomb_AngularIntegrals(la,lc,a->m,c->m);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                loop_3(a->es_indices[ib-1]);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << Jab(ic,id) << endl;    
                        assert(false);
                    }
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    RVec Rkac=loop_4_direct(c->es_indices[id-1],la,lc);
                    Jab(ic,id)=Akac*Rkac*norm;
                }
            }
        }
    }
    return J;
};

ERI4 AtomIE::MakeExchange(const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient* >(_a);
    const AtomIrrepIEClient* c=dynamic_cast<const AtomIrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        double na=a->ns(ia);
        for (size_t ic:c->indices())
        {
            int la=a->l, lc=c->l;
            RVec Akac=ExchangeAngularIntegrals(la,lc,a->m,c->m);
            double nac=na*c->ns(ic);
            for (size_t ib:a->indices(ia))
            {
                SMat& Kab=K(ia,ib);
                loop_2(a->es_indices[ib-1]);
                loop_3(c->es_indices[ic-1]);
                double nacb=nac*a->ns(ib);
                for (size_t id:c->indices())
                {
                    double norm=nacb*c->ns(id);
                    RVec RKac=loop_4_exchange(c->es_indices[id-1],la,lc);
                    if (ic==id)
                        Kab(ic,id)=Akac*RKac*norm; 
                    else if (id<ic)
                        Kab(id,ic)+=0.5*Akac*RKac*norm; 
                    else
                        Kab(ic,id)+=0.5*Akac*RKac*norm; 

                }
            }
        }
    }

    return K;
};

