// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include "Imp/Containers/ERI4.H"
#include <Cluster.H>
#include "oml/smatrix.h"

double AtomIE::FourPi2=4*4*pi*pi;

void AtomIE::Append(IrrepIEClient* iec)
{
    AnalyticIE<double>::Append(iec);
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

ERIJ1 AtomIE::MakeDirect  (const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient* >(_a);
    const AtomIrrepIEClient* c=dynamic_cast<const AtomIrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERIJ1 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        size_t iea=a->es_indices[ia-1]; //Unique exponent index. zero based.
        loop_1(iea); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            size_t iec=c->es_indices[ic-1]; //Unique exponent index. zero based.
            loop_2(iec);
            int la=a->l, lc=c->l;
            int ma=a->m, mc=c->m;
            RVec Akac=Coulomb_AngularIntegrals(la,lc,ma,mc);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                SMat& Jab=J(ia,ib);
                size_t ieb=a->es_indices[ib-1]; //Unique exponent index. zero based.
                loop_3(ieb);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    assert(la==a->l);
                    assert(lc==c->l);
                    size_t ied=c->es_indices[id-1]; //Unique exponent index. zero based.
                    if (Jab(ic,id)!=0.0)
                    {
                        cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << Jab(ic,id) << endl;    
                        assert(false);
                    }
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    RVec Rkac=loop_4_direct(ied,la,lc);
                    Jab(ic,id)=FourPi2*Akac*Rkac*norm;
//                    const SphericalGaussianCD* cd=iec->loop_4(id);
//                    J(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lc+1)*Akac*cd->Coulomb_Rk(la,lc)*norm;
                }
            }
        }
    }
    return J;
};

ERIJ1 AtomIE::MakeExchange(const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient* >(_a);
    const AtomIrrepIEClient* c=dynamic_cast<const AtomIrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERIJ1 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            int la=a->l, lc=c->l;
            RVec Akac=ExchangeAngularIntegrals(la,lc,a->m,c->m);
            for (size_t ib:a->indices(ia))
            {
                SMat& Kab=K(ia,ib);
                loop_2(a->es_indices[ib-1]);
                loop_3(c->es_indices[ic-1]);
                for (size_t id:c->indices())
                {
                    //if (ia==ib && id<ic) continue;
//                    if (ic<id && Kab(ic,id)!=0.0)
//                    {
//                        cout << "overwriting Knew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
//                        cout << Kab(ic,id) << endl;    
//                        assert(false);
//                    }
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    RVec RKac=loop_4_exchange(c->es_indices[id-1],la,lc);
                    if (ic==id)
                        Kab(ic,id)=FourPi2*Akac*RKac*norm; 
                    else if (id<ic)
                        Kab(id,ic)+=0.5*FourPi2*Akac*RKac*norm; 
                    else
                        Kab(ic,id)+=0.5*FourPi2*Akac*RKac*norm; 
//                    if (ia==ib) K(ia,ib,id,ic)=K(ia,ib,ic,id); //ERIK container does support this symmetry yet.
//                    else
//                    {
//                        if (K(ia,ib,id,ic)!=0.0)
//                        {
//                            double Kavg=0.5*(K(ia,ib,id,ic)+K(ia,ib,ic,id));
//                            K(ia,ib,id,ic)=K(ia,ib,ic,id)=Kavg;
//                        }
//                    }
//                    cout << "Knew(" << ia << " " << ic << " " << ib << " " << id << ")="; 
//                    cout << std::setprecision(8) << K(ia,ic,ib,id) << endl;    

                }
            }
        }
    }

    return K;
};

