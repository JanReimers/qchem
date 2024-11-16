// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include "Imp/Containers/ERI4.H"
#include <Cluster.H>
#include "oml/smatrix.h"

double AtomIE::FourPi2=4*4*pi*pi;

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

ERIJ AtomIE::MakeDirect  (const iriec* a, const iriec* c,const AtomIEClient* aiec) const 
{
    size_t Na=a->size(), Nc=c->size();
    ERIJ J(Na,Nc);
    for (size_t ia:a->indices())
    {
        aiec->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (size_t ic:c->indices())
        {
            aiec->loop_2(ic);
            int la=a->Ls(ia), lc=c->Ls(ic);
            int ma=a->Ms(ia), mc=c->Ms(ic);
            RVec Akac=Coulomb_AngularIntegrals(la,lc,ma,mc);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                aiec->loop_3(ib);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    assert(la==a->Ls(ib));
                    assert(lc==c->Ls(id));
                    if (J(ia,ib,ic,id)!=0.0)
                    {
                        cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        cout << J(ia,ib,ic,id) << endl;    
                        assert(false);
                    }
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    RVec Rkac=aiec->loop_4_direct(id,la,lc);
                    J(ia,ib,ic,id)=FourPi2*Akac*Rkac*norm;
//                    const SphericalGaussianCD* cd=iec->loop_4(id);
//                    J(ia,ib,ic,id)=FourPi2*(2*la+1)*(2*lc+1)*Akac*cd->Coulomb_Rk(la,lc)*norm;
                }
            }
        }
    }
    return J;
};

ERIK AtomIE::MakeExchange(const iriec* a, const iriec* b,const AtomIEClient* aiec) const 
{
    size_t Na=a->size(), Nb=b->size();
    ERIK K(Na,Nb);
    for (size_t ia:a->indices())
    {
        aiec->loop_1(ia); //Start a cache for SphericalGaussianCD*
        for (size_t ib:b->indices())
        {
            int la=a->Ls(ia), lb=b->Ls(ib);
            int ma=a->Ms(ia), mb=b->Ms(ib);
            RVec Akab=ExchangeAngularIntegrals(la,lb,ma,mb);

            for (size_t ic:a->indices())
            {
                if (ic<ia) continue;
                aiec->loop_2(ic);
                aiec->loop_3(ib);
                
                for (size_t id:b->indices())
                {
//                    if (id<ic) continue;
                    if (ia==ic && id<ib) continue;
                    //if (id<ib) continue;
                    assert(la==a->Ls(ic));
                    assert(lb==b->Ls(id));
                    if (K(ia,ic,ib,id)!=0.0)
                    {
                        cout << "overwriting Knew(" << ia << " " << ic << " " << ib << " " << id << ")="; 
                        cout << K(ia,ic,ib,id) << endl;    
                        assert(false);
                    }
                    double norm=a->ns(ia)*b->ns(ib)*a->ns(ic)*b->ns(id);
                    RVec RKab=aiec->loop_4_exchange(id,la,lb);
                    K(ia,ic,ib,id)=FourPi2*Akab*RKab*norm; 
//                    const SphericalGaussianCD* cd=iec->loop_4(id);
//                    K(ia,ic,ib,id)=FourPi2*(2*la+1)*(2*lb+1)*Akab*cd->ExchangeRk(la,lb)*norm; 
                    if (ia==ic) K(ia,ic,id,ib)=K(ia,ic,ib,id); //ERIK container does support this symmetry yet.
//                    cout << "Knew(" << ia << " " << ic << " " << ib << " " << id << ")="; 
//                    cout << std::setprecision(8) << K(ia,ic,ib,id) << endl;    

                }
            }
        }
    }

    return K;
};

void AtomIE::MakeDirect(erij_t& Jac, const ::IEClient* iec) const
{
    Jac.clear();
    const AtomIEClient& sg=*dynamic_cast<const AtomIEClient*>(iec);
    size_t NIrrep=sg.GetNumIrreps();
    for (size_t ia=1;ia<=NIrrep;ia++)
        for (size_t ic=1;ic<=NIrrep;ic++) //TODO run from ia n
        {
            const iriec* a=sg[ia];
            const iriec* c=sg[ic];
            Jac[ia][ic]=MakeDirect(a,c,&sg);
        }

}

void AtomIE::MakeExchange(erik_t& Kab, const ::IEClient* iec) const
{
    Kab.clear();
    const AtomIEClient& sg=*dynamic_cast<const AtomIEClient*>(iec);
    size_t NIrrep=sg.GetNumIrreps();
    for (size_t ia=1;ia<=NIrrep;ia++)
        for (size_t ib=1;ib<=NIrrep;ib++) //TODO run from ib 
        {
            const iriec* a=sg[ia];
            const iriec* b=sg[ib];
            Kab[ia][ib]=MakeExchange(a,b,&sg);
        }
    
}


