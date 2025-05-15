// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/Atom/IE.H"
#include "Imp/BasisSet/Atom/IEC.H"
#include "Imp/Containers/ERI4.H"
#include <Cluster.H>
#include <Irrep_BS.H>

template <class T> typename Integrals_Base<T>::SMat AtomIE_Overlap <T>::MakeOverlap() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= this->Overlap(a->es(i),a->es(j),2*l)*a->ns(i)*a->ns(j);

    return H;
}
template <class T> typename Integrals_Base<T>::SMat AtomIE_Grad2 <T>::MakeGrad2() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= this->Grad2(a->es(i),a->es(j),l,l)*a->ns(i)*a->ns(j);

    return H;
}
template <class T> typename Integrals_Base<T>::SMat AtomIE_Nuclear <T>::MakeNuclear(const Cluster* cl) const
{
    assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Z*Nuclear(a->es(i),a->es(j),2*l)*a->ns(i)*a->ns(j);

    return H;
}

#include <DHF_IBS.H>

template <class T> typename Integrals_Base<T>::Mat  AtomIE_XGrad2<T>::MakeGrad2(const Orbital_RKBS_IBS<T>* rkbs) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    const AtomIrrepIEClient* b=dynamic_cast<const AtomIrrepIEClient*>(rkbs);
    assert(a->l==b->l);
    size_t Na=a->size();
    size_t Nb=b->size();
    Matrix<double> Hk(Na,Nb);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols())
            Hk(i,j)=Grad2(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*b->ns(j);

    return Hk;
}

#include "Imp/BasisSet/Atom/IE_DFT.H"
template <class T> typename AtomIE_DFT<T>::ERI3 AtomIE_DFT<T>::MakeOverlap3C  (const fbs_t& _c) const
{
    const AtomIrrepIEClient& c=dynamic_cast<const AtomIrrepIEClient&>(_c);
    ERI3 s3;
    for (auto i:c.indices()) s3.push_back(MakeOverlap(c(i)));
    return s3;
}
template <class T> typename AtomIE_DFT<T>::ERI3 AtomIE_DFT<T>::MakeRepulsion3C(const fbs_t& _c) const
{
    const AtomIrrepIEClient& c=dynamic_cast<const AtomIrrepIEClient&>(_c);
    ERI3 s3;
    for (auto i:c.indices()) s3.push_back(MakeRepulsion(c(i)));
    return s3;
}
template <class T> typename AtomIE_DFT<T>::SMat AtomIE_DFT<T>::MakeOverlap  (const bf_tuple& c) const
{    
    const AtomIrrepIEClient* ab=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(ab);
    size_t N=ab->size();
    int Nc,Lc,Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=this->Overlap(ab->es(i)+ab->es(j),ec,ab->l+ab->l+Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}
template <class T> typename AtomIE_DFT<T>::SMat AtomIE_DFT<T>::MakeRepulsion(const bf_tuple& c) const
{    
    const AtomIrrepIEClient* ab=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(ab);
    size_t N=ab->size();
    int Nc,Lc,Mc;
    double ec,nc;
    std::tie(Nc,Lc,Mc,ec,nc)=c;
    SMat s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=this->Repulsion(ab->es(i)+ab->es(j),ec,ab->l,Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}

template class AtomIE_DFT<double>;


#include "Imp/BasisSet/Atom/IE_HF.H"

template <class T> void AtomIE_BS_2E<T>::Append(const IrrepIEClient* ciec)
{
    assert(ciec);
    DB_BS_2E<T>::Append(ciec);
    IrrepIEClient* iec=const_cast<IrrepIEClient*>(ciec);
    AtomIrrepIEClient* aiec=dynamic_cast<AtomIrrepIEClient*>(iec);
    assert(aiec);
    BFGrouper::Append(aiec);
}
template <class T> ERI4 AtomIE_BS_2E<T>::MakeDirect  (const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    typedef SMatrix<T> SMat;
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient* >(_a);
    const AtomIrrepIEClient* c=dynamic_cast<const AtomIrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
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
                        std::cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        std::cout << Jab(ic,id) << std::endl;    
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
template <class T> ERI4 AtomIE_BS_2E<T>::MakeExchange(const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    typedef SMatrix<T> SMat;
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient* >(_a);
    const AtomIrrepIEClient* c=dynamic_cast<const AtomIrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
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
template class AtomIE_BS_2E<double>;

#include "Imp/BasisSet/Atom/IE_Fit.H"
#include <Fit_IBS.H>

AtomIE_Fit::Vec  AtomIE_Fit::MakeCharge() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);
    Vec c(a->size());
    for (auto i:a->es.indices())  c(i)=Charge(a->es(i),a->l)*a->ns(i);
    return c;
}
AtomIE_Fit::SMat AtomIE_Fit::MakeRepulsion() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Repulsion(a->es(i),a->es(j),l,l)*a->ns(i)*a->ns(j);

    return H;
}
AtomIE_Fit::Mat  AtomIE_Fit::MakeRepulsion(const fbs_t& _b) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    const AtomIrrepIEClient* b=dynamic_cast<const AtomIrrepIEClient*>(&_b);
    assert(a);
    assert(b);
    size_t Na=a->es.size(), Nb=b->es.size();
    Mat s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=this->Repulsion(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*a->ns(j);

    return s;
}

#include "Imp/BasisSet/Atom/IE_DHF.H"

template class AtomIE_RKBL<double>;
template class AtomIE_RKBS<double>;
