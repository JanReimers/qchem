// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include "Imp/Containers/ERI4.H"
#include <Cluster.H>
#include <BasisSet.H>

using std::cout;
using std::endl;


template <class T> typename AtomIE_1E<T>::SMat AtomIE_1E<T>::MakeOverlap() const
{
    return MakeIntegrals(qchem::Overlap1);
}

template <class T> typename AtomIE_1E<T>::SMat AtomIE_1E<T>::MakeKinetic() const
{
    return MakeIntegrals(qchem::Kinetic1);
}

template <class T> typename AtomIE_1E<T>::SMat AtomIE_1E<T>::MakeNuclear(const Cluster* cl) const
{
    assert(cl);
    assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
    int Z=-cl->GetNuclearCharge(); 
    return Z*MakeIntegrals(qchem::Nuclear1,cl);
}

template <class T> typename AtomIE_1E<T>::SMat AtomIE_1E<T>::MakeIntegrals(qchem::IType t,const Cluster* cl) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Integral(t,a->es(i),a->es(j),l)*a->ns(i)*a->ns(j);

    return H;
}

template class AtomIE_1E<double>;

template <class T> typename AtomIE_DFT<T>::ERI3 AtomIE_DFT<T>::MakeOverlap3C(const bs_t& _c) const
{
    const AtomIrrepIEClient& c=dynamic_cast<const AtomIrrepIEClient&>(_c);
    ERI3 s3;
    for (auto i:c.indices()) s3.push_back(MakeOverlap(c(i)));
    return s3;
}
template <class T> typename AtomIE_DFT<T>::ERI3 AtomIE_DFT<T>::MakeRepulsion3C(const bs_t& _c) const
{
    const AtomIrrepIEClient& c=dynamic_cast<const AtomIrrepIEClient&>(_c);
    ERI3 s3;
    for (auto i:c.indices()) s3.push_back(MakeRepulsion(c(i)));
    return s3;
}
template <class T> typename AtomIE_DFT<T>::SMat AtomIE_DFT<T>::MakeOverlap(const bf_tuple& c) const
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
            s(i,j)=Overlap(ab->es(i)+ab->es(j),ec,ab->l+ab->l+Lc)*ab->ns(i)*ab->ns(j)*nc;            

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
            s(i,j)=Repulsion(ab->es(i)+ab->es(j),ec,ab->l,Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}

template class AtomIE_DFT<double>;



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
template class AtomIE_BS_2E<double>;

AtomIE_Fit::Vec AtomIE_Fit::MakeCharge() const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);
    Vec c(a->size());
    for (auto i:a->es.indices())  c(i)=Charge(a->es(i),a->l)*a->ns(i);
    return c;
}

AtomIE_Fit::SMat AtomIE_Fit::MakeOverlap() const
{
    return MakeIntegrals(qchem::Overlap1);
}
AtomIE_Fit::SMat AtomIE_Fit::MakeRepulsion() const
{
    return MakeIntegrals(qchem::Repulsion1);
}
AtomIE_Fit::Mat AtomIE_Fit::MakeRepulsion(const bs_t& _b) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    const AtomIrrepIEClient* b=dynamic_cast<const AtomIrrepIEClient*>(&_b);
    assert(a);
    assert(b);
    size_t Na=a->es.size(), Nb=b->es.size();
    Mat s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=Repulsion(a->es(i),b->es(j),a->l,b->l)*a->ns(i)*a->ns(j);

    return s;
}

AtomIE_Fit::SMat AtomIE_Fit::MakeIntegrals(qchem::IType t,const Cluster* cl) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Integral(t,a->es(i),a->es(j),l)*a->ns(i)*a->ns(j);

    return H;
}


#include <BasisSet.H>

// template <class T> typename AtomIE_RKBS<T>::Mat AtomIE_RKBS<T>::MakeKinetic(const IrrepBasisSet* L) const
// {
//     const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
//     const AtomIrrepIEClient* b=dynamic_cast<const AtomIrrepIEClient*>(L);
//     assert(a->l==b->l);
//     size_t Na=a->size();
//     size_t Nb=b->size();
//     Matrix<double> Hk(Na,Nb);
//     for (auto i:Hk.rows())
//         for (auto j:Hk.cols())
//             Hk(i,j)=Integral(qchem::Kinetic1,b->es(j),a->es(i),a->l)*a->ns(j)*b->ns(i);

//     return Hk;
// }

// template class AtomIE_RKBS<double>;