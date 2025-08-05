// File: AtomIE.C Common HF IE code for all atom basis sets.
module;
#include <iostream>
#include <memory>
#include <cassert>
module qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IEClient;

template <class T> void AtomIE_BS_2E<T>::Append(const IrrepIEClient* ciec)
{
    assert(ciec);
    DB_BS_2E<T>::Append(ciec);
    IrrepIEClient* iec=const_cast<IrrepIEClient*>(ciec);
    AtomIrrepIEClient* aiec=AtomIrrepIEClient::dcast(iec);
    BFGrouper::Append(aiec);
}
template <class T> ERI4 AtomIE_BS_2E<T>::MakeDirect  (const IrrepIEClient* _a, const IrrepIEClient* _c) const 
{
    typedef SMatrix<T> SMat;
    auto a=AtomIrrepIEClient::dcast(_a);
    auto c=AtomIrrepIEClient::dcast(_c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c->indices())
        {
            loop_2(c->es_indices[ic-1]);
            int la=a->l, lc=c->l;
            RVec Akac=itsAngular->Coulomb_AngularIntegrals(a,c);
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
    auto a=AtomIrrepIEClient::dcast(_a);
    auto c=AtomIrrepIEClient::dcast(_c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 K(Na,Nc);
    for (size_t ia:a->indices())
    {
        loop_1(a->es_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
        double na=a->ns(ia);
        for (size_t ic:c->indices())
        {
            int la=a->l, lc=c->l;
            RVec Akac=itsAngular->ExchangeAngularIntegrals(a,c);
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

