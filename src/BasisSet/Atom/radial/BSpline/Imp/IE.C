// File: BSpline/IE.C Common IE code for BSpline basis sets.
module;
#include <cassert>
#include <tuple>
#include <iostream>

module qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Internal.ERI4;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;

namespace BSpline
{
template <class T,size_t K> void IE_BS_2E<T,K>::Append(const ::IrrepIEClient* ciec)
{
    assert(ciec);
    DB_BS_2E<T>::Append(ciec);
    ::IrrepIEClient* iec=const_cast<::IrrepIEClient*>(ciec);
    IrrepIEClient<K>* bsiec=dynamic_cast<IrrepIEClient<K>*>(iec);
    assert(bsiec);
    BFGrouper<K>::Append(bsiec);
}
template <class T,size_t K> ERI4 IE_BS_2E<T,K>::MakeDirect  (const ::IrrepIEClient* _a, const ::IrrepIEClient* _c) const 
{
    typedef SMatrix<T> SMat;
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>* >(_a);
    const IrrepIEClient<K>* c=dynamic_cast<const IrrepIEClient<K>* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    for (size_t ia:a->indices())
    {
        size_t iau=a->sp_indices[ia-1]; //Absolute unique radial function index.
        loop_1(a->sp_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c->indices())
        {
            size_t icu=c->sp_indices[ic-1]; //Absolute unique radial function index.
            loop_2(c->sp_indices[ic-1]);
            int la=a->l, lc=c->l;
            RVec Akac=itsAngular->Coulomb_AngularIntegrals(a,c);
            for (size_t ib:a->indices())
            {
                if (ib<ia) continue; 
                size_t ibu=a->sp_indices[ib-1]; //Absolute unique radial function index.
                if (ibu>iau+K) continue;
                SMat& Jab=J(ia,ib);
                loop_3(a->sp_indices[ib-1]);
                for (size_t id:c->indices())
                {
                    if (id<ic) continue;
                    size_t idu=c->sp_indices[id-1]; //Absolute unique radial function index.

                    if (idu>icu+K) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        std::cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        std::cout << Jab(ic,id) << std::endl;    
                        assert(false);
                    }
                    // std::cout << "direct ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
                    double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                    RVec Rkac=loop_4_direct(c->sp_indices[id-1],la,lc);
                    // cout << la << " " << lc << " " << Akac << " " << Rkac << endl;
                    assert(Akac.size()==Rkac.size());
                    assert(Akac.GetLimits()==Rkac.GetLimits());
                    Jab(ic,id)=(Akac*Rkac)*norm;
                }
            }
        }
    }
    return J;
};
template <class T,size_t K> ERI4 IE_BS_2E<T,K>::MakeExchange(const ::IrrepIEClient* _a, const ::IrrepIEClient* _c) const 
{
    typedef SMatrix<T> SMat;
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>* >(_a);
    const IrrepIEClient<K>* c=dynamic_cast<const IrrepIEClient<K>* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 Kex(Na,Nc);
    for (size_t ia:a->indices())
    {
        size_t iau=a->sp_indices[ia-1]; //Absolute unique radial function index.
        loop_1(a->sp_indices[ia-1]); //Start a cache for Gaussian::RkEngine*
        double na=a->ns(ia);
        for (size_t ic:c->indices())
        {
            size_t icu=c->sp_indices[ic-1]; //Absolute unique radial function index.
            if (iau>icu+K || icu>iau+K) continue;
            int la=a->l, lc=c->l;
            RVec Akac=itsAngular->ExchangeAngularIntegrals(a,c);
            double nac=na*c->ns(ic);
            for (size_t ib:a->indices(ia))
            {
                size_t ibu=a->sp_indices[ib-1]; //Absolute unique radial function index.
                SMat& Kab=Kex(ia,ib);
                loop_2(a->sp_indices[ib-1]);
                loop_3(c->sp_indices[ic-1]);
                double nacb=nac*a->ns(ib);
                for (size_t id:c->indices())
                {
                    size_t idu=c->sp_indices[id-1]; //Absolute unique radial function index.
                    if (idu>ibu+K || ibu>idu+K) continue;
                    // std::cout << "exchange ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
                    double norm=nacb*c->ns(id);
                    RVec RKac=loop_4_exchange(c->sp_indices[id-1],la,lc);
                    assert(Akac.size()==RKac.size());
                    assert(Akac.GetLimits()==RKac.GetLimits());
                    
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

    return Kex;
};

#define INSTANCEk(k) template class IE_BS_2E<double,k>;
#include "../Instance.hpp"

} //namespace
