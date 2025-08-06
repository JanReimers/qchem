// File: BSpline/IE.C Common IE code for BSpline basis sets.
module;
#include <cassert>
#include <tuple>
#include <iostream>

module qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.Fit_IBS;
import qchem.IrrepBasisSet;
import qchem.BasisSet.Internal.ERI4;
import qchem.Orbital_DHF_IBS;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;

namespace BSpline
{
template <class T,size_t K> SMatrix<T> IE_Overlap<T,K>::MakeOverlap() const
{
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= pie->Overlap((*a)(i),(*a)(j),2*l)*a->ns(i)*a->ns(j);

    return H;
}
template <class T,size_t K> SMatrix<T> IE_Kinetic  <T,K>::MakeKinetic() const
{
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= (pie->Grad2((*a)(i),(*a)(j),l,l) + l*(l+1)*pie->Inv_r2((*a)(i),(*a)(j),l))*a->ns(i)*a->ns(j);

    return H;
}
template <class T,size_t K> SMatrix<T> IE_Inv_r1<T,K>::MakeNuclear(const Cluster* cl) const
{
    assert(cl);
    assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
    int Z=-cl->GetNuclearCharge(); 
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Z*pie->Inv_r1((*a)(i),(*a)(j),2*l)*a->ns(i)*a->ns(j);

    return H;
}
template <class T,size_t K> ERI3<T> IE_DFT<T,K>::MakeOverlap3C  (const Fit_IBS& _c) const
{
    const IrrepIEClient<K>& c=dynamic_cast<const IrrepIEClient<K>&>(_c);
    ERI3<T> s3;
    for (auto i:c.indices()) s3.push_back(MakeOverlap(c.tuple(i)));
    return s3;
}
template <class T,size_t K> ERI3<T> IE_DFT<T,K>::MakeRepulsion3C(const Fit_IBS& _c) const
{
    const IrrepIEClient<K>& c=dynamic_cast<const IrrepIEClient<K>&>(_c);
    ERI3<T> s3;
    for (auto i:c.indices()) s3.push_back(MakeRepulsion(c.tuple(i)));
    return s3;
}
template <class T,size_t K> SMatrix<T> IE_DFT<T,K>::MakeOverlap  (const bf_tuple& c) const
{    
    const IrrepIEClient<K>* ab=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(ab);
    size_t N=ab->size();
    int Lc;
    const spline_t* sc;
    double nc;
    std::tie(Lc,sc,nc)=c;
    SMatrix<T> s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=pie->Overlap((*ab)(i)+(*ab)(j),*sc,ab->l+ab->l+Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}
template <class T,size_t K> SMatrix<T> IE_DFT<T,K>::MakeRepulsion(const bf_tuple& c) const
{    
    const IrrepIEClient<K>* ab=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(ab);
    size_t N=ab->size();
    int Lc;
    const spline_t* sc;
    double nc;
    std::tie(Lc,sc,nc)=c;
    SMatrix<T> s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=pie->Repulsion((*ab)(i)+(*ab)(j),*sc,ab->l,Lc)*ab->ns(i)*ab->ns(j)*nc;            

    return s;
}
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
template <size_t K> Vector<double>  IE_Fit<K>::MakeCharge() const
{
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(a);
    Vector<double> c(a->size());
    for (auto i:a->indices())  c(i)=pie->Charge((*a)(i),a->l)*a->ns(i);
    return c;
}
template <size_t K> SMatrix<double>  IE_Fit<K>::MakeRepulsion() const
{
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>*>(this);
    assert(a);

    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= pie->Repulsion((*a)(i),(*a)(j),l,l)*a->ns(i)*a->ns(j);

    return H;
}
template <size_t K> Matrix<double> IE_Fit<K>::MakeRepulsion(const Fit_IBS& _b) const
{
    const IrrepIEClient<K>* a=dynamic_cast<const IrrepIEClient<K>*>(this);
    const IrrepIEClient<K>* b=dynamic_cast<const IrrepIEClient<K>*>(&_b);
    assert(a);
    assert(b);
    size_t Na=a->size(), Nb=b->size();
    Matrix<double> s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=pie->Repulsion((*a)(i),(*b)(j),a->l,b->l)*a->ns(i)*a->ns(j);

    return s;
}

#define INSTANCEk(k) template class IE_Overlap<double,k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class IE_Kinetic<double,k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class IE_Inv_r1<double,k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class IE_DFT<double,k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class IE_BS_2E<double,k>;
#include "../Instance.hpp"
#define INSTANCEk(k) template class IE_Fit<k>;
#include "../Instance.hpp"

} //namespace
