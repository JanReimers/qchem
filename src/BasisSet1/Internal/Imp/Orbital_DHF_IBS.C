// File: BasisSet1/Imp/Orbital_DHF_IBS.C
module;
#include <cassert>
#include <blaze/Math.h>
module qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace BasisSet
{

template <class T> const mat_t<T>& Orbital_RKBL_IBS<T>::Kinetic(const Orbital_RKBS_IBS<T>& rkbs) const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Has(IntegralsCache_Base::I2x::Kinetic,
            IntegralsCache_Base::IBS_ID_t(     RadialID(),     AngularID()),
            IntegralsCache_Base::IBS_ID_t(rkbs.RadialID(),rkbs.AngularID())
        )
        ? cache->GetMat() : cache->Set(MakeKinetic(rkbs));
    
}

template class Orbital_RKBL_IBS<double>;
template class Orbital_RKBS_IBS<double>;


template <class T> Orbital_RKB_IBS_Imp<T>::Orbital_RKB_IBS_Imp(Orbital_RKBL_IBS<T>* rkbl,Orbital_RKBS_IBS<T>* rkbs) 
    : itsRKBL(rkbl)
    , itsRKBS(rkbs)
{
    assert(itsRKBL);
    assert(itsRKBS);
}

template <class T> smat_t<T> Orbital_RKB_IBS_Imp<T>::merge_diag(const smat_t<T>& l,const smat_t<T>& s)
{
    size_t Nl=l.rows();
    size_t Ns=s.rows();
    smat_t<T> ls=zero<T>(Nl+Ns);
    for (auto i:iv_t(0,Nl))
        for (auto j:iv_t(i,Nl))
            ls(i,j)=l(i,j);
    for (auto i:iv_t(0,Ns))
        for (auto j:iv_t(i,Ns))
            ls(Nl+i,Nl+j)=s(i,j);
    return ls;
}
template <class T> smat_t<T> Orbital_RKB_IBS_Imp<T>::merge_off_diag(const mat_t<T>& ls)
{
    size_t Nl=ls.rows();
    size_t Ns=ls.columns();
    assert(Nl==Ns);
    smat_t<T> k=zero<T>(Nl+Ns);
    for (auto i:iv_t(0,Nl))
        for (auto j:iv_t(0,Nl))
            k(i,Ns+j)=ls(i,j);
   
    return k;
}

// Do not call the cached versions itsRKBL->Overlap() from here.
template <class T> smat_t<T> Orbital_RKB_IBS_Imp<T>::MakeOverlap() const
{
    smat_t<T> ol=itsRKBL->MakeOverlap();
    smat_t<T> os=itsRKBS->MakeOverlap();
    return merge_diag(ol,os);
}
template <class T> smat_t<T> Orbital_RKB_IBS_Imp<T>::MakeKinetic() const
{
    mat_t<T> kls=2*c_light*itsRKBL->MakeKinetic(*itsRKBS);
    return merge_off_diag(kls);
}
template <class T> smat_t<T> Orbital_RKB_IBS_Imp<T>::MakeNuclear(const Cluster* c) const
{
    smat_t<T> nl=itsRKBL->MakeNuclear(c);
    smat_t<T> ns=itsRKBS->MakeNuclear(c);
    return merge_diag(nl,ns);
}
template <class T> smat_t<T> Orbital_RKB_IBS_Imp<T>::MakeRestMass() const
{
    smat_t<T> rs=itsRKBS->MakeOverlap();
    smat_t<T> rl=zero<T>(rs.rows());
    return merge_diag(rl,rs);
}

template class Orbital_RKB_IBS_Imp<double>;

}