// File: Imp/DHF_IBS_Common.C  Common implementation for all Dirac-Hartree-Fock (HF) Irrep Basis Sets.
module;
#include <vector>
#include <memory>
#include <iostream>
#include <cassert>
#include <ranges>
#include "blaze/Math.h"
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.HeapDB;
import qchem.Conversions;

template <class T> Orbital_RKB_IBS_Common<T>::Orbital_RKB_IBS_Common
(const DB_cache<T>* db, int kappa,::Orbital_RKBL_IBS<T>* l,::Orbital_RKBS_IBS<T>* s)
    : Orbital_IBS_Common<T>()
    , DB_RKB<T>(db)
    , itsRKBL(l)
    , itsRKBS(s)
{
    assert(itsRKBL);
    assert(itsRKBS);
    s->Insert(itsRKBL);
}

template <class T> vec_t<T> Orbital_RKB_IBS_Common<T>::operator() (const rvec3_t&) const
{
    assert(false);
    return vec_t<T>();
}
template <class T> vec3vec_t<T> Orbital_RKB_IBS_Common<T>::Gradient   (const rvec3_t&) const
{
    assert(false);
    return vec3vec_t<T>();
}

template <class T> smat_t<T> Orbital_RKB_IBS_Common<T>::merge_diag(const smat_t<T>& l,const smat_t<T>& s)
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
template <class T> smat_t<T> Orbital_RKB_IBS_Common<T>::merge_off_diag(const mat_t<T>& ls)
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
template <class T> smat_t<T> Orbital_RKB_IBS_Common<T>::MakeOverlap() const
{
    smat_t<T> ol=itsRKBL->Overlap();
    smat_t<T> os=itsRKBS->Kinetic();
    return merge_diag(ol,os);
}
template <class T> smat_t<T> Orbital_RKB_IBS_Common<T>::MakeKinetic() const
{
    mat_t<T> kls=-itsRKBL->Kinetic(itsRKBS);
    return merge_off_diag(kls);
}
template <class T> smat_t<T> Orbital_RKB_IBS_Common<T>::MakeNuclear(const Cluster* c) const
{
    smat_t<T> nl=itsRKBL->Nuclear(c);
    smat_t<T> ns=itsRKBS->Nuclear(c);
    return merge_diag(nl,ns);
}
template <class T> smat_t<T> Orbital_RKB_IBS_Common<T>::MakeRestMass() const
{
    smat_t<T> rl=zero<T>(itsRKBL->GetNumFunctions());
    smat_t<T> rs=itsRKBS->Kinetic();
    return merge_diag(rl,rs);
}

template <class T> Orbital_RKBL_IBS_Common<T>::Orbital_RKBL_IBS_Common(int _kappa)
    : kappa(_kappa)
{
    assert(kappa!=0);
}

template <class T> Orbital_RKBS_IBS_Common<T>::Orbital_RKBS_IBS_Common(int _kappa)
    : kappa(_kappa)
    , large(0)
{
    assert(kappa!=0);
}

template <class T> void Orbital_RKBS_IBS_Common<T>::Insert(const Orbital_RKBL_IBS<T>* l)
{
    assert(l);
    large=l;
}
template class Orbital_RKB_IBS_Common<double>;
template class Orbital_RKBL_IBS_Common<double>;
template class Orbital_RKBS_IBS_Common<double>;