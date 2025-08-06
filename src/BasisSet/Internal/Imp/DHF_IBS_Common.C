// File: Imp/DHF_IBS_Common.C  Common implementation for all Dirac-Hartree-Fock (HF) Irrep Basis Sets.
module;
#include <vector>
#include <memory>
#include <iostream>
#include <cassert>
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.HeapDB;
import oml;

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

template <class T> Orbital_RKB_IBS_Common<T>::Vec     Orbital_RKB_IBS_Common<T>::operator() (const RVec3&) const
{
    assert(false);
    return Vec();
}
template <class T> Orbital_RKB_IBS_Common<T>::Vec3Vec Orbital_RKB_IBS_Common<T>::Gradient   (const RVec3&) const
{
    assert(false);
    return Vec3Vec();
}

template <class T> SMatrix<T> Orbital_RKB_IBS_Common<T>::merge_diag(const SMatrix<T>& l,const SMatrix<T>& s)
{
    size_t Nl=l.GetNumRows();
    size_t Ns=s.GetNumRows();
    SMatrix<T> ls(Nl+Ns);
    Fill(ls,0.0);
    for (auto i:l.rows())
        for (auto j:l.cols(i))
            ls(i,j)=l(i,j);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            ls(Nl+i,Nl+j)=s(i,j);
    return ls;
}
template <class T> SMatrix<T> Orbital_RKB_IBS_Common<T>::merge_off_diag(const Matrix<T>& ls)
{
    size_t Nl=ls.GetNumRows();
    size_t Ns=ls.GetNumCols();
    assert(Nl==Ns);
    SMatrix<T> k(Nl+Ns);
    Fill(k,0.0);
    for (auto i:ls.rows())
        for (auto j:ls.cols())
            k(i,Ns+j)=ls(i,j);
   
    return k;
}    
template <class T> SMatrix<T> Orbital_RKB_IBS_Common<T>::MakeOverlap() const
{
    SMatrix<T> ol=itsRKBL->Overlap();
    SMatrix<T> os=itsRKBS->Kinetic();
    return merge_diag(ol,os);
}
template <class T> SMatrix<T> Orbital_RKB_IBS_Common<T>::MakeKinetic() const
{
    Matrix<T> kls=-itsRKBL->Kinetic(itsRKBS);
    return merge_off_diag(kls);
}
template <class T> SMatrix<T> Orbital_RKB_IBS_Common<T>::MakeNuclear(const Cluster* c) const
{
    SMatrix<T> nl=itsRKBL->Nuclear(c);
    SMatrix<T> ns=itsRKBS->Nuclear(c);
    return merge_diag(nl,ns);
}
template <class T> SMatrix<T> Orbital_RKB_IBS_Common<T>::MakeRestMass() const
{
    SMatrix<T> rl(itsRKBL->GetNumFunctions());
    Fill(rl,0.0);
    SMatrix<T> rs=itsRKBS->Kinetic();
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