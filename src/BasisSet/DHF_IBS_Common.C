// File: DHF_IBS_Common.H  Common implementation for all Dirac-Hartree-Fock (HF) Irrep Basis Sets.

#include <vector>
#include <iostream>
#include "DHF_IBS_Common.H"


template <class T> Orbital_RKB_IBS_Common<T>::Orbital_RKB_IBS_Common
(const DB_cache<T>* db, Symmetry* sym,int kappa,::Orbital_RKBL_IBS<T>* l,::Orbital_RKBS_IBS<T>* s)
    : IBS_Common(sym)
    , Orbital_IBS_Common<T>()
    , DB_RKB<T>(db)
    , itsRKBL(l)
    , itsRKBS(s)
{
    assert(itsRKBL);
    assert(itsRKBS);
    s->InsertBasisFunctions(itsRKBL);
}

template <class T> typename Integrals_Base<T>::SMat Orbital_RKB_IBS_Common<T>::merge_diag(const SMat& l,const SMat& s)
{
    size_t Nl=l.GetNumRows();
    size_t Ns=s.GetNumRows();
    SMat ls(Nl+Ns);
    Fill(ls,0.0);
    for (auto i:l.rows())
        for (auto j:l.cols(i))
            ls(i,j)=l(i,j);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            ls(Nl+i,Nl+j)=s(i,j);
    return ls;
}
template <class T> typename Integrals_Base<T>::SMat Orbital_RKB_IBS_Common<T>::merge_off_diag(const Mat& ls)
{
    size_t Nl=ls.GetNumRows();
    size_t Ns=ls.GetNumCols();
    assert(Nl==Ns);
    SMat k(Nl+Ns);
    Fill(k,0.0);
    for (auto i:ls.rows())
        for (auto j:ls.cols())
            k(i,Ns+j)=ls(i,j);
   
    return k;
}    
template <class T> typename Integrals_Base<T>::SMat Orbital_RKB_IBS_Common<T>::MakeOverlap() const
{
    SMat ol=itsRKBL->Overlap();
    SMat os=itsRKBS->Kinetic();
    return merge_diag(ol,os);
}
template <class T> typename Integrals_Base<T>::SMat Orbital_RKB_IBS_Common<T>::MakeKinetic() const
{
    Mat kls=-itsRKBL->Kinetic(itsRKBS);
    return merge_off_diag(kls);
}
template <class T> typename Integrals_Base<T>::SMat Orbital_RKB_IBS_Common<T>::MakeNuclear(const Cluster* c) const
{
    SMat nl=itsRKBL->Nuclear(c);
    SMat ns=itsRKBS->Nuclear(c);
    return merge_diag(nl,ns);
}
template <class T> typename Integrals_Base<T>::SMat Orbital_RKB_IBS_Common<T>::MakeRestMass() const
{
    SMat rl(itsRKBL->size());
    Fill(rl,0.0);
    SMat rs=itsRKBS->Kinetic();
    return merge_diag(rl,rs);
}

template <class T> Orbital_RKBL_IBS_Common<T>::Orbital_RKBL_IBS_Common
(Symmetry* sym,int _kappa)
    : IBS_Common(sym)
    , TIBS_Common<T>()
    , kappa(_kappa)
{
    assert(kappa!=0);
}

template <class T> Orbital_RKBS_IBS_Common<T>::Orbital_RKBS_IBS_Common
(Symmetry* sym,int _kappa)
    : IBS_Common(sym)
    , TIBS_Common<T>()
    , kappa(_kappa)
{
    assert(kappa!=0);
}

template class Orbital_RKB_IBS_Common<double>;
template class Orbital_RKBL_IBS_Common<double>;
template class Orbital_RKBS_IBS_Common<double>;