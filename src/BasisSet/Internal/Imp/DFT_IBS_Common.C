// File: Imp/DFT_IBS_Common.C  Common implementation for all DFT Irrep Basis Sets.
module;
#include <vector>
#include <blaze/Math.h>
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.Conversions;
import oml;

template <class T> vec_t<T> Orbital_DFT_IBS_Common<T>::
Overlap3C(const smat_t<T>& Dcd, const Fit_IBS* ff) const
{
    vec_t<T> ret(ff->GetNumFunctions());
    auto& S=this->Overlap3C(*ff);
    for(auto i:iv_t(0,S.size()))
        ret[i]=sum(Dcd%S[i]);
    return ret;
}

template <class T> vec_t<T> Orbital_DFT_IBS_Common<T>::
Repulsion3C(const smat_t<T>& Dcd, const Fit_IBS* ff) const
{
    vec_t<T> ret(ff->GetNumFunctions());
    auto& R=this->Repulsion3C(*ff);
    for(auto i:iv_t(0,R.size()))
        ret[i]=sum(Dcd%R[i]);
    return ret;
}

template class Orbital_DFT_IBS_Common<double>;

