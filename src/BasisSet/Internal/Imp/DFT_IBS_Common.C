// File: Imp/DFT_IBS_Common.C  Common implementation for all DFT Irrep Basis Sets.
module;
#include <vector>
module qchem.BasisSet.Internal.IBS_Common;
import oml;

template <class T> Vector<T> Orbital_DFT_IBS_Common<T>::
Overlap3C(const SMatrix<T>& Dcd, const Fit_IBS* ff) const
{
    Vector<T> ret(ff->size());
    auto& S=this->Overlap3C(*ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,S[i-1]);
    return ret;
}

template <class T> Vector<T> Orbital_DFT_IBS_Common<T>::
Repulsion3C(const SMatrix<T>& Dcd, const Fit_IBS* ff) const
{
    Vector<T> ret(ff->size());
    auto& repulsion=this->Repulsion3C(*ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,repulsion[i-1]);
    return ret;
}

template class Orbital_DFT_IBS_Common<double>;

