// File: Imp/HF_IBS_Common.H  Common implementation for all Hartree-Fock (HF) Irrep Basis Sets.
module;
#include <cassert>
#include <blaze/math/SymmetricMatrix.h>
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.ERI4;

template <class T> smat_t<T> Orbital_HF_IBS_Common<T>::
Direct(const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(max(abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    return MatMul(ab->Direct(*cd),Dcd); //ERI4 Jabcd=ab->Direct(*cd);
}

template <class T> smat_t<T> Orbital_HF_IBS_Common<T>::
Exchange(const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(max(abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    return MatMul(ab->Exchange(*cd),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
}

template class Orbital_HF_IBS_Common<double>;

