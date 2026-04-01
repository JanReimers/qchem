// File: Imp/HF_IBS_Common.H  Common implementation for all Hartree-Fock (HF) Irrep Basis Sets.
module;
#include <cassert>
#include <blaze/math/SymmetricMatrix.h>
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.ERI4;

template <class T> void Orbital_HF_IBS_Common<T>::
AccumulateDirect(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(max(abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    MatMul(Sab,ab->Direct(*cd),Dcd); //ERI4 Jabcd=ab->Direct(*cd);
}

template <class T> void Orbital_HF_IBS_Common<T>::
AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(max(abs(Dcd))>0.0);  //Dcd should be pre-screened for zero.
    const Orbital_HF_IBS<T>* ab=this;
    MatMul(Sab,ab->Exchange(*cd),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
}

template class Orbital_HF_IBS_Common<double>;

