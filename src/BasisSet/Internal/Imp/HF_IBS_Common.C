// File: Imp/HF_IBS_Common.H  Common implementation for all Hartree-Fock (HF) Irrep Basis Sets.
module;
#include <cassert>
#include <vector>
#include <memory>
#include <iomanip>
module qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.ERI4;

template <class T> SMatrix<T> Orbital_HF_IBS_Common<T>::
Direct(const SMatrix<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const Orbital_HF_IBS<T>* ab=this;
    
    if (ab->GetID()<=cd->GetID())
    {
         ERI4 Jabcd=ab->Direct(*cd);
        return MatMul(Jabcd,Dcd);
    }
    else
    {
        ERI4 Jcdab=cd->Direct(*ab);
        return MatMul(Dcd,Jcdab);        
    }
}

#include <iomanip>
template <class T> SMatrix<T> Orbital_HF_IBS_Common<T>::
Exchange(const SMatrix<T>& Dcd, const Orbital_HF_IBS<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const Orbital_HF_IBS<T>* ab=this;

    if (ab->GetID()<=cd->GetID())
        return MatMul(ab->Exchange(*cd),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
    else
        return MatMul(Dcd,cd->Exchange(*ab)); // ERI4 Kcdab=cd->Exchange(*ab);    

}

template class Orbital_HF_IBS_Common<double>;

