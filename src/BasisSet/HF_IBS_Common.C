// File: TBasisSetImplementation.C

#include <cassert>
#include <vector>
#include <memory>
#include <iomanip>
#include "HF_IBS_Common.H"
#include "ERI4.H"


template <class T> typename Orbital_HF_IBS_Common<T>::SMat Orbital_HF_IBS_Common<T>::
Direct(const SMat& Dcd, const obs_t* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const TOrbital_HF_IBS<T>* ab=this;
    
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
template <class T> typename Orbital_HF_IBS_Common<T>::SMat Orbital_HF_IBS_Common<T>::
Exchange(const SMat& Dcd, const obs_t* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const TOrbital_HF_IBS<T>* ab=this;

    if (ab->GetID()<=cd->GetID())
        return MatMul(ab->Exchange(*cd),Dcd); // ERI4 Kabcd=ab->Exchange(*cd);
    else
        return MatMul(Dcd,cd->Exchange(*ab)); // ERI4 Kcdab=cd->Exchange(*ab);    

}

template class Orbital_HF_IBS_Common<double>;

