// File: TBasisSetImplementation.C


#include "Imp/BasisSet/HF_IBS_Common.H"
#include "Imp/Containers/ERI4.H"


template <class T> typename Orbital_HF_IBS_Common<T>::SMat Orbital_HF_IBS_Common<T>::
Direct(const SMat& Dcd, const obs_t* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const TOrbital_HF_IBS<T>* ab=this;
    
    if (ab->GetID()<=cd->GetID())
    {
         ERI4 Jabcd=ab->Direct(*cd);
        return Jabcd*Dcd;
    }
    else
    {
        //ERI4 Jcdab=GetDataBase()->GetDirect__4C(*cd,*ab);
        const TOrbital_HF_IBS<T>* cdhf=dynamic_cast<const TOrbital_HF_IBS<T>*>(cd);
        assert(cdhf);
        ERI4 Jcdab=cdhf->Direct(*ab);
        return Dcd*Jcdab;        
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
    {
        ERI4 Kabcd=ab->Exchange(*cd);
        return Kabcd*Dcd;
    }
    else
    {
        //ERI4 Kcdab=GetDataBase()->GetExchange4C(*cd,*ab);
        const TOrbital_HF_IBS<T>* cdhf=dynamic_cast<const TOrbital_HF_IBS<T>*>(cd);
        assert(cdhf);
        ERI4 Kcdab=cdhf->Exchange(*ab);
        return Dcd*Kcdab;        
    }

}

template class Orbital_HF_IBS_Common<double>;

