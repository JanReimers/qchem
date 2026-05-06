// File: BasisSet1/Imp/Orbital_DHF_IBS.C
module;
#include <cassert>
module qchem.BasisSet1.Internal.Orbital_DHF_IBS;
import qchem.BasisSet1.DB_Cache;

namespace BasisSet1
{

template <class T> const mat_t<T>& Orbital_RKBL_IBS<T>::Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Has(IntegralsCache_Base::I2x::Kinetic,
            IntegralsCache_Base::IBS_ID_t(      RadialID(),      AngularID()),
            IntegralsCache_Base::IBS_ID_t(rkbs->RadialID(),rkbs->AngularID())
        )
        ? cache->GetMat() : cache->Set(MakeKinetic(rkbs));
    
}

template class Orbital_RKBL_IBS<double>;
template class Orbital_RKBS_IBS<double>;

}