// File: BasisSet1/Imp/Orbital_DHF_IBS.C
module;
#include <cassert>
module qchem.BasisSet1.Orbital_DHF_IBS;
import qchem.BasisSet1.DB_Cache;

namespace BasisSet1
{

    template <class T> const smat_t<T>& Orbital_RKB_IBS<T>::RestMass() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Has(IntegralsCache_Base::I2C::RestMass,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
        ? cache->GetSMat() : cache->Set(MakeRestMass());
    
}

template class Orbital_RKB_IBS<double>;

}