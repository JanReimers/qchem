// File: BasisSet/Imp/Orbital_DHF_IBS.C
module;
#include <cassert>
module qchem.BasisSet.Orbital_DHF_IBS;
import qchem.BasisSet.Internal.DB_Cache;

namespace BasisSet
{

    template <class T> const smat_t<T>& Orbital_RKB_IBS<T>::RestMass() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::RestMass,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),
        [this]{ return MakeRestMass(); });
}

template class Orbital_RKB_IBS<double>;

}