// File: BasisSet/Imp/Orbital_DHF_IBS.C
module;
#include <memory>
#include <cassert>
module qchem.BasisSet.Orbital_DHF_IBS;
import qchem.BasisSet.BasisCache;

namespace qchem::BasisSet
{

    template <class T> const smat_t<T>& Orbital_RKB_IBS<T>::RestMass() const
{
    return this->Cache().template Get<T>(BasisCache::Own::RestMass, this->CacheIndex(), [this]{ return MakeRestMass(); });
}

template class Orbital_RKB_IBS<double>;

}