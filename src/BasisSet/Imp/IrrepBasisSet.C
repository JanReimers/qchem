// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <cassert>
module qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.Internal.DB_Cache;

namespace BasisSet
{

template <class T> const smat_t<T>& Integrals_Overlap<T>::Overlap() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::Overlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()),
        [this]{ return MakeOverlap(); });
}

template class Integrals_Overlap<double>;

} //namespace