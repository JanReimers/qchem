// File: BasisSet1/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <cassert>
module qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.DB_Cache;

namespace BasisSet
{

template <class T> const smat_t<T>& Integrals_Overlap<T>::Overlap() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Has(IntegralsCache_Base::I2C::Overlap,IntegralsCache_Base::IBS_ID_t(RadialID(),AngularID()))
        ? cache->GetSMat() : cache->Set(MakeOverlap());
}

template class Integrals_Overlap<double>;

} //namespace