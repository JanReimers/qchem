// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <cassert>
#include <map>
module qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.Internal.DB_Cache;

namespace BasisSet
{

template <class T> const hmat_t<T>& Integrals_Overlap<T>::Overlap() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::Overlap,this,
        [this]{ return MakeOverlap(); });
}

// theGlobalCache is double-only (the integral cache is a real Gaussian-basis facility); the
// plane-wave (dcmplx) lineage has no complex cache yet, so the complex Overlap() lazily buffers
// MakeOverlap() per instance.  Keeps the no-data-in-interface diamond rule (storage is external).
// A proper complex integral cache is pinned for later (see the plane-wave plan).
template <> const hmat_t<dcmplx>& Integrals_Overlap<dcmplx>::Overlap() const
{
    static std::map<const Integrals_Overlap<dcmplx>*,hmat_t<dcmplx>> buf;
    auto i=buf.find(this);
    if (i==buf.end()) i=buf.emplace(this,MakeOverlap()).first;
    return i->second;
}

template class Integrals_Overlap<double>;
template class Integrals_Overlap<dcmplx>;

} //namespace