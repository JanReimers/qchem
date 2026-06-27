// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <cassert>
#include <map>
#include <string>
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
// MakeOverlap().  Key by the geometry-aware BasisSetID() -- NOT the `this` pointer: a basis freed and
// another reallocated at the same address would otherwise collide (e.g. NaF's 113x113 overlap served to
// a later CsI 251x251 basis -> "Matrix sizes do not match").  Same geometry key the double cache uses.
// A proper complex integral cache is pinned for later (see the plane-wave plan).
template <> const hmat_t<dcmplx>& Integrals_Overlap<dcmplx>::Overlap() const
{
    static std::map<std::string,hmat_t<dcmplx>> buf;
    const std::string id=this->BasisSetID();
    auto i=buf.find(id);
    if (i==buf.end()) i=buf.emplace(id,MakeOverlap()).first;
    return i->second;
}

template class Integrals_Overlap<double>;
template class Integrals_Overlap<dcmplx>;

} //namespace