// File: BasisSet/Internal/Imp/CacheSlot.C  A block's slot in its parent's cache, and the cached overlap.
module;
#include <cstddef>
module qchem.BasisSet.Internal.CacheSlot;

namespace qchem::BasisSet
{

bool CacheSlot::SharesCacheWith(const CacheSlot& other) const
{
    return &other.Cache()==&Cache();
}

template <class T> const hmat_t<T>& CachedOverlap<T>::Overlap() const
{
    return this->Cache().template Get<T>(BasisCache::Own::Overlap, this->CacheIndex(), [this]{ return this->MakeOverlap(); });
}

template class CachedOverlap<double>;
template class CachedOverlap<dcmplx>;

} //namespace
