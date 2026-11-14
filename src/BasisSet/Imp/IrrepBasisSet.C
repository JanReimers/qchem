// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
#include <memory>
module qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.BasisCache;

namespace qchem::BasisSet
{

// Hermitian-correct for both lineages: hmat_t<double> is symmetric (rsmat_t), hmat_t<dcmplx> is chmat_t.  The table lives in the
// PARENT's cache under this block's index (doc/IntegralCacheOwnershipPlan.md S3a): no string key, so no collision class (the NaF 113x113
// overlap once served to a later CsI basis at the same address cannot recur: the cache dies with its basis).
template <class T> const hmat_t<T>& Integrals_Overlap<T>::Overlap() const
{
    return this->Cache().template Get<T>(BasisCache::Own::Overlap, this->CacheIndex(), [this]{ return MakeOverlap(); });
}

template class Integrals_Overlap<double>;
template class Integrals_Overlap<dcmplx>;

} //namespace