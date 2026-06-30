// File: BasisSet/IrrepBasisSet.C  Interface for an Irrep Basis Set (IBS)
module;
module qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.Internal.DB_Cache;

namespace qchem::BasisSet
{

// Hermitian-correct for both lineages now: hmat_t<double> is symmetric (rsmat_t), hmat_t<dcmplx> is
// chmat_t.  The cache keys by the geometry-aware BasisSetID() (NOT the `this` pointer -- a freed basis
// and another reallocated at the same address would otherwise collide, e.g. NaF's 113x113 overlap served
// to a later CsI 251x251 basis), so the one template serves the real and complex paths alike.
template <class T> const hmat_t<T>& Integrals_Overlap<T>::Overlap() const
{
    return theCache<T>().Get(IntegralsCache_Base::I2C::Overlap,this,
        [this]{ return MakeOverlap(); });
}

template class Integrals_Overlap<double>;
template class Integrals_Overlap<dcmplx>;

} //namespace