// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;

namespace BasisSet
{
template <class T> FIT_CD_ABS* BasisSet<T>::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateCDFitBasisSet(cl,mp);
}
template <class T> FIT_SF_ABS* BasisSet<T>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateVxcFitBasisSet(cl,mp);
}

// Auxiliary DFT-fit basis sets are a real Gaussian-basis facility; the plane-wave (dcmplx) lineage
// does not fit densities (it expands them in the plane-wave basis directly), so these are NA there.
template <> FIT_CD_ABS* BasisSet<dcmplx>::CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const { assert(false); return 0; }
template <> FIT_SF_ABS* BasisSet<dcmplx>::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const { assert(false); return 0; }

template class BasisSet<double>;
template class BasisSet<dcmplx>;

} //namespace