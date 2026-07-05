// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;

namespace qchem::BasisSet
{
template <class T> rFIT_CD_ABS* tBasisSet<T>::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateCDFitBasisSet(cl,mp);
}
template <class T> FIT_SF_ABS* tBasisSet<T>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateVxcFitBasisSet(cl,mp);
}

// Auxiliary DFT-fit basis sets are a real Gaussian-basis facility; the plane-wave (dcmplx) lineage
// does not fit densities (it expands them in the plane-wave basis directly), so these are NA there.
template <> rFIT_CD_ABS* tBasisSet<dcmplx>::CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const { assert(false); return 0; }
template <> FIT_SF_ABS* tBasisSet<dcmplx>::CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const { assert(false); return 0; }

template class tBasisSet<double>;
template class tBasisSet<dcmplx>;

} //namespace