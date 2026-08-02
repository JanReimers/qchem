// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
#include <memory>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Band_FT_IBS;   // the dcmplx (plane-wave) density-fit factory delegate

namespace qchem::BasisSet
{
template <class T> FIT_CD_ABS<T>* tBasisSet<T>::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateCDFitBasisSet(cl,mp);
}
template <class T> FIT_SF_ABS<T>* tBasisSet<T>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto dft=*Iterate<Orbital_DFT_IBS<double>>().begin();
    return dft->CreateVxcFitBasisSet(cl,mp);
}

// The plane-wave (dcmplx) density-fit basis is created THROUGH the orbital basis's own factory, exactly
// as the double path delegates to Orbital_DFT_IBS: iterate to the Band_FT_IBS (the reciprocal-space DFT
// capability, realized by the plane-wave basis) and let it build its auxiliary cFIT_CD_ABS.
template <> FIT_CD_ABS<dcmplx>* tBasisSet<dcmplx>::CreateCDFitBasisSet (const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto bft=*Iterate<Band_FT_IBS>().begin();
    return bft->CreateCDFitBasisSet(cl,mp);
}
// The plane-wave (dcmplx) Vxc fit basis is created THROUGH the orbital basis's own factory, exactly as the
// CD one: iterate to the Band_FT_IBS and let it build its auxiliary cFIT_SF_ABS.
template <> FIT_SF_ABS<dcmplx>* tBasisSet<dcmplx>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto bft=*Iterate<Band_FT_IBS>().begin();
    return bft->CreateVxcFitBasisSet(cl,mp);
}

// The XC quadrature (delta-fit) factory: generic T = the plain path (the Structure's own integration
// mesh, no fold -- molecules / any basis without an imposed-symmetry override).
template <class T> XCQuadrature tBasisSet<T>::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return {std::make_shared<const qcMesh::Mesh>(cl->CreateIntegrationMesh(mp)), {}};
}
// The plane-wave (dcmplx) path delegates THROUGH the orbital basis's factory, exactly as the fit bases
// do: the Band_FT_IBS block owns the cell + the imposed ops, so IT assembles the (invariant) quadrature.
template <> XCQuadrature tBasisSet<dcmplx>::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    auto bft=*Iterate<Band_FT_IBS>().begin();
    return bft->CreateXCQuadrature(cl,mp);
}

template class tBasisSet<double>;
template class tBasisSet<dcmplx>;

} //namespace