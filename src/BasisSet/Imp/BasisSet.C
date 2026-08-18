// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
#include <memory>
#include <stdexcept>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;   // the dcmplx (plane-wave) density-fit factory delegate

namespace qchem::BasisSet
{

// MIXED-AWARE first-block walk (doc/RealComplexPlan.md 3c-3): the run's fit basis / XC quadrature is a
// RUN property (V1.1) every block shares, so the FIRST block able to serve does.  After the factory
// flip a real TRIM block may sit at index 0 (Γ-first), where the same-scalar Iterate view throws -- so
// probe the cross-scalar view first: a real block IS-A Orbital_DFT_IBS<double,dcmplx>, whose FIT side
// is already the run's complex G-space basis, so both alternatives serve the same complex factories.
template <class F> static auto FirstPeriodicDFT(const tBasisSet<dcmplx>& bs, F&& serve)
{
    for (size_t i=0;i<bs.GetNumIBS();++i)
        if (const auto* rb=bs.GetRealIBS(i))
        {
            if (const auto* dft=dynamic_cast<const Orbital_DFT_IBS<double,dcmplx>*>(rb)) return serve(dft);
        }
        else if (const auto* dft=dynamic_cast<const Orbital_DFT_IBS<dcmplx>*>(bs[i])) return serve(dft);
    throw std::logic_error("tBasisSet<dcmplx>: no block carries the periodic DFT (fit-factory) face");
}
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
// as the double path delegates to Orbital_DFT_IBS: iterate to the Orbital_DFT_IBS<dcmplx> (the reciprocal-space DFT
// capability, realized by the plane-wave basis) and let it build its auxiliary cFIT_CD_ABS.
template <> FIT_CD_ABS<dcmplx>* tBasisSet<dcmplx>::CreateCDFitBasisSet (const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return FirstPeriodicDFT(*this, [&](const auto* dft){return dft->CreateCDFitBasisSet(cl,mp);});
}
// The plane-wave (dcmplx) Vxc fit basis is created THROUGH the orbital basis's own factory, exactly as the
// CD one: iterate to the Orbital_DFT_IBS<dcmplx> and let it build its auxiliary cFIT_SF_ABS.
template <> FIT_SF_ABS<dcmplx>* tBasisSet<dcmplx>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return FirstPeriodicDFT(*this, [&](const auto* dft){return dft->CreateVxcFitBasisSet(cl,mp);});
}

// The XC quadrature (delta-fit) factory: generic T = the plain path (the Structure's own integration
// mesh, no fold -- molecules / any basis without an imposed-symmetry override).
template <class T> XCQuadrature tBasisSet<T>::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return {std::make_shared<const qcMesh::Mesh>(cl->CreateIntegrationMesh(mp)), {}};
}
// The plane-wave (dcmplx) path delegates THROUGH the orbital basis's factory, exactly as the fit bases
// do: the Orbital_DFT_IBS<dcmplx> block owns the cell + the imposed ops, so IT assembles the (invariant) quadrature.
template <> XCQuadrature tBasisSet<dcmplx>::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return FirstPeriodicDFT(*this, [&](const auto* dft){return dft->CreateXCQuadrature(cl,mp);});
}

template class tBasisSet<double>;
template class tBasisSet<dcmplx>;

} //namespace