// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
#include <memory>
#include <stdexcept>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;   // the per-block fit factories this whole-set layer delegates to
import qchem.BasisSet.DeltaFit_IBS;      // DeltaFit_IBS -- the delta representation this layer builds itself
import qchem.Symmetry.Factory;           // BlochFactory (the delta basis's Gamma irrep)

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
template <class T> FIT_SF_ABS<T>* tBasisSet<T>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp,
                                                                    VxcFit fit) const
{
    // The REAL (molecular) path has exactly one representation -- the Gaussian auxiliary basis -- so a
    // delta request has nothing to build here: DeltaFit_IBS carries a Bloch irrep and the molecular XC
    // quadrature route is not wired to it.  Loud, not silent, so "it quietly fitted Gaussians" cannot happen.
    if (fit==VxcFit::Delta)
        throw std::logic_error("tBasisSet: VxcFit::Delta has no real (molecular) realization -- there is no "
            "delta fit basis on the double path; the molecular Vxc route fits the Gaussian auxiliary basis.");
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
// The XC quadrature (delta-fit) factory: generic T = the plain path (the Structure's own integration
// mesh, no fold -- molecules / any basis without an imposed-symmetry override).
template <class T> FitQuadrature tBasisSet<T>::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return {std::make_shared<const qcMesh::Mesh>(cl->CreateIntegrationMesh(mp)), {}};
}
// The plane-wave (dcmplx) path delegates THROUGH the orbital basis's factory, exactly as the fit bases
// do: the Orbital_DFT_IBS<dcmplx> block owns the cell + the imposed ops, so IT assembles the (invariant) quadrature.
template <> FitQuadrature tBasisSet<dcmplx>::CreateXCQuadrature(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return FirstPeriodicDFT(*this, [&](const auto* dft){return dft->CreateXCQuadrature(cl,mp);});
}

// The periodic (dcmplx) Vxc fit basis.  THE ONE SITE THAT CHOOSES A REPRESENTATION (2026-08-22):
//   Delta -- the delta basis over this run's XC quadrature.  The mesh work (grid build, group-averaging it
//            invariant under imposed ops, fold + Shubnikov tags) stays where it was, in the basis's own
//            CreateXCQuadrature; what changed is that its answer now comes back ATTACHED to a fit basis
//            instead of through a second factory the Hamiltonian had to call itself.  Gamma Bloch irrep:
//            the mesh is cell-periodic and carries no crystal momentum of its own.
//   else  -- the lineage's own fitted representation, from the first periodic DFT block's factory (an
//            unresolved Auto lands here, the historical pairing).
template <> FIT_SF_ABS<dcmplx>* tBasisSet<dcmplx>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp,
                                                                       VxcFit fit) const
{
    if (fit==VxcFit::Delta)
        return new DeltaFit_IBS(CreateXCQuadrature(cl,mp),
                                Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
    return FirstPeriodicDFT(*this, [&](const auto* dft){return dft->CreateVxcFitBasisSet(cl,mp);});
}


template class tBasisSet<double>;
template class tBasisSet<dcmplx>;

} //namespace