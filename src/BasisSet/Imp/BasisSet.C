// File: BasisSetImp.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <cassert>
#include <memory>
#include <stdexcept>
module qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;   // the per-block fit factories this whole-set layer delegates to
import qchem.BasisSet.DeltaFit_IBS;      // DeltaFit_IBS -- the delta representation this layer builds itself
import qchem.Symmetry.Factory;           // BlochFactory (the delta basis's Gamma irrep)
import qchem.Reporting;                  // report::Timed -- splitting the delta-fit build (ParallelAndOraclePlan 1.1(b))

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
// THE REAL (MOLECULAR) PATH'S FIRST DFT-CAPABLE BLOCK -- the double sibling of FirstPeriodicDFT above.
//
// ⛔ V1.15 (2026-09-09): THIS USED TO BE `*Iterate<Orbital_DFT_IBS<double>>().begin()`, WHICH IS A NULL
// DEREFERENCE IN RELEASE ON ANY BASIS WITHOUT THE FACE.  `D_IndexIterator::operator*` does the
// dynamic_cast and then `assert(d)` -- compiled out under NDEBUG, i.e. in every production run and every
// benchmark -- so a 1E/HF-only basis got a null back and called through it.  Worse, `begin()` on an EMPTY
// set equals `end()`, so the same expression also indexed block 0 of a basis that has none.  Same ruling
// as R1.0i and RequireSiteBlocks: a composition error no caller can act on THROWS, and it says what the
// basis was asked for.
//
// ⚠ AND THE `<double>` IN THE GENERIC BODY IS DELIBERATE, NOT AN OVERSIGHT.  `tBasisSet<dcmplx>`
// specializes every one of these factories (just below), so the template body IS the real/molecular path
// and there is no T for which asking for the dcmplx face would be right here.  Naming it once, here,
// beats repeating the reasoning at each site.
template <class T> static const Orbital_DFT_IBS<double>* FirstRealDFT(const tBasisSet<T>& bs, const char* who)
{
    for (size_t i=0;i<bs.GetNumIBS();++i)
        if (const auto* dft=dynamic_cast<const Orbital_DFT_IBS<double>*>(bs[i])) return dft;
    throw std::logic_error(std::string(who)+": no block of this basis carries the DFT fit-factory face "
        "(Orbital_DFT_IBS<double>), so it cannot build a fit basis.  A 1E- or HF-only basis reaches here "
        "when something asked it to fit a density or a potential -- that is a composition error, not a "
        "recoverable condition.");
}

template <class T> FIT_CD_ABS<T>* tBasisSet<T>::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    return FirstRealDFT(*this, "tBasisSet::CreateCDFitBasisSet")->CreateCDFitBasisSet(cl,mp);
}
template <class T> FIT_SF_ABS<T>* tBasisSet<T>::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp,
                                                                    VxcFit fit,
                                                                    FitQuadrature* quad) const
{
    // No quadrature on the real path: the Gaussian auxiliary basis keeps its own private quadrature mesh
    // (Fit_IBS::itsMesh, by value), and nothing molecular asks for atom-partitioned observables yet.
    if (quad) *quad=FitQuadrature();
    // The REAL (molecular) path has exactly one representation -- the Gaussian auxiliary basis -- so a
    // delta request has nothing to build here: DeltaFit_IBS carries a Bloch irrep and the molecular XC
    // quadrature route is not wired to it.  Loud, not silent, so "it quietly fitted Gaussians" cannot happen.
    if (fit==VxcFit::Delta)
        throw std::logic_error("tBasisSet: VxcFit::Delta has no real (molecular) realization -- there is no "
            "delta fit basis on the double path; the molecular Vxc route fits the Gaussian auxiliary basis.");
    return FirstRealDFT(*this, "tBasisSet::CreateVxcFitBasisSet")->CreateVxcFitBasisSet(cl,mp);
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
// ⏸ V1.15 asked for this neutral default to be HOISTED, since Orbital_DFT_IBS::CreateXCQuadrature carries a
// byte-identical body.  DECLINED 2026-09-09, and the reason is where the shared home would have to live:
// the only module below BOTH declarers is qchem.BasisSet.Fit_Types, which states in its own header that it
// "needs the mesh and the symmetry fold and nothing else" -- and this body needs `Structure`.  Hoisting two
// lines would drag Structure into a deliberate leaf to save a duplication that cannot drift silently (both
// are `the Structure's own mesh, no fold`, and a change to either shows up as a failing quadrature).
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
                                                                       VxcFit fit,
                                                                       FitQuadrature* quad) const
{
    if (quad) *quad=FitQuadrature();
    if (fit==VxcFit::Delta)
    {
        // ONE quadrature, built ONCE, INJECTED into both collaborators that need it (user, 2026-08-23):
        // the delta basis (for which it is constitutive) and -- optionally -- the XC strategy, which needs
        // the ATOMIC PARTITION the mesh carries (a general-purpose observable, no part of a fit basis's
        // job) and, since 2026-08-24, the ORBIT FOLD and the Shubnikov tags for the rho star-average --
        // which were reaching it as two MEMBERS on the fit face, where they never belonged.  Both hold the
        // same immutable Mesh through the same shared_ptr, so their orderings agree by construction rather
        // than by convention.  A raster (PlaneWave) fit basis has no such quadrature and leaves *quad
        // empty -- no Becke build is paid for on a route that would not use it.
        // ★ 1.1(b): three buckets, because this function is 23.75 s/call on MnO and the mesh build is only
        // a third of it.  Timed is exclusive, so the mesh build / orbit fold / Shubnikov buckets inside
        // CreateXCQuadrature stay children of the first one.
        FitQuadrature q;
        {
            qchem::report::Timed timed("setup: XC quadrature build (mesh + fold are its children)");
            q=CreateXCQuadrature(cl,mp);
        }
        if (quad)
        {
            // A full COPY of the quadrature -- mesh handle (cheap), but also the orbit fold, the σ tags
            // and the flip-fixed flags over every mesh point.  Priced because "hand the same bundle to
            // both collaborators" is an ownership decision, and nobody had ever asked what it costs.
            qchem::report::Timed timed("setup: XC quadrature copy (to the term stack)");
            *quad=q;
        }
        qchem::report::Timed timed("setup: DeltaFit_IBS ctor");
        return new DeltaFit_IBS(std::move(q),
                                Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
    }
    return FirstPeriodicDFT(*this, [&](const auto* dft){return dft->CreateVxcFitBasisSet(cl,mp);});
}


template class tBasisSet<double>;
template class tBasisSet<dcmplx>;

} //namespace