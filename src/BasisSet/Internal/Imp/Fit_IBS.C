// File: BasisSet/Internal/Imp/Fit_IBS.C  Implement a Gaussian-family fit basis's cached metric accessors and a Vxc fit's projection.
//
// Every integral here is ANALYTIC except a Vxc fit's projection of a pointwise field, Overlap(f), which runs over the projection mesh
// that fit owns (built from its Structure at construction).  A Coulomb fit owns no mesh (user 2026-11-11).
module;
#include <cassert>
#include <cmath>
#include <stdexcept>
module qchem.BasisSet.Internal.Fit_IBS;
import qchem.Mesh.Quadrature;          // qcMesh::Mesh, Overlap (over qcMath Vector/ScalarFunction)
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;
import qchem.Reporting;   // report::Scope (FitMesh)

namespace qchem::BasisSet
{

//! A Vxc fit's projection mesh, built under a report scope naming the fit (the builder reports one level in).
static qcMesh::Mesh FitMesh(const Structure& st, const qcMesh::MeshParams& mp)
{
    qchem::report::Scope here("auxiliaryFit");
    return st.CreateIntegrationMesh(mp);
}

Fit_IBS::Fit_IBS() = default;   // a Coulomb fit: no mesh

Fit_IBS::Fit_IBS(const Structure& st, const qcMesh::MeshParams& mp)
    : itsProjectionMesh(FitMesh(st,mp))         // the geometry's own mesh
{
    // Established ONCE, here, instead of re-asserted in every numerical accessor (R2.10).
    assert(itsProjectionMesh->size()>0 && "Fit_IBS: the structure produced an empty quadrature mesh");
}

rvec_t         Fit_IBS::Charge      () const {return itsCharge      .GetOrFill([this]{return MakeCharge      ();});}
const rsmat_t& Fit_IBS::Repulsion   () const {return itsRepulsion   .GetOrFill([this]{return MakeRepulsion   ();});}
const rsmat_t& Fit_IBS::InvOverlap  () const {return itsInvOverlap  .GetOrFill([this]{return MakeInvOverlap  ();});}
const rsmat_t& Fit_IBS::InvRepulsion() const {return itsInvRepulsion.GetOrFill([this]{return MakeInvRepulsion();});}
// Norm() is ANALYTIC (2026-11-11): a property of the functions alone.  It was a mesh quadrature, keyed by the mesh as well after the same
// fit basis built on a coarse seed mesh served its Norm to a later run on a finer one (a 585 ppm drift).
const rvec_t&  Fit_IBS::Norm        () const {return itsNorm        .GetOrFill([this]{return MakeNorm        ();});}

// The cross Coulomb against a fit basis \a b.  With itself (the only call today, ConstrainedFF) it is one of this basis's own
// integrals.  Against ANOTHER fit it is keyed by that partner, which needs the partner's instance token (S3, plan §3b), so it stays
// in the global cache until then.
const rmat_t& Fit_IBS::Repulsion(const rFIT_CD_ABS& b) const
{
    if (&b==static_cast<const rFIT_CD_ABS*>(this))
        return itsSelfRepulsion.GetOrFill([this,&b]{return MakeRepulsion(b);});
    return theCache<double>().Get(IntegralsCache_Base::I2x::Repulsion,this,&b
            ,[this,&b]{ return MakeRepulsion(b); });
}

rvec_t Fit_IBS::OverlapDiagonal() const
{
    throw std::logic_error("Fit_IBS::OverlapDiagonal: a non-orthogonal fit basis ('"+Name()+"') has no diagonal-only metric to give; "
                           "it takes the full S^-1 solve (doc/CleanupCandidates.md R1.0: delete, don't correct)");
}

rvec_t Fit_IBS::MakeNorm() const
{
    const rsmat_t& S=Integrals_Overlap<double>::Overlap();   // analytic, in op(r)'s normalisation
    rvec_t n(S.rows());
    for (size_t a=0; a<n.size(); a++) n[a]=1.0/std::sqrt(S(a,a));
    return n;
}

rvec_t Fit_IBS::Overlap(const Sf& f) const
{
    if (!itsProjectionMesh)
        throw std::logic_error("Fit_IBS::Overlap(f): the Coulomb fit basis '"+Name()+"' is analytic and owns no mesh to project a field on; "
                               "a field that is only a function on space integrates itself on its own points (BasisSet::SelfIntegratingField)");
    // <f_a|f> projection, normalised: component-wise multiply by the fit-basis norms.
    return qcMesh::Overlap(*itsProjectionMesh, *this, f) * Norm();
}

rsmat_t Fit_IBS::MakeInvOverlap  () const
{
    return blazem::inv(MakeOverlap());
}
rsmat_t Fit_IBS::MakeInvRepulsion() const
{
    return blazem::inv(MakeRepulsion());
}

} //namespace
