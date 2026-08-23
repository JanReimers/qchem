// File: BasisSet/Imp/Fit_IBS.C  Implement the numerical (mesh-quadrature) parts of a fit basis set.
//
// A fit basis OWNS its quadrature mesh (built from its Structure AT CONSTRUCTION).  The numerical integrals run over that mesh via the qcMesh free-function quadrature;
// Fit_IBS is-a pointwise VectorFunction, so it is passed straight to the quadrature (no adapter).
module;
#include <cassert>
module qchem.BasisSet.Fit_IBS;
import qchem.Mesh.Quadrature;          // qcMesh::Mesh, Normalize, Overlap (over qcMath Vector/ScalarFunction)
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace qchem::BasisSet
{

Fit_IBS::Fit_IBS(const Structure& st, const qcMesh::MeshParams& mp)
    : itsMesh  (st.CreateIntegrationMesh(mp))   // the geometry's own mesh
    , itsMeshID(mp.ID())                        // cache key axis for the mesh-dependent Norm() (see Norm below)
{
    // Established ONCE, here, instead of re-asserted in every numerical accessor (R2.10).
    assert(itsMesh.size()>0   && "Fit_IBS: the structure produced an empty quadrature mesh");
    assert(!itsMeshID.empty() && "Fit_IBS: MeshParams::ID() is the Norm cache key -- it must not be empty");
}

const  rvec_t& Fit_IBS::Charge   () const
{
    return theCache<double>().Get(IntegralsCache_Base::I1C::Charge,this,
        [this]{ return MakeCharge(); });
}
const rsmat_t& Fit_IBS::Repulsion() const
{
    return theCache<double>().Get(IntegralsCache_Base::I2C::Repulsion,this,
        [this]{ return MakeRepulsion(); });
}
const  rmat_t& Fit_IBS::Repulsion(const rFIT_CD_ABS& b) const
{
    return theCache<double>().Get(IntegralsCache_Base::I2x::Repulsion,this,&b
            ,[this,&b]{ return MakeRepulsion(b); });
}
const rsmat_t& Fit_IBS::InvOverlap() const
{
    return theCache<double>().Get(IntegralsCache_Base::I2C::InvOverlap,this,
        [this]{ return MakeInvOverlap(); });
}
const rsmat_t& Fit_IBS::InvRepulsion() const
{
    return theCache<double>().Get(IntegralsCache_Base::I2C::InvRepulsion,this,
        [this]{ return MakeInvRepulsion(); });
}

// Norm() is a MESH QUADRATURE (qcMesh::Normalize over itsMesh), so it MUST be keyed by the mesh as well
// as the basis: the same fit basis (same BasisSetID) built with a different mesh has a different Norm.
// We therefore use the mesh-keyed I1C cache variant (Mesh_ID = MeshParams::ID(), stamped in the ctor).
// Keying on BasisSetID alone silently served, e.g., the HF SAD bootstrap's coarse-seed-mesh Norm to a
// later production DFT run on a finer mesh -> a ~585 ppm energy drift (the analytic Charge/Repulsion/
// Inv* below are mesh-independent and correctly stay keyed on the basis alone).
const rvec_t& Fit_IBS::Norm() const
{
    return theCache<double>().Get(IntegralsCache_Base::I1C::Normalization,this,itsMeshID,
        [this]{ return MakeNorm(); });
}

// The expansion over this basis, evaluated: the two bodies moved here VERBATIM from the AO fitter, which
// used to reach in for op(r)/Gradient and sum them itself (Fitting/Internal/Imp/FunctionFitterImp.C).
// Same arithmetic, same order -- what changed is who owns it.
double Fit_IBS::EvalField(const rvec_t& c, const rvec3_t& r) const
{
    return blazem::trans(c) * (*this)(r);
}

rvec3_t Fit_IBS::EvalFieldGradient(const rvec_t& c, const rvec3_t& r) const
{
    vec_t<rvec3_t> br = this->Gradient(r);
    rvec3_t ret(0,0,0);
    auto ci(c.begin());
    auto b(br.begin());
    for (; b!=br.end() && ci!=c.end(); b++,ci++) ret+=(*ci) * (*b);
    return ret;
}

// <f_a|1> (FIT_SF_ABS::Integrals) IS <f_a|1> (FIT_CD_NonOrtho::Charge) -- one quantity, two faces of the
// one object, so the SF sibling forwards instead of re-deriving it on the mesh.  Both are in the
// NORMALISED convention (the per-function norm folded in), which is also the convention Overlap(f) and the
// InvOverlap() metric use, so c.Integrals() is the integral of the fitted field.
rvec_t Fit_IBS::Integrals() const {return Charge();}

// <f_a|f_a> = 1/Norm_a^2 -- the cached normalisation, un-inverted.  A non-orthogonal basis's fit does not
// use it (it solves with the full S^-1); it is here because the metric DIAGONAL is a question every scalar
// fit basis can answer, and answering it is what makes the orthogonal-vs-orthonormal distinction usable.
//
// WARNING, FOUND 2026-08-23: this is the ONLY member of the fit face in the UN-normalised convention.
// Overlap(f) below multiplies by Norm(), Charge()/Integrals() fold the norm in, and MakeOverlap()'s metric
// has a unit diagonal -- so a general orthogonal fitter pairing THIS diagonal with THIS projection would be
// wrong by Norm_a^2.  Nothing does: isOrtho()==false sends every Gaussian fit to the S^-1 solve, which
// never reads the diagonal.
//
// AND THE FIX IS NOT TO CORRECT IT (user ruling, 2026-08-23) -- it is to DELETE it, by moving
// OverlapDiagonal off the metric-NEUTRAL face onto a FIT_SF_Ortho refinement, so a non-orthogonal basis is
// never asked a question it has no honest answer to.  A corrected value here would be all ones, which is
// indistinguishable from the plane-wave answer and would quietly retire the orthogonal-vs-orthoNORMAL
// distinction the periodic gate exists to keep load-bearing.  Specced as the next increment (both fit
// faces together) in doc/CleanupCandidates.md R1.0; until it lands, THIS COMMENT is the only guard.
rvec_t Fit_IBS::OverlapDiagonal() const
{
    const rvec_t& nrm=Norm();
    rvec_t d(nrm.size());
    for (size_t a=0; a<nrm.size(); a++) d[a]=1.0/(nrm[a]*nrm[a]);
    return d;
}

rvec_t Fit_IBS::MakeNorm() const
{
    return qcMesh::Normalize(itsMesh, *this);
}

rvec_t Fit_IBS::Overlap(const Sf& f) const
{
    assert(itsMesh.size()>0);
    // <f_a|f> projection, normalised: component-wise multiply by the fit-basis norms.
    return qcMesh::Overlap(itsMesh, *this, f) * Norm();
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
