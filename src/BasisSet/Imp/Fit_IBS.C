// File: BasisSet/Imp/Fit_IBS.C  Implement the numerical (mesh-quadrature) parts of a fit basis set.
//
// A fit basis OWNS its quadrature mesh (the Becke molecular mesh, built from its Structure in
// SetMesh).  The numerical integrals run over that mesh via the qcMesh free-function quadrature;
// Fit_IBS is already a pointwise VectorFunction, exposed to qcMesh through tiny view adapters.
module;
#include <cassert>
module qchem.BasisSet.Fit_IBS;
import qchem.Mesh.Quadrature;          // qcMesh::Mesh, BasisField, ScalarField, Normalize, Overlap
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace qchem::BasisSet
{

namespace
{
// View a Fit_IBS (a pointwise VectorFunction) as a qcMesh::BasisField for the quadrature.
class FitBasisView : public qcMesh::BasisField<double>
{
    const Fit_IBS& its;
public:
    explicit FitBasisView(const Fit_IBS& b) : its(b) {}
    size_t     size()                       const override {return its.GetVectorSize();}
    rvec_t     operator()(const rvec3_t& r) const override {return its(r);}
    rvec3vec_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

// View an old ScalarFunction<double> as a qcMesh::ScalarField for the quadrature.
class ScalarFnView : public qcMesh::ScalarField<double>
{
    const ScalarFunction<double>& its;
public:
    explicit ScalarFnView(const ScalarFunction<double>& f) : its(f) {}
    double  operator()(const rvec3_t& r) const override {return its(r);}
    rvec3_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

} //anon

void Fit_IBS::SetMesh(const Structure& st, const qcMesh::MeshParams& mp)
{
    itsMesh   = st.CreateIntegrationMesh(mp);   // the geometry's own mesh (was MakeMolecularMesh(st,mp))
    itsMeshID = mp.ID();   // cache key axis for the mesh-dependent Norm() (see Norm below)
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
// We therefore use the mesh-keyed I1C cache variant (Mesh_ID = MeshParams::ID(), stamped in SetMesh).
// Keying on BasisSetID alone silently served, e.g., the HF SAD bootstrap's coarse-seed-mesh Norm to a
// later production DFT run on a finer mesh -> a ~585 ppm energy drift (the analytic Charge/Repulsion/
// Inv* below are mesh-independent and correctly stay keyed on the basis alone).
const rvec_t& Fit_IBS::Norm() const
{
    assert(!itsMeshID.empty());   // SetMesh must run before any numerical integral
    return theCache<double>().Get(IntegralsCache_Base::I1C::Normalization,this,itsMeshID,
        [this]{ return MakeNorm(); });
}

rvec_t Fit_IBS::MakeNorm() const
{
    assert(itsMesh.size()>0);   // SetMesh must run before any numerical integral
    return qcMesh::Normalize(itsMesh, FitBasisView(*this));
}

rvec_t Fit_IBS::Overlap(const Sf& f) const
{
    assert(itsMesh.size()>0);
    // <f_a|f> projection, normalised: component-wise multiply by the fit-basis norms.
    return qcMesh::Overlap(itsMesh, FitBasisView(*this), ScalarFnView(f)) * Norm();
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
