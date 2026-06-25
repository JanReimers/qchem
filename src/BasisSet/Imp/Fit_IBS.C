// File: BasisSet/Imp/Fit_IBS.C  Implement the numerical (mesh-quadrature) parts of a fit basis set.
//
// The fit quadrature now runs over a qcMesh1::Mesh via the qcMesh1 free-function quadrature.  Fit_IBS
// is already a pointwise VectorFunction (operator()(r), Gradient(r)); we expose it to qcMesh1 through
// tiny view adapters rather than changing its class hierarchy.
module;
#include <cassert>
#include <sstream>
#include <string>
module qchem.BasisSet.Fit_IBS;
import qchem.Mesh1.Quadrature;        // qcMesh1::Mesh, BasisField, ScalarField, Normalize, Overlap
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Blaze;

namespace BasisSet
{

namespace
{
// View a Fit_IBS (a pointwise VectorFunction) as a qcMesh1::BasisField for the quadrature.
class FitBasisView : public qcMesh1::BasisField<double>
{
    const Fit_IBS& its;
public:
    explicit FitBasisView(const Fit_IBS& b) : its(b) {}
    size_t     size()                       const override {return its.GetVectorSize();}
    rvec_t     operator()(const rvec3_t& r) const override {return its(r);}
    rvec3vec_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

// View an old ScalarFunction<double> as a qcMesh1::ScalarField for the quadrature.
class ScalarFnView : public qcMesh1::ScalarField<double>
{
    const ScalarFunction<double>& its;
public:
    explicit ScalarFnView(const ScalarFunction<double>& f) : its(f) {}
    double  operator()(const rvec3_t& r) const override {return its(r);}
    rvec3_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

// Cache key for a mesh (the geometry-free Mesh has no ID; mirror the old Mesh::ID summary).
std::string MeshKey(const qcMesh1::Mesh* m)
{
    std::ostringstream os;
    os << "N=" << m->size();
    if (m->size()>0) os << " " << m->Weights()[0] << " " << m->Points()[0];
    if (m->size()>1) { size_t l=m->size()-1; os << " " << m->Weights()[l] << " " << m->Points()[l]; }
    return os.str();
}
} //anon

const  rvec_t& Fit_IBS::Charge   () const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I1C::Charge,this,
        [this]{ return MakeCharge(); });
}
const rsmat_t& Fit_IBS::Repulsion() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::Repulsion,this,
        [this]{ return MakeRepulsion(); });
}
const  rmat_t& Fit_IBS::Repulsion(const Fit_IBS& b) const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2x::Repulsion,this,&b
            ,[this,&b]{ return MakeRepulsion(b); });
}
const rsmat_t& Fit_IBS::InvOverlap() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::InvOverlap,this,
        [this]{ return MakeInvOverlap(); });
}
const rsmat_t& Fit_IBS::InvRepulsion() const
{
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2C::InvRepulsion,this,
        [this]{ return MakeInvRepulsion(); });
}

const rvec_t& Fit_IBS::Norm(const qcMesh1::Mesh* m) const
{
    assert(m);
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I1C::Normalization,this,MeshKey(m),
        [this,m]{ return MakeNorm(m); });
}

const rmat_t& Fit_IBS::Overlap(const qcMesh1::Mesh* m,const Fit_IBS& b) const
{
    assert(m);
    auto cache=theGlobalCache;
    assert(cache);
    return cache->Get(IntegralsCache_Base::I2x::Overlap,this,&b,MeshKey(m),
        [this,m,&b]{ return MakeOverlap(m,b); });
}

rvec_t Fit_IBS::MakeNorm(const qcMesh1::Mesh* m) const
{
    return qcMesh1::Normalize(*m, FitBasisView(*this));
}

rmat_t Fit_IBS::MakeOverlap(const qcMesh1::Mesh* m,const Fit_IBS& b) const
{
    return qcMesh1::Overlap(*m, FitBasisView(*this), FitBasisView(b));
}

rvec_t Fit_IBS::Overlap(const qcMesh1::Mesh* m,const Sf& f) const
{
    const rvec_t& n=Norm(m);
    // <f_a|f> projection, normalised: component-wise multiply by the fit-basis norms.
    return qcMesh1::Overlap(*m, FitBasisView(*this), ScalarFnView(f)) * n;
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
