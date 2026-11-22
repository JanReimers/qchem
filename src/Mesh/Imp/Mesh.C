// File: Mesh/Imp/Mesh.C  The qcMesh::Mesh value type's bodies and the free helpers declared in Mesh.C (D-BUILD-BLAST, lane B
// 2026-11-18: an edit here recompiles this TU, not the ~490 importers of qchem.Mesh).
module;
#include <utility>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
module qchem.Mesh;

namespace qchem::qcMesh
{

Mesh::Mesh(rvec3vec_t r, rvec_t w) : itsR(std::move(r)), itsW(std::move(w)) {}

Mesh::Mesh(rvec3vec_t r, rvec_t w, std::vector<size_t> siteStart, rvec3vec_t ownerOffset)
    : itsR(std::move(r)), itsW(std::move(w)), itsSiteStart(std::move(siteStart)), itsU(std::move(ownerOffset))
{
    if (itsU.size()!=0 && itsU.size()!=itsW.size())
        throw std::invalid_argument("qcMesh::Mesh: the owner offsets are "+std::to_string(itsU.size())+" for "+std::to_string(itsW.size())+" points");
}

const rvec3vec_t& Mesh::Points() const {return itsR;}

std::unique_ptr<FourierSampler> Mesh::MakeFourierSampler(const FourierBox& b) const
{
    return itsSampling ? itsSampling->Make(itsR, b) : PointSampling::MakeFourierSampler(b);
}

void Mesh::InstallFourierSampling(std::shared_ptr<const FourierSamplingStrategy> s) {itsSampling=std::move(s);}
void Mesh::InstallWeightResponse(std::shared_ptr<const MeshWeightResponse> r) {itsWeightResponse=std::move(r);}
void Mesh::InstallForces(Symmetry::GroupPtr g)
{
    if (!g) throw std::logic_error("qcMesh::Mesh::InstallForces: a null group (a plain point cloud installs nothing: it forces {e})");
    itsForces=std::move(g);
}
Symmetry::GroupPtr Mesh::Forces() const {return itsForces ? itsForces : Symmetry::TrivialGroup();}

Outcome<Symmetry::GroupPtr,std::string> Forces(const MeshParams& mp, const Symmetry::Frame& frame)
{
    using R=Outcome<Symmetry::GroupPtr,std::string>;
    if (mp.cellKind==UnitCellKind::Uniform)
    {
        if (!frame.lattice) return R::Fail("a uniform cell mesh imposes H_R = T, and there is no lattice here to build it on -- To fix this, please state the XC "
                                           "mesh cellKind Becke or Auto, the atom-centred mesh");
        return R::Ok(Symmetry::TranslationGroup(*frame.lattice, Symmetry::kExactTolerance));
    }
    if (mp.angular==AngularKind::None)
    {
        if (!frame.centre) return R::Fail("the radial-only point set (angular none: one direction at 4 pi, exact only for a spherical integrand) "
                                          "imposes H_R = O(3) about one centre, and there is no such centre here -- To fix this, please state an angular rule "
                                          "(Lebedev or GaussLegendre)");
        return R::Ok(Symmetry::PointRotationGroup(*frame.centre, Symmetry::kExactTolerance));
    }
    return R::Ok(Symmetry::TrivialGroup());
}

SamplingTarget* Mesh::FourierTarget() const {return itsTarget.get();}

rvec_t SiteIntegrals(const Mesh& m, const rvec_t& f)
{
    assert(f.size()==m.size() && "SiteIntegrals: the field must be sampled on this mesh's points");
    assert(m.NSites()>0 && "SiteIntegrals: this mesh carries no site blocks (not an atom-centred build)");
    const rvec_t& w=m.Weights();
    rvec_t out(m.NSites(), 0.0);
    for (size_t a=0; a<m.NSites(); a++)
    {
        double s=0.0;
        for (size_t q=m.SiteBegin(a), e=m.SiteEnd(a); q<e; q++) s+=w[q]*f[q];
        out[a]=s;
    }
    return out;
}

double Integrate(const Mesh& m, const rvec_t& f)
{
    assert(f.size()==m.size() && "qcMesh::Integrate: the field must be sampled on this mesh's points");
    const rvec_t& w=m.Weights();
    double s=0.0;
    for (size_t g=0; g<f.size(); g++) s+=w[g]*f[g];
    return s;
}

// Append grows the SoA arrays by one.  blaze resize copies, so this is a setup-path convenience
// (the incremental Becke builder), NOT for hot loops -- when the size is known, use the
// from-arrays constructor instead.
void Mesh::Append(const rvec3_t& r, double w)
{
    itsSampling.reset();   // a point added: no longer the structure a builder recorded
    itsWeightResponse.reset();
    itsForces.reset();     // ... nor the symmetry it stamped
    itsU.resize(0);        // ... nor the owners it recorded
    size_t n=itsW.size();
    assert(itsR.size()==n);
    itsR.resize(n+1,true);
    itsW.resize(n+1,true);
    itsR[n]=r;
    itsW[n]=w;
}

void Mesh::ShiftOrigin(const rvec3_t& o)
{
    itsSampling.reset();   // the raster it described has moved
    itsWeightResponse.reset();
    itsForces.reset();     // a centre or lattice that moved
    for (size_t i=0; i<itsR.size(); i++) itsR[i]+=o;
}

int UniformDivisions(double cellEdge, double eCut)
{
    assert(cellEdge>0.0 && "UniformDivisions: a cell edge is a positive length");
    assert(eCut>0.0 && "UniformDivisions: a cutoff of 0 does not select a resolution (use nUniform)");
    return max(1, int(ceil(2.0*cellEdge*sqrt(2.0*eCut)/Pi)));
}

double UniformCutoff(double cellEdge, int n)
{
    assert(cellEdge>0.0 && "UniformCutoff: a cell edge is a positive length");
    assert(n>0 && "UniformCutoff: a mesh has at least one division per axis");
    const double G=n*Pi/(2.0*cellEdge);   // invert n = 2 a sqrt(2E)/pi for sqrt(2E)
    return 0.5*G*G;
}

} //namespace qchem::qcMesh
