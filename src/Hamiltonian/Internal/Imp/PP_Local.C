// File: Hamiltonian/Internal/Imp/PP_Local.C  Local-pseudopotential term (mesh quadrature of V_loc(r)).
//
// The term owns only the MODEL (V_loc as a real-space field).  It asks the STRUCTURE for its integration
// mesh -- CreateIntegrationMesh is a polymorphic geometry capability (Atom -> radial x angular, Molecule ->
// Becke, lattice -> uniform), so the mesh TYPE follows the geometry, not this term -- and quadratures
// <i|V_loc|j> with the generic qcMesh::WeightedOverlap (the same routine the XC path uses).  The term is
// thus geometry-neutral without any basis-specific dispatch.
module;
#include <cassert>
#include <iostream>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Mesh.Quadrature;           // qcMesh::Mesh, WeightedOverlap, ScalarField, BasisField
import qchem.VectorFunction;            // the orbital basis IS-A VectorFunction (the integrand source)
import qchem.Math;                      // norm(Vector3D)

namespace qchem::Hamiltonian
{

namespace {

// Adapt the orbital basis (a VectorFunction: r -> [chi_i(r)]) to the mesh-quadrature BasisField.
class BFView : public qcMesh::BasisField<double>
{
    const VectorFunction<double>& its;
public:
    explicit BFView(const VectorFunction<double>& v) : its(v) {}
    size_t     size()                       const override {return its.GetVectorSize();}
    rvec_t     operator()(const rvec3_t& r) const override {return its(r);}
    rvec3vec_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

// The real-space local pseudopotential as a scalar field: V_loc(r) = Sum_atoms V_loc(Z_a, |r - R_a|).
class VlocField : public qcMesh::ScalarField<double>
{
    const Structure& cl;
    const Pseudopotential::LocalPotential_R& v;
public:
    VlocField(const Structure& c, const Pseudopotential::LocalPotential_R& vl) : cl(c), v(vl) {}
    double operator()(const rvec3_t& r) const override
    {
        double s=0;
        for (size_t i=0;i<cl.GetNumAtoms();i++) { const Atom* a=cl[i]; s+=v.Vloc(a->itsZ, norm(r-a->itsR)); }
        return s;
    }
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}   // unused by WeightedOverlap
};

} //anon

PP_Local::PP_Local(const st_t& st, vloc_t vloc, const qcMesh::MeshParams& mp)
    : theStructure(st), itsVloc(std::move(vloc)), itsMeshParams(mp)
{
    assert(theStructure);
    assert(itsVloc);
}

rsmat_t PP_Local::CalculateMatrix(const robs_t* bs, const Spin&) const
{
    qcMesh::Mesh mesh = theStructure->CreateIntegrationMesh(itsMeshParams);   // the geometry's own mesh
    return qcMesh::WeightedOverlap(mesh, BFView(*bs), VlocField(*theStructure, *itsVloc));
}

void PP_Local::GetEnergy(EnergyBreakdown& te, const rDM_CD* cd) const
{
    te.Een += cd->DM_Contract(this);   // electron-ion (local PP) energy = Tr(D V_loc)
}

std::ostream& PP_Local::Write(std::ostream& os) const
{
    return os << "    Local pseudopotential V_loc(r) (mesh quadrature), " << theStructure->GetNumAtoms()
              << " atom(s)." << std::endl;
}

} //namespace
