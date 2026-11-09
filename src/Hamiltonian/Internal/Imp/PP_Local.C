// File: Hamiltonian/Internal/Imp/PP_Local.C  Local-pseudopotential term.
//
// The term owns only the MODEL (V_loc as a real-space field).  A basis that can assemble the field
// (Orbital_PP_IBS, a capability cross-cast) answers in closed form.  Otherwise the term asks the STRUCTURE for its integration
// mesh -- CreateIntegrationMesh is a polymorphic geometry capability (Atom -> radial x angular, Molecule ->
// Becke, lattice -> uniform), so the mesh TYPE follows the geometry, not this term -- and quadratures
// <i|V_loc|j> with the generic qcMesh::MatrixOverlap (the same routine the XC path uses).  The term is
// thus geometry-neutral without any basis-specific dispatch.
module;
#include <cassert>
#include <iostream>
#include <stdexcept>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Mesh.Quadrature;           // qcMesh::Mesh, MatrixOverlap (over qcMath Scalar/VectorFunction)
import qchem.VectorFunction;            // the orbital basis IS-A VectorFunction (the integrand source, passed direct)

namespace qchem::Hamiltonian
{

namespace {

// The real-space local pseudopotential as a scalar field: V_loc(r) = Sum_atoms V_loc(Z_a, |r - R_a|).
class VlocField : public ScalarFunction<double>
{
    const Structure& cl;
    const BasisSet::SpeciesRadialField& v;
public:
    VlocField(const Structure& c, const BasisSet::SpeciesRadialField& vl) : cl(c), v(vl) {}
    double operator()(const rvec3_t& r) const override
    {
        double s=0;
        for (size_t i=0;i<cl.GetNumAtoms();i++) { const Site& a=cl.Sites()[i]; s+=v.ValueR(a.Z, norm(r-a.R), BasisSet::FieldRange::Full); }
        return s;
    }
    //! No gradient: MatrixOverlap samples values only.  It returned a silent 0 until the hand-gradient audit (2026-11-09).
    rvec3_t Gradient(const rvec3_t&) const override {throw std::logic_error("VlocField::Gradient: the sampled local-PP field has no gradient (MatrixOverlap samples values only)");}
};

} //anon

PP_Local::PP_Local(const st_t& st, vloc_t vloc, mesh_t mesh)
    : theStructure(st), itsVloc(std::move(vloc)), itsMesh(std::move(mesh))
{
    assert(theStructure);
    assert(itsVloc);
}

rsmat_t PP_Local::MakeMatrix(const robs_t* bs, const Spin&) const
{
    if (auto* pp=dynamic_cast<const BasisSet::Orbital_PP_IBS<double>*>(bs); pp && pp->CanAssembleField(*itsVloc))
        return pp->MakeSpeciesFieldMatrix(theStructure.get(), *itsVloc, BasisSet::FieldRange::Full);
    if (!itsMesh) throw std::logic_error("PP_Local: this basis needs the PP quadrature, and the Hamiltonian built none");
    return qcMesh::MatrixOverlap(*itsMesh, *bs, VlocField(*theStructure, *itsVloc));
}

void PP_Local::GetEnergy(EnergyBreakdown& te, const rDM_CD* cd) const
{
    const double e=cd->DM_Contract(this);   // electron-ion (local PP) energy = Tr(D V_loc)
    te.Add("Een", e, EnergyRole::Potential, e);   // linear
}

std::ostream& PP_Local::Write(std::ostream& os) const
{
    return os << "    Local pseudopotential V_loc(r) (" << (itsMesh ? "mesh quadrature" : "closed form") << "), " << theStructure->GetNumAtoms()
              << " atom(s)." << std::endl;
}

} //namespace
