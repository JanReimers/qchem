// File: Hamiltonian/Internal/Imp/PP_Local.C  Local-pseudopotential term (mesh quadrature of V_loc(r)).
//
// The term owns only the MODEL (V_loc as a real-space field) and asks the BASIS to assemble the matrix:
// it obtains a mesh-integrator from the orbital basis (CreateMeshIntegrator -- the field-operator analog of
// CreateVxcFitBasisSet) and calls <i|V_loc|j>.  The mesh + quadrature are the basis's business, exactly as
// the XC path's are (doc/MolecularPP_HarmonizationFindings.md §3a): the term is now mesh-agnostic, so the
// same term would serve a plane-wave basis (which realizes the field-operator in G-space) unchanged.
module;
#include <cassert>
#include <iostream>
#include <memory>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.ScalarFunction;                 // V_loc(r) presented as a ScalarFunction<double>
import qchem.BasisSet.Mesh_Integrated_IBS;   // MeshIntegratorSource -> Mesh_Integrated_IBS (the field-operator)
import qchem.Math;                            // norm(Vector3D)

namespace qchem::Hamiltonian
{

namespace {

// The real-space local pseudopotential as a scalar field: V_loc(r) = Sum_atoms V_loc(Z_a, |r - R_a|).
class VlocField : public ScalarFunction<double>
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
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}   // unused by the field-operator
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
    // Ask the orbital basis to make a mesh-integrator, then assemble <i|V_loc|j> on it.  The basis owns the
    // mesh (Becke grid for a molecule); the term owns only the model.  A plane-wave basis would answer this
    // cast in G-space -- same term.
    auto src=dynamic_cast<const BasisSet::MeshIntegratorSource<double>*>(bs);
    assert(src && "PP_Local requires a mesh-integrating (MeshIntegratorSource) orbital basis");
    std::unique_ptr<const BasisSet::Mesh_Integrated_IBS<double>>
        mi(src->CreateMeshIntegrator(theStructure.get(), itsMeshParams));
    return mi->Overlap(VlocField(*theStructure, *itsVloc));
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
