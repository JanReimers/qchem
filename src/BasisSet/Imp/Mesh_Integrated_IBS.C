// File: BasisSet/Imp/Mesh_Integrated_IBS.C  Concrete Becke-mesh field-operator (molecular realization).
module;
#include <memory>
module qchem.BasisSet.Mesh_Integrated_IBS;
import qchem.Mesh.Quadrature;          // qcMesh::Mesh, WeightedOverlap, Integrate, ScalarField, BasisField
import qchem.Structure.MolecularMesh;  // MakeMolecularMesh (the atom-centred Becke mesh)
import qchem.Types;                    // hmat_t<double> (== rsmat_t)

namespace qchem::BasisSet
{

namespace {

// Adapt the orbital basis (a VectorFunction: r -> [phi_i(r)]) to the mesh-quadrature BasisField.
class BFView : public qcMesh::BasisField<double>
{
    const VectorFunction<double>& its;
public:
    explicit BFView(const VectorFunction<double>& v) : its(v) {}
    size_t     size()                       const override {return its.GetVectorSize();}
    rvec_t     operator()(const rvec3_t& r) const override {return its(r);}
    rvec3vec_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

// Adapt a real-space ScalarFunction f(r) to the mesh-quadrature ScalarField (the integration weight V).
class SFView : public qcMesh::ScalarField<double>
{
    const ScalarFunction<double>& its;
public:
    explicit SFView(const ScalarFunction<double>& f) : its(f) {}
    double  operator()(const rvec3_t& r) const override {return its(r);}
    rvec3_t Gradient  (const rvec3_t& r) const override {return its.Gradient(r);}
};

// Molecular field-operator: owns an atom-centred Becke mesh, references the orbital basis.  <i|f|j> and
// integral f are the two mesh quadratures the XC path already runs -- reused here for any scalar field.
class BeckeMeshIntegrator : public Mesh_Integrated_IBS<double>
{
    const VectorFunction<double>& itsBasis;
    qcMesh::Mesh                  itsMesh;
public:
    BeckeMeshIntegrator(const VectorFunction<double>& basis, const Structure* cl, const qcMesh::MeshParams& mp)
        : itsBasis(basis), itsMesh(MakeMolecularMesh(*cl, mp)) {}
    hmat_t<double> Overlap (const ScalarFunction<double>& f) const override
        {return qcMesh::WeightedOverlap(itsMesh, BFView(itsBasis), SFView(f));}
    double         Integral(const ScalarFunction<double>& f) const override
        {return qcMesh::Integrate(itsMesh, SFView(f));}
};

} //anon

Mesh_Integrated_IBS<double>* MakeBeckeMeshIntegrator(const VectorFunction<double>& basis,
                                                     const Structure* cl, const qcMesh::MeshParams& mp)
{
    return new BeckeMeshIntegrator(basis, cl, mp);
}

} //namespace
