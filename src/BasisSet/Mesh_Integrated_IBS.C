// File: BasisSet/Mesh_Integrated_IBS.C  Neutral "operator matrix from a real-space scalar field" capability.
//
// The field-operator dual of the fit-based potential path.  A basis that can integrate a real-space scalar
// field f(r) against its orbitals -- <i|f|j> = integral phi_i* f phi_j, and the plain integral integral f --
// on its OWN mesh (a Becke grid for a molecule, a uniform G-grid/FFT for plane waves).  "Tell, don't ask":
// the only input is a real-space ScalarFunction, there are no getters, and the integration scheme is
// entirely the realizer's business (cf. the Band_DFT_IBS design note, which now derives from this base).
//
// This is what a STATIC, smooth external potential (a local pseudopotential V_loc(r)) wants: a raw
// quadrature <i|V_loc|j>, NOT a least-squares fit onto an auxiliary basis (V_loc is density-independent, so
// there is no per-iteration efficiency argument and no reason to accept the fitting error).  It is the
// molecule/solid-agnostic home for that assembly -- the first of the two PP integral types being routed
// through the BasisSet/Fit_ABS network (see doc/MolecularPP_HarmonizationFindings.md).
module;
export module qchem.BasisSet.Mesh_Integrated_IBS;
export import qchem.ScalarFunction;    // ScalarFunction<double> -- the real-space field to integrate
export import qchem.VectorFunction;    // the orbital basis as [phi_i(r)] -- the integrand source
import qchem.Types;                    // hmat_t<T>
import qchem.Structure;                // Structure (the mesh-integrator factory's geometry)
import qchem.Mesh;                     // qcMesh::MeshParams (the mesh resolution)

export namespace qchem::BasisSet
{

//! \brief A basis that integrates a real-space scalar field against its orbitals on its own mesh.  The
//! field-operator capability the KS/PP potential terms require.  \tparam T the matrix scalar (double for a
//! molecule/atom, dcmplx for plane waves).
template <class T> class Mesh_Integrated_IBS
{
public:
    virtual ~Mesh_Integrated_IBS() = default;
    //! Weighted overlap \f$\langle i|f|j\rangle=\int\phi_i^* f\,\phi_j\,d^3r\f$ of the real-space field \a f.
    //! Direct -- \a f is the only weight, no \f$1/r\f$ kernel.  Not cached (a ScalarFunction has no cache ID).
    virtual hmat_t<T> Overlap (const ScalarFunction<double>& f) const=0;
    //! Scalar integral \f$\int f\,d^3r\f$ over the basis's own mesh.
    virtual double    Integral(const ScalarFunction<double>& f) const=0;
};

//! \brief A basis that can PRODUCE a mesh-integrator for a given structure + mesh resolution -- the
//! field-operator analog of \c CreateVxcFitBasisSet / \c CreateCDFitBasisSet.  The orbital basis is the
//! FACTORY of the integrator, exactly as it is the factory of its fit basis: a basis that owns its mesh
//! intrinsically (plane waves) may return itself; a mesh-free orbital (Gaussian) basis builds a fresh Becke
//! integrator over itself.  The caller owns the returned integrator.  A capability MIXIN -- only bases that
//! support real-space mesh integration inherit it (a term \c dynamic_casts ACROSS to it), so atoms and
//! plane waves are not forced to implement it before they need it.
template <class T> class MeshIntegratorSource
{
public:
    virtual ~MeshIntegratorSource() = default;
    virtual Mesh_Integrated_IBS<T>* CreateMeshIntegrator(const Structure*, const qcMesh::MeshParams&) const=0;
};

//! \brief Build a mesh-integrator for the orbital \a basis (a VectorFunction \f$[\phi_i(r)]\f$) over the
//! integration mesh of \a structure at resolution \a mp -- the concrete real-space realization of
//! \c Mesh_Integrated_IBS<double>.  The MESH TYPE is the structure's choice (Structure::CreateIntegrationMesh:
//! Becke for a molecule/atom, uniform/unit-cell-Becke for a lattice), so this helper is geometry-agnostic.
//! The returned object holds a REFERENCE to \a basis (which must outlive it) and owns its mesh.  Caller owns
//! the result.  (Free helper so the orbital IBS mixins need no knowledge of the concrete integrator.)
Mesh_Integrated_IBS<double>* MakeMeshIntegrator(const VectorFunction<double>& basis,
                                                const Structure*, const qcMesh::MeshParams&);

} //namespace
