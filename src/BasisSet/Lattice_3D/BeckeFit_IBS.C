// File: BasisSet/Lattice_3D/BeckeFit_IBS.C  Becke-mesh quadrature fit basis (the Becke XC route's mesh owner).
//
// doc/SymmetryUpgradePlan.md §6a (user diagnosis, 2026-08-01): the Becke XC route conflated the term
// with a raw mesh, where every other fitted term holds a FIT BASIS that owns its integration mesh --
// which left the T2 symmetrization homeless.  This class restores the separation: the QUADRATURE face
// of the Becke route, created by GPW_IBS::CreateVxcFitBasisSet (the grid-policy home) exactly like its
// sibling PlaneWaveFit_IBS, and OWNING
//   - the periodic atom-centred Becke mesh,
//   - the §3-imposed symmetry: when the ctor-injected crystal ops are non-empty (reduceBZ asserted),
//     the mesh is group-averaged INVARIANT (W1's generic route, MakeInvariant) and the orbit fold is
//     prepared, so SymmetrizeRaster becomes the EXACT pointwise star-average (orbit mean) -- the Becke
//     sibling of PlaneWaveFit_IBS's voxel-shift raster average.  Free run: plain mesh, no-op hook.
// The E/H adjoint needs NO term-side change: on the invariant, orbit-symmetric-weight mesh the
// projector is self-adjoint in the w-inner product and v(rho_sym) is already symmetric, so the
// existing Phi^dag diag(w v) Phi quadrature IS the exact derivative of the symmetrized energy (§6a).
//
// A quadrature carries NO fit functions: the IrrepBasisSet face answers honestly empty (zero
// functions, empty op()/Gradient) -- nothing on the Becke route ever asks for fit functions.
module;
#include <iosfwd>
#include <string>
#include <vector>
export module qchem.BasisSet.Lattice_3D.BeckeFit_IBS;
export import qchem.BasisSet.Fit_IBS;                    // cFIT_SF_Quadrature (the mesh-quadrature face)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.SymmetrizeMesh;                             // MakeInvariant / FoldMesh (qcStructure: where Mesh meets Symmetry)
import qchem.Symmetry.Lattice_3D.SpaceGroup;             // DirectOp {W|τ} (ctor-injected crystal ops)
import qchem.Symmetry.Lattice_3D.Fold;                   // Fold + SymmetrizeValues (the orbit-mean star-average)
import qchem.Symmetry;                                   // sym_t (the Bloch irrep label)
import qchem.Types;                                      // dcmplx, cvec_t
import qchem.Matrix3D;                                   // Matrix3D (the cell matrix A)

export namespace qchem::BasisSet::Lattice_3D
{

class BeckeFit_IBS
    : public virtual BasisSet::cFIT_SF_Quadrature       // FIT_SF_Quadrature<dcmplx> : cFIT_SF_ABS
    , public         BasisSet::IrrepBasisSetImp<dcmplx> // GetSymmetry/GetSymt/GetIrrep
{
public:
    //! \param mesh   the periodic Becke mesh (from Structure::CreateIntegrationMesh).
    //! \param sym    the Bloch irrep label (Γ; the mesh is k-independent).
    //! \param A      the cell matrix (columns = lattice vectors; the ops act on fractional coords).
    //! \param directOps  the §3-IMPOSED crystal ops \f$\{W|\tau\}\f$ -- non-empty ONLY when the caller
    //!               asserted imposition (reduceBZ); then the mesh is group-averaged invariant and
    //!               \c SymmetrizeRaster becomes the exact pointwise star-average.  Empty = free run.
    //! \param meshID the MeshParams identity string (the BasisSetID axis).
    BeckeFit_IBS(qcMesh::Mesh mesh, const sym_t& sym, const Matrix3D<double>& A,
                 std::vector<Symmetry::Lattice_3D::DirectOp> directOps, const std::string& meshID);

    //! No fit ever runs on a quadrature face (there is no metric to be non-orthonormal IN); declaring
    //! orthonormal keeps the \c FIT_SF_NonOrtho contract un-implied.
    bool isOrtho() const override {return true;}

    const qcMesh::Mesh& IntegrationMesh() const override {return itsMesh;}

    //! The pointwise star-average (orbit mean) over the invariant mesh's fold -- exact projector under
    //! the imposed group (§6a W1).  No-op on a free run (empty ops at construction).
    void SymmetrizeRaster(rvec_t& rho) const override
    {
        if (itsImpose) Symmetry::Lattice_3D::SymmetrizeValues(itsFold, rho);
    }

    // --- IrrepBasisSet face: a quadrature carries NO fit functions (honest empties, no stubs). ---
    size_t  GetNumFunctions() const override {return 0;}
    cvec_t  operator()(const rvec3_t&) const override {return cvec_t{};}
    cvec3vec_t Gradient(const rvec3_t&) const override {return cvec3vec_t{};}
    std::string Name      () const override {return "BeckeFit";}
    std::string BasisSetID() const override;
    std::ostream& Write(std::ostream&) const override;

private:
    qcMesh::Mesh                  itsMesh;      //!< the (invariant, when imposed) Becke quadrature
    Symmetry::Lattice_3D::Fold    itsFold;      //!< its orbit partition (imposed runs only)
    bool                          itsImpose=false;
    std::string                   itsMeshID;    //!< MeshParams::ID() -- the identity axis
};

} // namespace
