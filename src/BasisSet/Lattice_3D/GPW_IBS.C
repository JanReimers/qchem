// File: BasisSet/Lattice_3D/GPW_IBS.C  Gaussian-And-Plane-Waves irrep basis set for one k-point.
//
// The GPW sibling of PlaneWave_IBS: a complex (dcmplx) Orbital_1E_IBS whose functions are periodic GAUSSIANS
// (Bloch sums of contracted Gaussians standing at the cell's atoms) rather than plane waves.  As with the PW
// basis, all the geometry lives in a shared evaluator (GPW_Evaluator) and the evaluator-templated
// EPW_Orbital1E_IBS<E> mixin forwards op()/Gradient/GetNumFunctions/MakeOverlap/MakeKinetic/MakeNuclear to it
// -- "GPW is a new evaluator, not a new IBS" (doc/MolecularPP_HarmonizationRound2.md section 2.5), so this
// class is almost empty (ctors + identity).  The periodic Gaussian 1E integrals themselves are computed by the
// molecular Gaussian basis GPW_Evaluator owns (via Molecule::LatticeSum1E), never here.
//
// SCOPE (first increment): the 1E tier at the GAMMA point.  The DFT tier (Hartree/XC by collocation) and the
// external pseudopotential are NOT mixed in yet -- they are later increments.
module;
#include <iosfwd>
#include <memory>
#include <string>

export module qchem.BasisSet.Lattice_3D.GPW_IBS;
import qchem.BasisSet.Lattice_3D.Evaluators.GPW;  // GPW_Evaluator (base subobject) -- NOT re-exported (internal)
import qchem.BasisSet.Lattice_3D.IBS;             // EPW_Orbital1E_IBS<E> + EPW_Orbital_DFT_IBS<E> (mixins)
import qchem.BasisSet.Lattice_3D.PlaneWaveFit_IBS; // the auxiliary PW fit basis the DFT factory returns
import qchem.BasisSet.Internal.IrrepBasisSetImp;  // IrrepBasisSetImp<dcmplx>: GetSymmetry/GetSymt/GetIrrep
export import qchem.BasisSet.Band_FT_IBS;          // Band_FT_IBS (the DFT capability; Create*FitBasisSet)
export import qchem.BasisSet.Fit_IBS;              // cFIT_CD_ABS / cFIT_SF_ABS + qcMesh::MeshParams
export import qchem.Pseudopotential.Integrals_Pseudo; // Integrals_Pseudo<dcmplx> (external-PP capability) + the models
export import qchem.BasisSet;                      // Real_BS (the molecular Gaussian basis handed to the ctor)
export import qchem.UnitCell;                      // UnitCell (the direct lattice handed to the ctor)
import qchem.Symmetry;                            // sym_t (the Bloch irrep)
import qchem.Structure;                           // Structure (Create*FitBasisSet arg)
import qchem.Types;

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The lattice-image MODE.  \c Periodic (the default): every lattice sum is an
//! \f$\varepsilon\f$-CONVERGED SERIES enumerated internally per shell pair -- THERE IS NO CUT in the R
//! direction and no radius parameter exists (user pin, doc/GPWPlan.md).  \c HomeCellOnly: the
//! FINITE-molecule configuration (no images anywhere; 1E == the finite matrices, KB bra = the raw home
//! orbital) -- what the finite==lattice box gates compare against.  An \c enum \c class so a stray numeric
//! argument can never silently select a mode.
enum class CellImages { Periodic, HomeCellOnly };

//! \brief GPW basis for a single k-point: periodic Gaussians on a lattice.  This increment: the 1E tier at
//! \f$\Gamma\f$.  Built from a molecular Gaussian basis (over the cell's atoms) + the cell.
class GPW_IBS
    : public EPW_Orbital1E_IBS<GPW_Evaluator>       // op()/Gradient/GetNumFunctions/MakeOverlap/MakeKinetic/MakeNuclear
    , public EPW_Orbital_DFT_IBS<GPW_Evaluator>     // DFT tier (IS-A Band_FT_IBS): MakeOverlap/MakeRepulsion3C/MakeOverlap3C
    , public BasisSet::IrrepBasisSetImp<dcmplx>     // supplies GetSymmetry/GetSymt/GetIrrep + itsSymmetry
    , public Pseudopotential::Integrals_Pseudo<dcmplx> // external-PP assembly (real-space); PW_Pseudo casts ACROSS to this
    , public GPW_Evaluator                          // the shared Gaussian evaluator (Cast() target for the mixins)
{
public:
    //! \c MakeOverlap is in BOTH the 1E mixin (no-arg \f$\langle i|j\rangle\f$) and the DFT mixin (the field
    //! bridge \f$\langle i|f|j\rangle\f$); merge them into one overload set so a concrete-class call is not an
    //! ambiguous multi-base lookup.
    using EPW_Orbital1E_IBS<GPW_Evaluator>::MakeOverlap;
    using EPW_Orbital_DFT_IBS<GPW_Evaluator>::MakeOverlap;

    //! \brief Primary constructor: the Bloch symmetry IS the k-label (\f$k=\f$ Symmetry::Lattice_3D::Getk).
    //! \param cell  the direct lattice (its atoms carry the Gaussian centres; source of the cell geometry).
    //! \param irrep the Bloch irrep (a BlochQN); this increment requires \f$k=\Gamma\f$.
    //! \param mol   the molecular Gaussian orbital basis built over \a cell's atoms (kept alive by the evaluator).
    //! \param densityEcut  the DENSITY-grid selector (GPW's only cutoff): \f$<0\f$ = AUTOMATIC floor
    //!                   \a cutoffFactor\f$\cdot\alpha_{\max}\f$ (recommended -- no Hartree value needed);
    //!                   \f$=0\f$ = DFT tier OFF (1E-only); \f$>0\f$ = explicit Hartree cutoff (warned on \c cerr
    //!                   if below the floor).  See \c GPW_Evaluator.
    //! \param images  the lattice-image MODE (see \c CellImages; default \c Periodic).
    //! \param cutoffFactor  \f$C\ge4\f$ in the density-grid floor \f$C\cdot\alpha_{\max}\f$ (default 4).
    GPW_IBS(const UnitCell& cell, const sym_t& irrep,
            std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut = 0.0,
            CellImages images = CellImages::Periodic, double cutoffFactor = 2.0,
            RasterPolicy raster = RasterPolicy::AliasFree);

    //! \brief Convenience constructor in BZ-grid indices: builds the Bloch irrep \c BlochFactory(N,kIndex).
    GPW_IBS(const UnitCell& cell, const ivec3_t& N, const ivec3_t& kIndex,
            std::shared_ptr<const BasisSet::Real_BS> mol, double densityEcut = 0.0,
            CellImages images = CellImages::Periodic, double cutoffFactor = 2.0,
            RasterPolicy raster = RasterPolicy::AliasFree);

    //! \brief The DFT factory seam (Band_FT_IBS): the auxiliary density/potential fit basis is a plane-wave grid
    //! over GPW's OWN density grid -- so the collocated \f$\tilde\rho\f$'s \f$\{G\}\f$ matches the fitter's.  A
    //! GPW density lives on a plane-wave grid whatever the orbitals are (never orbital==fit).
    virtual BasisSet::cFIT_CD_ABS* CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const override;
    virtual BasisSet::cFIT_SF_ABS* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const override;

    //! \brief The external-PP capability (Integrals_Pseudo<dcmplx>): assemble \f$\langle i|V_{loc}|j\rangle\f$ /
    //! \f$\langle i|V_{NL}|j\rangle\f$ in REAL SPACE (cross-cast the model to its \c *_R face, delegate to the
    //! evaluator's mesh quadrature).  This is what lets the plane-wave \c PW_Pseudo term (and thus the whole
    //! \c Ham_PW_DFT) drive a GPW basis unchanged -- GPW answers the same abstract cross-cast a PW basis does.
    virtual hmat_t<dcmplx> MakeLocalPotential   (const Structure* cl, const Pseudopotential::LocalPotential&     loc) const override;
    //! \brief The long/short split of \c MakeLocalPotential (doc/GPWPlan.md 0e-PP): the LONG (softened-Coulomb)
    //! part is assembled analytically on the smooth density-grid integrate-back (\c GPW_Evaluator::MakeLocalPPLong,
    //! folded into \c PW_Hartree); the SHORT (compact poly-Gaussian) part rides the sharp-field local-PP sweep.
    virtual hmat_t<dcmplx> MakeLocalPotentialLong (const Structure* cl, const Pseudopotential::LocalPotential& loc) const override;
    virtual hmat_t<dcmplx> MakeLocalPotentialShort(const Structure* cl, const Pseudopotential::LocalPotential& loc) const override;
    virtual hmat_t<dcmplx> MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& nl ) const override;

    virtual std::string Name      () const override {return "GPW";}
    virtual std::string BasisSetID() const override; // geometry-aware cache key (Name + molecular ID + k + nR)

    virtual std::ostream& Write(std::ostream&) const override;

protected:
    //! \brief The DFT 3-centre tables (Band_FT_IBS) built over the REQUESTED fit basis's OWN grid -- NOT the
    //! block's \c DensityGrid.  \a c is the fit basis \c CreateCD/VxcFitBasisSet produced (it IS-A
    //! PW_Grid_Evaluator carrying the density-fit \f${G}\f$/grid policy), so we hand its grid to the evaluator:
    //! the table returned is the one REQUESTED, honouring the factory's grid choice rather than silently
    //! overriding it with the block's own (doc/GPWPlan §0e).  Overrides the shared \c EPW_Orbital_DFT_IBS mixin
    //! (which dropped \a c); the block's own no-arg tensors remain the convenience/test path.
    virtual G_ERI3 MakeRepulsion3C(const cFIT_CD_ABS& c) const override;
    virtual G_ERI3 MakeOverlap3C  (const cFIT_SF_ABS& c) const override;
};

} //namespace
