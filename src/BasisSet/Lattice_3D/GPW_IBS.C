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
import qchem.BasisSet.Lattice_3D.IBS;             // EPW_Orbital1E_IBS<E> (the evaluator-templated mixin)
import qchem.BasisSet.Internal.IrrepBasisSetImp;  // IrrepBasisSetImp<dcmplx>: GetSymmetry/GetSymt/GetIrrep
export import qchem.BasisSet;                      // Real_BS (the molecular Gaussian basis handed to the ctor)
export import qchem.UnitCell;                      // UnitCell (the direct lattice handed to the ctor)
import qchem.Symmetry;                            // sym_t (the Bloch irrep)
import qchem.Types;

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief GPW basis for a single k-point: periodic Gaussians on a lattice.  This increment: the 1E tier at
//! \f$\Gamma\f$.  Built from a molecular Gaussian basis (over the cell's atoms) + the cell.
class GPW_IBS
    : public EPW_Orbital1E_IBS<GPW_Evaluator>       // op()/Gradient/GetNumFunctions/MakeOverlap/MakeKinetic/MakeNuclear
    , public BasisSet::IrrepBasisSetImp<dcmplx>     // supplies GetSymmetry/GetSymt/GetIrrep + itsSymmetry
    , public GPW_Evaluator                          // the shared Gaussian evaluator (Cast() target for the mixin)
{
public:
    //! \brief Primary constructor: the Bloch symmetry IS the k-label (\f$k=\f$ Symmetry::Lattice_3D::Getk).
    //! \param cell  the direct lattice (its atoms carry the Gaussian centres; source of the translation set).
    //! \param irrep the Bloch irrep (a BlochQN); this increment requires \f$k=\Gamma\f$.
    //! \param mol   the molecular Gaussian orbital basis built over \a cell's atoms (kept alive by the evaluator).
    //! \param Rcut  lattice-translation sphere radius (a.u.); \f$\le 0\f$ = home cell only (the finite limit).
    GPW_IBS(const UnitCell& cell, const sym_t& irrep,
            std::shared_ptr<const BasisSet::Real_BS> mol, double Rcut = 0.0);

    //! \brief Convenience constructor in BZ-grid indices: builds the Bloch irrep \c BlochFactory(N,kIndex).
    GPW_IBS(const UnitCell& cell, const ivec3_t& N, const ivec3_t& kIndex,
            std::shared_ptr<const BasisSet::Real_BS> mol, double Rcut = 0.0);

    virtual std::string Name      () const override {return "GPW";}
    virtual std::string BasisSetID() const override; // geometry-aware cache key (Name + molecular ID + k + nR)

    virtual std::ostream& Write(std::ostream&) const override;
};

} //namespace
