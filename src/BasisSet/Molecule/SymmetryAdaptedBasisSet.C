//! \file
//! \brief A molecular basis re-expressed in symmetry-adapted (SALC) blocks: one \c SymmetryAdapted_IBS per
//! point-group irrep, each carrying its Mulliken label.  The SCF iterator / accelerators iterate these IBSs
//! exactly as they do the \f$l\f$-channels of an atom -- the molecular case is now "atoms with
//! \f$O\neq\mathbb{1}\f$".
module;
#include <memory>
export module qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;
export import qchem.BasisSet;                         // tBasisSet<double>, Orbital_1E_IBS
export import qchem.Structure;                          // Structure (for the factory hook)
import qchem.BasisSet.Internal.BasisSetImp;           // BasisSetImp (iteration/Insert)
import qchem.BasisSet.SymmetryAdapted_IBS;            // the per-irrep decorator
import qchem.Symmetry.Molecule.SALC;                           // SALCs (the transform O + labels)

export namespace qchem::BasisSet::Molecule
{

//! \brief The whole-molecule basis presented as per-irrep SALC blocks (built by \c PG::SymmetryAdapt from a
//! raw basis + its \c SALCs transform \f$O\f$).  Each block is a \c SymmetryAdapted_IBS decorating the raw
//! basis with its irrep's slice of \f$O\f$.
class SymmetryAdaptedBasisSet
    : public virtual ::qchem::BasisSet::tBasisSet<double>
    , public ::qchem::BasisSet::BasisSetImp<double>
{
public:
    //! \param raw the whole-molecule AO basis -- REFERENCED by the per-irrep decorators (not owned here), so
    //! it must outlive this object.  \param salc the SALC transform \f$O\f$ + irrep labels from \c BuildSALCs.
    SymmetryAdaptedBasisSet(const ::qchem::BasisSet::Orbital_1E_IBS<double>* raw, const Symmetry::Molecule::SALCs& salc);

    //! Optionally hold the raw basis alive (used by \c SymmetryAdapt so the returned object is self-contained
    //! and the caller need not manage the raw basis lifetime separately).
    void KeepAlive(std::shared_ptr<const ::qchem::BasisSet::tBasisSet<double>> raw) {itsRawBasis = raw;}

private:
    std::shared_ptr<const ::qchem::BasisSet::tBasisSet<double>> itsRawBasis;  //!< raw AO basis (kept alive)
};

// The SymmetryAdapt(rawBasis, cl) factory is basis-specific (it extracts AO shells from the concrete
// PGData) and lives in the basis tree: PG_Cart::SymmetryAdapt builds this same
// Molecule-general SymmetryAdaptedBasisSet.

} //namespace
