// File: BasisSet/Molecule/PG_Cart/SymmetryAdapt.C
// PG's factory hook for building a (Molecule-general) SymmetryAdaptedBasisSet from a raw PG basis.
// The SymmetryAdaptedBasisSet class itself is basis-agnostic (qchem.BasisSet.Molecule.*); only this
// glue is PG-specific (it extracts AO shells from PG's PGData).  Lives in the PG namespace so it
// does not collide with the old PG_Cart::SymmetryAdapt used by the legacy SALC test.
module;
#include <memory>
export module qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;
export import qchem.BasisSet;                              // BasisSet<double>
export import qchem.Structure;                               // Structure
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;    // the general SALC basis (return type)

export namespace qchem::BasisSet::Molecule::PG_Cart
{

//! \brief Build a symmetry-adapted basis from a raw molecular AO basis + its structure:
//! \c Structure + \c RawBasisSet \f$\rightarrow\f$ \c SALC_IBS.  Pipeline: the orbital IBS yields its own AO
//! shells (\c Molecule::Orbital_1E_IBS::GetAoShells) \f$\rightarrow\f$ detect the point group
//! \f$\rightarrow\f$ \c BuildSALCs \f$\rightarrow\f$ wrap.  The returned object owns the raw basis (KeepAlive).
//! \image html salc_call_flow.svg "Structure + RawBasisSet -> SALC_IBS dataflow" width=640
//! (source: doc/diagrams/salc_call_flow.svg; add doc/diagrams to the Doxyfile IMAGE_PATH to render.)
::qchem::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::qchem::BasisSet::BasisSet<double>> rawBasis, const Structure&, double tol=1e-4);

} //namespace
