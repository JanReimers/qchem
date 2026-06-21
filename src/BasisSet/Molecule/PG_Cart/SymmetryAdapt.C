// File: BasisSet/Molecule/PG_Cart/SymmetryAdapt.C
// PG's factory hook for building a (Molecule-general) SymmetryAdaptedBasisSet from a raw PG basis.
// The SymmetryAdaptedBasisSet class itself is basis-agnostic (qchem.BasisSet.Molecule.*); only this
// glue is PG-specific (it extracts AO shells from PG's PGData).  Lives in the PG namespace so it
// does not collide with the old PG_Cart::SymmetryAdapt used by the legacy SALC test.
module;
#include <memory>
export module qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;
export import qchem.BasisSet;                              // BasisSet<double>
export import qchem.Cluster;                               // Cluster
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;    // the general SALC basis (return type)

export namespace BasisSet::Molecule::PG_Cart
{

// Build the symmetry-adapted basis from a raw PG molecular AO basis + its cluster: extract shells ->
// detect point group -> BuildSALCs -> wrap.  The returned object owns the raw basis (KeepAlive).
::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis, const Cluster& cl, double tol=1e-4);

} //namespace
