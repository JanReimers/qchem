// File: BasisSet/Molecule/PolarizedGaussian1/SymmetryAdapt.C
// PG1's factory hook for building a (Molecule-general) SymmetryAdaptedBasisSet from a raw PG1 basis.
// The SymmetryAdaptedBasisSet class itself is basis-agnostic (qchem.BasisSet.Molecule.*); only this
// glue is PG-specific (it extracts AO shells from PG1's PGData).  Lives in the PG1 namespace so it
// does not collide with the old PolarizedGaussian::SymmetryAdapt used by the legacy SALC test.
module;
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian1.SymmetryAdapt;
export import qchem.BasisSet;                              // BasisSet<double>
export import qchem.Cluster;                               // Cluster
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;    // the general SALC basis (return type)

export namespace BasisSet::Molecule::PolarizedGaussian1
{

// Build the symmetry-adapted basis from a raw PG1 molecular AO basis + its cluster: extract shells ->
// detect point group -> BuildSALCs -> wrap.  The returned object owns the raw basis (KeepAlive).
::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis, const Cluster& cl, double tol=1e-4);

} //namespace
