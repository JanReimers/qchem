// File: BasisSet/Molecule/PolarizedGaussian/SymmetryAdapt.C
// Legacy PG's factory hook for building a (Molecule-general) SymmetryAdaptedBasisSet from a raw PG
// basis.  The SymmetryAdaptedBasisSet class is basis-agnostic (qchem.BasisSet.Molecule.*); only this
// glue is PG-specific (it extracts AO shells from PG's PGData).  In the PG namespace so it does not
// collide with PolarizedGaussian1::SymmetryAdapt.
module;
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian.SymmetryAdapt;
export import qchem.BasisSet;                              // BasisSet<double>
export import qchem.Cluster;                               // Cluster
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;    // the general SALC basis (return type)

export namespace BasisSet::Molecule::PolarizedGaussian
{

::BasisSet::Molecule::SymmetryAdaptedBasisSet*
SymmetryAdapt(std::shared_ptr<const ::BasisSet::BasisSet<double>> rawBasis, const Cluster& cl, double tol=1e-4);

} //namespace
