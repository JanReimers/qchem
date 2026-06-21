// File: BasisSet/Molecule/PolarizedGaussian/Symmetry.C
// Bridge from the Cartesian-Gaussian molecular basis to the symmetry machinery: extract the
// AoShell layout (centers, monomials, normalization, shell types) the representation/SALC
// builders consume, and the nuclear point set for point-group detection.  Stage 5 wiring of
// the molecular-symmetry plan (feeds Symmetry::BuildAbelianGroup / BuildSALCs).
module;
#include <vector>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Symmetry;
export import qchem.Symmetry.SALC;        // AoShell, SymPoint, BuildSALCs (transitively)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.Cluster;

export namespace BasisSet::Molecule::PolarizedGaussian
{
using namespace ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

// The AO shell layout of a (flattened) PG basis: one AoShell per Gaussian block, with its
// center, Cartesian monomials, per-component normalization, a center-independent shellType
// (so symmetry-equivalent shells match), and its offset in the global AO ordering.
std::vector<Symmetry::AoShell> ExtractAoShells(const PGData& pg);

// The nuclear point set (species = Z, position) for point-group detection.
std::vector<Symmetry::SymPoint> ClusterToSymPoints(const Cluster& cl);

} //namespace
