//! \file
//! \brief Bridge from the Cartesian-Gaussian molecular basis to the symmetry machinery: extract the
//! \c AoShell layout (centres, monomials, normalization, shell types) the representation/SALC builders
//! consume, plus the nuclear point set for point-group detection.  Feeds \c Symmetry::Molecule::BuildAbelianGroup /
//! \c BuildSALCs.
module;
#include <vector>
export module qchem.BasisSet.Molecule.PG_Cart.Symmetry;
export import qchem.Symmetry.Molecule.SALC;        // AoShell, SymPoint, BuildSALCs (transitively)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.Structure;

export namespace qchem::BasisSet::Molecule::PG_Cart
{
using namespace ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

//! \brief The AO-shell layout of a (flattened) Cartesian PG basis: one \c AoShell per Gaussian block,
//! carrying its centre, a \c CartesianShellRep over its monomials, per-component normalization, a
//! centre-independent \c shellType (so symmetry-equivalent shells match), and its \c offset in the global
//! AO ordering.
std::vector<Symmetry::Molecule::AoShell> ExtractAoShells(const PGData&);

//! \brief The nuclear point set (species \f$=Z\f$, position) for point-group detection.
std::vector<Symmetry::Molecule::SymPoint> StructureToSymPoints(const Structure&);

} //namespace
