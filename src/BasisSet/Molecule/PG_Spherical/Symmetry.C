//! \file
//! \brief Bridge from the in-house spherical-Gaussian molecular basis (PG_Spherical_MnD) to the symmetry
//! machinery: extract the \c AoShell layout the representation/SALC builders consume -- the spherical
//! counterpart of \c PG_Cart::ExtractAoShells.  Each shell's components are the real solid harmonics, and
//! the rep is built from the basis's OWN Cartesian expansion (\c SphData::comps[].terms) in the exact
//! \f$m\f$-ordering and coefficients the basis uses -- self-consistent, no foreign convention to match
//! (see doc/SphericalSALCPlan.md S3a).
module;
#include <vector>
export module qchem.BasisSet.Molecule.PG_Spherical.Symmetry;
export import qchem.Symmetry.Molecule.SALC;        // AoShell, BuildSALCs (transitively)
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;   // SphData

export namespace qchem::BasisSet::Molecule::PG_Spherical
{
using Evaluators::PG_Spherical_MnD::SphData;

//! \brief The AO-shell layout of a (flattened) spherical-Gaussian basis: one \c AoShell per radial shell,
//! its components the real solid harmonics.  Its \c SphericalShellRep is built from each harmonic's Cartesian
//! expansion (in the basis's own \f$m\f$-ordering) with per-component normalization, so \c BuildOperationRep
//! produces the correct spherical operation rep.  (Point set reuses \c PG_Cart::StructureToSymPoints --
//! basis-independent.)
std::vector<Symmetry::Molecule::AoShell> ExtractAoShells(const SphData&);

} //namespace
