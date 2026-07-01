// File: BasisSet/Molecule/PG_Spherical/Symmetry.C
// Bridge from the in-house spherical-Gaussian molecular basis (PG_Spherical_MnD) to the symmetry machinery:
// extract the AoShell layout the representation/SALC builders consume.  The spherical counterpart of
// PG_Cart::ExtractAoShells -- each shell's components are the real solid harmonics, and the AoShell.c2s is
// the basis's OWN Cartesian expansion (SphData::comps[].terms), so the rep is built in the exact m-ordering
// and coefficients the basis uses (self-consistent -- no foreign convention to match; see doc/SphericalSALCPlan.md S3a).
module;
#include <vector>
export module qchem.BasisSet.Molecule.PG_Spherical.Symmetry;
export import qchem.Symmetry.SALC;        // AoShell, BuildSALCs (transitively)
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;   // SphData

export namespace qchem::BasisSet::Molecule::PG_Spherical
{
using Evaluators::PG_Spherical_MnD::SphData;

// The AO shell layout of a (flattened) spherical-Gaussian basis: one AoShell per radial shell, its
// components the real solid harmonics.  c2s carries each harmonic's Cartesian expansion (in the basis's own
// m-ordering) and norm the per-component normalization, so BuildOperationRep builds the correct spherical
// operation rep.  Reuses the nuclear point set helper from PG_Cart::StructureToSymPoints (basis-independent).
std::vector<Symmetry::AoShell> ExtractAoShells(const SphData&);

} //namespace
