// File: BasisSet/Molecule/PG_Spherical/LatticeView.C
//
// The SPHERICAL LATTICE VIEW (doc/SphericalLatticePlan.md I1): a Real_BS whose single orbital block
// answers Real_OIBS + Molecule::LatticeSum1E in the SPHERICAL span -- the peer implementation the GPW
// factory consumes by plain virtual dispatch (it never learns which family answered).
//
// STRATEGY (private, invisible to consumers): compose the wrapped CARTESIAN basis's own faces and
// transform at the boundary with the fixed per-shell cart->sphere matrix T (nCart x nSph; s,p blocks
// identity, each Cartesian-d 6-block -> 5 solid-harmonic columns, f 10 -> 7):
//     matrices  M_sph = T^T M_cart T          vectors  b_sph = T^T b_cart
//     densities D_cart = T D_sph T^T (once per consumer call; the collocation streams, ladder and
//     T3 fold then run UNCHANGED on the Cartesian side)
// No new integrals, no spherical lattice sums -- and ENGINE-BLIND by construction: the inner basis is
// held ONLY through the abstract faces (Molecule::Orbital_1E_IBS + LatticeSum1E), so whichever engine
// (qchem-MnD today, libCint later) computes the Cartesian primitives, the view neither knows nor cares.
//
// The harmonics come from Math::SphericalShell (RAW relative weights); each column is normalised
// against the inner shell's own overlap block, so the m-ordering and normalisation conventions are
// taken from the code's single source of truth rather than re-derived (the SphericalSALCPlan S3
// lesson).  The shell composition is read through the ShellRep::Monomials() soft capability.
module;
#include <memory>
export module qchem.BasisSet.Molecule.PG_Spherical.LatticeView;
import qchem.BasisSet;   // Real_BS (the only type in the exported surface)

export namespace qchem::BasisSet::Molecule::PG_Spherical
{

//! \brief Wrap a CARTESIAN-monomial molecular basis in its spherical (contaminant-free) lattice view.
//! The returned Real_BS has one orbital block implementing Real_OIBS + LatticeSum1E in the spherical
//! span (nSph = nCart - #contaminants); hand it to \c L3::GPWFactory exactly like the Cartesian one.
//! Throws \c std::runtime_error if the wrapped basis's orbital block does not expose the faces the
//! view composes (LatticeSum1E + GetAoShells with Cartesian-monomial shells).
std::shared_ptr<const Real_BS> MakeSphericalLatticeView(std::shared_ptr<const Real_BS> cart);

} //namespace
