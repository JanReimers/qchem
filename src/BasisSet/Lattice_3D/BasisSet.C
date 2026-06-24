// File: BasisSet/Lattice_3D/BasisSet.C  Factory + container for 3D-periodic (plane-wave) basis sets.
//
// The public entry point for building a crystal basis set: hand it a Lattice_3D (the cell + its
// Brillouin-zone grid) and a cutoff, get back an abstract BasisSet<dcmplx>.  The concrete container
// (PW_BasisSet) and the per-k PlaneWave_IBS list it owns are an implementation detail (Imp/BasisSet.C);
// callers are forced through the polymorphic BasisSet interface, exactly as for molecular bases.
module;
#include <memory>
export module qchem.BasisSet.Lattice_3D.BasisSet;
export import qchem.BasisSet;                          // Complex_BS (= BasisSet<dcmplx>)
export import qchem.Lattice_3D;                        // Lattice_3D (the crystal structure + BZ grid)
export import qchem.BasisSet.Lattice_3D.PlaneWave_IBS; // PlaneWave_IBS + LocalPotential/SeparablePotential

export namespace BasisSet::Lattice_3D
{

//! \brief Which 3D-periodic basis to build.  PW = plane waves (lineage A); APW/LAPW (lineage B) follow.
enum class Type { PW };

//! \brief Build a 3D-periodic basis set for \a lat at energy cutoff \a Ecut.  Returns an abstract
//! BasisSet<dcmplx> on the heap (caller owns), so callers must use the polymorphic interface.
//! \param loc,nl  optional external pseudopotential (local + KB nonlocal) configured onto every
//!   plane-wave block; if \a loc is null the external term falls back to the bare nuclear \f$-Z/r\f$.
//! \note Single-k for now: the returned basis holds ONE Bloch block at \f$\Gamma\f$.  Phase 2
//!   generalises this to the full BZ k-list (intended to be the only k-loop in the framework).
Complex_BS* Factory(Type type, const ::Lattice_3D& lat, double Ecut,
                    const LocalPotential* loc=nullptr, const SeparablePotential* nl=nullptr);

} //namespace
