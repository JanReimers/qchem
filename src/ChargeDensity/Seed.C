// File: ChargeDensity/Seed.C  SCF seed-density strategy + factory (see doc/SCFSeedingPlan.md).
module;
#include <vector>
#include <utility>

export module qchem.ChargeDensity.Seed;
export import qchem.ChargeDensity;   // tDM_CD<T>
import qchem.BasisSet;               // BasisSet<T>, and (re-exported) Structure / ElectronConfiguration

export namespace qchem::ChargeDensity
{

//! Structure-wide electronegativity heuristic for IonicSAD: given each atom's atomic number and valence
//! electron count (\c {Z, Nval}, e.g. NaF = {{11,1},{9,7}}), assign a charge-CONSERVING integer formal
//! charge per atom.  Electrons flow from the least- to the most-electronegative atoms (Pauling, from the
//! periodic table), each donor capped at its valence (down to the noble core) and each acceptor at the
//! next closed shell {2,8,18,32}; only across a real electronegativity gap.  Returns the formal charges in
//! the SAME order as \a atoms (Na -> +1, F -> -1); a single species / no EN gap gives all-zero (neutral).
std::vector<int> IonicFormalCharges(const std::vector<std::pair<int,int>>& atoms);


//! How the SCF loop is seeded with an initial charge density.
//!   - \c CoreGuess  : null density (today's \c cd=0) -- the density-independent core Hamiltonian
//!                     (kinetic + external only), free-electron / core guess.
//!   - \c Uniform    : \f$\rho(r)=N/V\f$, i.e. \f$D=(N/n)\,I\f$ on the first block.  The plane-wave
//!                     default (Hartree+XC active from iteration 0); centralizes the old per-test
//!                     boilerplate.
//!   - \c SAD        : superposition of neutral atomic densities (Phases 1-2, not yet implemented).
//!   - \c IonicSAD   : superposition of ionic atomic densities (Phase 3, not yet implemented).
//!   - \c Default    : resolved per matrix-element type -- molecular (\c double) -> \c CoreGuess,
//!                     plane-wave (\c dcmplx) -> \c Uniform (the behaviour each path has today).
enum class SeedStrategy { Default, CoreGuess, Uniform, SAD, IonicSAD };

//! Build the initial SCF density for basis \a bs / configuration \a ec under strategy \a s.  Returns a
//! heap density the caller owns (the SCFIterator wraps it in a \c shared_ptr), or \c nullptr for
//! \c CoreGuess.  The seed is a \c tChargeDensity (the DFT Fock-build face): \c Uniform/\c CoreGuess give
//! a matrix-backed \c tDM_CD, \c SAD a fit-backed \c NumericCD.  \a st is the molecular/crystal
//! structure -- unused for \c CoreGuess/\c Uniform, threaded for the SAD seeds (Phases 1-3).
template <class T> tChargeDensity<T>* MakeSeedDensity(SeedStrategy s, const BasisSet::BasisSet<T>* bs,
                                                      const Structure* st, const ElectronConfiguration* ec);

} //namespace
