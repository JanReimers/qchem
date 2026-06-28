// File: ChargeDensity/Seed.C  SCF seed-density strategy + factory (see doc/SCFSeedingPlan.md).
module;

export module qchem.ChargeDensity.Seed;
export import qchem.ChargeDensity;   // tDM_CD<T>
import qchem.BasisSet;               // BasisSet<T>, and (re-exported) Structure / ElectronConfiguration

export namespace qchem::ChargeDensity
{

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
//! \c CoreGuess.  \a st is the molecular/crystal structure -- unused for \c CoreGuess/\c Uniform,
//! threaded for the SAD seeds (Phases 1-3).
template <class T> tDM_CD<T>* MakeSeedDensity(SeedStrategy s, const BasisSet::BasisSet<T>* bs,
                                              const Structure* st, const ElectronConfiguration* ec);

} //namespace
