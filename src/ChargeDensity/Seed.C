// File: ChargeDensity/Seed.C  SCF seed-density strategy + factory (see doc/SCFSeedingPlan.md).
module;
#include <vector>
#include <utility>
#include <map>
#include <string>

export module qchem.ChargeDensity.Seed;
export import qchem.ChargeDensity;   // tDM_CD<T>
import qchem.BasisSet;               // tBasisSet<T>, and (re-exported) Structure / ElectronConfiguration

export namespace qchem::ChargeDensity
{

//! Structure-wide electronegativity heuristic for IonicSAD: given each atom's atomic number and valence
//! electron count (\c {Z, Nval}, e.g. NaF = {{11,1},{9,7}}), assign a charge-CONSERVING integer formal
//! charge per atom.  Electrons flow from the least- to the most-electronegative atoms (Pauling, from the
//! periodic table), each donor capped at its valence (down to the noble core) and each acceptor at the
//! next closed shell {2,8,18,32}; only across a real electronegativity gap.  Returns the formal charges in
//! the SAME order as \a atoms (Na -> +1, F -> -1); a single species / no EN gap gives all-zero (neutral).
std::vector<int> IonicFormalCharges(const std::vector<std::pair<int,int>>& atoms);

//! \brief The COLLINEAR magnetic decoration of a structure -- Shubnikov S3, doc/SymmetryUpgradePlan.md
//! §7 step 7: one label per atom, \f$\pm1\f$ = the site's majority-channel sign (the atom's
//! \c itsSpinFlip bit), \f$0\f$ = non-magnetic.  MAGNETIC-OR-NOT follows the SEED's OWN resolution rule
//! (\c SeedCD): a species is magnetic iff the library carries a spin-resolved pair at its TARGET
//! electron count (\a ionicNvalByZ entry, else the neutral valence count) -- so the decoration and the
//! seed that plants the moments cannot drift apart.  Assembled HERE, beside that rule, and threaded to
//! the basis factory's symmetry-policy resolution (\c GPWParams::siteSpins): the factory then imposes
//! the SHUBNIKOV group of this decoration instead of the grey group, which would erase the order.
std::vector<int> MagneticDecoration(const Structure* st, const std::string& functional = "LDA",
                                    const std::map<size_t,int>& ionicNvalByZ = {});

//! \brief The IonicSAD per-species TARGET valence counts, derived from the structure exactly as the
//! IonicSAD seed derives them (formal charges from \c IonicFormalCharges, \f$N_{val}\f$ from the
//! neutral valence density) -- the ONE resolution the seed and the S3 \c MagneticDecoration share, so
//! a driver decorating an IonicSAD run passes \c IonicSADTargets(st) and cannot drift from its seed.
std::map<size_t,int> IonicSADTargets(const Structure* st, const std::string& functional = "LDA");


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
//!
//! \a polarized (the Hamiltonian's \c IsPolarized -- the SCFIterator passes it through) upgrades the
//! plane-wave SAD/IonicSAD seeds to the TWO-CHANNEL \c PolarizedSeedCD (doc/SCFSeedingPlan.md §10): the
//! library's per-species Hund pairs + the structure's per-atom \c itsSpinFlip bits assemble a genuinely
//! spin-resolved seed (an AFM staggering, a doublet's majority channel), instead of the spin-agnostic
//! total that collapses to rho/2 in the polarized terms.  Ignored (harmless) on unpolarized runs.
template <class T> tChargeDensity<T>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<T>* bs,
                                                      const Structure* st, const ElectronConfiguration* ec,
                                                      bool polarized=false);

} //namespace
