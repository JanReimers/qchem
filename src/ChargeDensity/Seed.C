// File: ChargeDensity/Seed.C  SCF seed-density strategy + factory (see doc/SCFSeedingPlan.md).
module;
#include <vector>
#include <utility>
#include <map>
#include <string>

export module qchem.ChargeDensity.Seed;
export import qchem.ChargeDensity;   // tDM_CD<T>
export import qchem.ChargeDensity.AtomicDensity;   // SeedKey (the functional + per-Z PP variant the seeds are keyed by)
import qchem.BasisSet;               // tBasisSet<T>, and (re-exported) Structure
import qchem.ElectronConfiguration;
import qchem.UnitCell;               // UnitCell (the plane-wave seeds' cell)
import qchem.Mesh;                   // qcMesh::MeshParams (the molecular SAD seed's own quadrature)

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
//! \c Site::spinFlip bit), \f$0\f$ = non-magnetic.  MAGNETIC-OR-NOT follows the SEED's OWN resolution rule
//! (\c SeedCD): a SITE is magnetic iff the library carries a spin-resolved pair at its TARGET
//! electron count (\a ionicNvalBySite entry, else the neutral valence count) -- so the decoration and the
//! seed that plants the moments cannot drift apart.  Assembled HERE, beside that rule, and threaded to
//! the basis factory's symmetry-policy resolution (\c GPWParams::siteSpins): the factory then imposes
//! the SHUBNIKOV group of this decoration instead of the grey group, which would erase the order.
std::vector<int> MagneticDecoration(const Structure* st, const SeedKey& key = SeedKey{},
                                    const std::vector<int>& ionicNvalBySite = {});

//! \brief The IonicSAD per-SITE TARGET valence counts (one per atom, in atom-index order), derived from the
//! structure exactly as the IonicSAD seed derives them (formal charges from \c IonicFormalCharges,
//! \f$N_{val}\f$ from the neutral valence density) -- the ONE resolution the seed and the S3
//! \c MagneticDecoration share, so a driver decorating an IonicSAD run passes \c IonicSADTargets(st) and
//! cannot drift from its seed.  PER SITE, not per species (OpenWork §3, 2026-10-26): a mixed-valence cell
//! (LiMn2O4: Mn3+ and Mn4+) keeps each atom's own charge state, and the targets sum to the cell's charge.
std::vector<int> IonicSADTargets(const Structure* st, const SeedKey& key = SeedKey{});


//! How the SCF loop is seeded with an initial charge density.
//!   - \c CoreGuess  : null density (today's \c cd=0) -- the density-independent core Hamiltonian
//!                     (kinetic + external only), free-electron / core guess.
//!   - \c Uniform    : \f$\rho(r)=N/V\f$, i.e. \f$D=(N/n)\,I\f$ on the first block.  The plane-wave
//!                     default (Hartree+XC active from iteration 0); centralizes the old per-test
//!                     boilerplate.
//!   - \c SAD        : superposition of NEUTRAL atomic densities (molecular: the all-electron library into
//!                     a NumericCD; periodic: the valence library's form factors, spin-resolved when the
//!                     run is polarized).
//!   - \c IonicSAD   : superposition of IONIC valence densities at the structure's formal charges
//!                     (\c IonicSADTargets).  Periodic only; the facade default there (V2.2), because
//!                     \c Uniform has a stable wrong basin for electron-sparse cells.
//!   - \c Default    : resolved per matrix-element type -- molecular (\c double) -> \c CoreGuess,
//!                     plane-wave (\c dcmplx) -> \c Uniform.  This is the RAW iterator's fallback (a
//!                     plane-wave run whose species may have no library entry); the facades choose above it.
enum class SeedStrategy { Default, CoreGuess, Uniform, SAD, IonicSAD };

//! Build the initial SCF density for basis \a bs / configuration \a ec under strategy \a s.  Returns a
//! heap density the caller owns (the SCFIterator wraps it in a \c shared_ptr), or \c nullptr for
//! \c CoreGuess.  The seed is a \c tChargeDensity (the DFT Fock-build face): \c Uniform/\c CoreGuess give
//! a matrix-backed \c tDM_CD, \c SAD a fit-backed \c NumericCD.  \a st is the molecular structure --
//! unused for \c CoreGuess/\c Uniform, read by the SAD seed.  The periodic overload demands its \a cell.
//!
//! \a polarized (the Hamiltonian's \c IsPolarized -- the SCFIterator passes it through) upgrades the
//! plane-wave SAD/IonicSAD seeds to the TWO-CHANNEL \c PolarizedSeedCD (doc/SCFSeedingPlan.md §10): the
//! library's per-species Hund pairs + the structure's per-atom \c Site::spinFlip bits assemble a genuinely
//! spin-resolved seed (an AFM staggering, a doublet's majority channel), instead of the spin-agnostic
//! total that collapses to rho/2 in the polarized terms.  Ignored (harmless) on unpolarized runs.
//!
//! \a key names the functional and, per Z, the pseudopotential variant q whose seed densities the
//! plane-wave seeds read (the run setup builds it from its \c PseudopotentialSet); the default is unkeyed LDA.
//! \a quadrature (molecular): the mesh a matrix-free SAD seed integrates itself on -- a Coulomb fit has no grid of its own
//! (doc/Pins.md), so the field that is only a function on space brings one; the run passes its XC mesh.
tChargeDensity<double>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<double>* bs,
                                        const Structure* st, const ElectronConfiguration* ec, const qcMesh::MeshParams& quadrature);
tChargeDensity<dcmplx>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<dcmplx>* bs,
                                        const UnitCell& cell, const ElectronConfiguration* ec,
                                        bool polarized=false, const SeedKey& key=SeedKey{});

} //namespace
