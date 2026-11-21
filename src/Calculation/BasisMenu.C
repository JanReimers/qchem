// File: Calculation/BasisMenu.C  Which orbital bases a structure admits: the menu a front end offers, as facts it can sort.
//
// The user, 2026-11-21: "for a single atom structure the user should automatically see only atom friendly orbital basis sets in the
// drop down menu.  The GUI cannot be structure neutral, but this particular aspect is on track to be that."  So a front end ASKS
// the library which bases a structure admits instead of keeping its own list.  Design: side-session note laneG-basis-menu-design.md
// (lane G, 2026-11-21; Q1-Q3 ruled by the user "I agree with your recommendations", Q4 by the controller).
//
// A POLICY-LAYER composition, never a library sort (lane H's Group-face note §9.2): it answers per entry (1) is it admitted, and
// if not, why not, and (2) which group its blocks carry natively (G_blk).  The front end SORTS ("the lib does not need think about
// the ordering, it just needs to supply info needed for the GUI to sort the list" -- user, quoted in lane H's note §0).
//
// TODAY the answers are composed from what the policy layer already knows (the route rule, the reader's coverage test, the blocking
// each facade does).  Group-face migration STEP 4 (`BasisSet::Admits(G) -> Admission`, `BasisSet::BlockGroup()`, a tier-1 hub
// change with its own ruling) REFILLS `admitted` / `reason` / `blockGroup` from the basis itself; the record does not change.
//
// It files no `policy` row: it refuses nothing.  The run's own refusals (BasisCoverage, BasisUnstated) still fire at Admit.
module;
#include <string>
#include <vector>
export module qchem.BasisMenu;
import qchem.Structure;
import qchem.RunOptions;            // RunOptions::symmetryTol: the molecule's detection tolerance, as its setup reads it

export namespace qchem
{

//! One orbital-basis choice for a run on a structure: what selects it, whether it is admitted (and why not), and the facts a
//! front end sorts by.  Every choice the structure's route has is listed; an inadmissible one is listed GREYED with the reason,
//! never dropped, so the user sees why it is missing.
struct BasisOffer
{
    std::string route;      //!< the facade that runs it: "atom" (the radial atom) | "molecule" | "crystal"
    std::string key;        //!< the deck spelling that selects it: "atom.type=Slater", "basis.data=dzvp", "basis.family=planeWave"; "" for a route entry
    std::string name;       //!< the menu's word: "Slater", "dzvp", "plane waves"
    bool        admitted;   //!< false = greyed (\c reason says why)
    std::string reason;     //!< why not, the library's sentence; "" when admitted
    std::string blockGroup; //!< G_blk described: the group whose irreps label the solve's blocks ("O(3) (by l, m analytic)", "T (by k)", ...)
    std::string note;       //!< a fact the user should see beside it ("four-component (RKB): relativistic", "Cartesian d only for MN")
    bool        optIn;      //!< offered only when the user asks (user ruling Q1, 2026-11-21): shown greyed-out-unless-asked
};

//! Every orbital-basis choice for a run on \a st, in the library's order (unsorted).  The ROUTE is the deck's rule (`deck::KindOf`):
//! a periodic structure is a crystal, one site is the radial atom, else a molecule.
//! - the radial ATOM: its basis families (deck `atom.type`), all admitted, the RKB types with their relativistic note (Q3); then ONE
//!   opt-in route entry, "3-D Gaussians (one-atom molecule)" (Q1): selecting it is the molecule route, whose menu is the Gaussian files.
//! - a MOLECULE: every registered Gaussian basis file, admitted when it holds a block for every element (the reader's own test, as
//!   the BasisCoverage refusal), G_blk the detected point group's abelian subgroup when imposed.
//! - a CRYSTAL: the same Gaussian files on the GPW route (G_blk T, by k), then the plane-wave family.
//! \a symmetryTol: the molecule's point-group detection tolerance (its setup's).
std::vector<BasisOffer> OrbitalBasisOffers(const Structure& st, double symmetryTol = RunOptions{}.symmetryTol);

//! The MOLECULE route's offers for a non-periodic  st, whatever its site count: what the opt-in one-atom molecule entry leads to
//! (its menu is the molecule's), and what \c OrbitalBasisOffers lists for two sites or more.
std::vector<BasisOffer> MoleculeBasisOffers(const Structure& st, double symmetryTol = RunOptions{}.symmetryTol);

//! The radial atom's basis families in the deck's spelling (`atom.type`), the menu's order -- the names \c OrbitalBasisOffers keys.
//! A unit test reads each back through the deck (the deck's table owns the spelling; this one must agree with it).
std::vector<std::string> RadialAtomBasisNames();

} // namespace qchem
