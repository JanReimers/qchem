// File: Calculation/Setup.C  qchem::Setup -- what EVERY run setup holds, whatever the structure (V-FACADE-PARITY P1, rulings R3/R4/R6).
//
// The three setups (AtomSetup, MoleculeSetup, SolidSetup) grew apart, each naming the same decision in its own words: a solid's
// `polarized` was a molecule's promoted `opts.spin`, a molecule's `seed` was a solid's `opts.seed`, `mesh` and `xcMesh`, `reference` and
// `cell`.  The user (2026-11-17, R4): the rows are nearly all ONE concept -- "accidents of how each setup grew ... with this refactor the
// compiler will then be able to help us keep things harmonized".  So the shared decisions live here, under ONE name each, and a
// derived setup adds only what is truly its own (the basis and its recipe, the lattice and its k mesh, the radial atom's grids).
//
// WHERE A DECISION LIVES: a decision that HAS an options field is resolved IN PLACE in `opts` (the spin group, the seed, the model);
// the members below hold what the options cannot say (the set the pp name resolves to, the electron counts, the run's XC mesh,
// whether the imposition survived the run policy).  One home per decision, so two readers cannot drift.
module;
#include <memory>
#include <optional>
#include <vector>
export module qchem.Setup;

export import qchem.RunOptions;               // RunOptions + AcceleratorOptions: the one options tree
export import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type, AcceleratorParams
export import qchem.Pseudopotential.Set;      // PseudopotentialSet
export import qchem.Structure;                // Structure
export import qchem.ChargeDensity.AtomicDensity;   // SeedKey
export import qchem.ChargeDensity.Factory;         // DensityMaker
export import qchem.Mesh;                     // qcMesh::MeshParams
export import qchem.Symmetry.Spin;            // SpinGroup
export import qchem.SetupRefusal;             // SetupRefusal + Outcome: every Admit door's failure half

export namespace qchem
{

//! \brief THE COMMON RUN SETUP: the options as resolved, plus the decisions every structure makes ONCE from them.
//!
//! A value record, immutable once its derived setup's \c Resolve returns, held as \c shared_ptr<const ...>.  Not built on its own: an
//! \c AtomSetup, \c MoleculeSetup or \c SolidSetup is (the constructors are protected).
struct Setup
{
    //! THE KIND'S VIEW of the options asked for (\c AtomOptions / \c MoleculeOptions / \c SolidOptions: the kind's defaults filled, a
    //! foreign option refused), with every decision that has an options field RESOLVED IN PLACE: \c spin (\c ResolveSpinGroup: a
    //! multiplicity > 1 forces Polarized, R6) and \c seed (\c Default resolved per model and structure).
    RunOptions         opts;
    AcceleratorOptions acc;   //!< the accelerator block as asked
    //! \c acc resolved over the kind's defaults (\c ResolveAcceleratorParams): every knob concrete, the DIIS engage gate included.
    SCFAccelerators::AcceleratorParams accParams;
    //! The accelerator TYPE the run states, by name (\c AcceleratorTypeByName: an unknown name is refused, never a silent DIIS).  A
    //! solid may leave it \c Auto: its stages decide per regime.
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::DIIS;
    //! The run's pseudopotential set (\c opts.pp with \c opts.species applied); EMPTY = all-electron.  A periodic run always has one
    //! today.  (The radial atom's pseudo-atom is still stated in its \c atom block, the GTH-LDA at \c atom.valence: not a set yet.)
    std::optional<Pseudopotential::PseudopotentialSet> pp;
    //! The REFERENCE geometry, the TRUE atomic charges: the structure the setup was resolved at (a periodic UnitCell stays periodic).
    std::shared_ptr<const Structure> structure;
    int Ne    = 0;            //!< the electrons the SCF carries (the VALENCE count under a pseudopotential)
    int nUp   = 0, nDown = 0; //!< per spin channel, from \c opts.multiplicity (\c ResolveSpinChannels, the rule every structure shares)
    //! \brief THE RUN'S DENSITY MAKER (D-MUTABLE-STATICS-RUNPOLICY group E, the user's ruling 2026-11-23, E2'): the rho route every block
    //! density of this run is made on, handed to the wave function, the seed and a restart.  A molecule and an atom take the qchem default
    //! (the factored route: they have no policy block); a solid's setup resolves it from its policy (stated, cp2kCompat, or the default)
    //! and reports the decision (policy DensityRoute).
    ChargeDensity::DensityMaker densityMaker{ChargeDensity::RhoRoute::PivotedCholesky};
    ChargeDensity::SeedKey seedKey;   //!< the seed library's key: the run's own pseudo-atoms (\c SeedKeyOf), unkeyed all-electron
    qcMesh::MeshParams     xcMesh;    //!< THE RUN'S XC MESH, every unstated axis resolved (an entry rung, \c Auto, the radial-only set)
    //! The run IMPOSES its symmetry: the caller asked AND the policy allowed it (a solid's run policy can veto it).  Nothing reads
    //! \c opts.imposeSymmetry after \c Resolve, or the two would drift.
    bool imposed = false;
    //! The magnetic decoration the run IMPOSES (one +1/-1/0 per site; was a solid's \c siteSpins): EMPTY when nothing magnetic is imposed.
    std::vector<int> spins;

    bool Polarized() const {return opts.spin==SpinGroup::Polarized;}   //!< the spin-native two-channel run

protected:
    Setup() = default;
    Setup(const Setup&) = default;
    Setup& operator=(const Setup&) = default;
    ~Setup() = default;   //!< never deleted through the base: a derived setup is built by make_shared, which keeps its own deleter
};

} //export namespace qchem
