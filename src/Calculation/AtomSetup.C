// File: Calculation/AtomSetup.C
//
// qchem::AtomSetup -- THE RUN SETUP of a radial-atom calculation (V-FACADE-PARITY P1, user ruling R3 2026-11-17: "OK Yes we need
// AtomSetup even if it seems trivial"), the twin of MoleculeSetup and SolidSetup: the options as asked, and every decision RESOLVED
// from them ONCE -- the atom's view of the options, the accelerator by name (an unknown one REFUSED, where the atom used to run DIIS
// silently), its engage gate, the seed and the seed library's key.  An immutable value, shared: an AtomCalculation is built FROM it.
// The atom's loose (Z, charge, opts, acc) constructor arguments are this setup now.
module;
#include <memory>
export module qchem.AtomSetup;

export import qchem.Setup;   // Setup: what every run setup holds
export import qchem.Mesh.Maker;   // qcMesh::MeshMaker (meshMaker)
import qchem.Hamiltonian.Factory;   // Hamiltonian::Model (the view's default)

export namespace qchem
{

//! The atom's DEFAULT XC mesh (unstated \c RunOptions::xcMesh): the RADIAL-ONLY set (\c AngularKind::None, one direction at 4 pi), exact
//! because the atom's density is spherical by construction; \c nRadial 0 = unstated: the run enters at 50 MHL points and climbs to eps.
inline const qcMesh::MeshParams kAtomMeshDefault = {.radial = qcMesh::RadialKind::MHL, .nRadial = 0, .mhl_m = 3, .mhl_alpha = 2.0,
                                                    .angular = qcMesh::AngularKind::None, .angularDegree = 0, .beckeOrder = 2};
//! THE ATOM'S VIEW of the options tree: \a asked with the atom's defaults filled (model HF, \c kAtomMeshDefault, ortho Auto) and what an
//! atom does not consume REFUSED (policy ForeignOption): the molecule and solid blocks, Nelec / charge / multiplicity (the atom's Z and
//! charge are its structure's, its spin from the configuration), pp / species (the atom block's pseudopotential + valence), the
//! imposition, spins, orthoTol, label.  (\c saveStateTo is every kind's checkpoint, V-FACADE-PARITY P3.)
RunOptions AtomOptions(const RunOptions& asked);
//! A caller that edits ONE field of an atom's XC mesh starts from the atom's default (\c kAtomMeshDefault), filled on first use.
inline qcMesh::MeshParams& AtomXCMesh(RunOptions& o) { if (!o.xcMesh) o.xcMesh=kAtomMeshDefault; return *o.xcMesh; }

//! \brief THE RUN SETUP of a radial atom: \c RunOptions as asked, plus every decision resolved from them ONCE.
//!
//! The shared decisions are the base's (\c Setup): \c opts (the atom's view; \c seed resolved: an atom seeds from the core guess unless
//! told otherwise; the old radial-only spelling, Lebedev at degree 0, read as \c AngularKind::None), \c acc, \c accParams (over the molecular
//! defaults at the atom's Z-scaled DIIS engage gate), \c accelerator, \c structure (the lone atom at the origin, its net charge), \c Ne,
//! \c seedKey (the pseudo-atom's LDA valence density under a pseudopotential, else unkeyed) and \c xcMesh (\c *opts.xcMesh; the radial count the
//! run climbs to is the calculation's grid policy).  NOT decided here: \c pp (the pseudo-atom is still stated in the \c atom block: the
//! GTH-LDA at \c PseudoValence(), not a set yet), \c nUp / \c nDown (the electron configuration's, built by \c AtomCalculation), \c imposed
//! and \c spins (an atom's symmetry is spherical by construction).
struct AtomSetup : Setup
{
    //! \brief THE ONE RESOLVE PHASE, as an OUTCOME: the door for USER INPUT (a deck, a GUI form, a command line).  \a atom must hold ONE
    //! atom; its Z and net charge are the run's.  An option an atom does not consume, an unknown accelerator name, a gradient consumer
    //! on a route that is not dE/dD come back as a \c SetupRefusal naming the policy that refused.  A broken invariant still throws.
    static Outcome<std::shared_ptr<const AtomSetup>, SetupRefusal> Admit(const Structure& atom, const RunOptions& opts = {},
                                                                     const AcceleratorOptions& acc = {});
    //! \brief \c Admit, THROWING its refusal.  For CODE that states its own options (tests, probes, the valence-basis generator), where a
    //! refusal is the caller's programming error (CLAUDE.md).
    static std::shared_ptr<const AtomSetup> Resolve(const Structure& atom, const RunOptions& opts = {}, const AcceleratorOptions& acc = {});
    //! \c Resolve for the atom of nuclear charge \a Z and net charge \a charge (electrons = Z - charge), at the origin.
    static std::shared_ptr<const AtomSetup> Resolve(int Z, int charge = 0, const RunOptions& opts = {}, const AcceleratorOptions& acc = {});

    int  Z = 0;               //!< the nuclear charge
    //! The \c directmin accelerator name (an alias kept from the atom's json escape hatch): GDM stepping from the first iteration,
    //! seeded by a diagonalization (\c accelerator is GDM).
    bool directMin = false;
    //! The mesh maker over \c structure (a \c StructureMeshMaker<Atom>, V-MESH-FACTORY §3): the atom never moves, so its one geometry
    //! point is the reference and the setup can hold the maker.  Every mesh the atom run builds at a recipe comes from it.
    std::shared_ptr<const qcMesh::MeshMaker> meshMaker;

    //! The pseudo-atom's valence (Zion): \c opts.atom.valence when stated, else \c Ne (a neutral pseudo-atom).  Meaningful only under
    //! \c opts.atom.pseudopotential.
    int PseudoValence() const {return opts.atom.valence>0 ? opts.atom.valence : Ne;}
};

} // namespace qchem
