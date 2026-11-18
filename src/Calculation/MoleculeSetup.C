// File: Calculation/MoleculeSetup.C
//
// qchem::MoleculeSetup -- THE RUN SETUP of a molecular calculation (doc/CalcNetworkReview.md §2 and §7, step 7), the twin of
// qchem::SolidSetup: the options as asked, and every decision RESOLVED from them ONCE -- the pseudopotential set (or
// all-electron), the spin group (an open shell promoted to polarized), the electron counts, the seed, the accelerator.  An
// immutable value, shared: a geometry point (a MoleculeCalculation) is built FROM it and changes none of it.  RunOptions lives here
// because the setup is what reads it; qchem.MoleculeCalculation re-exports it.
module;
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
export module qchem.MoleculeSetup;

import qchem.Structure;            // Structure
import qchem.Hamiltonian.Factory;  // Hamiltonian::Model
export import qchem.Symmetry.Spin; // SpinGroup -- the imposed spin subgroup, a RunOptions field
import qchem.Mesh;                 // qcMesh::MeshParams (the DFT integration grid)
import qchem.ChargeDensity.Seed;   // SeedStrategy
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSet::Gaussian::Angular -- the ONE angular option, shared with SolidSetup
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type (the resolved accelerator)
export import qchem.RunOptions;        // RunOptions: the one options tree (V-FACADE-PARITY P0b-2)
export import qchem.SetupDecisions;   // AcceleratorOptions (the accelerator block, every structure), the spin rule
export import qchem.SCFRecipe;         // SCFRecipe (AcceleratorOptions::recipe)
export import qchem.Pseudopotential.Set;   // PseudopotentialSet -- the setup's resolved choice
export import qchem.Accuracy;            // AccuracyTarget -- a RunOptions field
export import qchem.SetupRefusal;        // SetupRefusal + Outcome: the Admit door's failure half

export namespace qchem
{

using Hamiltonian::Model;          // {E1, HF, DE1, DHF}




//! \brief THE RUN SETUP of a molecular calculation: \c RunOptions as asked, plus every decision resolved from them ONCE.
//!
//! Only geometry-INDEPENDENT decisions live here (doc/CalcNetworkReview.md §2): the pseudopotential set, the spin group, the
//! electron counts, the seed, the accelerator and its EMax gate.  The basis, the Hamiltonian and the electron configuration
//! belong to a geometry point and are built from this by \c MoleculeCalculation.  Built only by \c Resolve, held as
//! \c shared_ptr<const MoleculeSetup>.
//! The Becke XC mesh's ENTRY RUNG on a molecule (pin 31; doc/GeometryOptimizationPlan.md §14.1): an unstated mesh starts here and
//! climbs only up, so an anchor moves only where this rung misses eps.  Today's former default, kept so nothing moves that met eps.
//! The molecule's DEFAULT XC mesh (unstated \c RunOptions::xcMesh): MHL m3 alpha 2, Lebedev, Becke order 2, both counts UNSTATED (0) so
//! the run climbs the Becke ladder from the entry rung below (pin 31 rule 1).  Was \c CalcOptions::mesh's member default.
inline const qcMesh::MeshParams kMolecularMeshDefault = {.radial = qcMesh::RadialKind::MHL, .nRadial = 0, .mhl_m = 3, .mhl_alpha = 2.0,
                                                         .angular = qcMesh::AngularKind::Lebedev, .angularDegree = 0, .beckeOrder = 2};
//! THE MOLECULE'S VIEW of the options tree: \a asked with the molecule's defaults filled (model HF, all-electron, \c kMolecularMeshDefault)
//! and what a molecule does not consume REFUSED (policy ForeignOption): the `solid` block, Nelec / charge, ortho, spins, the top-level
//! accelerator / recipe.
RunOptions MoleculeOptions(const RunOptions& asked);
//! A caller that edits ONE field of a molecule's XC mesh starts from the molecule's default (\c kMolecularMeshDefault), filled on first use.
inline qcMesh::MeshParams& MoleculeXCMesh(RunOptions& o) { if (!o.xcMesh) o.xcMesh=kMolecularMeshDefault; return *o.xcMesh; }
inline constexpr int kMolecularMeshEntryRadial = 80;
inline constexpr int kMolecularMeshEntryDegree = 29;

struct MoleculeSetup
{
    //! \brief THE ONE RESOLVE PHASE, as an OUTCOME: the door for USER INPUT (the GUI binding, a command line).  Options that contradict
    //! each other or the structure -- a multiplicity whose parity disagrees with the electron count, an unknown accelerator name, an
    //! unstated basis or one that does not cover an element, a route the code does not realize -- come back as a \c SetupRefusal naming
    //! the policy that refused.  A broken invariant still throws.
    static Outcome<std::shared_ptr<const MoleculeSetup>, SetupRefusal> Admit(const Structure& st, const RunOptions& opts, const AcceleratorOptions& acc = {});
    //! \brief \c Admit, THROWING its refusal (a \c policy::Refused, a \c std::invalid_argument, with the refusal's sentence).  For CODE that
    //! states its own options -- tests, probes, a campaign driver -- where a refusal is the caller's programming error (CLAUDE.md).  Never for
    //! user input: the GUI binding and \c scfrun go through \c Admit.
    static std::shared_ptr<const MoleculeSetup> Resolve(const Structure& st, const RunOptions& opts, const AcceleratorOptions& acc = {});

    RunOptions        opts;   //!< as asked, except \c spin: an open shell (multiplicity > 1) is promoted to Polarized
    AcceleratorOptions acc;    //!< as asked
    //! The resolved choice (\c *opts.pp with \c opts.species); EMPTY = all-electron.
    std::optional<Pseudopotential::PseudopotentialSet> pp;
    //! The REFERENCE geometry, deep-copied (a periodic UnitCell stays periodic), with the TRUE atomic charges.
    std::shared_ptr<const Structure> reference;
    int    Ne = 0;             //!< electrons: the VALENCE count under a pseudopotential
    int    nUp = 0, nDown = 0; //!< the electron count per spin channel, from \c opts.multiplicity (\c ResolveSpinChannels, the rule every structure shares)
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;   //!< \c Default resolved per model
    //! The seed library's key (\c SeedKeyOf): per element, the pseudopotential q of the run's pseudo-atom; unkeyed all-electron.
    qchem::ChargeDensity::SeedKey seedKey;
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::DIIS;   //!< \c acc.type, by name
    double eMax = 0.0;         //!< \c acc.diis.engageBelowFD, or the model's default when it was left unstated
    SCFAccelerators::AcceleratorParams accParams;   //!< \c acc resolved over the molecular defaults (\c ResolveAcceleratorParams)
    //! THE RUN'S XC MESH: \c *opts.xcMesh with an unstated axis at the entry rung (\c kMolecularMeshEntryRadial / \c kMolecularMeshEntryDegree), or the rung a pilot settled
    //! on (\c MoleculeCalculation::SizedSetup).  Every Hamiltonian of a point is built on this, never on \c opts.xcMesh->
    qcMesh::MeshParams mesh;
    bool meshStated   = false;   //!< the caller stated \c *opts.xcMesh (either number): audited against eps, never raised (pin 31 rule 3)
    bool meshResolved = false;   //!< a pilot settled \c mesh already (a session's later points): no climb, no audit

    //! The structure a geometry point is built over: \a st's atoms, each carrying its VALENCE charge Z - Zion under a
    //! pseudopotential (so the electron count is the valence count), or \a st as is when all-electron.  Deep-copied.
    std::shared_ptr<const Structure> PointStructure(const Structure& st) const;
    //! THE PRE-RUN PP / FUNCTIONAL CHECK (user 2026-11-12: track, warn and guide, mostly in the GUI): every element's pseudopotential
    //! against this setup's model: its fit (CP2K label + libxc ids), match yes/no, and the model to state to match it.  What the run will warn about,
    //! answered without running, so a GUI reads it here and never recomputes it.  Empty for an all-electron or non-DFT setup.  THROWS on a fit with no libxc mapping.
    std::vector<Hamiltonian::PPFunctionalCheck> PPFunctionalChecks() const;
};

} // namespace qchem
