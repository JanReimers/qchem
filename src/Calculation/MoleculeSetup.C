// File: Calculation/MoleculeSetup.C
//
// qchem::MoleculeSetup -- THE RUN SETUP of a molecular calculation (doc/CalcNetworkReview.md §2 and §7, step 7), the twin of
// qchem::SolidSetup: the options as asked, and every decision RESOLVED from them ONCE -- the pseudopotential set (or
// all-electron), the spin group (an open shell promoted to polarized), the electron counts, the seed, the accelerator.  An
// immutable value, shared: a geometry point (a MoleculeCalculation) is built FROM it and changes none of it.  RunOptions lives here
// because the setup is what reads it; qchem.MoleculeCalculation re-exports it.
module;
#include <functional>
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
export import qchem.Mesh.Maker;    // qcMesh::MeshMaker: a geometry point carries the maker over its own structure (V-MESH-FACTORY)
import qchem.ChargeDensity.Seed;   // SeedStrategy
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSet::Gaussian::Angular -- the ONE angular option, shared with SolidSetup
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type (the resolved accelerator)
export import qchem.RunOptions;        // RunOptions: the one options tree (V-FACADE-PARITY P0b-2)
export import qchem.SetupDecisions;   // AcceleratorOptions (the accelerator block, every structure), the spin rule
export import qchem.SCFRecipe;         // SCFRecipe (AcceleratorOptions::recipe)
export import qchem.Pseudopotential.Set;   // PseudopotentialSet -- the setup's resolved choice
export import qchem.Accuracy;            // AccuracyTarget -- a RunOptions field
export import qchem.SetupRefusal;        // SetupRefusal + Outcome: the Admit door's failure half
export import qchem.Setup;               // Setup: what every run setup holds (V-FACADE-PARITY P1)

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

struct MoleculeSetup : Setup
{
    //! \brief THE ONE RESOLVE PHASE, as an OUTCOME: the door for USER INPUT (the GUI binding, a command line).  Options that contradict
    //! each other or the structure -- a multiplicity whose parity disagrees with the electron count, an unknown accelerator name, an
    //! unstated basis or one that does not cover an element, a route the code does not realize -- come back as a \c SetupRefusal naming
    //! the policy that refused.  A broken invariant still throws.
    //! TYPED (V-MESH-FACTORY §8 step 2): the door takes the CONCRETE structure, so the setup can build every geometry point's mesh maker
    //! without a cast (CLAUDE.md, "constructors pass concrete objects").  A caller holding only a \c Structure& spells out its concrete type.
    static Outcome<std::shared_ptr<const MoleculeSetup>, SetupRefusal> Admit(const Molecule& st, const RunOptions& opts, const AcceleratorOptions& acc = {});
    //! \brief \c Admit, THROWING its refusal (a \c policy::Refused, a \c std::invalid_argument, with the refusal's sentence).  For CODE that
    //! states its own options -- tests, probes, a campaign driver -- where a refusal is the caller's programming error (CLAUDE.md).  Never for
    //! user input: the GUI binding and \c scfrun go through \c Admit.
    static std::shared_ptr<const MoleculeSetup> Resolve(const Molecule& st, const RunOptions& opts, const AcceleratorOptions& acc = {});

    //! The shared decisions are the base's (\c Setup): \c opts (\c spin promoted for an open shell, \c seed resolved per model), \c acc,
    //! \c accParams (the DIIS engage gate: the model's default when unstated), \c accelerator, \c pp, \c structure (deep-copied, TRUE charges),
    //! \c Ne / \c nUp / \c nDown, \c seedKey, and \c xcMesh: \c *opts.xcMesh with an unstated axis at the entry rung
    //! (\c kMolecularMeshEntryRadial / \c kMolecularMeshEntryDegree), or the rung a pilot settled on (\c MoleculeCalculation::SizedSetup).
    //! Every Hamiltonian of a point is built on \c xcMesh, never on \c opts.xcMesh.
    bool meshStated   = false;   //!< the caller stated \c *opts.xcMesh (either number): audited against eps, never raised (pin 31 rule 3)
    bool meshResolved = false;   //!< a pilot settled \c xcMesh already (a session's later points): no climb, no audit

    //! \brief ONE GEOMETRY POINT: the structure it is built over and the mesh maker over THAT structure (V-MESH-FACTORY §8 step 2).
    //! The structure is the reference's atoms (moved, for \c MakePoint(R)), each ADDING its core charge Z - Zion to the net charge it
    //! states under a pseudopotential (so the electron count is the valence count, an ion's charge kept), or as is when all-electron.
    //! WHY THE SETUP BUILDS THE MAKER: a point's concrete structure type is known only where the reference was handed in (the typed
    //! \c Resolve); a point rebuilt later from the abstract face (\c Clone + \c MoveSites) has lost it, and getting it back would be a cast.
    //! So \c Resolve stores a TYPED point builder, and every point (the reference, or a geometry step) comes out of it with its maker.
    struct Point
    {
        std::shared_ptr<const Structure>         structure;
        std::shared_ptr<const qcMesh::MeshMaker> meshMaker;   //!< builds meshes over \c structure, never over the reference
    };
    Point MakePoint() const;                                //!< at the reference geometry
    Point MakePoint(const std::vector<rvec3_t>& R) const;   //!< the reference's sites moved to \a R (one position per atom, atom order)
    //! THE PRE-RUN PP / FUNCTIONAL CHECK (user 2026-11-12: track, warn and guide, mostly in the GUI): every element's pseudopotential
    //! against this setup's model: its fit (CP2K label + libxc ids), match yes/no, and the model to state to match it.  What the run will warn about,
    //! answered without running, so a GUI reads it here and never recomputes it.  Empty for an all-electron or non-DFT setup.  THROWS on a fit with no libxc mapping.
    std::vector<Hamiltonian::PPFunctionalCheck> PPFunctionalChecks() const;

private:
    //! The typed point builder \c Resolve stores: copies the reference as its concrete type, moves its sites (\a R non-null), adds the
    //! core charges (\a pp non-null), and pairs it with \c StructureMeshMaker of that type.
    using PointBuilder=std::function<Point(const std::vector<rvec3_t>* R, const Pseudopotential::PseudopotentialSet* pp)>;
    PointBuilder itsPointBuilder;
    //! The structure-neutral body of \c Resolve, once the typed door has stored \c structure and \c itsPointBuilder in \a s.
    static std::shared_ptr<const MoleculeSetup> Decide(std::shared_ptr<MoleculeSetup> s, const RunOptions& asked, const AcceleratorOptions& acc);
};

} // namespace qchem
