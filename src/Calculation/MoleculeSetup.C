// File: Calculation/MoleculeSetup.C
//
// qchem::MoleculeSetup -- THE RUN SETUP of a molecular calculation (doc/CalcNetworkReview.md §2 and §7, step 7), the twin of
// qchem::SolidSetup: the options as asked, and every decision RESOLVED from them ONCE -- the pseudopotential set (or
// all-electron), the spin group (an open shell promoted to polarized), the electron counts, the seed, the accelerator.  An
// immutable value, shared: a geometry point (a MoleculeCalculation) is built FROM it and changes none of it.  CalcOptions lives here
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
export import qchem.Symmetry.Spin; // SpinGroup -- the imposed spin subgroup, a CalcOptions field
import qchem.Mesh;                 // qcMesh::MeshParams (the DFT integration grid)
import qchem.ChargeDensity.Seed;   // SeedStrategy
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSet::Gaussian::Angular -- the ONE angular option, shared with SolidSetup
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type (the resolved accelerator)
export import qchem.SetupDecisions;   // AcceleratorOptions (the accelerator block, every structure), the spin rule
export import qchem.SCFRecipe;         // SCFRecipe (AcceleratorOptions::recipe)
export import qchem.Pseudopotential.Set;   // PseudopotentialSet -- the setup's resolved choice
export import qchem.Accuracy;            // AccuracyTarget -- a CalcOptions field
export import qchem.SetupRefusal;        // SetupRefusal + Outcome: the Admit door's failure half

export namespace qchem
{

using Hamiltonian::Model;          // {E1, HF, DE1, DHF}

//! Orbital-integral engine: the in-house MnD recursion (default) or the libcint foreign engine.
enum class Engine  { MnD, LibCint };
//! Angular basis convention: Cartesian Gaussians (default) or real solid-harmonic (spherical).  THE library's own axis, shared by the
//! molecular \c CalcOptions::angular and the periodic \c GaussianBasisRecipe::angular (user, 2026-11-12: "we should strive to make
//! seeding and spherical/Cartesian options structure neutral"), so one option means one thing on every structure.
using BasisSet::Gaussian::Angular;

//! How to set up the calculation.  Designated-initializer friendly:
//!     MoleculeCalculation calc(MoleculeSetup::Resolve(water, {.basis="dzvp"}));
//! `model` selects the Hamiltonian; pass 1 wires the non-DFT factory path (HF is the default).
struct CalcOptions
{
    //! The orbital basis: a file stem in the basis data directory ("dzvp", "valence_lowq_sr", ...; BasisNames()).  NO DEFAULT: an
    //! unstated basis is refused at \c Resolve (policy BasisUnstated), naming the registered bases that cover the structure -- the
    //! old default "sto-3g" had no file, so a run that relied on it could never build.
    std::string basis;
    //! A .bsd PATH (it contains a '/'; a bare name means the shipped directory) used INSTEAD of the shipped file \c basis names (\c basis is then only its label) -- the molecular twin of
    //! \c GaussianBasisRecipe::file (lane G 2026-11-10, controller ruling: parity with the solid side; \c basis stays a NAME, never a
    //! path).  What the GUI's Basis Set Exchange download points at (~/.solid2/bases/...).  The report states the path and the
    //! file's checksum, and a saved state refuses a different checksum, so an edited or missing file can never pass as the same
    //! basis.  Empty (the default) = the shipped file \c basis names, exactly as before.  (V-FACADE-OPTS: the options tree takes it
    //! beside its solid twin.)
    std::string basisFile;
    Model       model = Model::HF;   //!< HF (default) | Xalpha | LDA | E1/DE1/DHF (test-only)
    SpinGroup         spin  = SpinGroup::UnPolarized;   //!< the imposed spin subgroup (V1.37)
    //! Spin multiplicity 2S+1.  0 (default) = minimal spin: closed-shell singlet for even Ne, doublet for
    //! odd -- the historical behaviour.  Set explicitly for an open shell: 3 = triplet, 2 = doublet, ...
    //! The facade converts it to (nUp,nDown) [nUp-nDown = 2S = multiplicity-1, nUp+nDown = Ne] and PROMOTES
    //! the calculation to SpinGroup::Polarized when 2S>0 (unrestricted open shell needs distinct up/down densities).
    //! A multiplicity whose parity disagrees with Ne (e.g. a singlet for odd Ne) is rejected.
    int         multiplicity = 0;
    //! Basis construction variants (threaded into BasisSet::Gaussian::Factory).  Defaults reproduce
    //! today's behaviour (in-house MnD, Cartesian).  angular==Spherical + symmetry works with the MnD engine (the SALC builder
    //! handles spherical shells, doc/OldPlans/SphericalSALCPlan.md S3a/S4); only engine LibCint + Spherical + symmetry is
    //! refused (S3b: libcint's real-harmonic layout is unwired).
    Engine      engine  = Engine::MnD;
    Angular     angular = Angular::Cartesian;
    //! DFT-only knobs (ignored when model is HF/1-e/Dirac).  xalpha: the Slater exchange parameter, used
    //! only by model==Xalpha.
    double      xalpha = 0.7;
    //! \brief THE XC MESH, SIZED FROM EPS BY DEFAULT (pin 31 rule 1; doc/GeometryOptimizationPlan.md §14): \c nRadial 0 and
    //! \c angularDegree 0 say UNSTATED, and the run climbs the Becke ladder (\c qchem.Mesh.BeckeClimb) at its converged D from the entry
    //! rung (\c kMolecularMeshEntryRadial / \c kMolecularMeshEntryDegree: 80 / 29), only UP, to the rung that meets \c accuracy on what the client consumes.  A STATED mesh
    //! (either number set) is honoured: audited and reported against eps, never raised (rule 3); an axis left at 0 takes the entry
    //! rung's.  The other fields (MHL m3 alpha 2, Lebedev, Becke order 2) are the rule the climb walks on.
    //! The entry rung was the converged default of 2026-10-31.  MEASURED 2026-11-08 at the converged D against (300, 65), after the
    //! closed-form PP integrals: water LDA |ΔE| 2e-7, max|ΔH| 1e-6; F₂ GTH q7 |ΔE| 2.3e-6, max|ΔH| 2.4e-5 -- so F₂ meets an energy
    //! run's eps x N_atoms and a geometry run raises it (its 2026-10-31 "27 μHa/atom" was the PP's own mesh error, gone since).
    //! It was MHL 30 / Lebedev 5 (12 directions),
    //! 20 mHa/atom LOW on water LDA, 0.52 Ha/atom on F₂ -- the error is ANGULAR and NOT variational, so a coarse mesh reads as a
    //! better energy.
    qcMesh::MeshParams mesh = {.radial  = qcMesh::RadialKind::MHL,   .nRadial   = 0,
                              .mhl_m    = 3,                         .mhl_alpha = 2.0,
                              .angular  = qcMesh::AngularKind::Lebedev, .angularDegree = 0,   // unstated: climbed from (80, 29)
                              .beckeOrder = 2};
    //! The Vxc FIT BASIS (2026-10-31): \c Auto = \c Delta, the DEFAULT (user ruling "δ default, A1 kept") -- the δ basis on
    //! \c mesh's own points, v_xc and eps_xc sampled where the quadrature integrates, no fit; \c Auxiliary = the Gaussian auxiliary
    //! basis (DGauss A1).  Both reach the one XC term through the fit face (pin 30: Vxc = system x fit basis x grid).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
    //! WHICH \f$\rho\f$ XC SEES (user ruling 2026-11-09, \c Hamiltonian::XCρ): \c Exact = \f$\rho[D]\f$, the default; \c Fit = the Vxc fit's
    //! projection \f$P\rho\f$ (variational: \f$E=\int e(P\rho)\f$, \f$H=\partial E/\partial D\f$).  It matters only for \c VxcFit::Auxiliary, whose
    //! \f$P\f$ is inexact: there \c Exact is the Gaussian x δ cell, which runs WITH A WARNING (H is not dE/dD, pin 32).
    Hamiltonian::XCρ xcRho = Hamiltonian::XCρ::Exact;
    //! The Hartree route (2026-11-06): by default the Gaussian Coulomb fit; \c .fit=CoulombFit::Delta puts \f$\rho\f$ on a mesh, with a
    //! first-order energy and the mesh's charge error reported.  See \c Hamiltonian::CoulombRoute.
    Hamiltonian::CoulombRoute coulomb = {};
    //! Point-group SALC blocking + per-irrep aufbau.  GUARDED TO THE CARTESIAN PG BASIS: the SALC
    //! builder needs a PolarizedGaussian (PGData) orbital IBS and throws otherwise.  Since the facade
    //! only builds the default Cartesian basis today, this is always the supported path; the guard
    //! future-proofs the day spherical/libcint deliveries are exposed (see doc/OldPlans/SphericalSALCPlan.md).
    bool        symmetry    = false;
    double      symmetryTol = 1e-4;   //!< geometry tolerance for point-group detection
    //! Replace the bare nuclear attraction with the GTH/HGH pseudopotential (valence-only): the atoms keep
    //! only their valence electrons and the Hamiltonian gets V_loc(r) + the KB-separable projectors +
    //! the ion-ion term, in place of Ven.  \c pp NAMES the pseudopotential set (src/Pseudopotential/Data/pseudopotentials.json,
    //! e.g. "gth-lda-valence": each element's smallest-q GTH-LDA entry); EMPTY (default) = all-electron.  The set also fixes the
    //! functional the entries were generated for.  Multi-species is native (a per-Z router over the molecule's elements).
    //! Requires a valence basis in \c basis (an all-electron basis has spurious core shells under a PP) -- e.g.
    //! {.basis="sipp", .model=Model::LDA, .pp="gth-lda-valence"} for silicon.  The PP front door is LSDA, so \c model must be LDA:
    //! any other is refused at \c Resolve (it would be ignored -- a silent substitution).
    //! See doc/MolecularPseudopotentialPlan.md and doc/CalcNetworkReview.md §4.
    std::string pp;
    //! Per-element valence overrides on top of \c pp's set: `(element, q)` replaces the set's q for that element (or adds the element).
    //! E.g. a semicore Na q9 run: {.pp="gth-lda-valence", .ppOverride={{"Na",9}}}.
    std::vector<std::pair<std::string,int>> ppOverride;
    //! SCF seed.  Default == auto: a DFT run seeds from SAD (superposition of atomic densities) when the atomic-density library has
    //! every element, else from the core guess, as HF/1-e do (the choice files a `policy` row, AutoSeed).  Set explicitly to override -- e.g. {.seed=SAD} drives HF's matrix-free-seed
    //! bootstrap (the iterator manufactures a D0 via a one-step LDA sibling).
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;
    //! \brief The ONE accuracy target eps (doc/Pins.md pin 31; V-FACADE-ACCURACY, user ruling 2026-11-03: every structure type carries
    //! it).  CONSUMED BY NO TERM of this facade yet -- the eps inventory (doc/Records/EpsilonInventory.md) wires terms to it one at a
    //! time -- and the run says so (banner + report: "consumed by: none yet"), never silently.
    AccuracyTarget accuracy = {};
};


//! \brief THE RUN SETUP of a molecular calculation: \c CalcOptions as asked, plus every decision resolved from them ONCE.
//!
//! Only geometry-INDEPENDENT decisions live here (doc/CalcNetworkReview.md §2): the pseudopotential set, the spin group, the
//! electron counts, the seed, the accelerator and its EMax gate.  The basis, the Hamiltonian and the electron configuration
//! belong to a geometry point and are built from this by \c MoleculeCalculation.  Built only by \c Resolve, held as
//! \c shared_ptr<const MoleculeSetup>.
//! The Becke XC mesh's ENTRY RUNG on a molecule (pin 31; doc/GeometryOptimizationPlan.md §14.1): an unstated mesh starts here and
//! climbs only up, so an anchor moves only where this rung misses eps.  Today's former default, kept so nothing moves that met eps.
inline constexpr int kMolecularMeshEntryRadial = 80;
inline constexpr int kMolecularMeshEntryDegree = 29;

struct MoleculeSetup
{
    //! \brief THE ONE RESOLVE PHASE, as an OUTCOME: the door for USER INPUT (the GUI binding, a command line).  Options that contradict
    //! each other or the structure -- a multiplicity whose parity disagrees with the electron count, an unknown accelerator name, an
    //! unstated basis or one that does not cover an element, a route the code does not realize -- come back as a \c SetupRefusal naming
    //! the policy that refused.  A broken invariant still throws.
    static Outcome<std::shared_ptr<const MoleculeSetup>, SetupRefusal> Admit(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc = {});
    //! \brief \c Admit, THROWING its refusal (a \c policy::Refused, a \c std::invalid_argument, with the refusal's sentence).  For CODE that
    //! states its own options -- tests, probes, a campaign driver -- where a refusal is the caller's programming error (CLAUDE.md).  Never for
    //! user input: the GUI binding and \c scfrun go through \c Admit.
    static std::shared_ptr<const MoleculeSetup> Resolve(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc = {});

    CalcOptions        opts;   //!< as asked, except \c spin: an open shell (multiplicity > 1) is promoted to Polarized
    AcceleratorOptions acc;    //!< as asked
    //! The resolved choice (\c opts.pp with \c opts.ppOverride); EMPTY = all-electron.
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
    //! THE RUN'S XC MESH: \c opts.mesh with an unstated axis at the entry rung (\c kMolecularMeshEntryRadial / \c kMolecularMeshEntryDegree), or the rung a pilot settled
    //! on (\c MoleculeCalculation::SizedSetup).  Every Hamiltonian of a point is built on this, never on \c opts.mesh.
    qcMesh::MeshParams mesh;
    bool meshStated   = false;   //!< the caller stated \c opts.mesh (either number): audited against eps, never raised (pin 31 rule 3)
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
