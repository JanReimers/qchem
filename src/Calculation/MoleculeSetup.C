// File: Calculation/MoleculeSetup.C
//
// qchem::MoleculeSetup -- THE RUN SETUP of a molecular calculation (doc/CalcNetworkReview.md §2 and §7, step 7), the twin of
// qchem::SolidSetup: the options as asked, and every decision RESOLVED from them ONCE -- the pseudopotential set (or
// all-electron), the spin group (an open shell promoted to polarized), the electron counts, the seed, the accelerator.  An
// immutable value, shared: a geometry point (a Calculation) is built FROM it and changes none of it.  CalcOptions lives here
// because the setup is what reads it; qchem.Calculation re-exports it.
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
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type (the resolved accelerator)
export import qchem.Pseudopotential.Set;   // PseudopotentialSet -- the setup's resolved choice

export namespace qchem
{

using Hamiltonian::Model;          // {E1, HF, DE1, DHF}

//! Orbital-integral engine: the in-house MnD recursion (default) or the libcint foreign engine.
enum class Engine  { MnD, LibCint };
//! Angular basis convention: Cartesian Gaussians (default) or real solid-harmonic (spherical).
enum class Angular { Cartesian, Spherical };

//! How to set up the calculation.  Designated-initializer friendly:
//!     Calculation calc(MoleculeSetup::Resolve(water, {.basis="dzvp"}));
//! `model` selects the Hamiltonian; pass 1 wires the non-DFT factory path (HF is the default).
struct CalcOptions
{
    std::string basis = "sto-3g";
    Model       model = Model::HF;   //!< HF (default) | Xalpha | LDA | E1/DE1/DHF (test-only)
    SpinGroup         spin  = SpinGroup::UnPolarized;   //!< the imposed spin subgroup (V1.37)
    //! Spin multiplicity 2S+1.  0 (default) = minimal spin: closed-shell singlet for even Ne, doublet for
    //! odd -- the historical behaviour.  Set explicitly for an open shell: 3 = triplet, 2 = doublet, ...
    //! The facade converts it to (nUp,nDown) [nUp-nDown = 2S = multiplicity-1, nUp+nDown = Ne] and PROMOTES
    //! the calculation to SpinGroup::Polarized when 2S>0 (unrestricted open shell needs distinct up/down densities).
    //! A multiplicity whose parity disagrees with Ne (e.g. a singlet for odd Ne) is rejected.
    int         multiplicity = 0;
    //! Basis construction variants (threaded into BasisSet::Gaussian::Factory).  Defaults reproduce
    //! today's behaviour (in-house MnD, Cartesian).  angular==Spherical + symmetry is rejected until
    //! the Spherical SALC track (doc/OldPlans/SphericalSALCPlan.md) lands -- the SALC builder needs Cartesian PGData.
    Engine      engine  = Engine::MnD;
    Angular     angular = Angular::Cartesian;
    //! DFT-only knobs (ignored when model is HF/1-e/Dirac).  xalpha: the Slater exchange parameter, used
    //! only by model==Xalpha.  mesh: the numerical XC integration grid; a designated initializer overrides just
    //! the resolution you care about, e.g. {.nRadial=50}.
    //! The default is CONVERGED, MHL 80 / Lebedev 29 (302 directions), measured 2026-10-31 against MHL 150 / Lebedev 35:
    //! per atom, water LDA 1 μHa, N₂ Xalpha 3 μHa, F₂ GTH q7 27 μHa.  It was MHL 30 / Lebedev 5 (12 directions), proven on
    //! HF only (no mesh there) and 20 mHa/atom LOW on water LDA, 0.52 Ha/atom on F₂ -- the error is ANGULAR (30/29 is 0.7 mHa
    //! from converged on water, 80/5 still 65 mHa) and NOT variational, so a coarse mesh reads as a better energy.  Cost:
    //! ~20x the 30/5 SCF wall on water (0.14 s -> 2.8 s single-threaded); 100/35 buys ~1 μHa/atom more for 1.75x again.
    double      xalpha = 0.7;
    qcMesh::MeshParams mesh = {.radial  = qcMesh::RadialKind::MHL,   .nRadial   = 80,
                              .mhl_m    = 3,                         .mhl_alpha = 2.0,
                              .angular  = qcMesh::AngularKind::Lebedev, .angularDegree = 29,   // the 302-direction rule
                              .beckeOrder = 2};
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
    //! {.basis="sipp", .pp="gth-lda-valence"} for silicon.  Ignores \c model (the PP front door is always LSDA today).
    //! See doc/MolecularPseudopotentialPlan.md and doc/CalcNetworkReview.md §4.
    std::string pp;
    //! Per-element valence overrides on top of \c pp's set: `(element, q)` replaces the set's q for that element (or adds the element).
    //! E.g. a semicore Na q9 run: {.pp="gth-lda-valence", .ppOverride={{"Na",9}}}.
    std::vector<std::pair<std::string,int>> ppOverride;
    //! SCF seed.  Default == auto: DFT seeds from SAD (superposition of atomic densities), HF/1-e take
    //! the core guess.  Set explicitly to override -- e.g. {.seed=SAD} drives HF's matrix-free-seed
    //! bootstrap (the iterator manufactures a D0 via a one-step LDA sibling).
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;
};

//! DIIS accelerator knobs.  These were stringly-typed json keys grepped out of tests; here they
//! are documented fields with defaults.  eMax<=0 means "auto" -- the Z-scaled heuristic the
//! molecular recipe has always used.  The long tail still rides through the json factory directly.
struct AcceleratorOptions
{
    int    nProj = 4;
    double eMax  = 0.0;     //!< <=0 => Z*Z*0.1/32 (the proven molecular default)
    double eMin  = 1e-7;
    double svTol = 5e-9;
    //! SCF accelerator: "DIIS" (default) | "GDM" (direct minimisation, robust for hard/open-shell SCF) |
    //! "Ladder" (auto None->DIIS->GDM).  A hard pseudopotential (open-shell atom, diffuse valence) that
    //! limit-cycles under DIIS often converges under GDM or Ladder.
    std::string type = "DIIS";
};

//! \brief THE RUN SETUP of a molecular calculation: \c CalcOptions as asked, plus every decision resolved from them ONCE.
//!
//! Only geometry-INDEPENDENT decisions live here (doc/CalcNetworkReview.md §2): the pseudopotential set, the spin group, the
//! electron counts, the seed, the accelerator and its EMax gate.  The basis, the Hamiltonian and the electron configuration
//! belong to a geometry point and are built from this by \c Calculation.  Built only by \c Resolve, held as
//! \c shared_ptr<const MoleculeSetup>.
struct MoleculeSetup
{
    //! \brief THE ONE RESOLVE PHASE.  THROWS on options that contradict each other or the structure: a multiplicity whose parity
    //! disagrees with the electron count, an unknown accelerator name, a pseudopotential set that does not cover an element.
    static std::shared_ptr<const MoleculeSetup> Resolve(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc = {});

    CalcOptions        opts;   //!< as asked, except \c spin: an open shell (multiplicity > 1) is promoted to Polarized
    AcceleratorOptions acc;    //!< as asked
    //! The resolved choice (\c opts.pp with \c opts.ppOverride); EMPTY = all-electron.
    std::optional<Pseudopotential::PseudopotentialSet> pp;
    //! The REFERENCE geometry, deep-copied (a periodic UnitCell stays periodic), with the TRUE atomic charges.
    std::shared_ptr<const Structure> reference;
    int    Ne = 0;             //!< electrons: the VALENCE count under a pseudopotential
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;   //!< \c Default resolved per model
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::DIIS;   //!< \c acc.type, by name
    double eMax = 0.0;         //!< \c acc.eMax, or the model's default when it was left <= 0

    //! The structure a geometry point is built over: \a st's atoms, each carrying its VALENCE charge Z - Zion under a
    //! pseudopotential (so the electron count is the valence count), or \a st as is when all-electron.  Deep-copied.
    std::shared_ptr<const Structure> PointStructure(const Structure& st) const;
};

} // namespace qchem
