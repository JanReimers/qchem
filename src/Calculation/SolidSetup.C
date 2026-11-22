// File: Calculation/SolidSetup.C
//
// qchem::SolidSetup -- THE RUN SETUP of a periodic calculation (doc/CalcNetworkReview.md §2 and §7, step 6): the options as
// asked, and every decision RESOLVED from them ONCE -- the pseudopotential set, the imposition (the run policy's veto applied),
// the magnetic decoration, the XC quadrature (Auto resolved), the spin counts.  An immutable value, shared: a geometry point
// (a SolidCalculation) is built FROM it and changes none of it, which is what lets a geometry optimizer move the atoms under one
// fixed setup (step 8).  RunOptions lives here because the setup is what reads it; qchem.SolidCalculation re-exports it.
module;
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
export module qchem.SolidSetup;

import qchem.Lattice_3D;                      // Lattice_3D (crystal + BZ mesh)
import qchem.Structure;                       // Structure
import qchem.BasisSet;                        // Real_BS (the molecular orbital basis)
import qchem.BasisSet.Lattice.BasisSet;       // CellImages, GPWTolerances, ShellTrim
import qchem.BasisSet.Gaussian.Point.Factory; // BasisSetData (the registered basis files)
import qchem.Hamiltonian.Factory;             // HubbardManifold, VxcFit
import qchem.SCFAccelerator.Factory;          // Type, AcceleratorParams
import qchem.ChargeDensity.Seed;              // SeedStrategy
import qchem.Mesh;                            // qcMesh::MeshParams / UnitCellKind
import qchem.Mesh.XCPolicy;                   // qcMesh::XCMeshDefault (the BeckeXCParams recipe under Auto)
import qchem.RunPolicy;                       // RunPolicySpec
export import qchem.RunOptions;               // RunOptions: the one options tree (V-FACADE-PARITY P0b-2)
export import qchem.SetupDecisions;           // AcceleratorOptions (the accelerator block, every structure), the spin rule
export import qchem.SCFRecipe;                // SCFRecipe (the stage's material regime, beside its accelerator)
export import qchem.SetupRefusal;             // SetupRefusal + Outcome: the Admit door's failure half
export import qchem.Setup;                    // Setup: what every run setup holds (V-FACADE-PARITY P1)
export import qchem.Symmetry.Spin;            // SpinGroup -- the imposed spin subgroup, a RunOptions field (R6)
export import qchem.AuditFootprint;           // AuditFootprintCalibration (the audit RSS policy's expert knob)
export import qchem.Pseudopotential.Set;      // PseudopotentialSet -- the setup's resolved choice

export namespace qchem
{

//! \brief How to set up a periodic calculation.  Designated-initializer friendly:
//!     auto calc=SolidCalculation::Prepare(SolidSetup::Resolve(lattice, basis, {.Nelec=8, .species={{"Si",4}}})).TakeValue();
//!     auto answer=calc->Converge();
//! A Hubbard manifold spelled the way the literature quotes it: \a U in eV.  (0.147 Ha = 4.0 eV.)  The projector is
//! the DEFAULT one -- the pseudo-atom's own \f$l\f$ orbital, independent of the orbital basis (HubbardProjectorPlan
//! increment 5, 2026-10-25; atomic, i.e. QE's `atomic`).
inline Hamiltonian::HubbardManifold HubbardU(size_t site, int l, double U_eV)
{
    return Hamiltonian::HubbardManifold{site, l, U_eV/27.211386245988};
}
//! CP2K's LOWDIN manifold, NAMED: every shell of angular momentum \a l on the site, \f$T=S^{1/2}[:,M]\f$ (the block's
//! own columns).  The CP2K-parity route (the CP2K-gated +U tests and the CrossCode decks name it).
inline Hamiltonian::HubbardManifold HubbardU_SiteShells(size_t site, int l, double U_eV)
{
    Hamiltonian::HubbardManifold M{site, l, U_eV/27.211386245988};
    M.projector=Hamiltonian::HubbardProjectorKind::SiteShells;
    return M;
}
//! The pseudo-atomic manifold, Löwdin-orthogonalised in the TRUE overlap among every ortho-atomic manifold of the run
//! (QE's `ortho-atomic`); list the spectators (O 2p, O 2s, Mn 4s) at U=0 to reproduce QE's full set.
inline Hamiltonian::HubbardManifold HubbardU_OrthoAtomic(size_t site, int l, double U_eV)
{
    Hamiltonian::HubbardManifold M=HubbardU(site, l, U_eV);
    M.orthoAtomic=true;
    return M;
}


//! \brief REFUSES (\c policy::Refused, rows SiteSpins / StatedMagneticOrdering) unless \a o states its magnetic ordering where the run needs
//! one: a polarized run carrying a Hubbard manifold must give \c spins (\c greyImposition, the erasure control, is exempt), and a stated
//! \c spins must have one entry per site of \a st, each -1, 0 or +1.  User input, so the refusal reaches a caller as a \c SetupRefusal
//! through \c SolidSetup::Admit; \c deck::Vet calls it too, so \c scripts/vet-decks catches a deck before it runs.
//! THE SOLID'S VIEW of the options tree: \a asked with the solid's defaults filled (model LDA, pp gth-lda-valence, \c XCMeshDefault,
//! Cholesky; \c solid.grid.eps = \c accuracy.eps) and what a periodic run does not consume REFUSED (policy ForeignOption): the
//! `molecule` block, xalpha, a Hartree route, symmetryTol, a stated \c solid.grid.eps.
RunOptions SolidOptions(const RunOptions& asked);
//! A caller that edits ONE field of a periodic run's XC mesh starts from the solid's default (\c qcMesh::XCMeshDefault), filled on first use.
inline qcMesh::MeshParams& SolidXCMesh(RunOptions& o) { if (!o.xcMesh) o.xcMesh=qcMesh::XCMeshDefault(); return *o.xcMesh; }
void RequireStatedMagneticOrdering(const RunOptions& o, const Structure& st);
//! \brief REFUSES (\c policy::Refused, row SiteCharges) unless a stated \c o.solid.siteCharges has one entry per site of \a st, sums to zero
//! (formal charges conserve the cell's charge, as the derived ones do), and seeds an \c IonicSAD run.  Empty passes.  Called by the
//! resolve phase and by \c deck::Vet, as \c RequireStatedMagneticOrdering is.
void RequireConsistentSiteCharges(const RunOptions& o, const Structure& st);
//! \brief The electron count per spin channel, (nUp, nDown), for \a o.multiplicity and the RESOLVED electron count \a nel (decision 2).
//! REFUSES (\c policy::Refused, row SpinMultiplicity) an impossible request: a multiplicity whose parity disagrees with \a nel or that exceeds it -- in particular multiplicity 1
//! (the EXPLICIT two-channel singlet, nUp = nDown) on an odd \a nel, unless \c o.solid.spinsShareFermi lets one Fermi level set the channels
//! (a metal: the counts are only the split the occupations start from).  Called by the resolve phase and by \c deck::Resolve.
std::pair<int,int> SpinChannelCounts(const RunOptions& o, int nel);
//! \brief THE PERIODIC GGA ROUTE (GGA stage 3), asked once the XC mesh and the imposition are resolved (\a xcMesh, \a imposed): a GGA runs
//! on the FREE Becke δ point route, and every other periodic route REFUSES it by name (policy GGARoute).  A non-GGA passes silently.
//! \a choseBecke: the setup itself picked the Becke mesh for an unstated one (the row says so).  CLIENTS: the setup (\c Resolve, on what it
//! resolved) and \c deck::Vet (on the deck's stated mesh; an unstated one is what the setup will choose, Becke).
void RouteSolidGradient(const RunOptions& o, const qcMesh::MeshParams& xcMesh, bool imposed, bool choseBecke=false);

//! \brief HOW the molecular orbital basis of a GPW run is made at a geometry: the frozen RECIPE a geometry step rebuilds
//! (doc/CalcNetworkReview.md §3.2 item 2, §7).
//!
//! The trim (doc/Pins.md pin 22) is STATED, or decided ONCE by the vet-stage loop when the setup is resolved -- at the
//! REFERENCE geometry -- and then FROZEN into \c trim: re-vetting at every geometry would make the energy surface
//! discontinuous.  It used to be the deck's (\c deck::Run) and nothing else's, so a calculation built in code got no trim
//! at all; the recipe is its one home now.
struct GaussianBasisRecipe
{
    BasisSet::Gaussian::BasisSetData data = BasisSet::Gaussian::BasisSetData::VALENCE_LOWQ_SR;   //!< the registered basis file
    std::string                       file;   //!< a .bsd path (or a bare name in the shipped directory) used INSTEAD of \c data; empty = \c data
    BasisSet::Gaussian::ShellTrim     trim;   //!< shells read the file WITHOUT (element, l, exponent); empty = the file as written
    //! Decide \c trim by the vet-stage loop at the run's \c orthoTol on the reference geometry's full k-mesh -- EXCLUSIVE with a
    //! stated \c trim.  \c SolidSetup::Resolve runs it once and stores the frozen result (\c autoTrim false, \c trim filled).
    //! Needs an \c orthoTol > 0, which it trims at; \c orthoTol alone never trims (user ruling 2026-11-04).
    bool                              autoTrim = false;
    //! The angular form, the SAME option as the molecular \c CalcOptions::angular (user, 2026-11-12, structure-neutral options).
    //! \c Spherical wraps the MnD Cartesian lattice basis in the spherical lattice view (contaminant-free d/f).
    BasisSet::Gaussian::Angular       angular = BasisSet::Gaussian::Angular::Cartesian;
    //! The pair-loop tolerances \c Make BINDS on the basis it builds (D-ENV-TOL option A, 2026-10-31): set by \c SolidSetup::Resolve
    //! from \c RunOptions::tolerances (the one source -- a caller's value here is replaced), frozen with the trim so a geometry
    //! step's basis computes at the same setting.  The GPW evaluator checks them; it no longer sets them.
    BasisSet::Gaussian::GPWTolerances tolerances;

    //! The basis over \a cell's atoms, its pair-loop \c tolerances bound before anything can compute from it.  CHECKS the file's
    //! declared pseudopotential provenance against \a pp first (D-BASIS-PP): THROWS naming the file, the element and both q on a
    //! mismatch.  \a cell must outlive the basis.
    std::shared_ptr<const BasisSet::Real_BS> Make(const Structure& cell, const Pseudopotential::PseudopotentialSet& pp) const;
};

//! \brief THE RUN SETUP of a periodic calculation: \c RunOptions as asked, plus every decision resolved from them ONCE.
//!
//! WHAT IS HERE, AND WHAT IS NOT (doc/CalcNetworkReview.md §2).  A run splits into three lifetimes: the SETUP (the choices and
//! what is resolved from them), the GEOMETRY POINT (the Bloch basis, S, the Hamiltonian's geometry-keyed terms, the electron
//! configuration over the basis's k-blocks) and the SCF STATE.  Only the first lives here.  Every member is a decision that does
//! not change when the atoms move -- the set, the imposition, the decoration, the XC quadrature (sized from exponents and the cell,
//! never from positions), the spin counts -- so a geometry optimizer can rebuild the point at new positions under the SAME setup.
//! A field that would depend on positions does not belong here.
//!
//! Built only by \c Resolve, held as \c shared_ptr<const SolidSetup>: immutable once made.
struct SolidSetup : Setup
{
    //! \brief THE ONE RESOLVE PHASE (doc/Pins.md pin 11), as an OUTCOME: the door for USER INPUT (a deck, a GUI form, a command line).
    //! Installs \a opts.solid.policy as the process run policy FIRST (every factory after it consults the policy), then decides.  Options
    //! that contradict each other or the structure -- an unstated magnetic ordering, a multiplicity whose parity disagrees with
    //! \c Nelec, a stated Nelec the species disagree with, \c numBands too small, a route the code does not realize, a basis that does
    //! not cover the cell -- come back as a \c SetupRefusal naming the policy that refused (its row is in the run report too).  A
    //! broken invariant still throws.  \a mol: the molecular orbital basis at \a lat's geometry -- the XC-mesh sharpness reads its
    //! α_max, and a GPW point is built over it; a plane-wave run ignores it.
    static Outcome<std::shared_ptr<const SolidSetup>, SetupRefusal> Admit(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                                     const RunOptions& opts,
                                                     const AcceleratorOptions& acc = {});
    //! \c Admit from a basis RECIPE: the molecular basis is made at \a lat's geometry by \a recipe (its provenance checked against
    //! the run's set), after the vet-stage loop -- if the recipe asks for it -- has decided the trim, which is then frozen in
    //! \c recipe.  A setup with a recipe can rebuild its basis at moved atoms (a geometry step); one from a pre-built \a mol cannot.
    static Outcome<std::shared_ptr<const SolidSetup>, SetupRefusal> Admit(const Lattice_3D& lat, const GaussianBasisRecipe& recipe,
                                                     const RunOptions& opts,
                                                     const AcceleratorOptions& acc = {});
    //! \brief \c Admit, THROWING its refusal (a \c policy::Refused, a \c std::invalid_argument, with the refusal's sentence).  For CODE
    //! that states its own options -- the tests, the probes (\c gpwprobe), a campaign driver: there a refusal is the caller's programming
    //! error, which no one at run time can act on, so it throws (CLAUDE.md).  Never for user input: a deck, the GUI binding or a command
    //! line goes through \c Admit and reports the refusal.
    static std::shared_ptr<const SolidSetup> Resolve(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                                     const RunOptions& opts,
                                                     const AcceleratorOptions& acc = {});
    //! \c Resolve from a basis RECIPE (\c Admit's recipe door, throwing; for code that states its own options).
    static std::shared_ptr<const SolidSetup> Resolve(const Lattice_3D& lat, const GaussianBasisRecipe& recipe,
                                                     const RunOptions& opts,
                                                     const AcceleratorOptions& acc = {});
    //! THE PRE-RUN PP / FUNCTIONAL CHECK (user 2026-11-12: track, warn and guide, mostly in the GUI): every element of the reference
    //! cell, its pseudopotential's fit (CP2K label + libxc ids) against \c opts.model, match yes/no, and the model to state to match it.  What the
    //! run will warn about, answered without running; a GUI reads it here and never recomputes it.  THROWS on a fit with no libxc mapping.
    std::vector<Hamiltonian::PPFunctionalCheck> PPFunctionalChecks() const;

    //! The shared decisions are the base's (\c Setup): \c opts (\c ortho: a plane-wave run's default Cholesky becomes Orthonormal; \c seed
    //! and \c Nelec resolved; \c solid.tolerances' raster demands, see \c rastersResolved), \c acc and \c accParams (resolved over the GPW
    //! production recipe, the solid's per-kind defaults; a stage derives its own from it, a restored stage's Ladder rung), \c accelerator
    //! (\c opts.accelerator: \c Auto lets each stage decide), \c pp (always set: \c opts.pp with the \c opts.species overrides), \c structure
    //! (the reference cell, the one the reference basis is built over, held so it outlives that basis), \c Ne / \c nUp / \c nDown, \c seedKey
    //! (\c pp's functional and, per element, its q: the seed is the SAME pseudo-atom as the terms), \c xcMesh (\c Auto RESOLVED, the facade's
    //! decision 1), \c imposed (the caller's \c imposeSymmetry AND the run policy's permission: CP2K parity vetoes the capability) and \c spins
    //! (the imposition's MAGNETIC DECORATION, what makes the group Shubnikov: stated, or -- for an imposed polarized run without +U only --
    //! derived from the seed; EMPTY when the run imposes nothing or the grey group).
    //! The REFERENCE lattice the setup was resolved at (the trim vetted, the imposed group detected); \c structure is its cell.
    Lattice_3D                               lattice;
    std::shared_ptr<const BasisSet::Real_BS> mol;    //!< the molecular orbital basis at the REFERENCE geometry
    //! How \c mol was made, the trim FROZEN; empty when the caller handed a pre-built basis (that setup cannot follow moved atoms).
    std::optional<GaussianBasisRecipe>       recipe;
    bool planeWave = false;   //!< \c opts.solid.orbitalBasis is PlaneWave
    bool ggaChoseBecke = false;   //!< a GGA with the XC mesh unstated: the setup chose the Becke delta route (policy GGARoute, GGA stage 3)
    //! \brief THE ANSATZ, decided ONCE from k and q (user ruling 2026-11-18): a TRIM block may run REAL (\c GPWParams::hamPreservesReal;
    //! a non-TRIM block is complex regardless).  False exactly when \c opts.solid.responseQMeshes plans a response.  The run's ONE source
    //! of the decision; there is no user knob (\c forceComplex retired 2026-11-19).  A test comparing real against complex on purpose sets
    //! it on a COPY of a resolved setup (tests may cheat).
    bool realTRIMBlocks = true;
    std::vector<int> decoration;   //!< \c spins, or all zeros: the per-site spins the run was built with (the +U site groups read it)
    //! \brief THE GROUP A GEOMETRY STEP MUST KEEP, read at \a at: the space group detected there (\c GetSpaceGroup's \f$\{W|\tau\}\f$, or
    //! the Shubnikov group under \c spins) when the setup is \c imposed, P1 (the identity alone) when it is free.  ONE rule for its two
    //! clients, so they cannot drift: \c SolidSession::ConvergeAt refuses a geometry that keeps fewer ops than the reference, and
    //! \c MakeSolidGeometryRun builds the coordinate map that keeps exactly these.  It is the DETECTED group, not the k-mesh intersection
    //! the density is star-averaged under: the positions keep the whole group, which is what the session checks.
    std::vector<Symmetry::SymOp> KeptOps(const Lattice_3D& at) const;
    //! \brief THE SHELL WHAT-IF of this setup's basis (the user 2026-11-11, the GUI's conditioning table): its recipe with the trim
    //! CLEARED, built at the reference geometry and analysed on the full k-mesh (\c ShellConditioning).  THROWS for a plane-wave setup
    //! (no shells) and for one built from a pre-built basis (no recipe to clear the trim of).
    std::unique_ptr<BasisSet::Lattice::ShellConditioning> MakeShellConditioning() const;
    //! \brief What the vet-stage trim removes from the UNTRIMMED recipe at \a orthoTol -- the same loop \c Resolve runs for autoTrim,
    //! for any tolerance (the what-if's reference: "autoTrim would remove these").  Writes its passes one level in, as in \c Resolve.
    BasisSet::Gaussian::ShellTrim MakeVetTrim(double orthoTol) const;
    //! \brief The GPW raster demands in \c opts.solid.tolerances are FINAL (C2 step 3): stated (accuracy 0), derived from a stated calibration,
    //! or settled by a pilot audit (\c SolidCalculation::SizedSetup).  false = this setup's first SCF is the PILOT: it runs on the ball's
    //! own rasters and \c SolidCalculation sizes them from what it delivered.  A plane-wave run has no raster: true.
    bool rastersResolved = true;
    //! \brief WHO SET the raster demands this setup runs (the report's \c grids.rasterSource): the user (\c Stated, accuracy 0), a stated
    //! calibration at Resolve (\c Calibration), the pilot audit (\c Pilot: on the ball's own rasters until it raises them), a saved state
    //! sized for this eps (\c SavedState), or the start raised for an ill-conditioned overlap (\c Conditioning).  Every place that writes
    //! \c opts.solid.tolerances' raster demands writes this beside it; a reader of the grids can then tell a deck value from a measured one.
    enum class RasterSource { Stated, Calibration, Pilot, SavedState, Conditioning };
    RasterSource rasterSource = RasterSource::Stated;

private:
    SolidSetup(Pseudopotential::PseudopotentialSet p, const Lattice_3D& lat) : lattice(lat)
    {
        pp=std::move(p);
        structure=lat.GetStructure();
    }
    //! The decisions, shared by both \c Resolve overloads (the policy is already installed).
    static std::shared_ptr<SolidSetup> Decide(const Lattice_3D&, std::shared_ptr<const BasisSet::Real_BS> mol,
                                              const RunOptions&, const AcceleratorOptions&);
};

} //export namespace qchem
