// File: RunOptions.C  THE ONE OPTIONS TREE every structure kind takes (V-FACADE-PARITY P0b-2, 2026-11-17; was SolidCalcOptions,
// CalcOptions and AtomCalcOptions).  Flat at the top -- the schema-2 deck's top-level keys, one to one -- with one sub-block per
// geometry (`molecule`, `solid`; `atom` in step iii) and the accelerator block.  A kind REFUSES what it does not consume (a stated
// field or block another kind owns), never ignores it.  Defaults that differ by kind are UNSTATED here (std::optional) and each kind's
// setup fills its own, filing a `policy` row (OpenWork "Per-kind defaults: harmonize where the physics allows").
module;
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>
export module qchem.RunOptions;
export import qchem.Symmetry.Spin;              // SpinGroup
export import qchem.Hamiltonian.Factory;        // Model, VxcFit, XCρ, CoulombRoute, HubbardManifold
export import qchem.Mesh;                       // qcMesh::MeshParams
export import qchem.Mesh.XCPolicy;              // qcMesh::BeckeLadderCalibration, XCMeshDefault
export import qchem.BasisSet.Lattice.BasisSet;  // CellImages, GPWTolerances, GridAccuracy
export import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSet::Gaussian::Angular
export import qchem.RunPolicy;                  // RunPolicySpec
export import qchem.AuditFootprint;             // AuditFootprintCalibration
export import qchem.PhiTablesPolicy;            // policy::XCOrbitalTables (the δ fit's orbital tables)
export import qchem.ChargeDensity.Seed;         // SeedStrategy
export import qchem.LASolver;                   // qchem::Ortho
export import qchem.SCFAccelerator.Factory;     // SCFAccelerators::Type, ScheduleSignal
export import qchem.SCFRecipe;                  // SCFRecipe
export import qchem.SCFParams;                  // SCFParams (R10: the run's SCF settings live in the options)
export import qchem.Accuracy;                   // AccuracyTarget
export import qchem.BasisSet.Radial.Factory;    // BasisSet::Radial::Type, BasisSetAccuracy (the atom's exponent pool)
import qchem.Types;                             // rvec3_t

export namespace qchem
{
//! Orbital-integral engine: the in-house MnD recursion (default) or the libcint foreign engine.
using BasisSet::Gaussian::Engine;   // the library's own axis (one enum, as Angular)
//! Angular basis convention: Cartesian Gaussians (default) or real solid-harmonic (spherical).  THE library's own axis, shared by the
//! molecular \c CalcOptions::angular and the periodic \c GaussianBasisRecipe::angular (user, 2026-11-12: "we should strive to make
//! seeding and spherical/Cartesian options structure neutral"), so one option means one thing on every structure.
using BasisSet::Gaussian::Angular;

using AtomType = BasisSet::Radial::Type;    //!< the atom's basis family {Slater, Gaussian, BSpline6, BSpliner6, Gaussian_RKB, Slater_RKB}
using BasisSet::Radial::BasisSetAccuracy;   //!< the atom's preset exponent pool {Low, Medium, High}

//! \brief THE ACCELERATOR BLOCK of the options tree, every structure (V-FACADE-PARITY P0b-1, 2026-11-17; was the molecule's flat
//! \c AcceleratorOptions and the solid's \c SolidAcceleratorOptions).  Every knob is OPTIONAL: unstated takes the KIND's default, which
//! the kind's setup fills (\c ResolveAcceleratorParams) -- the molecule's and the atom's Z-scaled DIIS gate, Nproj 4 and a Ladder GDM rung
//! that always steps; the solid's GPW production recipe (\c SCFAccelerators::AcceleratorParams' member defaults).  Those differences
//! are per-kind DEFAULTS, logged for harmonizing (OpenWork "Per-kind defaults").  The knobs mirror \c DIISParams / \c GDMParams and
//! the Ladder's triggers one to one, so a knob added there has one spelling here.
struct AcceleratorOptions
{
    //! (WHICH accelerator, and the SCF recipe, are \c RunOptions::accelerator / \c recipe -- their ONE home, every kind, V-FACADE-PARITY
    //! P2; until then a molecule and an atom named theirs here, as text.  This block holds only how each one is TUNED.)
    struct DIIS
    {
        std::optional<size_t> nProj;          //!< history depth
        //! DIIS starts extrapolating once [F,D] < this.  A molecule's / atom's unstated gate is Z-scaled (\c DIISEngageGate).
        std::optional<double> engageBelowFD;
        //! The ABSOLUTE floor on svMin(B) (OFF for a molecule since 2026-11-11: B scales as |[F,D]|^2, so it purged the history once
        //! [F,D] < ~7e-5) and the RELATIVE prune svMin(B) < svTolRel * max_i B_ii, which guards a genuinely dependent history (doc/Records/
        //! OpenWork_History5.md, "The molecular DIIS switches itself off").
        std::optional<double> svTol, svTolRel;
    };
    struct GDM
    {
        //! Engage the geodesic step once [F,D] < this.  A molecule's / atom's UNSTATED standalone-GDM gate is the DIIS gate (stated or
        //! Z-scaled), as the flat \c eMax always was; the Ladder's GDM rung keeps its own (1e10 there: it always steps).
        std::optional<double> engageBelowFD;
        std::optional<double> trust;          //!< the largest geodesic rotation (radians) per step
        std::optional<double> trustBackoff, trustMin, pcFloor, aufbauTol;
        std::optional<int>    aufbauMaxRepicks;
    };
    DIIS diis;
    GDM  gdm;      //!< standalone GDM
    //! LADDER (DIIS -> GDM): its GDM rung and the hand-off triggers (\c SCFAcceleratorLadder): energy still moving by more than
    //! \c eThresh for \c stall non-improving steps (STALL, \c floor its low backstop), or the signal below \c switchAt (TAIL; 0 = none;
    //! 1e-6 by default, the CLAUDE.md rule).  \c signal: |ΔE/E| (\c EnergyChange) or [F,D] (\c Error).
    struct Ladder
    {
        GDM gdm;
        std::optional<double> eThresh;
        std::optional<int>    stall;
        std::optional<double> floor, switchAt;
        std::optional<SCFAccelerators::ScheduleSignal> signal;
        bool startOnGDM = false;   //!< the seed is a RESTORED state: start on the GDM rung (a restart's own decision, never a user's)
    } ladder;
};

struct RunOptions
{
    //! \name The system
    //!@{
    //! VALENCE electrons per cell (Zion, not Z).  0 = DERIVED: \f$\sum_{\rm atoms} q\f$ of the run's species minus \c charge.  A stated value
    //! must equal that, or \c SolidSetup::Resolve REFUSES the run naming both numbers (2026-10-31: a stated Na q9 once ran NaF with the
    //! material's 8 electrons instead of 16 -- a cell charged +8, silently).
    int Nelec = 0;
    //! The cell's stated NET CHARGE (electrons removed; -1 = one extra electron).  A charged cell is legal only when stated here.
    int charge = 0;
    //! Spin multiplicity 2S+1 -- S ONLY (R6, 2026-11-17: the same two knobs as the molecule, \c ResolveSpinGroup).  0 (default) =
    //! minimal spin: a closed shell for even \c Nelec, a doublet for odd.  1 = the singlet (nUp = nDn).  > 1 forces \c spin Polarized.
    int multiplicity = 0;
    //! The imposed spin subgroup: \c Polarized = the SPIN-NATIVE two-channel run, which an AFM or broken-symmetry singlet needs
    //! (multiplicity 1 with \c Polarized is the explicit two-channel singlet, the ζ=0 cross-check of the polarized machinery against the
    //! unpolarized anchors).  A multiplicity > 1 forces it.  Until R6 the solid had no such knob and read multiplicity >= 1 as polarized.
    SpinGroup spin = SpinGroup::UnPolarized;
    //! The Hamiltonian MODEL.  UNSTATED = the kind's default: a molecule's and an atom's HF, a solid's LDA (a \c policy row).  A solid
    //! takes LDA (Dirac + VWN5) or PADE (the GTH-PADE PPs' own, user ruling 2026-11-12), STATED, never switched for you: a PP fitted
    //! with another functional is a warn-and-run row (\c PPFittedFunctional).
    std::optional<Hamiltonian::Model> model;
    //! The Slater exchange parameter (model Xalpha only).
    double xalpha = 0.7;
    //! The pseudopotential SET (src/Pseudopotential/Data/pseudopotentials.json; doc/CalcNetworkReview.md §4): per element the GTH
    //! valence variant and the functional it was generated for, chosen ONCE.  UNSTATED = the kind's default: a molecule all-electron,
    //! a solid "gth-lda-valence" (a deck takes its material's set).  An EMPTY string states all-electron.
    std::optional<std::string> pp;
    //! Per-element valence overrides on top of \c pp's set, as `(element, valence)`: replaces the set's q for that element (a
    //! semicore Mn q15 run) or adds it.  Empty (default) = the set as is.  Multi-species is the norm for a compound; a per-Z
    //! router PP is built over the elements of the structure so each atom gets its own GTH entry.
    std::vector<std::pair<std::string,int>> species;
    //!@}

    //! \name XC and the Hartree route
    //!@{
    //! The XC quadrature mesh.  UNSTATED = the kind's default: a molecule's Becke mesh climbed to eps from its entry rung (80 / 29),
    //! an atom's radial-only set climbed from 50, a solid's \c qcMesh::XCMeshDefault (cellKind Auto: the cheaper of a uniform and a
    //! Becke mesh).  A stated mesh is honoured (audited, never raised).
    std::optional<qcMesh::MeshParams> xcMesh;
    //! WHICH fit basis represents v_xc -- ORTHOGONAL to the grid above.  \c Auto picks Delta whenever the
    //! plane-wave fit cannot serve (Becke grid, or a polarized run).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
    //! WHICH density XC sees (the deck's `xc.rho`, user ruling 2026-11-09): \c Exact = rho[D], \c Fit = the Vxc fit's projection
    //! P rho.  The two coincide on every cell a solid realizes Fit on: the delta fit, and the plane-wave ball over plane-wave orbitals
    //! (>= 4 Ecut holds every pair); a stated \c Fit is reported there.  Over GPW orbitals no ball holds rho[D], so Fit on the plane-wave
    //! fit is REFUSED by name, never run as Exact (\c qchem.XCRhoPolicy, \c ReportSolidXCRho).
    Hamiltonian::XCρ xcRho = Hamiltonian::XCρ::Exact;
    //! The Hartree route (2026-11-06): by default the Gaussian Coulomb fit; \c .fit=CoulombFit::Delta puts \f$\rho\f$ on a mesh, with a
    //! first-order energy and the mesh's charge error reported.  See \c Hamiltonian::CoulombRoute.
    Hamiltonian::CoulombRoute coulomb = {};
    //! \brief THE ONE accuracy target eps (doc/Pins.md pin 31; V-FACADE-ACCURACY): Ha per atom on an energy run, per matrix element on
    //! a geometry run.  A solid's grid calibrations ride \c solid.grid; its eps is THIS one.
    AccuracyTarget accuracy = {};
    //!@}

    //! \name Convergence machinery
    //!@{
    //! \brief THE ACCELERATOR TYPE, every kind: its ONE home (V-FACADE-PARITY P2; a molecule and an atom used to name theirs as text in
    //! \c AcceleratorOptions, the solid here).  \c Auto (default) is the KIND's policy.  A molecule and an atom run \c DIIS under it and
    //! REFUSE \c Null (no density mixer to pair it with); a text name (a GUI form, the atom's json escape hatch) becomes a type through
    //! \c AcceleratorTypeByName, which refuses an unknown one.  A SOLID resolves \c Auto per stage, PAIRED with an \c Auto mixer
    //! (OpenWork §3 "Kerker/Pulay as the DEFAULT mixer", measured 2026-10-22: Pulay on ρ beside DIIS on F fights -- free Si
    //! k222 stalls, Si Γ diverges -- while Kerker/Pulay with NO Fock extrapolator beat linear + DIIS on every probe; the
    //! DEFAULT since 2026-10-26): a seeded or carried stage -> \c Null + Pulay on D (no Kerker, no fit: pin 28); a RESTORED
    //! stage (Restart, same Hamiltonian) at kT = 0 -> the Ladder on its GDM rung, a restored smeared one -> \c DIIS, either
    //! beside a linear step capped at 0.3.
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::Auto;
    //! The SCF RECIPE by material regime (qchem.SCFRecipe; deck key `scf.recipe`, doc/SCFRobustnessPlan.md R4.1 as ruled for B5):
    //! \c Auto (default) = the policy above, named and printed; \c MagneticInsulator / \c Metal fill Kerker+Pulay on a solid's cold or
    //! carried stage and REFUSE kT = 0.  A molecule takes \c Auto or \c Insulator only (ruling (c)): the other two are REFUSED there (their
    //! Kerker+Pulay mixes in G space).  A stated accelerator or mixer knob beside it wins (a deviation, printed).
    SCFRecipe recipe = SCFRecipe::Auto;
    //! Impose the detected space group (IBZ k-fold + per-iteration density star-average + site-adapted
    //! Becke mesh).  OPT-IN (default \c false): it roughly doubles the XC grid, and see the warning below.
    //!
    //! \note DELIBERATELY DIVERGES from the test harness's \c GpwOptions, which defaults this to \c true --
    //! while its own comment says "OPT-IN per the §3 pin".  That comment states the intent and the value
    //! contradicts it (doc/CleanupCandidates.md V1.30).  This facade follows the stated intent, because
    //! default-ON is the more dangerous way to be wrong: a symmetry imposition you did not ask for is
    //! invisible in the result, whereas a missing one merely costs time.
    //! \warning Do NOT combine with a polarized ANTIFERROMAGNETIC density until the fold is spin-aware:
    //! the op set comes from the CHEMICAL space group, which contains sublattice-exchanging operations,
    //! so star-averaging each channel under it averages the magnetic sublattices together and silently
    //! demagnetises the run (doc/CleanupCandidates.md V1.28).
    bool imposeSymmetry = false;
    double symmetryTol = 1e-4;   //!< geometry tolerance for point-group detection (a molecule's imposition)
    //! The SCF seed.  \c Default = THE Auto rule both setups share (\c qchem.SeedPolicy, user 2026-11-12): on a lattice that is
    //! \c IonicSAD (V2.2, 2026-09-14), a superposition of the library's per-species
    //! valence densities at their formal charges (Na+ F-, Mn2+ O2-; an elemental solid is the neutral SAD).
    //! \c Uniform (\f$\rho=N/V\f$, the old default) is an explicit OPT-IN, because it has a STABLE WRONG
    //! BASIN for electron-sparse systems: the Na-doublet campaign converged a lone-electron doublet 72 mHa
    //! high from it with every health metric green (doc/CleanupCandidates.md V2.2).  The molecular facade
    //! already defaults DFT to SAD for the same reason.
    //! \note IonicSAD needs every species in \c atomic_valence_densities.json; under \c Default a missing species takes the core guess
    //! with a WARN (ruling 2026-11-12); a STATED IonicSAD with a missing species still THROWS at seed time.  Resolved into the setup's
    //! \c opts.seed, so everything downstream reads the seed that runs.
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;
    //! The orthogonalisation of the generalised eigenproblem.  UNSTATED = the kind's default: a solid's Cholesky, a molecule's and an
    //! atom's Auto.
    std::optional<qchem::Ortho> ortho;
    //! THE NEAR-DEPENDENCE FLOOR on S's eigenvalues: what the vet-stage trim trims at (basis autoTrim, and a geometry step's re-vet), and
    //! the floor a TRUNCATING ortho solver drops directions under -- one quantity (the vet's metric is the truncating solver's own drop
    //! count at it), two readers.  The solver reads it only when it truncates (\c SolverTruncation); the solid default, Cholesky, does
    //! not, so there the floor is the trim's alone (OpenWork §3 "solid.orthoTol has two meanings", user feedback 36, 2026-11-12).
    double       orthoTol = 0.0;
    //! \brief The tolerance handed to the ortho SOLVER: \c orthoTol when \c ortho truncates (Eigen, SVD, CholeskyPivoted), else 0 --
    //! Cholesky, Auto (it finds its own gap) and Orthonormal take none.
    double SolverTruncation() const
    { return (ortho==qchem::Eigen || ortho==qchem::SVD || ortho==qchem::CholeskyPivoted) ? orthoTol : 0.0; }
    //! \brief THE MAGNETIC DECORATION: one entry per cell site, +1 / -1 / 0, naming which sublattice each
    //! site belongs to.  It is what makes an imposition SHUBNIKOV rather than merely spatial -- without it
    //! an "imposed AFM" run star-averages under the spatial group and the order is ERASED.
    //!
    //! EMPTY (the default) = DERIVE it from the seed, which reproduces the historical driver behaviour and
    //! is right for the common case.  SET IT to state a specific ordering: a cell usually admits several
    //! (AFM-I vs AFM-II vs ferrimagnetic), the derived one is only the seed's guess, and comparing orderings
    //! is the whole point of a magnetic campaign.  It is the run's ONE source of the ordering (2026-11-12): the polarized seed
    //! reads it too (+1 majority up, -1 flipped, 0 no moment), overriding the structure's per-atom spins, which are only the
    //! default when it is absent.  The imposition reads it only when \c imposeSymmetry.
    //! \note \c greyImposition forces the empty decoration regardless -- that arm is the erasure control.
    std::vector<int> spins = {};
    //! Names this run in its own console output -- the fingerprint, the stage banners, the order trace.
    std::string label = "gpw";
    //! \brief SAVE THE STATE HERE after every \c Converge -- each anneal stage, each ACBN0 outer step -- converged
    //! or not (the file records which).  EMPTY (default) = never.  The unattended-run safety net (CK-1): a crash
    //! costs at most the stage in flight, and \c SolidCalculation::Restart picks up from the last one.  A write
    //! failure here WARNS and the run continues -- a checkpoint must never cost the computation it protects
    //! (\c SaveState called directly throws as usual).
    std::string saveStateTo;
    //! THE RUN'S SCF SETTINGS (R10, 2026-11-17): carried here since P0b-2; \c Converge() reads them from the setup in P3.
    SCFParams scf = {};
    //!@}

    //! \brief THE ATOM-ONLY OPTIONS: the radial exponent pool (the atom's basis), its functional override, the pseudo-atom, the stated
    //! configuration, the accelerator-sweep escape hatch.  Any other kind REFUSES a stated one.  (The pseudo-atom's bool + valence and the
    //! json escape hatch stay as they were: harmonizing them with the top-level \c pp / accelerator block is OpenWork row "Atom options".)
    struct Atom
    {
        AtomType         type     = AtomType::Slater;            //!< atomic basis family
        BasisSetAccuracy basisPool = BasisSetAccuracy::Medium;   //!< the preset exponent POOL (a model choice, not the accuracy target; used when N<=0)
        //! Explicit exponent pool -- when N>0 these override the basisPool preset (the {type,N,emin,emax} json
        //! the scaffold built by hand for the A_SG/A_SL exponent-sweep tests).
        int    N    = 0;
        double emin = 0.0;
        double emax = 0.0;
        //! Explicit exponent LIST -- the "bring your own exponents" path (used by the valence-basis generator).
        //! When non-empty it overrides N/emin/emax and the basisPool preset: every occupied angular irrep is built
        //! from these exponents, trimmed per-l by \c ltrim ((3-l)*ltrim smallest + l*ltrim largest dropped; 0 =
        //! every irrep gets the full list).  So the SAME exponents that are validated here can be emitted verbatim.
        std::vector<double> exponents;
        size_t              ltrim = 0;
        //! INDEPENDENT per-l exponent lists (l, exponents) -- highest precedence.  Each occupied angular irrep is
        //! built from its OWN list verbatim (no auto-trim), so s and p can carry different windows.  The right
        //! path for pseudopotential valence bases (the ExponentScaler N-4l trim assumes absent core exponents).
        std::vector<std::pair<int,std::vector<double>>> exponentsByL;
        //! DFT exchange-correlation override.  When set, the DFT Hamiltonian uses this functional (e.g. a
        //! libxc LDA) via the public XC selector, instead of the model's built-in Xalpha/LDA.  (model is then
        //! only used to decide HF/Dirac vs DFT -- any DFT model works.)
        std::optional<qchem::Hamiltonian::XCFunctional> xc;
        //! Pseudo-atom: replace the all-electron nuclear attraction with the GTH local + KB nonlocal PP (LDA
        //! XC), and use a valence-only PseudoAtom_EC.  The element is looked up from Z.  `valence` is the GTH
        //! Zion; 0 means "use the electron count Z-charge" (a neutral pseudo-atom).  For a pseudo-ion set both:
        //! electrons = Z - charge, Zion = valence, net charge = Zion - electrons (e.g. F⁻: valence=7, one extra e).
        bool   pseudopotential = false;
        int    valence         = 0;
        //! \brief The electron configuration, in NIST's notation ("[Ar] 3d^6 4s^1"; the "^" optional).  EMPTY = derived: the isoelectronic
        //! neutral's (all-electron) or the pseudo-ion rule (\c PseudoAtom_EC).  Stated, it must hold the run's electrons, and on a pseudo-atom its
        //! [X] core must be the pseudopotential's own (Z - q electrons; the text then names the valence, Mn2+ q7 = "[Ar] 3d^5"); a mismatch, a
        //! shape the radial atom cannot hold, or a Dirac model is REFUSED (policy AtomConfiguration).  Why it exists: a d- or f-block CATION
        //! is not its isoelectronic neutral (Fe+ is 3d^6 4s^1, not Mn's 3d^5 4s^2; V+ and Co+ move electrons from 4s INTO 3d), and no removal
        //! rule reproduces NIST's first cations (OpenWork row "All-electron ION seeds", 2026-11-15).
        std::string configuration;
        //! Accelerator escape hatch for the long tail (accelerator TUNING: GDM/Ladder/directmin + their knobs).
        //! Empty => the typed AcceleratorOptions DIIS path.  Non-empty => a full accelerator config merged over
        //! the Z-scaled DIIS defaults, with "type" selecting DIIS/GDM/Ladder/directmin (mirrors the long-serving
        //! scfrun knobs).  The common path stays typed; this is only for the accelerator-sweep driver.
        nlohmann::json accelerator;
        bool operator==(const Atom&) const = default;   //!< "is the block stated?" (another kind refuses it)
    } atom;

    //! \brief THE MOLECULE-ONLY OPTIONS: the orbital basis BY NAME and its construction.  A periodic kind REFUSES a stated one.
    struct Molecule
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
        //! Basis construction variants (threaded into BasisSet::Gaussian::Factory).  Defaults reproduce
        //! today's behaviour (in-house MnD, Cartesian).  angular==Spherical + symmetry works with the MnD engine (the SALC builder
        //! handles spherical shells, doc/OldPlans/SphericalSALCPlan.md S3a/S4); only engine LibCint + Spherical + symmetry is
        //! refused (S3b: libcint's real-harmonic layout is unwired).
        Engine      engine  = Engine::MnD;
        Angular     angular = Angular::Cartesian;
        bool operator==(const Molecule&) const = default;   //!< "is the block stated?" (a periodic kind refuses it)
    } molecule;

    //! \brief THE PERIODIC-ONLY OPTIONS (V-FACADE-PARITY P0b-2): the deck's `solid` block, one to one.  A non-periodic kind REFUSES
    //! a stated one (it is not ignored).
    struct Solid
    {
        //! \name The orbital basis FAMILY (doc/OpenWork.md §2 PW-1)
        //!@{
        //! \c GPW (default): periodic Gaussians over the ctor's molecular basis.  \c PlaneWave: plane waves at
        //! \c planeWaveEcut, the SAME Hamiltonian stack (the terms ask the basis, never its family) -- with QE's
        //! pseudopotential that makes a QE comparison a Hamiltonian-only one.  A PW run ignores the molecular
        //! basis (its +U projectors are pseudo-atomic, basis-independent); the grid knobs below are GPW's.
        //! First increment: Γ only, no imposition (the PW factory builds one Bloch block and no symmetry).
        enum class OrbitalBasis { GPW, PlaneWave };
        OrbitalBasis orbitalBasis = OrbitalBasis::GPW;
        double planeWaveEcut = 0.0;   //!< Ha, the orbital cutoff \f$\tfrac12|k+G|^2<E_{cut}\f$ (QE's ecutwfc/2); PlaneWave only
        //!@}

        //! \name Grids (the efficiency lever for ionic systems)
        //!@{
        double densityEcut  = -1.0;   //!< <0 AUTO: floor at \c cutoffFactor*alpha_max (recommended).
        double cutoffFactor = 2.0;    //!< C in that floor; 2 = the density's own product exponent.  The auto ball's START: the pilot judges it after.
        double ladderFactor = 4.0;    //!< multigrid ladder gradation (DEPTH is automatic).
        BasisSet::Gaussian::CellImages   images = BasisSet::Gaussian::CellImages::Periodic;
        rvec3_t kShift = rvec3_t(0,0,0);
        //! Where a geometry run's pilot enters the Becke ladder, and the rotation of its rules (expert deck `solid.beckeLadder`).
        qcMesh::BeckeLadderCalibration beckeLadder;
        //! The audit RSS policy's calibration: estimate = eta x byte count + kappa x raster points, and the band of measured / estimate it
        //! refuses and warns on (expert deck `solid.auditFootprint`; the default carries its measurement record, \c kAuditFootprintCalibration).
        AuditFootprintCalibration auditFootprint = kAuditFootprintCalibration;
        //! The GPW basis's numerical TOLERANCES (local-PP G-ball eps, local-PP kappa, the HartreeOnly floor fraction, an explicit
        //! multigrid cutoff list).  Tier-2 numerical policy, typed (D-ENV step 5): the defaults are today's behaviour, and the
        //! environment (GPW_VLOC_EPS, ...) is applied on top as an OVERRIDE LAYER when the run resolves its options -- every
        //! override is named on the run banner.
        qchem::BasisSet::Gaussian::GPWTolerances tolerances = {};
        //! The declared CP2K deviations (routes and accelerations CP2K does not run) -- the deck's `policy` block (D-ENV step 6b).  Unstated
        //! routes take the qchem default, or CP2K's value under \c cp2kCompat; a stated one wins over the umbrella.  The facade INSTALLS it
        //! (\c SetRunPolicy) as the first act of construction; the resolved table is printed on the run banner.  No environment variable sets it.
        qchem::RunPolicySpec policy = {};
        //! \brief DFT+U manifolds (programme step 5; doc/Pins.md pin 23): one entry per (site, l, U) -- the
        //! manifold is an INPUT, never the transition metal's d shell by assumption.  \c U in HARTREE here as
        //! everywhere inside the library; \c HubbardU(site, l, eV) below converts from the eV the literature
        //! quotes.  Empty (default) = no Hubbard term.  Increment 1: shell-averaged U, Löwdin projectors,
        //! occupations from \f$D_{out}\f$ (a declared CP2K deviation: they mix P).
        std::vector<Hamiltonian::HubbardManifold> hubbard;
        //! A SMEARED run fills each k-block to its OWN chemical potential -- the non-default (D-FERMI-DEFAULT, user 2026-10-05,
        //! landed 2026-10-21; was \c globalFermi, default false).  A k-point is a label, not a reservoir: per-k μ has no physical
        //! counterpart in a crystal ("a very obscure feature for anyone to want"), and it produced SrVO3's wrong-sign χ0(q≠0) and
        //! a ground state 15 mHa too high.  Default false = ONE μ across the k-mesh.  At kT = 0 every block keeps its prescribed
        //! integer count either way (Crystal_EC never ranks across k); per SPIN channel μ stays \c spinsShareFermi's business.
        bool perKFermiLevels = false;
        //! FREE run: build the XC mesh invariant under the space group DETECTED at the run's geometry (KAxisPlan ruling 1: a free run
        //! gets symmetric numerics, not a symmetric answer).  \c false = the plain atom-centred mesh, whatever the geometry's group --
        //! what a \c SolidSession REQUIRES of a free setup: a step re-detects the group, so the invariant mesh would SWITCH along a path
        //! (measured, NaF Γ, F along x, Becke(20,2,24): 48 ops / 16752 points on the site, 8 ops / 17536 off it; a 2.5e-5 Ha step at the
        //! site and a 9.4e-5 Ha/bohr odd term the mirror forbids).  The plain mesh moves with the atoms and is the same kind of mesh at
        //! every point; it is also ~2.2x fewer points there (7.6k).  Ignored on an imposed run (the imposed group's mesh is the run's).
        //! WHY TWO MESHES (for whoever decides whether one can serve both): the INVARIANT mesh keeps a single-geometry free run's
        //! INTEGRATOR group-equivariant, which ruling 1 consumes -- it reduces a free run by the ops whose measured
        //! \f$\|D_{gk}-MD_kM^\dagger\|\f$ passes ~1e-6, and on a non-invariant mesh no op ever passed (F 5e-4, D 2e-3; GPW_IBS
        //! itsMeshOps) -- so a session on the plain mesh gives that reduction up.  The PLAIN mesh keeps E continuous along a path that
        //! BREAKS the group, which the invariant mesh cannot: its group, and so its points, change with the geometry.  The open question
        //! (OpenWork, X): a session invariant under the group the WHOLE session keeps (a symmetry-preserving relaxation keeps its space
        //! group), so one mesh serves both.
        bool xcMeshSymmetric = true;
        //! \brief HOW THE δ FIT HOLDS ITS ORBITAL TABLES ON THE XC MESH -- a MEMORY-vs-SPEED TRADE-OFF the user may decide (lane S 2026-11-19;
        //! the user: "For some of these the user should be allowed to decide").  Both give the same numbers (to roundoff):
        //!   \c dense   one table \f$\Phi^k_{gi}\f$ per k-block: nk x points x functions x 16 B.  The fastest contractions, the largest RAM.
        //!   \c images  one k-independent record of the lattice-image values, shared by every k-block and phased when used: ~25-30 B per
        //!              (point, function) whatever nk -- smaller from about two complex k-blocks on, a few percent slower.
        //!   \c auto    (default) the realization with the lower PEAK memory, unless the other is faster by more than
        //!              \c policy::kPhiSpeedGainThreshold and still fits (qchem.PhiTablesPolicy).
        //! A stated value wins; one the policy finds peculiar (dense near the memory budget, images where they are larger) runs with a
        //! warning naming what \c auto would pick.  The [policy] line prints both realizations' peak, average and wall either way.
        policy::XCOrbitalTables xcOrbitalTables = policy::XCOrbitalTables::Auto;
        //! THE CELL'S STRAIN ALLOWANCE s (linear, a fraction; the cell half's frozen rule, doc/Records/CellHalf_Design.md §3, ruled
        //! 2026-11-10).  0 (default): the cell is fixed -- a geometry session refuses a moved cell (\c PESFailure::CellChange) and every size
        //! is decided on the cell itself, as always.  > 0: every discrete size (G sets, raster and FFT counts, uniform-mesh divisions) is
        //! decided ONCE on the SIZING BOX (1+s) A_ref and held as the cell strains (CP2K's CELL_REF), so E(cell) is smooth; a session then
        //! follows any cell whose largest principal stretch from the reference is <= 1+s and refuses one beyond.  The price: the
        //! reference point itself runs on sets sized for a cell s larger (a POLICY, reported).
        double       cellStrainAllowance = 0.0;
        //! Bands (eigenpairs) kept per Bloch block, QE's nbnd: 0 (default) = all of them.  > 0 = a PARTIAL eigensolve of the
        //! lowest \c numBands -- a plane-wave block is thousands wide and the SCF occupies a few dozen.  Orthonormal basis only
        //! (the plane-wave family); must exceed the larger spin channel's electron count; refused beside a direct minimiser
        //! (GDM / the Ladder's GDM rung), which rotates into the WHOLE virtual space.
        size_t       numBands = 0;
        //! ANSATZ POLICY (doc/RealComplexPlan.md §6): force every Bloch block COMPLEX even where the run
        //! admits real ones (a TRIM k-point under a realness-preserving Hamiltonian -- with the default
        //! LDA stack every Γ / zone-boundary block runs REAL: real S/T/V/KB, real eigensolve, real
        //! quadrature GEMMs).  The physics is identical either way (the 3c-3 acceptance gate pins ON==OFF
        //! to machine precision); complex is the DOWNGRADE direction (§1), kept as the escape hatch for
        //! complex-instability experiments and as the gate's A/B door.
        //!
        //! NOT needed for MOM any more (R2.21, 2026-08-17).  It used to be: a real block's MOM reference had
        //! nowhere to live in a run-typed reference home, so switching \c SCFParams::UseMOM on threw mid-SCF
        //! and the documented workaround was to come here and pay for an all-complex run.  The occupation
        //! state now keys each block's reference by the BLOCK's scalar, so MOM runs REAL on a mixed mesh --
        //! measured identical to the forced-complex twin (Si (3,1,1), ΔE = 1.8e-15) at the real run's cost.
        //! No knob, and no surprise: turning MOM on no longer changes what the run is made of.
        bool forceComplex = false;
        //! \brief THE LINEAR RESPONSES THIS RUN WILL RUN, one q-mesh each (empty = none).  The ansatz is decided ONCE per run, at setup,
        //! from k AND q (user ruling 2026-11-18, OpenWork "LR on REAL TRIM blocks"): a response's dD at q != 0 couples block k to block
        //! k+q through a complex phase, so a run that will run one builds EVERY block complex from the start -- never a real block
        //! re-typed, rebuilt or released mid-run.  q = 0 alone counts too: the response face is complex-only (route B, typing each object
        //! from its own quantum numbers, is the hub change that would lift it).  \c SolidSetup::realTRIMBlocks is the decision.  A deck
        //! fills it per run (its two-run split: none for the ground state, the planned responses for the response run); code that calls
        //! \c SolidCalculation::HubbardLinearResponse states it (or \c forceComplex).
        std::vector<ivec3_t> responseQMeshes = {};

        //! \name Knobs the GPW test driver carried that this struct did not (2026-08-25)
        //! Added so a periodic recipe can be stated in ONE place instead of being assembled across a test's
        //! options struct, its run driver and its accelerator helper (user: "some parameters are set in TESTF,
        //! some in RunGpw, some in MakeGpwAccelerator, all spread out by 100's of lines").
        //!@{
        //! One Fermi level shared by BOTH spin channels (the moment becomes an OUTPUT) rather than a fixed
        //! per-channel electron count.
        bool spinsShareFermi = false;
        //! Impose the DETECTED grey group with the magnetic decoration suppressed -- the erasure control.
        bool greyImposition  = false;
        //! Pin the MOM reference to the SEED's own occupied subspace before iteration 1, so a hot stage's
        //! CHARACTER continues, not just its density.
        bool momFromSeed     = false;
        //! \note MANDATORY for a polarized run that carries a Hubbard manifold (OpenWork §3 "The magnetic decoration is
        //! DERIVED FROM THE SEED", ruled 2026-10-23: the ordering is chemistry, stated like the Hubbard subspace, never
        //! inferred -- a seed change must never change the imposed symmetry or the manifolds' site groups).  A
        //! non-magnetic polarized run states all zeros.  See \c RequireStatedMagneticOrdering.
        //! \brief THE CHARGE ORDER of an \c IonicSAD seed: one formal charge per cell site (Li +1, Mn +3 / +4, O -2 ...), summing to
        //! zero.  EMPTY (the default) = DERIVED by the electronegativity heuristic (\c ChargeDensity::IonicFormalCharges), which is
        //! right where every site of an element takes the same charge, but on a MIXED-VALENCE cell puts the charge states of
        //! symmetry-equivalent sites by ATOM INDEX (LiMn₂O₄'s Mn³⁺/Mn⁴⁺ on the 16d sites): a charge order nobody chose.  State it, like
        //! \c siteSpins, when the order matters.  Only the IonicSAD seed reads it (a stated order on another seed THROWS: nothing would
        //! consume it).  See \c RequireConsistentSiteCharges.
        std::vector<int> siteCharges = {};
        //! The grid accuracy CALIBRATIONS (doc/GridEpsilonPlan.md §1; deck `solid.accuracyCalibration`): the eta of each raster, the
        //! V_loc routing, the ladder's rung gate, the sampler.  Its eps is the run's \c accuracy (set by \c SolidSetup::Resolve; a stated
        //! one here is refused).
        qchem::BasisSet::Gaussian::GridAccuracy grid = {};
        bool operator==(const Solid&) const = default;   //!< "is the block stated?" (a non-periodic kind refuses it)
    } solid;
};

} // namespace qchem
