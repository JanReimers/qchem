// File: Calculation/SolidSetup.C
//
// qchem::SolidSetup -- THE RUN SETUP of a periodic calculation (doc/CalcNetworkReview.md §2 and §7, step 6): the options as
// asked, and every decision RESOLVED from them ONCE -- the pseudopotential set, the imposition (the run policy's veto applied),
// the magnetic decoration, the XC quadrature (Auto resolved), the spin counts.  An immutable value, shared: a geometry point
// (a SolidCalculation) is built FROM it and changes none of it, which is what lets a geometry optimizer move the atoms under one
// fixed setup (step 8).  SolidCalcOptions lives here because the setup is what reads it; qchem.SolidCalculation re-exports it.
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
import qchem.SCFAccelerator.Factory;          // Type, SolidAcceleratorOptions
import qchem.ChargeDensity.Seed;              // SeedStrategy
import qchem.ChargeDensity.AtomicDensity;     // SeedKey
import qchem.Mesh;                            // qcMesh::MeshParams / UnitCellKind
import qchem.Mesh.XCPolicy;                   // qcMesh::XCMeshDefault (the BeckeXCParams recipe under Auto)
import qchem.RunPolicy;                       // RunPolicySpec
import qchem.LASolver;                        // qchem::Ortho
export import qchem.AuditFootprint;           // AuditFootprintCalibration (the audit RSS policy's expert knob)
export import qchem.Pseudopotential.Set;      // PseudopotentialSet -- the setup's resolved choice

export namespace qchem
{

//! \brief How to set up a periodic calculation.  Designated-initializer friendly:
//!     SolidCalculation calc(SolidSetup::Resolve(lattice, basis, {.Nelec=8, .species={{"Si",4}}}));
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

struct SolidCalcOptions
{
    //! \name The system
    //!@{
    //! VALENCE electrons per cell (Zion, not Z).  0 = DERIVED: \f$\sum_{\rm atoms} q\f$ of the run's species minus \c charge.  A stated value
    //! must equal that, or \c SolidSetup::Resolve REFUSES the run naming both numbers (2026-10-31: a stated Na q9 once ran NaF with the
    //! material's 8 electrons instead of 16 -- a cell charged +8, silently).
    int Nelec = 0;
    //! The cell's stated NET CHARGE (electrons removed; -1 = one extra electron).  A charged cell is legal only when stated here.
    int charge = 0;
    //! Spin multiplicity 2S+1.  0 (default) = minimal spin: the ζ=0 collapse for even \c Nelec, a doublet
    //! for odd.  >=1 promotes to the SPIN-NATIVE two-channel run (nUp-nDown = 2S).  NB unlike the molecular
    //! facade, 1 here means the EXPLICIT two-channel singlet (nUp=nDn) -- the ζ=0 cross-check of the
    //! polarized machinery against the unpolarized anchors, not a shortcut to it.
    int multiplicity = 0;
    //! The pseudopotential SET the run uses (src/Pseudopotential/Data/pseudopotentials.json; doc/CalcNetworkReview.md §4): per element
    //! the GTH valence variant, and the XC functional the entries were generated for.  Chosen ONCE here; the Hamiltonian's PP terms,
    //! the XC-grid sharpness gather, the +U projector and the saved-state fingerprint are all handed the same resolved value.
    std::string pp = "gth-lda-valence";
    //! Per-element valence overrides on top of \c pp's set, as `(element, valence)`: replaces the set's q for that element (a
    //! semicore Mn q15 run) or adds it.  Empty (default) = the set as is.  Multi-species is the norm for a compound; a per-Z
    //! router PP is built over the elements of the structure so each atom gets its own GTH entry.
    std::vector<std::pair<std::string,int>> species;
    //!@}

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
    //! \brief The XC real-space quadrature.  DEFAULT \c Auto = "you choose": the class costs a uniform
    //! mesh sized to this run's sharpness against the atom-centred Becke mesh and takes the cheaper, with
    //! a margin (V1.26/V2.4).  Pin \c cellKind to \c Uniform or \c Becke to overrule it and be obeyed --
    //! the grid is a user decision, and the library's job is only to warn about a dominated one.  The recipe fields default to
    //! \c BeckeXCParams()'s (\c qcMesh::XCMeshDefault), so stating ONLY \c cellKind Becke runs the calibrated mesh.
    qcMesh::MeshParams xcMesh=qcMesh::XCMeshDefault();
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
    //! \brief The ONE grid accuracy target eps (deck `solid.accuracy`; Ha per atom on an energy run, per matrix element on a geometry
    //! run -- \c SolidCalculation::Consumes) and its calibration eta (expert, `solid.accuracyCalibration`), doc/GridEpsilonPlan.md §1.
    //! C2 step 3: the GPW raster demands (\c tolerances.rasterX ...) are DERIVED from it -- from the stated calibration at Resolve, or
    //! by the run's own pilot audit after its SCF (the default).  eps = 0 is NO target: the rasters are \c tolerances as stated (0 =
    //! the ball's own raster, CP2K's arrangement), and stating a raster demand with eps > 0 is refused (the target derives them).
    qchem::BasisSet::Gaussian::GridAccuracy accuracy = {};
    //! The declared CP2K deviations (routes and accelerations CP2K does not run) -- the deck's `policy` block (D-ENV step 6b).  Unstated
    //! routes take the qchem default, or CP2K's value under \c cp2kCompat; a stated one wins over the umbrella.  The facade INSTALLS it
    //! (\c SetRunPolicy) as the first act of construction; the resolved table is printed on the run banner.  No environment variable sets it.
    qchem::RunPolicySpec policy = {};
    //! WHICH fit basis represents v_xc -- ORTHOGONAL to the grid above.  \c Auto picks Delta whenever the
    //! plane-wave fit cannot serve (Becke grid, or a polarized run).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
    //! WHICH density XC sees (the deck's `xc.rho`, user ruling 2026-11-09): \c Exact = rho[D], \c Fit = the Vxc fit's projection
    //! P rho.  The two coincide on every cell a solid realizes Fit on: the delta fit, and the plane-wave ball over plane-wave orbitals
    //! (>= 4 Ecut holds every pair); a stated \c Fit is reported there.  Over GPW orbitals no ball holds rho[D], so Fit on the plane-wave
    //! fit is REFUSED by name, never run as Exact (\c qchem.XCRhoPolicy, \c ReportSolidXCRho).
    Hamiltonian::XCρ xcRho = Hamiltonian::XCρ::Exact;
    //! \brief DFT+U manifolds (programme step 5; doc/Pins.md pin 23): one entry per (site, l, U) -- the
    //! manifold is an INPUT, never the transition metal's d shell by assumption.  \c U in HARTREE here as
    //! everywhere inside the library; \c HubbardU(site, l, eV) below converts from the eV the literature
    //! quotes.  Empty (default) = no Hubbard term.  Increment 1: shell-averaged U, Löwdin projectors,
    //! occupations from \f$D_{out}\f$ (a declared CP2K deviation: they mix P).
    std::vector<Hamiltonian::HubbardManifold> hubbard;
    //!@}

    //! \name Convergence machinery
    //!@{
    //! The Fock-side accelerator.  \c Auto (default) is the facade's policy, resolved per stage, PAIRED with an \c Auto mixer
    //! (OpenWork §3 "Kerker/Pulay as the DEFAULT mixer", measured 2026-10-22: Pulay on ρ beside DIIS on F fights -- free Si
    //! k222 stalls, Si Γ diverges -- while Kerker/Pulay with NO Fock extrapolator beat linear + DIIS on every probe; the
    //! DEFAULT since 2026-10-26): a seeded or carried stage -> \c Null + Pulay on D (no Kerker, no fit: pin 28); a RESTORED
    //! stage (Restart, same Hamiltonian) at kT = 0 -> the Ladder on its GDM rung, a restored smeared one -> \c DIIS, either
    //! beside a linear step capped at 0.3.
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::Auto;
    //! A SMEARED run fills each k-block to its OWN chemical potential -- the non-default (D-FERMI-DEFAULT, user 2026-10-05,
    //! landed 2026-10-21; was \c globalFermi, default false).  A k-point is a label, not a reservoir: per-k μ has no physical
    //! counterpart in a crystal ("a very obscure feature for anyone to want"), and it produced SrVO3's wrong-sign χ0(q≠0) and
    //! a ground state 15 mHa too high.  Default false = ONE μ across the k-mesh.  At kT = 0 every block keeps its prescribed
    //! integer count either way (Crystal_EC never ranks across k); per SPIN channel μ stays \c spinsShareFermi's business.
    bool perKFermiLevels = false;
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
    qchem::Ortho ortho    = qchem::Cholesky;
    //! THE NEAR-DEPENDENCE FLOOR on S's eigenvalues: what the vet-stage trim trims at (basis autoTrim, and a geometry step's re-vet), and
    //! the floor a TRUNCATING ortho solver drops directions under -- one quantity (the vet's metric is the truncating solver's own drop
    //! count at it), two readers.  The solver reads it only when it truncates (\c SolverTruncation); the solid default, Cholesky, does
    //! not, so there the floor is the trim's alone (OpenWork §3 "solid.orthoTol has two meanings", user feedback 36, 2026-11-12).
    double       orthoTol = 0.0;
    //! \brief The tolerance handed to the ortho SOLVER: \c orthoTol when \c ortho truncates (Eigen, SVD, CholeskyPivoted), else 0 --
    //! Cholesky, Auto (it finds its own gap) and Orthonormal take none.
    double SolverTruncation() const
    { return (ortho==qchem::Eigen || ortho==qchem::SVD || ortho==qchem::CholeskyPivoted) ? orthoTol : 0.0; }
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
    std::vector<int> siteSpins = {};
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
    // NB: the test driver's `realTRIMBlocks` is NOT added here -- it is \c forceComplex with the opposite
    // polarity (both feed GPWParams::hamPreservesReal).  One decision, one field.
    //! Names this run in its own console output -- the fingerprint, the stage banners, the order trace.
    std::string label = "gpw";
    //! \brief SAVE THE STATE HERE after every \c Converge -- each anneal stage, each ACBN0 outer step -- converged
    //! or not (the file records which).  EMPTY (default) = never.  The unattended-run safety net (CK-1): a crash
    //! costs at most the stage in flight, and \c SolidCalculation::Restart picks up from the last one.  A write
    //! failure here WARNS and the run continues -- a checkpoint must never cost the computation it protects
    //! (\c SaveState called directly throws as usual).
    std::string saveStateTo;
    // (No onIteration any more -- 2026-10-30.  The per-iteration observer is a SolidCalculation / Restart / SolidSession
    //  CONSTRUCTOR argument, as on the molecular MoleculeCalculation: an observer is not a resolved decision.)
    // (No orderName/orderProbe any more -- 2026-09-14.  The order parameter is the INTEGRATED site moment
    //  every polarized run carries on its EnergyBreakdown (the m_site cell of the iteration row); the
    //  caller-supplied POINT probe it replaced was a spin density, not a moment.)
    //!@}
    //!@}
};

//! \brief THROWS unless \a o states its magnetic ordering where the run needs one: a polarized run carrying a
//! Hubbard manifold must give \c siteSpins (\c greyImposition, the erasure control, is exempt), and a stated
//! \c siteSpins must have one entry per site of \a st, each -1, 0 or +1.  A broken invariant, not an Outcome: the
//! facade calls it first thing, and \c deck::Vet calls it so \c scripts/vet-decks catches a deck before it runs.
void RequireStatedMagneticOrdering(const SolidCalcOptions& o, const Structure& st);
//! \brief THROWS unless a stated \c o.siteCharges has one entry per site of \a st, sums to zero (formal charges conserve the
//! cell's charge, as the derived ones do), and seeds an \c IonicSAD run.  Empty passes.  Called by the facade and by \c deck::Vet,
//! as \c RequireStatedMagneticOrdering is.
void RequireConsistentSiteCharges(const SolidCalcOptions& o, const Structure& st);
//! \brief The electron count per spin channel, (nUp, nDown), for \a o.multiplicity and the RESOLVED electron count \a nel (decision 2).
//! THROWS on an impossible request: a multiplicity whose parity disagrees with \a nel or that exceeds it -- in particular multiplicity 1
//! (the EXPLICIT two-channel singlet, nUp = nDown) on an odd \a nel, unless \c o.spinsShareFermi lets one Fermi level set the channels
//! (a metal: the counts are only the split the occupations start from).  Called by the resolve phase and by \c deck::Resolve.
std::pair<int,int> SpinChannelCounts(const SolidCalcOptions& o, int nel);

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
    //! from \c SolidCalcOptions::tolerances (the one source -- a caller's value here is replaced), frozen with the trim so a geometry
    //! step's basis computes at the same setting.  The GPW evaluator checks them; it no longer sets them.
    BasisSet::Gaussian::GPWTolerances tolerances;

    //! The basis over \a cell's atoms, its pair-loop \c tolerances bound before anything can compute from it.  CHECKS the file's
    //! declared pseudopotential provenance against \a pp first (D-BASIS-PP): THROWS naming the file, the element and both q on a
    //! mismatch.  \a cell must outlive the basis.
    std::shared_ptr<const BasisSet::Real_BS> Make(const Structure& cell, const Pseudopotential::PseudopotentialSet& pp) const;
};

//! \brief THE RUN SETUP of a periodic calculation: \c SolidCalcOptions as asked, plus every decision resolved from them ONCE.
//!
//! WHAT IS HERE, AND WHAT IS NOT (doc/CalcNetworkReview.md §2).  A run splits into three lifetimes: the SETUP (the choices and
//! what is resolved from them), the GEOMETRY POINT (the Bloch basis, S, the Hamiltonian's geometry-keyed terms, the electron
//! configuration over the basis's k-blocks) and the SCF STATE.  Only the first lives here.  Every member is a decision that does
//! not change when the atoms move -- the set, the imposition, the decoration, the XC quadrature (sized from exponents and the cell,
//! never from positions), the spin counts -- so a geometry optimizer can rebuild the point at new positions under the SAME setup.
//! A field that would depend on positions does not belong here.
//!
//! Built only by \c Resolve, held as \c shared_ptr<const SolidSetup>: immutable once made.
struct SolidSetup
{
    //! \brief THE ONE RESOLVE PHASE (doc/Pins.md pin 11).  Installs \a opts.policy as the process run policy FIRST (every factory
    //! after it consults the policy), then decides.  THROWS on options that contradict each other or the structure -- an unstated
    //! magnetic ordering (\c RequireStatedMagneticOrdering), a multiplicity whose parity disagrees with \c Nelec, \c numBands too
    //! small, a plane-wave run asked for a k-mesh, a shift or an imposition: broken invariants of the request, not outcomes.
    //! \a mol: the molecular orbital basis at \a lat's geometry -- the XC-mesh sharpness reads its α_max, and a GPW point is built
    //! over it; a plane-wave run ignores it.
    static std::shared_ptr<const SolidSetup> Resolve(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                                     const SolidCalcOptions& opts,
                                                     const SCFAccelerators::SolidAcceleratorOptions& acc = {});
    //! The same from a basis RECIPE: the molecular basis is made at \a lat's geometry by \a recipe (its provenance checked against
    //! the run's set), after the vet-stage loop -- if the recipe asks for it -- has decided the trim, which is then frozen in
    //! \c recipe.  A setup with a recipe can rebuild its basis at moved atoms (a geometry step); one from a pre-built \a mol cannot.
    static std::shared_ptr<const SolidSetup> Resolve(const Lattice_3D& lat, const GaussianBasisRecipe& recipe,
                                                     const SolidCalcOptions& opts,
                                                     const SCFAccelerators::SolidAcceleratorOptions& acc = {});

    SolidCalcOptions                         opts;   //!< as asked, except \c ortho (a plane-wave run's default Cholesky becomes Orthonormal) and a \c Default seed (resolved)
    SCFAccelerators::SolidAcceleratorOptions acc;    //!< as asked; a stage derives its own from it (a restored stage's Ladder rung)
    Pseudopotential::PseudopotentialSet      pp;     //!< \c opts.pp with the \c opts.species overrides applied
    //! The REFERENCE geometry: the lattice the setup was resolved at (the trim vetted, the imposed group detected), and its
    //! structure -- the one the reference basis is built over, held so it outlives that basis.
    Lattice_3D                               reference;
    std::shared_ptr<const Structure>         cell;
    //! The seed library's key: \c pp's functional and, per element of the structure, its q (the seed is the SAME pseudo-atom as the terms).
    ChargeDensity::SeedKey                   seedKey;
    std::shared_ptr<const BasisSet::Real_BS> mol;    //!< the molecular orbital basis at the REFERENCE geometry
    //! How \c mol was made, the trim FROZEN; empty when the caller handed a pre-built basis (that setup cannot follow moved atoms).
    std::optional<GaussianBasisRecipe>       recipe;
    bool planeWave = false;   //!< \c opts.orbitalBasis is PlaneWave
    bool polarized = false;   //!< the spin-native two-channel run (\c opts.multiplicity >= 1)
    //! The caller's \c imposeSymmetry AND the run policy's permission (CP2K parity vetoes the capability).  Nothing reads
    //! \c opts.imposeSymmetry after Resolve, or the two would drift.
    bool imposed   = false;
    //! The imposition's MAGNETIC DECORATION (one +1/-1/0 per site; what makes the group Shubnikov): stated, or -- for an imposed
    //! polarized run without +U only -- derived from the seed.  EMPTY when the run imposes nothing or imposes the grey group.
    std::vector<int> siteSpins;
    std::vector<int> decoration;   //!< \c siteSpins, or all zeros: the per-site spins the run was built with (the +U site groups read it)
    //! \brief THE GROUP A GEOMETRY STEP MUST KEEP, read at \a at: the space group detected there (\c GetSpaceGroup's \f$\{W|\tau\}\f$, or
    //! the Shubnikov group under \c siteSpins) when the setup is \c imposed, P1 (the identity alone) when it is free.  ONE rule for its two
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
    qcMesh::MeshParams xcMesh;     //!< the XC quadrature, \c Auto RESOLVED (the facade's decision 1)
    int nUp = 0, nDown = 0;        //!< the electron counts per spin channel (decision 2: multiplicity -> (nUp, nDown))
    //! \brief The GPW raster demands in \c opts.tolerances are FINAL (C2 step 3): stated (accuracy 0), derived from a stated calibration,
    //! or settled by a pilot audit (\c SolidCalculation::SizedSetup).  false = this setup's first SCF is the PILOT: it runs on the ball's
    //! own rasters and \c SolidCalculation sizes them from what it delivered.  A plane-wave run has no raster: true.
    bool rastersResolved = true;
    //! \brief WHO SET the raster demands this setup runs (the report's \c grids.rasterSource): the user (\c Stated, accuracy 0), a stated
    //! calibration at Resolve (\c Calibration), the pilot audit (\c Pilot: on the ball's own rasters until it raises them), a saved state
    //! sized for this eps (\c SavedState), or the start raised for an ill-conditioned overlap (\c Conditioning).  Every place that writes
    //! \c opts.tolerances' raster demands writes this beside it; a reader of the grids can then tell a deck value from a measured one.
    enum class RasterSource { Stated, Calibration, Pilot, SavedState, Conditioning };
    RasterSource rasterSource = RasterSource::Stated;

private:
    SolidSetup(Pseudopotential::PseudopotentialSet p, const Lattice_3D& lat)
        : pp(std::move(p)), reference(lat), cell(lat.GetStructure()) {}
    //! The decisions, shared by both \c Resolve overloads (the policy is already installed).
    static std::shared_ptr<SolidSetup> Decide(const Lattice_3D&, std::shared_ptr<const BasisSet::Real_BS> mol,
                                              const SolidCalcOptions&, const SCFAccelerators::SolidAcceleratorOptions&);
};

} //export namespace qchem
