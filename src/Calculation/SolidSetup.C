// File: Calculation/SolidSetup.C
//
// qchem::SolidSetup -- THE RUN SETUP of a periodic calculation (doc/CalcNetworkReview.md §2 and §7, step 6): the options as
// asked, and every decision RESOLVED from them ONCE -- the pseudopotential set, the imposition (the run policy's veto applied),
// the magnetic decoration, the XC quadrature (Auto resolved), the spin counts.  An immutable value, shared: a geometry point
// (a SolidCalculation) is built FROM it and changes none of it, which is what lets a geometry optimizer move the atoms under one
// fixed setup (step 8).  SolidCalcOptions lives here because the setup is what reads it; qchem.SolidCalculation re-exports it.
module;
#include <memory>
#include <string>
#include <utility>
#include <vector>
export module qchem.SolidSetup;

import qchem.Lattice_3D;                      // Lattice_3D (crystal + BZ mesh)
import qchem.Structure;                       // Structure
import qchem.BasisSet;                        // Real_BS (the molecular orbital basis)
import qchem.BasisSet.Lattice.BasisSet;       // RasterPolicy, CellImages, GPWTolerances
import qchem.Hamiltonian.Factory;             // HubbardManifold, VxcFit
import qchem.SCFAccelerator.Factory;          // Type, SolidAcceleratorOptions
import qchem.SCFIterator;                     // SolidSCFIterator::Observer
import qchem.ChargeDensity.Seed;              // SeedStrategy
import qchem.ChargeDensity.AtomicDensity;     // SeedKey
import qchem.Mesh;                            // qcMesh::MeshParams / UnitCellKind
import qchem.RunPolicy;                       // RunPolicySpec
import qchem.LASolver;                        // qchem::Ortho
export import qchem.Pseudopotential.Set;      // PseudopotentialSet -- the setup's resolved choice

export namespace qchem
{

//! \brief How to set up a periodic calculation.  Designated-initializer friendly:
//!     SolidCalculation calc(lattice, basis, {.Nelec=8, .species={{"Si",4}}});
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
    int Nelec = 0;                                        //!< VALENCE electrons per cell (Zion, not Z).
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
    double cutoffFactor = 2.0;    //!< C in that floor; 2 = the density's own product exponent.
    double ladderFactor = 4.0;    //!< multigrid ladder gradation (DEPTH is automatic).
    BasisSet::PlaneWave::RasterPolicy raster = BasisSet::PlaneWave::RasterPolicy::BallOnly;
    BasisSet::Gaussian::CellImages   images = BasisSet::Gaussian::CellImages::Periodic;
    rvec3_t kShift = rvec3_t(0,0,0);
    //! \brief The XC real-space quadrature.  DEFAULT \c Auto = "you choose": the class costs a uniform
    //! mesh sized to this run's sharpness against the atom-centred Becke mesh and takes the cheaper, with
    //! a margin (V1.26/V2.4).  Pin \c cellKind to \c Uniform or \c Becke to overrule it and be obeyed --
    //! the grid is a user decision, and the library's job is only to warn about a dominated one.
    qcMesh::MeshParams xcMesh{.cellKind=qcMesh::UnitCellKind::Auto};
    //! The GPW basis's numerical TOLERANCES (local-PP G-ball eps, local-PP kappa, the HartreeOnly floor fraction, an explicit
    //! multigrid cutoff list).  Tier-2 numerical policy, typed (D-ENV step 5): the defaults are today's behaviour, and the
    //! environment (GPW_VLOC_EPS, ...) is applied on top as an OVERRIDE LAYER when the run resolves its options -- every
    //! override is named on the run banner.
    qchem::BasisSet::Gaussian::GPWTolerances tolerances = {};
    //! The declared CP2K deviations (routes and accelerations CP2K does not run) -- the deck's `policy` block (D-ENV step 6b).  Unstated
    //! routes take the qchem default, or CP2K's value under \c cp2kCompat; a stated one wins over the umbrella.  The facade INSTALLS it
    //! (\c SetRunPolicy) as the first act of construction; the resolved table is printed on the run banner.  No environment variable sets it.
    qchem::RunPolicySpec policy = {};
    //! WHICH fit basis represents v_xc -- ORTHOGONAL to the grid above.  \c Auto picks Delta whenever the
    //! plane-wave fit cannot serve (Becke grid, or a polarized run).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
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
    //! The SCF seed.  Default \c IonicSAD (V2.2, 2026-09-14): a superposition of the library's per-species
    //! valence densities at their formal charges (Na+ F-, Mn2+ O2-; an elemental solid is the neutral SAD).
    //! \c Uniform (\f$\rho=N/V\f$, the old default) is an explicit OPT-IN, because it has a STABLE WRONG
    //! BASIN for electron-sparse systems: the Na-doublet campaign converged a lone-electron doublet 72 mHa
    //! high from it with every health metric green (doc/CleanupCandidates.md V2.2).  The molecular facade
    //! already defaults DFT to SAD for the same reason.
    //! \note IonicSAD needs every species in \c atomic_valence_densities.json (today: O, F, Na q1, Si, Mn);
    //! a missing species THROWS at seed time -- state \c Uniform, or generate the entry (CLIapps/valgen).
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::IonicSAD;
    qchem::Ortho ortho    = qchem::Cholesky;
    double       orthoTol = 0.0;
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
    //! is the whole point of a magnetic campaign.  Ignored unless \c imposeSymmetry.
    //! \note \c greyImposition forces the empty decoration regardless -- that arm is the erasure control.
    std::vector<int> siteSpins = {};
    //! \note MANDATORY for a polarized run that carries a Hubbard manifold (OpenWork §3 "The magnetic decoration is
    //! DERIVED FROM THE SEED", ruled 2026-10-23: the ordering is chemistry, stated like the Hubbard subspace, never
    //! inferred -- a seed change must never change the imposed symmetry or the manifolds' site groups).  A
    //! non-magnetic polarized run states all zeros.  See \c RequireStatedMagneticOrdering.
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
    //! \brief Per-iteration telemetry, live FROM CONSTRUCTION.
    //!
    //! It belongs in the options block rather than only on \c OnIteration because the constructor
    //! CONVERGES: an observer attached afterwards has already missed stage 0 -- which, on an annealed
    //! recipe, is the stage whose trajectory the fingerprint most wants.  Put it here and the run is
    //! observed from its first iteration, with the rest of the recipe, in one place.
    qchem::SCFIterator::SolidSCFIterator::Observer onIteration = nullptr;
    // (No orderName/orderProbe any more -- 2026-09-14.  The order parameter is the INTEGRATED site moment
    //  every polarized run carries on its EnergyBreakdown (SCFProgress::order, the m_site column); the
    //  caller-supplied POINT probe it replaced was a spin density, not a moment.)
    //!@}
    //!@}
};

//! \brief THROWS unless \a o states its magnetic ordering where the run needs one: a polarized run carrying a
//! Hubbard manifold must give \c siteSpins (\c greyImposition, the erasure control, is exempt), and a stated
//! \c siteSpins must have one entry per site of \a st, each -1, 0 or +1.  A broken invariant, not an Outcome: the
//! facade calls it first thing, and \c deck::Vet calls it so \c scripts/vet-decks catches a deck before it runs.
void RequireStatedMagneticOrdering(const SolidCalcOptions& o, const Structure& st);

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

    SolidCalcOptions                         opts;   //!< as asked, except \c ortho: a plane-wave run's default Cholesky becomes Orthonormal
    SCFAccelerators::SolidAcceleratorOptions acc;    //!< as asked; a stage derives its own from it (a restored stage's Ladder rung)
    Pseudopotential::PseudopotentialSet      pp;     //!< \c opts.pp with the \c opts.species overrides applied
    //! The seed library's key: \c pp's functional and, per element of the structure, its q (the seed is the SAME pseudo-atom as the terms).
    ChargeDensity::SeedKey                   seedKey;
    std::shared_ptr<const BasisSet::Real_BS> mol;    //!< the molecular orbital basis at the REFERENCE geometry
    bool planeWave = false;   //!< \c opts.orbitalBasis is PlaneWave
    bool polarized = false;   //!< the spin-native two-channel run (\c opts.multiplicity >= 1)
    //! The caller's \c imposeSymmetry AND the run policy's permission (CP2K parity vetoes the capability).  Nothing reads
    //! \c opts.imposeSymmetry after Resolve, or the two would drift.
    bool imposed   = false;
    //! The imposition's MAGNETIC DECORATION (one +1/-1/0 per site; what makes the group Shubnikov): stated, or -- for an imposed
    //! polarized run without +U only -- derived from the seed.  EMPTY when the run imposes nothing or imposes the grey group.
    std::vector<int> siteSpins;
    std::vector<int> decoration;   //!< \c siteSpins, or all zeros: the per-site spins the run was built with (the +U site groups read it)
    qcMesh::MeshParams xcMesh;     //!< the XC quadrature, \c Auto RESOLVED (the facade's decision 1)
    int nUp = 0, nDown = 0;        //!< the electron counts per spin channel (decision 2: multiplicity -> (nUp, nDown))

private:
    explicit SolidSetup(Pseudopotential::PseudopotentialSet p) : pp(std::move(p)) {}
};

} //export namespace qchem
