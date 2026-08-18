// File: Calculation/SolidCalculation.C
//
// qchem::SolidCalculation -- the production front door for a PERIODIC (GPW / plane-wave) SCF.  The
// solid peer of qchem::Calculation (molecules) and qchem::AtomCalculation (atoms): it owns the whole
// object graph -- basis, electron configuration, Hamiltonian, accelerator, iterator -- runs the
// canonical assemble->converge recipe, and exposes high-level questions rather than the graph.
//
// WHY IT EXISTS, beyond symmetry with the other two (doc/CleanupCandidates.md Step 4).  The user's
// reason is the operative one: there needs to be something SPECIFIC AND NAMED to point at when a
// decision belongs ABOVE the SCF iterator.  Those decisions were previously spread through a driver
// function in an integration test, which meant that in practice they had no home at all -- D6 moved
// the Becke recipe out of that file, V1.26 moved the grid selector, and each time the question "where
// does this belong?" had no better answer than "not there".  This class is that answer.
//
// DECISIONS THIS CLASS OWNS (the point of it):
//   - WHICH XC QUADRATURE the run uses.  It gathers the run's sharpness (basis alpha_max + local-PP
//     alpha_pp, both off ABSTRACT capability faces) and hands it to qcMesh::ResolveXCMesh, which costs
//     the two grids and picks (V1.26/V2.4).  A caller may still pin cellKind and be obeyed.
//   - WHICH ACCELERATOR, by policy name rather than by constructing internals.
//   - The SPIN BOOKKEEPING: multiplicity -> (nUp,nDown) -> Crystal_EC, with the parity check that
//     catches "singlet with an odd electron count" before it silently empties a channel.
//
// WHY IT COULD NOT BE WRITTEN UNTIL NOW.  Both the solid Hamiltonian and the complex accelerator were
// reachable ONLY through .Internal. modules, which per CLAUDE.md may be imported across a library
// boundary by unit tests and nothing else.  That -- not carelessness -- is why the driver lived in a
// test file: it was the one place the cheat was legal.  Step 4 opened the two public doors first
// (qchem.Hamiltonian.Factory's cHamiltonian overload; qchem.SCFAccelerator.Factory's typed-options
// overload), so this file imports ZERO internals, exactly as the molecular facades already manage.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
export module qchem.SolidCalculation;

import qchem.Lattice_3D;                      // Lattice_3D (crystal + BZ mesh)
import qchem.Structure;                       // Structure
import qchem.ScalarFunction;                  // ScalarFunction<double> (the rho(r) face)
import qchem.BasisSet;                        // Real_BS (orbital source), Complex_BS (the Bloch basis)
import qchem.BasisSet.Lattice_3D.BasisSet;    // GPWFactory, GPWParams, RasterPolicy, CellImages
import qchem.Hamiltonian.Factory;             // Pol, VxcFit, the cHamiltonian solid door
import qchem.SCFAccelerator.Factory;          // Type, SolidAcceleratorOptions, the typed solid door
import qchem.SCFIterator;                     // SolidSCFIterator, SCFParams, SCFProgress, EnergyBreakdown
import qchem.ChargeDensity;                   // cDM_CD
import qchem.ChargeDensity.Seed;              // SeedStrategy
import qchem.Mesh;                            // qcMesh::MeshParams / UnitCellKind
import qchem.LASolver;                        // qchem::Ortho

export namespace qchem
{

//! \brief How to set up a periodic calculation.  Designated-initializer friendly:
//!     SolidCalculation calc(lattice, basis, {.Nelec=8, .species={{"Si",4}}});
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
    //! Pseudopotential species as `(element, valence)`.  Multi-species is the norm for a compound; a
    //! per-Z router PP is built so each atom gets its own GTH entry.
    std::vector<std::pair<std::string,int>> species;
    //!@}

    //! \name Grids (the efficiency lever for ionic systems)
    //!@{
    double densityEcut  = -1.0;   //!< <0 AUTO: floor at \c cutoffFactor*alpha_max (recommended).
    double cutoffFactor = 2.0;    //!< C in that floor; 2 = the density's own product exponent.
    double ladderFactor = 4.0;    //!< multigrid ladder gradation (DEPTH is automatic).
    BasisSet::Lattice_3D::RasterPolicy raster = BasisSet::Lattice_3D::RasterPolicy::BallOnly;
    BasisSet::Lattice_3D::CellImages   images = BasisSet::Lattice_3D::CellImages::Periodic;
    rvec3_t kShift = rvec3_t(0,0,0);
    //! \brief The XC real-space quadrature.  DEFAULT \c Auto = "you choose": the class costs a uniform
    //! mesh sized to this run's sharpness against the atom-centred Becke mesh and takes the cheaper, with
    //! a margin (V1.26/V2.4).  Pin \c cellKind to \c Uniform or \c Becke to overrule it and be obeyed --
    //! the grid is a user decision, and the library's job is only to warn about a dominated one.
    qcMesh::MeshParams xcMesh{.cellKind=qcMesh::UnitCellKind::Auto};
    //! WHICH fit basis represents v_xc -- ORTHOGONAL to the grid above.  \c Auto picks Delta whenever the
    //! plane-wave fit cannot serve (Becke grid, or a polarized run).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
    //!@}

    //! \name Convergence machinery
    //!@{
    SCFAccelerators::Type accelerator = SCFAccelerators::Type::DIIS;
    bool globalFermi    = false;  //!< metal: one μ across the k-mesh.
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
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Uniform;
    qchem::Ortho ortho    = qchem::Cholesky;
    double       orthoTol = 0.0;
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
    //!@}
};

//! \brief A converged (or attempted) periodic SCF.
class SolidCalculation
{
public:
    using sf_t     = ScalarFunction<double>;
    using Observer = qchem::SCFIterator::SolidSCFIterator::Observer;

    //! Build the whole graph and converge once with \a params.  \a mol is the MOLECULAR orbital basis the
    //! GPW basis is built over (shared: one per k-block).  \a acc tunes the chosen accelerator; its
    //! defaults are the GPW production recipe.
    SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                     const SolidCalcOptions& opts = {},
                     const SCFParams& params = {},
                     const SCFAccelerators::SolidAcceleratorOptions& acc = {});
    ~SolidCalculation();

    SolidCalculation(const SolidCalculation&)            = delete;   // owns raw resources
    SolidCalculation& operator=(const SolidCalculation&) = delete;

    //! Re-run the SCF (e.g. with tighter tolerances).  Invalidates references handed out by Density().
    bool Converge(const SCFParams& params = {});
    //! Live per-iteration telemetry.  Set BEFORE a Converge() call to observe it (the ctor's initial
    //! convergence runs before any observer can be attached).
    void OnIteration(Observer obs);

    bool                   Converged()      const;
    size_t                 IterationCount() const;
    double                 Energy()         const;   //!< total energy (hartree)
    qchem::EnergyBreakdown EnergyTerms()    const;
    double                 TotalCharge()    const;   //!< the converged electron count (a charge-conservation check)

    //! rho(r) of the converged state.  A periodic density is dcmplx-backed but rho is real, so this is the
    //! same \c ScalarFunction<double> face the molecular facade returns -- sample it, or take its Gradient.
    const sf_t& Density() const;

    //! \brief WHICH XC quadrature this run actually used, after \c Auto was resolved.
    //!
    //! Exposed because the resolution is a real decision made here on the run's behalf, and a caller that
    //! wants to reproduce, pin, or report the run needs to see what was chosen rather than re-deriving it.
    const qcMesh::MeshParams& ResolvedXCMesh() const;

private:
    struct Imp;
    std::unique_ptr<Imp> itsImp;
};

} //export namespace qchem
