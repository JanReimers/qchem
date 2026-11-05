// File: Calculation/Imp/SolidCalculation.C  The periodic SCF front door.
//
// Imports ZERO .Internal. modules -- the same discipline Imp/Calculation.C and Imp/AtomCalculation.C
// already keep, and the reason Step 4 had to open two public doors before this file could exist.
module;
#include <algorithm>   // std::max (the T2 site-moment scan)
#include <cmath>       // std::fabs
#include <complex>     // the response U table (chi0^-1 - chi^-1)
#include <limits>      // quiet_NaN (lastCommutator before any iteration)
#include <map>         // the magnetic decoration's IonicSAD targets
#include <iomanip>     // the stage summary's stated precision
#include <iostream>    // the per-stage anneal banner
#include <cassert>
#include <cstdlib>   // std::getenv -- the banner's thread state
#include <memory>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
module qchem.SolidCalculation;

import qchem.BasisSet.AoShellSource;           // FunctionName: the audit's largest element, by atom, shell and l
import qchem.BasisSet.Conditioning;            // ReportConditioning, MinOverlapEigenvalue
import qchem.BasisSet.Orbital_1E_IBS;          // Real_OIBS / Complex_OIBS (the per-irrep bases to iterate)
import qchem.Pseudopotential.Set;             // PseudopotentialSet -- the run's chosen pseudopotential, resolved once by SolidSetup::Resolve
import qchem.Pseudopotential.GTH_Potentials;  // MakeMultiSpeciesLocal -- the local the ladder is sized for (the Hamiltonian's own builder)
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.ElectronConfiguration.Crystal;    // Crystal_EC
import qchem.Symmetry.Irrep;                   // Irrep, Spin
import qchem.Blaze;                            // blazem::conjs (the u_bands weights)
import qchem.PseudoAtomOrbital;                // the pseudo-atom l orbital: the +U projector (HubbardProjectorPlan D5)
import qchem.Orbitals;                         // TOrbitals/TOrbital -- the ACBN0 feed (EstimateHubbardU)
import qchem.Diagnostics;                      // u_bands: the per-band +U slot weights (QCHEM_DIAGNOSTICS=u_bands[=N])
import qchem.BasisSet.Orbital_DFT_IBS;         // the block an orbital set lives on (EstimateHubbardU)
import qchem.Environment;                   // WarnRetiredEnvironment (D-ENV step 6a)
import qchem.Mesh.XCPolicy;                     // BeckeXCParams / BeckeLadderAbove (a geometry run's XC route and its nested ladder)
import qchem.Parallel;                         // WorkerThreads() -- half of the thread state a row must state
import qchem.SCFCap;                           // SetupSCF -- the seed's atomic solves are not the run's iterations
import qchem.Reporting;                       // report::Timed -- the facade's own setup buckets (1.1, 2026-09-06)
import qchem.SCFRecipeReport;                        // ReportSCFRecipe: the `scf` header every facade files (V-FACADE-PARITY)

namespace qchem
{


//---------------------------------------------------------------------------------------------------
struct SolidCalculation::Imp
{
    //! THE RUN SETUP (doc/CalcNetworkReview.md §2, step 6): the options and every decision resolved from them ONCE -- the
    //! pseudopotential set, the imposition, the decoration, the XC quadrature, the spin counts.  Shared and immutable: everything
    //! below it is this GEOMETRY POINT (structure, Bloch basis, EC, Hamiltonian) or its SCF STATE.
    std::shared_ptr<const SolidSetup>           setup;
    //! The accelerator options a STAGE starts from: \c setup->acc, except where a caller's stage overrides the Ladder's starting
    //! rung for its own duration (the FD chi's restart-like stages).  Per-run state, not a decision of the setup.
    SCFAccelerators::SolidAcceleratorOptions    accOpts;
    //! The caller's per-iteration telemetry: the constructor's (or \c Restart's) \a onIteration, replaced by \c OnIteration.
    Observer                                    onIteration;
    std::shared_ptr<const Structure>            st;
    std::unique_ptr<BasisSet::Complex_BS>       bs;
    std::unique_ptr<Crystal_EC>                 ec;
    //! OWNED HERE (R2.22).  It used to be a bare pointer handed to the iterator, which deleted it -- so
    //! every anneal stage had to build a fresh one (15.5 s/call).  The Hamiltonian is a pure function of
    //! (structure, basis, species, functional, xcMesh, vxcFit), and a stage changes NONE of those, so ONE
    //! serves the whole schedule.  Declared BEFORE `scf` so it outlives every iterator built over it.
    std::unique_ptr<qchem::Hamiltonian::cHamiltonian>     ham;
    //! OWNED HERE, and deliberately REPLACED per stage: a stage invalidates the Pulay/DIIS history and may
    //! change the accelerator TYPE outright (anneal on Ladder, finish on GDM).  Cheap to build, so that is
    //! correct -- unlike the Hamiltonian above, which was only ever rebuilt because of the `delete`.
    std::unique_ptr<SCFAccelerators::SCFAccelerator>      accel;
    std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scf;
    std::unique_ptr<qchem::ChargeDensity::cDM_CD>         cd;    // the converged density (outlives the WF)
    //! m(r) of the converged state.  OWNED: WaveFunction::GetSpinDensity() BUILDS it and hands over the
    //! unique_ptr (V1.25).  EMPTY on an unpolarized run -- that WF does not implement the face at all (V1.17).
    //! m(r), built on the FIRST \c Converged::SpinDensity() call (reset by each Converge): built eagerly it re-collocated
    //! both spin channels' D at the end of every stage, two n×n copies beside the iterator's own state -- the run's peak.
    mutable std::unique_ptr<SolidCalculation::sf_t>       spin;
    SCFAccelerators::Type                       stageAccel = SCFAccelerators::Type::DIIS;  //!< the CURRENT stage's, RESOLVED (banner, [F,D])
    SCFAccelerators::Type                       stageAccelAsked = SCFAccelerators::Type::Auto;  //!< ...as ASKED: what an FD/U-loop stage re-resolves
    bool                                        stageRestored = false;  //!< the CURRENT stage starts from a RESTORED state (Restart: same Hamiltonian)
    //! An EXACT RESUME's saved energy, until the first Converge reports the resumed one beside it (then cleared): the line
    //! belongs to the restored stage's convergence, whoever runs it (\c Restart, or a \c PrepareRestart caller).
    std::optional<double>                       resumedFrom;
    bool   converged = false;
    //! What the LAST converge stage filled with and how well it converged its eigenvalues -- the linear
    //! response's reference is built over the SAME occupancy rule, gated against the SAME measured noise
    //! (doc/LinearResponsePlan.md E1).  The commutator is the final iteration's [F,D].
    qchem::OccupationConfig lastOccupation;
    double lastCommutator = std::numeric_limits<double>::quiet_NaN();
    double lastDensityChange = std::numeric_limits<double>::quiet_NaN();   //!< the last iteration's Δρ (the saved state's residual)
    double charge    = 0.0;
    //! The OUTCOME DETECTORS' record (N1/T3-T4).  Accumulates across every attempt this object makes --
    //! an annealed schedule and a re-Converge are CONTINUATIONS of one SCF, and an order that died in
    //! stage 1 must not become invisible because stage 2 started from the corpse.
    RunDiagnostics diag;
    //! The +U manifolds AS THE HAMILTONIAN CARRIES THEM: resolved by BuildRun, and kept in step with every
    //! ACBN0 outer step's Apply -- so a saved state records the U its density was converged at, not the
    //! options' starting U.
    std::vector<qchem::Hamiltonian::HubbardManifold> hubbard;
    //! The last Converge's recipe (kT, and the banner facts, for the saved state's summary).
    double      lastkT = 0.0;
    std::string lastRecipe;
    //! THE LATTICE, HELD (2026-10-27, user: "SolidCalculation lacks the lattice -- this seems very strange").  It used to be
    //! taken by const& at construction and kept only in FRAGMENTS (the cell matrix for the fingerprint, then a group op list
    //! and site positions for the q-star fold), each later need re-deriving a piece of what the lattice answers.  A copy:
    //! Lattice_3D is an immutable value (copies share its space-group cache).  optional only because the Imp is built first.
    std::optional<Lattice_3D> lat;
    //! The MOLECULAR orbital basis the Bloch basis was built over: held for the grid audit, which rebuilds the Bloch basis on another raster.
    std::shared_ptr<const BasisSet::Real_BS> mol;
    //! THIS GEOMETRY POINT's imposition and raster ladder (V-BS-RASTER), resolved here and handed to the GPW factory: the facade holds
    //! them, so the grid report and the raster pilot ASK them (never a basis block).  Null on a plane-wave run.
    std::shared_ptr<const CrystalImposition>                  imposition;
    std::shared_ptr<const BasisSet::Gaussian::GpwRasterLadder> ladder;
    //! What the run's client consumes (\c SolidCalculation::Consumes): an energy, unless a SolidSession built this point.
    BasisSet::Gaussian::ClientConsumes consumes = BasisSet::Gaussian::ClientConsumes::Energy;
    //! The pilot audit has SETTLED this run's rasters (\c CalibrateRasters): \c setup->opts.tolerances are the final demands.
    bool rastersSettled = false;
    SolidCalculation::RasterPilotCost pilotCost;   //!< what the pilot cost (\c PilotCost)
    //! The pilot's start rasters were already raised once for the overlap's conditioning (\c Converge): never a second time.
    bool startRaisedForConditioning = false;
    //! \brief THE rho-tilde SAMPLER'S ONE TIGHTENING (step 2, user rulings 2026-11-04).  The sampler's tau is fixed from the seed and never
    //! changes inside an SCF; what it delivers is the samp/τ column.  When an SCF of this run ends with a row over target (> 1), the NEXT
    //! SCF BOUNDARY that rebuilds the Hamiltonian (\c BuildRun: the pilot's next round, a re-seat) scales tau by 1/(1.5 x worst) --
    //! ONCE PER RUN and TIGHTEN-ONLY (the scale never rises back), so it can neither oscillate nor creep.  An annealed stage that keeps
    //! its Hamiltonian, or the last SCF of the run, cannot be re-sized: its over-target is REPORTED (console + run report
    //! `xcSamplingOverTarget`: the SCF's first over iteration, the worst delivered/tau, tau), never silent.
    double samplingTauScale   = 1.0;    //!< applied to tau at every Hamiltonian build (the recipe's \c SamplingAccuracy::tauScale)
    bool   samplingTightened  = false;  //!< the one tightening has happened
    double samplingScfWorst   = 0.0;    //!< this SCF's worst samp/τ
    size_t samplingScfFirstOver = 0;    //!< this SCF's first iteration over target (0 = none)
    double samplingPending    = 0.0;    //!< a finished SCF's over-target worst awaiting the next boundary (0 = none)
};

//---------------------------------------------------------------------------------------------------
//  The INTEGRATED order parameter: max_A |mu_A| over the XC quadrature's atom-centred partition.  Empty
//  basins (a uniform XC mesh, or an unpolarized run) give an empty vector, which the caller reads as "not
//  measurable" -- never as "zero", because those are different facts and only one of them is a failure.
//  Every SCF iterate's moments arrive on the EnergyBreakdown the observer already receives
//  (ChargeBreakdown::siteMoments, R1.0h); only the RAW SEED -- which has no energy pass -- is ASKED for
//  through the Hamiltonian's face, once, before Init consumes it.
static double MaxSiteMoment(const rvec_t& m, bool& hasBasins)
{
    hasBasins = hasBasins || m.size()>0;
    double mx=0.0;
    for (size_t a=0;a<m.size();a++) mx=std::max(mx, std::fabs(m[a]));
    return mx;
}

//====================================================================================================
//  THE SELF-DESCRIBING RUN BANNER (doc/OpenWork.md N5/T5, and index item 2's "benchmark protocol").
//
//  WHY IT IS UNCONDITIONAL, and why it is HERE.  doc/Benchmark.md's rows are only comparable if each one
//  declares its thread state and its qchem-only accelerations -- and until now "declare them" was a
//  matter of DISCIPLINE, which is why the low-rank rho route has been ON BY DEFAULT since 07d13bf6 with
//  no row since saying so.  A row that describes itself needs no discipline.  Two measured defects this
//  closes directly: nothing stated the thread counts, and with KerkerG0=0 the fall back to linear
//  D-mixing was ENTIRELY SILENT (the mixer identity appeared only in the Verbose per-iteration column).
//
//  It is printed by the FACADE because the facade is what RESOLVES these choices -- the Auto XC mesh, the
//  accelerator, the spin bookkeeping.  A banner assembled anywhere else would be re-deriving them and
//  hoping the two agree, which is the failure mode it exists to remove.
//====================================================================================================
static const char* MeshName(qcMesh::UnitCellKind k)
{
    switch (k)
    {
        case qcMesh::UnitCellKind::Becke:   return "Becke";
        case qcMesh::UnitCellKind::Uniform: return "Uniform";
        default:                            return "Auto(unresolved)";
    }
}
static const char* SeedName(qchem::ChargeDensity::SeedStrategy s)
{
    using S=qchem::ChargeDensity::SeedStrategy;
    switch (s)
    {
        case S::CoreGuess: return "CoreGuess";
        case S::Uniform:   return "Uniform";
        case S::SAD:       return "SAD";
        case S::IonicSAD:  return "IonicSAD";
        default:           return "Default";
    }
}
static const char* AccelName(SCFAccelerators::Type t)
{
    using T=SCFAccelerators::Type;
    switch (t)
    {
        case T::DIIS:   return "DIIS";
        case T::GDM:    return "GDM";
        case T::Ladder: return "Ladder";
        case T::Null:   return "Null";
        case T::Auto:   return "Auto";
        default:        return "?";
    }
}
// WHAT THE RUN IS MADE OF: system, grids, symmetry, threads, deviations.  Emitted once, from the ctor.
static void EmitRunBanner(const SolidCalcOptions& o, const qcMesh::MeshParams& xc, size_t nAtoms,
                          bool imposed, const std::vector<qchem::Hamiltonian::HubbardManifold>& hubbard)
{
    const char* omp=std::getenv("OMP_NUM_THREADS");
    std::cout<<"["<<o.label<<" run] system: "<<nAtoms<<" atoms, "<<o.Nelec<<" valence e, multiplicity "
             <<o.multiplicity<<(o.multiplicity>=1 ? " (POLARIZED)" : " (unpolarized)")
             <<", seed="<<SeedName(o.seed)<<std::endl;
    std::cout<<"["<<o.label<<" run] grids: densityEcut="
             <<(o.densityEcut<0 ? std::string("auto") : std::to_string(o.densityEcut))
             <<" C="<<o.cutoffFactor<<" raster x(density, V_long, XC L0)=("<<o.tolerances.rasterX<<", "<<o.tolerances.vlocRasterX
             <<", "<<o.tolerances.xcRasterX<<")"
             <<" xcMesh="<<MeshName(xc.cellKind);
    // The radial/angular pair describes an ATOM-CENTRED mesh; printing it beside "Uniform" would be
    // stating a number that had no effect on the run, which is the opposite of self-description.
    if (xc.cellKind==qcMesh::UnitCellKind::Becke) std::cout<<" (nR="<<xc.nRadial<<" L="<<xc.angularDegree<<")";
    else                                          std::cout<<" (eCut="<<xc.eCut<<" Ha)";
    std::cout<<std::endl;
    // The ONE accuracy target and what it implies (doc/GridEpsilonPlan.md §1): the rigorous demand beside the knobs the sizing uses.
    if (o.planeWaveEcut<=0.0) std::cout<<"["<<o.label<<" run] accuracy: "<<o.accuracy.Describe(o.tolerances)<<std::endl;
    std::cout<<"["<<o.label<<" run] symmetry: "
             <<(imposed ? (o.greyImposition ? "IMPOSED (grey -- the erasure control)"
                                            : "IMPOSED (Shubnikov from the decoration)")
                        : "FREE")
             // A VETOED imposition must never be silent -- a caller who asked for one and did not get it
             // would otherwise read the row as though it had, which is the mirror of the hazard that made
             // imposeSymmetry opt-in in the first place.
             <<(o.imposeSymmetry && !imposed ? "  [asked for, VETOED by policy (cp2kCompat / imposeSymmetry=false)]" : "")
             <<";  threads: OMP_NUM_THREADS="<<(omp?omp:"unset")
             <<" ["<<qchem::ThreadSummary()<<"]"
             // THE WAIT POLICY BELONGS ON THIS LINE, not in a footnote: with libomp's 200 ms default
             // spin, 65% of a 12-thread run's billed CPU was the barrier, so a threaded row is
             // uninterpretable unless it says which policy produced it (see StopOmpThreadsBusyWaiting).
             <<";  OMP wait="<<[]{const char* b=std::getenv("KMP_BLOCKTIME"); return b?b:"?";}()
             <<" ms spin"<<std::endl;
    std::cout<<"["<<o.label<<" run] "<<theRunPolicy().Banner()<<std::endl;
    // +U's manifolds are a RUN OPTION, so they are stated here beside the knob banner (the FORM --
    // eigenvalues vs CP2K's diagonal populations -- is the knob QCHEM_U_EIGEN on that banner).  What is
    // declared here beyond the manifolds is the one deviation that is neither: occupations from D_out, the
    // DM-backed source of the mixed density, where CP2K mixes P (same fixed point, different trajectory).
    if (!o.hubbard.empty())
    {
        // ORBITAL RESOLUTION (increment 2): a manifold with Uirrep carries one U per (site irrep, grey
        // parent) slot; the slot table itself -- which irrep is which, its dimension -- is the TERM's to say
        // (it prints the "[+U] site s l=..: U slots" line when it builds the table, pin 17), so the banner
        // states what the RUN asked for: the U vector and the two groups it is resolved under.
        bool resolved=false; for (const auto& M : hubbard) resolved|=!M.Uirrep.empty();
        std::cout<<"["<<o.label<<" run] +U: LOWDIN (CP2K's all-l-shells manifold), "
                 <<(resolved ? "ORBITAL-RESOLVED* (one U per site-irrep slot, see the [+U] slot table)" : "shell-averaged")<<", "
                 <<(theRunPolicy().HubbardEigen() ? "eigenvalue form (Dudarev)" : "DIAGONAL populations (CP2K's form)")
                 <<", occupations from D_out*;  manifolds:";
        for (const auto& M : hubbard)             // the list AS HANDED to the Hamiltonian (site groups filled)
        {
            std::cout<<" (site "<<M.site<<", l="<<M.l;
            if (M.Uirrep.empty()) std::cout<<", U="<<M.U*27.211386245988<<" eV";
            else { std::cout<<", Uirrep="; for (size_t k=0;k<M.Uirrep.size();k++) std::cout<<(k?",":"")<<M.Uirrep[k]*27.211386245988; std::cout<<" eV"; }
            if (M.alpha!=0.0 || M.alphaUp!=0.0 || M.alphaDown!=0.0)   // the static alpha*P perturbation: say it, per spin
                std::cout<<", alpha*P="<<M.alpha*27.211386245988<<" eV (up +"<<M.alphaUp*27.211386245988<<", down +"<<M.alphaDown*27.211386245988<<")";
            if (!M.alphaSlot.empty())   // the slot-resolved part, in the [+U] slot-table order
            { std::cout<<", alphaSlot=+"; for (size_t k=0;k<M.alphaSlot.size();k++) std::cout<<(k?",":"")<<M.alphaSlot[k]*27.211386245988; std::cout<<" eV"; }
            if (M.projector==qchem::Hamiltonian::HubbardProjectorKind::PseudoAtomic)
                std::cout<<", projector: PSEUDO-ATOMIC* ("<<(M.orthoAtomic ? "ORTHO-atomic" : "atomic")<<", independent of the orbital basis)";
            else
                std::cout<<", projector: site shells (CP2K's LOWDIN, every shell)";
            std::cout
                     <<", site group "<<(M.siteOps.empty() ? 1 : M.siteOps.size())<<" ops / grey "
                     <<(M.greyOps.empty() ? (M.siteOps.empty() ? 1 : M.siteOps.size()) : M.greyOps.size())<<")";
        }
        std::cout<<"   [* = differs from CP2K]"<<std::endl;
    }
}
// WHAT THE SCF IS DOING: the mixer, the accelerator, the occupation machinery.  Emitted per Converge,
// because that is where these take effect -- and an anneal changes them stage by stage.
//! \a p is RESOLVED (\c tSCFIterator::ResolveMixing); \a autoMixer = the caller left \c Mixer at \c Auto.
static void EmitSCFBanner(const std::string& label, const SCFParams& p, SCFAccelerators::Type acc, bool autoMixer, bool autoAccel)
{
    using M = SCFParams::Mixing;
    // ⚠ THE PRECONDITIONER AND THE HISTORY ARE INDEPENDENT, AND THE BANNER USED TO HIDE IT (2026-08-28).
    // It printed "Pulay(depth 8, start 5)" whenever PulayDepth>0 and never mentioned Kerker, so a reader
    // of a benchmark row would conclude the G-space preconditioner had been swapped OUT.  It has not:
    // PulayParams carries G0 into PulayMixer, so the two COMPOSE -- which matters
    // because that composition is exactly CP2K's own recipe (its &MIXING BETA is the Kerker damping
    // denominator, applied alongside BROYDEN_MIXING).  This line is the one doc/Benchmark.md tells people
    // to copy beside a row, so it has to say both.
    // A D-MIXER (linear, or Pulay at KerkerG0 = 0 -- PulayDMixer) hands the Fock the mixed density MATRIX: XC sees rho[D].
    const bool dMixer = p.Mixer==M::Linear || (p.Mixer==M::Pulay && p.KerkerG0<=0.0);
    std::cout<<"["<<label<<" scf] mixer: "
             <<(p.Mixer==M::Linear ? std::string("LINEAR D-mixing (no G-space preconditioner)")
              : dMixer             ? std::string("D-mixing (no G-space preconditioner)")
              :                      "Kerker(G0="+std::to_string(p.KerkerG0)+")")
             <<(p.Mixer==M::Pulay ? " + Pulay history(depth "+std::to_string(p.PulayDepth)+", start "
                                    +std::to_string(p.PulayStart)+")"
                                  : std::string())
             <<(autoMixer ? " [mixer: Auto]" : "")
             <<" alpha="<<p.StartingRelaxRo
             <<";  XC rho source: "<<(dMixer                   ? "rho[D] (D-mixing)"
                                     : theRunPolicy().XCFromDM() ? "rho[D] WHOLESALE"
                                     : p.XCCuspDeficit         ? "rho_mix + cusp deficit"
                                     :                           "rho_mix")
             <<";  accel: "<<AccelName(acc)<<(autoAccel ? " [accel: Auto]" : "")
             <<";  kT="<<p.SmearingkT<<" MOM="<<(p.UseMOM?"on":"off")
             <<" NMaxIter="<<p.NMaxIter<<std::endl;
}

//! The per-stage summary the campaign reads its numbers from.  AT A STATED PRECISION: these are the
//! figures doc/Benchmark.md's rows are transcribed from, and 6 s.f. cannot express a sub-mHa delta on a
//! -61 Ha crystal ("A=E-TS=-61.4" was once the whole energy this line reported).
void EmitStageSummary(const std::string& label, size_t s, size_t n, double kT,
                      bool converged, size_t iters, const qchem::EnergyBreakdown& E)
{
    const std::streamsize prec0=std::cout.precision();
    std::cout << "["<<label<<" stage "<<s+1<<"/"<<n<<"] kT="<<kT<<" conv="<<converged<<" iters="<<iters
              << std::setprecision(10)
              << " A=E-TS="<<E.GetTotalEnergy()<<" -TS="<<E["MinusTS"]
              << " E(internal)="<<(E.GetTotalEnergy()-E["MinusTS"])
              << std::setprecision(prec0) << std::endl;
}


//---------------------------------------------------------------------------------------------------
static std::vector<Math::Gaussian> PseudoAtomicFunction(const Structure& st, size_t site, int l,
                                                        const Pseudopotential::PseudopotentialSet& pp);

SolidCalculation::SolidCalculation(BuildOnly) : itsImp(std::make_unique<Imp>()) {}

SolidCalculation::SolidCalculation(std::shared_ptr<const SolidSetup> setup, const SCFParams& params, Observer onIteration)
    : itsImp(std::make_unique<Imp>())
{
    itsImp->onIteration=std::move(onIteration);   // before Start: BuildRun attaches it ahead of iteration 1
    // ★ THE CTOR'S RESIDUE BUCKET (doc/ParallelAndOraclePlan.md 1.1(a)).  Timed is EXCLUSIVE, so this
    // outer scope charges itself exactly what the named setup buckets below do NOT -- the decisions
    // between them (the irrep/EC build, the banner).  With the same bracket on Converge and on Iterate, the ledger
    // PARTITIONS the run and "everything not in a bucket" stops being a place time can hide.
    qchem::report::Timed residue("setup: facade ctor (residue -- decisions between the named buckets)");
    const SCFAccelerators::Type accel=setup->opts.accelerator;
    Start(std::move(setup), params, accel);
    CalibrateRasters(params);   // C2 step 3: the pilot, after the (only) stage
}

void SolidCalculation::Start(std::shared_ptr<const SolidSetup> setup, const SCFParams& params, SCFAccelerators::Type accel)
{
    BuildFirstStage(std::move(setup), params, accel);
    (void)Converge(params);   // the ctor ATTEMPTS; the caller faces the result via Result()
}

void SolidCalculation::BuildFirstStage(std::shared_ptr<const SolidSetup> setup, const SCFParams& params, SCFAccelerators::Type accel)
{
    const Lattice_3D& reference=setup->reference;   // held by the setup, which the point keeps alive
    auto mol=setup->mol;
    BuildPoint(std::move(setup), reference, std::move(mol));
    BuildRun(nullptr, params, accel, /*restored*/false, /*exactResume*/false);
}

std::unique_ptr<SolidCalculation> SolidCalculation::Prepare(std::shared_ptr<const SolidSetup> setup, const std::vector<SCFStage>& schedule,
                                                            Observer onIteration)
{
    if (schedule.empty())
        throw std::runtime_error("SolidCalculation::Prepare: an EMPTY schedule has no meaning -- pass at least one stage.");
    std::unique_ptr<SolidCalculation> c(new SolidCalculation(BuildOnly{}));
    c->itsImp->onIteration=std::move(onIteration);   // before BuildRun: it attaches the observer ahead of iteration 1
    qchem::report::Timed residue("setup: facade ctor (residue -- decisions between the named buckets)");
    c->BuildFirstStage(std::move(setup), schedule.front().params, schedule.front().accelerator);
    return c;
}

namespace {
//! A GPW point's imposition and raster ladder, resolved from the params' own statements of them and the run's PP, and written into
//! \a p so the factory uses (and checks) them.  The ladder is sized for the local PP the Hamiltonian will build: the SAME builder
//! (MakeMultiSpeciesLocal over the PP set's species), so V_loc long's audit row is a pure function of what the facade holds.
void ResolveGpwPoint(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol, const Pseudopotential::PseudopotentialSet& pp,
                     BasisSet::Lattice::GPWParams& p)
{
    namespace L3 = BasisSet::Lattice;
    const std::shared_ptr<const Structure> held=lat.GetStructure();   // held: GetStructure may hand back a fresh object
    const Structure& st=*held;
    p.imposition=std::make_shared<const CrystalImposition>(ResolveImposition(lat, ImpositionOptions{
        .imposeSymmetry=p.imposeSymmetry, .siteSpins=p.siteSpins, .freeMeshInvariant=p.freeMeshInvariant, .kShift=p.kShift}));
    const double ppSharpness=BasisSet::Gaussian::LocalPPSharpness(st, *Pseudopotential::MakeMultiSpeciesLocal(pp.Species(st), pp.Functional()));
    p.ladder=L3::MakeRasterLadder(lat, std::move(mol), p, ppSharpness);
}
} //anon

// THE FIRST HALF OF THE BUILD: the GEOMETRY POINT over a resolved setup -- everything the saved-state fingerprint is judged
// against (the Bloch blocks), and nothing that depends on the seed.  Restart stops here to compare before paying for the
// Hamiltonian.  Every DECISION is the setup's (SolidSetup::Resolve); this builds.
void SolidCalculation::BuildPoint(std::shared_ptr<const SolidSetup> setup, const Lattice_3D& lat,
                                  std::shared_ptr<const BasisSet::Real_BS> mol)
{
    itsImp->setup       = std::move(setup);
    itsImp->accOpts     = itsImp->setup->acc;
    itsImp->st          = lat.GetStructure();
    itsImp->lat.emplace(lat);
    itsImp->mol         = mol;
    const SolidSetup&       S    = *itsImp->setup;
    const SolidCalcOptions& opts = S.opts;

    namespace L3 = BasisSet::Lattice;
    // THE WORKING-TYPE DECISION (doc/RealComplexPlan.md §3, Step 3c-3): a block is real ⇔ its irrep is
    // (TRIM) ∧ every term preserves realness.  This composition root builds the LDA GPW stack --
    // kinetic, the PP trio, Hartree, XC -- every member of which PreservesReal(), so the term half is
    // TRUE here; it is asserted against the constructed Hamiltonian below, so a future term that
    // breaks realness must also flip this forecast.  forceComplex is the §6 ansatz-policy downgrade.
    const bool hamPreservesReal = !opts.forceComplex;
    // ★ THE FACADE'S SETUP IS BUCKETED (doc/ParallelAndOraclePlan.md 1.1, 2026-09-06).  Every MnO row drives
    // SolidCalculation directly, and until now only the phases INSIDE these calls were timed -- so the basis
    // build, the Hamiltonian ctor's non-mesh half, the seed and the ortho landed in no bucket at all.  That
    // was ~48 s of a 119 s threaded MnO run, i.e. the largest non-scaling block in it, and unnamed.  Labels
    // match the GPW_SCF harness's so the two paths' ledgers read the same.
    if (S.planeWave)
    {
        qchem::report::Timed timed("setup: PW basis build");
        itsImp->bs.reset(L3::Factory(L3::Type::PW, lat, opts.planeWaveEcut, S.siteSpins));   // {} unless a magnetic imposition
        std::cout<<"["<<opts.label<<" run] plane-wave orbital basis: Ecut "<<opts.planeWaveEcut<<" Ha, "
                 <<itsImp->bs->GetNumFunctions()<<" PWs at Γ"<<std::endl;
    }
    else
    {
        qchem::report::Timed timed("setup: GPW basis build");
        // THE TOLERANCES are the typed value and nothing else (D-ENV step 6a: the environment is no longer an override layer; the retired
        // GPW_* variables are reported once and ignored).  A non-default value is named on the banner.
        qchem::WarnRetiredEnvironment();
        const auto& tol=opts.tolerances;
        if (!tol.Describe().empty())
            std::cout<<"["<<opts.label<<" run] GPW tolerances (non-default):"<<tol.Describe()<<std::endl;
        L3::GPWParams p{
            .densityEcut = opts.densityEcut, .cutoffFactor = opts.cutoffFactor,
            .images = opts.images, .kShift = opts.kShift, .ladderFactor = opts.ladderFactor,
            .imposeSymmetry = S.imposed, .siteSpins = S.siteSpins, .freeMeshInvariant = opts.xcMeshSymmetric, .tol = tol,
            .hamPreservesReal = hamPreservesReal};
        ResolveGpwPoint(lat, mol, S.pp, p);
        itsImp->imposition=p.imposition;
        itsImp->ladder    =p.ladder;
        itsImp->bs.reset(L3::GPWFactory(lat, mol, p));
    }
    // THE RUN REPORTS ITS OWN BASIS AND GRIDS (TE phase 2, 2026-09-15; the reporting rule: each class emits
    // contemporaneously with its own activity).  These two sections used to be emitted by the integration
    // test driver around this same construction; the facade is the orchestrator now.  Order matters and is
    // the one the driver established: the conditioning pre-flight on the ANALYTIC overlap first (no grids
    // needed), then the grid ladder's report (the point's ladder, ASKED: V-BS-RASTER; only with a run report open).
    // Never an abort: a dependent basis is the caller's to
    // handle through a pivoted ortho, and an unhandled one fails loudly in the iterator.
    // (Plane waves are orthonormal and grid-free: neither section has anything to say about them.)
    if (!S.planeWave)
    {
        qchem::report::Timed timed("setup: vet basis (analytic S) + grids report");
        size_t nRemoved=0;
        {
            qchem::report::Log("vetting basis conditioning");
            qchem::report::Section basis("basis");
            qchem::report::Set("nFunctions", (long)itsImp->bs->GetNumFunctions());
            nRemoved=BasisSet::ReportConditioning(*itsImp->bs);
        }
        if (nRemoved>0)
            qchem::report::Log("basis is rank-deficient ("+std::to_string(nRemoved)+" redundant functions, see basis.removed)");
        if (qchem::report::Depth()>0)
        {
            qchem::report::Log("the grid ladder");
            qchem::report::Section grids("grids");   // the point's grids: the ladder and the lattice sums report one level in
            L3::EmitGpwGrids(*itsImp->ladder, *itsImp->mol, *itsImp->imposition);
        }
    }

    // The electron configuration over THIS basis's k-blocks, from the setup's spin counts (decision 2).
    auto irreps = itsImp->bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry Sum_k)
    itsImp->ec = std::make_unique<Crystal_EC>(irreps, S.nUp, S.nDown, !opts.perKFermiLevels, opts.spinsShareFermi);
    itsImp->hubbard = opts.hubbard;   // the U the run will carry (BuildRun fills the groups and radials)
    // A run that will SAVE its state must be able to, and must say so NOW -- not after the SCF it was asked to protect (OpenWork
    // §3 "PW states cannot be saved", 2026-10-26: a 52-minute converged run banked nothing).  The fingerprint is exactly what a
    // save writes first, and everything it reads is built by this point.
    if (!opts.saveStateTo.empty())
        try { (void)Fingerprint(0.0); }
        catch (const std::logic_error& e)
        {
            throw std::invalid_argument("SolidCalcOptions.saveStateTo: this run cannot save its state ("+std::string(e.what())+")");
        }
}


// THE SECOND HALF: the Hamiltonian, the seed (the caller's when given -- a restored state), the iterator.
void SolidCalculation::BuildRun(std::unique_ptr<qchem::ChargeDensity::cChargeDensity> seed, const SCFParams& params,
                                SCFAccelerators::Type accel, bool restored, bool exactResume)
{
    const SolidSetup&       S    = *itsImp->setup;
    const SolidCalcOptions& opts = S.opts;
    const Lattice_3D&       lat  = *itsImp->lat;
    const std::vector<int>& siteSpins = S.siteSpins;
    const bool polarized        = S.polarized;
    const bool hamPreservesReal = !opts.forceComplex;   // the BuildPoint forecast, restated for the assert below
    const bool imposed          = S.imposed;

    // THE HUBBARD MANIFOLDS' SITE GROUPS (step 5 increment 2): each manifold is labelled by the stabiliser of
    // its site in the STATED DECORATION's Shubnikov group (sigma=None ops), grey parentage beside it -- the
    // same decoration the imposition uses, whether or not the run imposes: the site symmetry of the ordered
    // state is the physical question either way.  No longer derived from the seed (OpenWork §3, 2026-10-26):
    // a polarized +U run was refused at SolidSetup::Resolve without siteSpins.  A caller who filled siteOps by hand keeps them.
    std::vector<qchem::Hamiltonian::HubbardManifold> hubbard=opts.hubbard;
    if (!hubbard.empty())
    {
        std::vector<int> decoration=siteSpins;
        if (decoration.empty() && polarized && !opts.greyImposition) decoration=opts.siteSpins;
        assert(!polarized || opts.greyImposition || !decoration.empty());
        for (auto& M : hubbard)
        {
            if (M.projector==qchem::Hamiltonian::HubbardProjectorKind::PseudoAtomic && M.atomicFunction.empty())
                M.atomicFunction=PseudoAtomicFunction(*itsImp->st, M.site, M.l, itsImp->setup->pp);
            if (M.siteOps.empty()) M.siteOps=lat.SiteRotations(M.site, decoration);
            // PARENTAGE is chemistry: the point group of the site's coordination polyhedron, not the
            // cell's grey stabiliser (which on the rhombohedral AFM-II supercell is D_3d too, and would
            // name nothing -- measured 2026-09-21).  Every site op is among these, so the slots are a
            // clean refinement.
            if (M.greyOps.empty()) M.greyOps=lat.SiteEnvironmentRotations(M.site);
        }
    }
    {
        qchem::report::Timed timed("setup: hamiltonian ctor (fit bases + becke mesh)");
        qchem::report::Scope hamiltonian("hamiltonian");   // its meshes and terms report one level in (user rulings 2026-11-02, 2026-11-04)
        // THE ACCURACY a Fourier series sampled on the XC quadrature must meet (the rho-tilde sampler, step 2; pin 31): the run's
        // target for what its client consumes, the density gate its sampled residual is judged on, and the two GPW calibrations.
        // The XC term derives the one tolerance tau from it on the seed.  No target (eps 0) leaves it unset: the step-1 fixed width.
        qcMesh::MeshParams xc=itsImp->setup->xcMesh;
        if (opts.accuracy.eps>0.0)
        {
            // THE ONE TIGHTENING, at this SCF boundary, when a finished SCF of the run measured the sampling over target (see the Impl).
            if (itsImp->samplingPending>1.0 && !itsImp->samplingTightened)
            {
                const double old=itsImp->samplingTauScale;
                itsImp->samplingTauScale=std::min(old, 1.0/(1.5*itsImp->samplingPending));
                itsImp->samplingTightened=true;
                std::cout<<"["<<opts.label<<" XC sampling] TIGHTENED once: the last SCF delivered "<<itsImp->samplingPending
                         <<" x tau, so tau is scaled by "<<itsImp->samplingTauScale<<" (was "<<old<<") from this SCF on"<<std::endl;
                qchem::report::Scope here("xcSamplingTightened");
                qchem::report::Set("worstDeliveredOverTau", itsImp->samplingPending);
                qchem::report::Set("tauScaleOld", old);
                qchem::report::Set("tauScaleNew", itsImp->samplingTauScale);
            }
            itsImp->samplingPending=0.0;
            xc.sampling.tauScale=itsImp->samplingTauScale;
            const bool onH=itsImp->consumes==BasisSet::Gaussian::ClientConsumes::MatrixElements;
            xc.sampling.target        =opts.accuracy.Target(itsImp->consumes, itsImp->st->GetNumAtoms());
            xc.sampling.perElement    =onH;
            xc.sampling.residualGate  =(params.Δρmeasure==SCFParams::Measure::L1PerElectron && params.MinΔρPerElectron)
                                       ? *params.MinΔρPerElectron : 0.0;   // a MaxΔD gate is not an L1 of the sampled pair
            xc.sampling.residualMargin=opts.accuracy.samplerResidualMargin;
            xc.sampling.eta           =opts.accuracy.samplerEta;
        }
        itsImp->ham.reset(qchem::Hamiltonian::Factory(
            polarized ? qchem::SpinGroup::Polarized : qchem::SpinGroup::UnPolarized,
            itsImp->st, itsImp->bs.get(), itsImp->setup->pp, xc, opts.vxcFit, hubbard));
    }
    // The forecast crosscheck: the basis was built on the promise that every term preserves realness
    // (the AND's term half, above); the constructed Hamiltonian must agree, or real blocks were built
    // that its terms cannot serve.
    assert((!hamPreservesReal || itsImp->ham->PreservesReal()) &&
           "SolidCalculation: the term stack no longer preserves realness -- update the forecast above");

    // DECISION 3 -- the accelerator, by policy (Auto pairs it with the mixer), through the public typed door.
    MakeAccelerator(accel, params, restored);

    // THE SEED IS BUILT HERE, not inside the iterator -- the same factory call the SeedStrategy ctor
    // would have made (ChargeDensity::MakeSeedDensity with the Hamiltonian's own polarization), handed
    // to the explicit-seed ctor instead.  WHY THE FACADE TAKES IT OVER (N1/T2, 2026-08-26): the
    // postcondition needs the order the run STARTED with, and the iterator consumes the seed inside
    // Init -- which builds a Fock from it, diagonalizes and FILLS, so the earliest density a caller can
    // reach is already one aufbau fill downstream.  For a run whose order dies in that very first fill
    // the difference is the whole measurement: a Na2 seed staggered at +/-1 e reads +/-0.07 e by the
    // time Init hands its density back, which is below any honest floor and made the postcondition
    // silently skip.  Measured before it is consumed, the baseline is the seed's own.
    const bool polarizedHam = itsImp->ham->GetSpinGroup()==SpinGroup::Polarized;
    const auto seedRole = !seed       ? qchem::SCFIterator::SeedRole::Guess
                        : exactResume ? qchem::SCFIterator::SeedRole::ExactResume
                        :               qchem::SCFIterator::SeedRole::WorkingDensity;   // restored = carried
    if (!seed)   // a RESTORED state arrives pre-built (Restart); otherwise the strategy the options name
    {
        qchem::report::Timed timed("setup: seed density (SAD/IonicSAD atomic solves)");
        qchem::SetupSCF setup;   // the atomic solves are setup, not the run's convergence (a panel budget does not count them)
        seed.reset(qchem::ChargeDensity::MakeSeedDensity<dcmplx>(opts.seed, itsImp->bs.get(), itsImp->st.get(),
                                                                 itsImp->ec.get(), polarizedHam,
                                                                 itsImp->setup->seedKey));
    }
    // MEASURED ONLY WHERE IT CAN MEAN SOMETHING.  SiteMoments rasters BOTH spin channels before it can
    // discover it has no basins to integrate over, and unlike the per-iteration probe -- which rides a
    // raster the Fock build has already made for this density serial -- nothing has been built yet at
    // seed time, so this one would be paid in full.  An unpolarized run has m == 0 identically and a
    // uniform XC mesh has no basins at all, and the facade knows both facts here without asking.
    if (seed && polarizedHam && itsImp->setup->xcMesh.cellKind==qcMesh::UnitCellKind::Becke)
    {
        // "paid in full" (see above) -- so it gets its own bucket rather than hiding inside the seed's.
        qchem::report::Timed timed("setup: seed order probe (site moments)");
        itsImp->diag.itsSeedOrder = MaxSiteMoment(itsImp->ham->SiteMoments(seed.get()), itsImp->diag.itsHasBasins);
    }

    {
        // The iterator's Init builds a Fock from the seed, diagonalizes and fills -- so the FIRST Fock's
        // lazy heavy builds (collocation task list, local-PP sweep, KB, analytic 1E) are children here.
        qchem::report::Timed timed("setup: seed + ortho (iterator ctor)");
        itsImp->scf = std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
            itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
            seed.release(), itsImp->st.get(), opts.ortho, opts.orthoTol, seedRole, opts.numBands);   // the iterator consumes it in Init
    }

    // Observe from iteration ONE: the ctor converges, so an observer attached afterwards has already
    // missed stage 0 (hence the constructors' onIteration argument).
    AttachProbes();

    // MOM CONTINUATION FROM THE SEED (S0e): pin the reference to the seed's OWN freshly-filled occupied
    // subspace before iteration 1, so the CONFIGURATION the seed chose survives, not merely its density.
    // No-op unless SCFParams::UseMOM is also set.
    if (opts.momFromSeed) itsImp->scf->AdoptMOMReference(*itsImp->scf->GetWaveFunction());

    // (T2's gate -- itsImp->setup->imposed -- was resolved by SolidSetup::Resolve: an imposition is an ASSERTION about the answer,
    // and only then is losing the order a contradiction rather than physics.)

    EmitRunBanner(opts, itsImp->setup->xcMesh, itsImp->st->GetNumAtoms(), imposed, hubbard);
    itsImp->hubbard = hubbard;        // RESOLVED (groups, radials): what the Hamiltonian carries and a save records
}

// THE ANNEALED CTOR.  Delegates to the single-stage one with the FIRST stage's parameters/accelerator --
// so the graph is built exactly once and stage 0 runs as the plain ctor's convergence -- then continues
// through the rest of the schedule.  (Delegating rather than duplicating the build is what keeps the two
// ctors from drifting; every DECISION in the build is made in one place.)
SolidCalculation::SolidCalculation(std::shared_ptr<const SolidSetup> setup, const std::vector<SCFStage>& schedule, Observer onIteration)
    : itsImp(std::make_unique<Imp>())
{
    itsImp->onIteration=std::move(onIteration);
    if (schedule.empty())
        throw std::runtime_error("SolidCalculation: an EMPTY anneal schedule has no meaning -- pass at "
                                 "least one stage, or use the single-SCFParams constructor.");
    {
        qchem::report::Timed residue("setup: facade ctor (residue -- decisions between the named buckets)");
        Start(std::move(setup), schedule.front().params, schedule.front().accelerator);
    }
    RunRestOf(schedule);
    CalibrateRasters(schedule.back().params);   // C2 step 3: the pilot, after the FINAL stage
}

// Stage 0 has already run (Start converged it) -- report it, then continue through the rest, each stage seeded from the one
// before it.
void SolidCalculation::RunRestOf(const std::vector<SCFStage>& schedule)
{
    EmitStageSummary(itsImp->setup->opts.label, 0, schedule.size(), schedule.front().params.SmearingkT,
                     itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->scf->GetEnergy());
    for (size_t s=1; s<schedule.size(); ++s)
    {
        BuildStage("annealed schedule stage "+std::to_string(s), schedule[s].accelerator, std::move(itsImp->cd), schedule[s].params);
        (void)Converge(schedule[s].params);
        EmitStageSummary(itsImp->setup->opts.label, s, schedule.size(), schedule[s].params.SmearingkT,
                         itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->scf->GetEnergy());
    }
}

SolidCalculation::~SolidCalculation() = default;

//---------------------------------------------------------------------------------------------------
// THE ANSWERS LIVE ON THE PROOF (doc/OpenWork.md N1/T1).  Each of these used to sit on SolidCalculation
// itself with no precondition, so a run that never converged served a plausible number to anyone who did
// not think to ask Converged() first.  Reaching them now REQUIRES having been handed a Converged, and the
// only source of one is a successful attempt.
double SolidCalculation::Converged::Energy()      const {return itsImp->scf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown SolidCalculation::Converged::EnergyTerms() const {return itsImp->scf->GetEnergy();}
double SolidCalculation::Converged::TotalCharge() const {return itsImp->charge;}
size_t SolidCalculation::Converged::IterationCount() const {return itsImp->scf->GetIterationCount();}

const ScalarFunction<double>* SolidCalculation::Converged::SpinDensity() const
{
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!itsImp->spin && wf->GetSpinGroup()==SpinGroup::Polarized) itsImp->spin=wf->GetSpinDensity();
    return itsImp->spin.get();
}
const qchem::ChargeDensity::cDM_CD& SolidCalculation::Converged::DensityMatrix() const {return *itsImp->cd;}

const ScalarFunction<double>& SolidCalculation::Converged::Density() const
{
    // Still an assert, and legitimately so: reaching here without a density would be a BROKEN INVARIANT
    // (a Converged is only ever minted beside one), not a user error.  The user-error case is what the
    // type now prevents outright.
    assert(itsImp->cd && "SolidCalculation::Converged::Density: converged without a density");
    return *itsImp->cd;
}

//====================================================================================================
//  THE OUTCOME DETECTORS (doc/OpenWork.md N1/T3-T4).  Every rule below judges the run against ITS OWN
//  trajectory, never against a constant: the quantities are extensive, so an absolute threshold would
//  be one cell's number wearing a library's clothes.
//====================================================================================================
double RunDiagnostics::OrderPeak() const
{
    double mx=itsSeedOrder;                       // the RAW seed counts: see the header on why both ends mislead
    for (size_t i=0;i<itsOrder.size();i++) mx=std::max(mx, itsOrder[i]);
    return mx;
}
double RunDiagnostics::OrderFinal() const {return itsOrder.empty() ? 0.0 : itsOrder.back();}
bool   RunDiagnostics::HasOrder()   const {return itsHasBasins && OrderPeak() > kOrderFloor;}

// DEAD == below 1% of the high-water mark AND STAYING there.  The "staying" half is what separates a
// collapse from a zero crossing: an order parameter that dips through zero on its way to the other sign
// is not dead, and a single small iterate proves nothing.
size_t RunDiagnostics::OrderDiedAt() const
{
    if (!HasOrder() || itsOrder.empty()) return 0;
    const double dead=kCollapseFraction*OrderPeak();
    size_t died=itsOrder.size();                  // first index from which |order| stays below dead
    while (died>0 && std::fabs(itsOrder[died-1])<=dead) --died;
    return died<itsOrder.size() ? died+1 : 0;     // 1-based STEP number (see the header); 0 == never died
}
bool RunDiagnostics::OrderCollapsed() const {return OrderDiedAt()>0;}

bool   RunDiagnostics::HasHartree()   const {return itsEee.size()>1;}
double RunDiagnostics::HartreeFloor() const
{
    double mn=0.0;
    for (size_t i=0;i<itsEee.size();i++) mn = (i==0) ? itsEee[i] : std::min(mn, itsEee[i]);
    return mn;
}
double RunDiagnostics::HartreePeak() const
{
    double mx=0.0;
    for (size_t i=0;i<itsEee.size();i++) mx=std::max(mx, itsEee[i]);
    return mx;
}
// THE RATIO IS TAKEN AT THE END, NOT AT THE PEAK, and that is the whole design of this detector.  Every
// SCF passes through a transient on its way out of the seed -- MEASURED on the healthy MnO baseline,
// 2026-08-26: Eee runs 14.17, 15.23, 12.51 and then settles at 13.18, so a peak-based ratio would read
// 1.22 on a run that did nothing wrong, and on a slower system it would read far more.  A run that
// overshoots and COMES BACK is a healthy run doing its job.  What the measured collapses have in common
// is not a spike but a LEVEL: they finish high and stay there (13.48 healthy against 29.0 and 35.1),
// because the low-G charge mode is no longer damped and the density has genuinely piled up.  So: where
// the run ENDED, against the lowest it ever managed.
double RunDiagnostics::SloshRatio() const
{
    const double floor=HartreeFloor();
    return (!HasHartree() || floor<=0.0) ? 1.0 : itsEee.back()/floor;
}
// !itsConverged is part of the PREDICATE, not of the caller's discipline (see the header): a converged
// density is stationary, so "it sloshed" is not a thing that can be true of it.
bool RunDiagnostics::ChargeSloshed() const
{return !itsConverged && HasHartree() && SloshRatio() > kSloshFactor;}

std::string RunDiagnostics::Summary() const
{
    std::ostringstream os;
    os<<std::setprecision(4);
    if (HasOrder())
    {
        os<<"order(integrated site moment): seed "<<itsSeedOrder<<" e, peak "<<OrderPeak()
          <<" e, final "<<OrderFinal()<<" e";
        if (const size_t d=OrderDiedAt(); d>0) os<<" -- DIED at step "<<d;
        else                                   os<<" -- SURVIVED";
    }
    else if (itsHasBasins) os<<"order: none to lose (peak "<<OrderPeak()<<" e, below the "
                             <<RunDiagnostics::kOrderFloor<<" e floor)";
    else                   os<<"order: not measurable (no atom-centred basins on this XC mesh)";
    if (HasHartree())
        os<<"; Eee "<<HartreeFloor()<<" -> "<<itsEee.back()<<" Ha (peak "<<HartreePeak()
          <<", end/floor "<<SloshRatio()<<")";
    return os.str();
}

//---------------------------------------------------------------------------------------------------
//! Mint the outcome from the state the last attempt left behind -- one place, so Converge() and Result()
//! cannot drift apart in what they call a success.
//!
//! THE ORDER OF THE CHECKS IS THE DESIGN.  A collapsed run usually trips more than one of them, and the
//! caller reads exactly one \c why, so the list runs most-mechanistic first: "the Hartree term ran away"
//! tells you WHAT to fix, "it lost the order it imposed" tells you the answer is to a different question,
//! and "it hit the iteration cap" -- the only one of the three that names no mechanism -- is last of the
//! three failures a collapse can trip.  Whichever fires, \c details carries the whole post-mortem, so
//! nothing that was measured is lost to the ranking.
Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Outcome_() const
{
    using O = Outcome<Converged, SCFFailure>;
    const RunDiagnostics& d = itsImp->diag;
    auto fail=[&](SCFFailure::Why why, std::string what)
    {
        SCFFailure f;
        f.why        = why;
        f.iterations = itsImp->scf->GetIterationCount();
        f.lastEnergy = itsImp->scf->GetEnergy().GetTotalEnergy();
        f.details    = std::move(what);
        if (const std::string s=d.Summary(); !s.empty()) f.details += "  [" + s + "]";
        return O::Fail(std::move(f));
    };

    // ★ T3 -- THE CHARGE-SLOSH DETECTOR, as a REFINEMENT OF NON-CONVERGENCE.  "Ran out of iterations" is
    // true of every collapse and explains none of them; Eee is the term the low-G charge failure MOVES,
    // so it is what turns the useless half of the story into a mechanism.  Measured 2026-08-25 it ordered
    // four MnO collapses correctly without having been designed for any of them -- INCLUDING the one it
    // must stay silent on (the linear-D-mix arm, whose moment died with Eee unmoved at 13.6: same
    // symptom, different mechanism, and conflating the two is the trap this detector must not fall into).
    // ⚠ ChargeSloshed() is false on a converged run BY CONSTRUCTION: see RunDiagnostics::ChargeSloshed.  A converged
    // density is stationary, so there is nothing sloshing -- while a healthy run that RESTRUCTURES on its
    // way to the answer can raise Eee a long way (Na2: 1.71x, converged and correct).  Confining the
    // detector to runs that already failed means a mistake here costs a LABEL, never a good answer.
    if (!itsImp->converged && !itsImp->scf->UnphysicalStop().empty())
        return fail(SCFFailure::Why::UnphysicalEnergy, itsImp->scf->UnphysicalStop());
    if (!itsImp->converged)
        return d.ChargeSloshed()
             ? fail(SCFFailure::Why::ChargeSlosh,
                    "the Hartree term ran away and the SCF never converged: Eee ended at "
                    +std::to_string(d.SloshRatio())+"x the lowest value this run reached, which is the "
                    "signature of an UNDAMPED low-G charge mode -- the density has piled up, so the last "
                    "iterate's energy is not comparable with a healthy run's")
             : fail(SCFFailure::Why::NotConverged,
                    "the SCF reached its iteration limit with the residual still above tolerance");

    // ★ T2/T4 -- THE POSTCONDITION ON AN IMPOSITION.  Imposing a MAGNETIC (Shubnikov) group asserts that
    // the solution carries that order.  A run that imposes it and then loses the order has not found a
    // worse answer -- it has answered a DIFFERENT QUESTION than it was asked, and its energy is not
    // comparable with the one that was wanted.  So it is a FAILURE, not a diagnostic.
    // Gated three ways so it cannot false-positive: the run must have IMPOSED, the order must have been
    // MEASURABLE (Becke basins; a uniform mesh has none and the check correctly skips), and the run must
    // have CARRIED order at some point.  A FREE run that finds m=0 is PHYSICS and is never touched.
    // T4 is the same detector doing the second half of its job: the trajectory rule ("rose, then stayed
    // below 1% of its peak") is what the MnO test used to compute for itself, so a moment that dies
    // MID-RUN is caught even when the final iterate is not the evidence.
    if (itsImp->setup->imposed && d.OrderCollapsed())
        return fail(SCFFailure::Why::OrderLost,
                    "the run IMPOSED a magnetic symmetry and did not keep it: the integrated site moment "
                    "peaked at "+std::to_string(d.OrderPeak())+" e and was dead ( < "
                    +std::to_string(RunDiagnostics::kCollapseFraction)+" of that) from step "
                    +std::to_string(d.OrderDiedAt())+" onward, so this energy answers a different "
                    "question than the one asked");

    return O::Ok(Converged(itsImp.get()));
}

// THE ACCELERATOR: the caller's, or the Auto policy's (SolidCalcOptions::accelerator).  The iterator then resolves an Auto mixer
// against the accelerator built here (tSCFIterator::ResolveMixing asks ExtrapolatesFock()), so the two halves cannot disagree.
void SolidCalculation::MakeAccelerator(SCFAccelerators::Type asked, const SCFParams& params, bool restored)
{
    using T=SCFAccelerators::Type;
    // THE Auto POLICY (OpenWork §3 "Kerker/Pulay as the DEFAULT mixer"; DEFAULT 2026-10-26, user: D-space Pulay):
    //   a SEEDED or CARRIED stage  -> Null, so an Auto mixer resolves to Pulay on D (ResolveMixing: no Fock history beside a
    //                                 density history -- they extrapolate at once and fight, measured 2026-10-22).  A CARRIED
    //                                 stage (anneal, U loop, +-alpha) may have a NEW Hamiltonian or kT: its density is a good
    //                                 start, not the answer -- GDM there false-converges on its own small steps (OpenWork §3).
    //   a RESTORED stage at kT = 0 -> the Ladder on its GDM rung: SAME Hamiltonian, the density IS the answer (NaF 9 -> 2);
    //   a RESTORED smeared stage   -> DIIS (GDM would minimise an integer D' and never apply the smearing -- see Restart);
    //   a RESTORED partial-eigensolve stage (numBands > 0, a plane-wave run) -> DIIS: GDM rotates into virtuals never built.
    // Beside either Fock extrapolator an Auto mixer stays LINEAR, its step capped at 0.3 (Converge: a restart's linear step).
    SCFAccelerators::SolidAcceleratorOptions ao=itsImp->accOpts;
    T t = asked;
    if (asked==T::Auto)
    {
        if (!restored)                    t=T::Null;
        else if (params.SmearingkT==0.0 && itsImp->setup->opts.numBands==0) { t=T::Ladder; ao.ladderStartOnGDM=true; }
        else                              t=T::DIIS;   // smeared, or a PARTIAL eigensolve (GDM needs every virtual)
    }
    // A direct minimiser rotates the occupied block into the WHOLE virtual space, which a partial eigensolve never built.
    if (itsImp->setup->opts.numBands>0 && (t==T::GDM || t==T::Ladder))
        throw std::invalid_argument("SolidCalcOptions: numBands (a partial eigensolve) beside the GDM / Ladder accelerator -- "
                                    "they need every virtual; use Null or DIIS, or numBands 0");
    itsImp->stageAccelAsked = asked;
    itsImp->stageAccel      = t;
    itsImp->stageRestored   = restored;
    // The facade sets the context, the accelerator files what it was built with (THE RUN REPORT): one `accelerators` row per
    // accelerator built in this run (a solid rebuilds one per stage), its Ladder rungs one level further in.
    qchem::report::Row built("accelerators");
    itsImp->accel.reset(SCFAccelerators::Factory(t, ao));
}

// ONE STAGE: a FRESH Hamiltonian + accelerator (a kT change must not carry stale DIIS history across the
// re-seed), an iterator seeded from \a carried when there is one, and MOM continuation from \a prev.
// Ownership: the iterator OWNS ham and accel and deletes them, so \a prev must outlive the adoption and is
// released immediately after it -- the same ordering the driver this replaces used.
void SolidCalculation::BuildStage(const std::string& why, SCFAccelerators::Type accType,
                                  std::unique_ptr<qchem::ChargeDensity::cDM_CD> carried, const SCFParams& params)
{
    // ★ R2.22: A STAGE NO LONGER REBUILDS THE HAMILTONIAN.  It used to, and NOT for a physics reason --
    // the iterator deleted the one it had been handed, so the facade had no choice but to build another.
    // The Hamiltonian is a pure function of (structure, basis, species, functional, xcMesh, vxcFit) and a
    // stage changes none of them; at 15.5 s/call that rebuild was 19% of a threaded MnO run
    // (doc/ParallelAndOraclePlan.md 1.1(b), doc/CleanupCandidates.md R2.22).  `itsImp->ham` now simply
    // persists, and the ledger's `setup: hamiltonian ctor` reads [x1] for any schedule length.
    //
    // What a stage DOES rebuild, and why each one is right: the ACCELERATOR (its Pulay/DIIS history is
    // invalid across a re-seed, and the TYPE itself changes -- anneal on Ladder, finish on GDM) and the
    // ITERATOR (whose ctor Inits, so the stage's first Fock is a child of this bucket).  Both are cheap;
    // the residue bucket below measures exactly that and has stayed ~0.36 s of a 121 s run.
    qchem::report::Timed residue("setup: anneal stage rebuild (residue -- accel + iterator + MOM adopt)");
    // R1.3: a stage boundary is a hand-off -- every history (accelerator, mixer) is dropped and rebuilt over the carried D.
    PrintHandOff(0, {"stage "+std::string(AccelName(itsImp->stageAccel)), "fresh stage", why},
                 " [accelerator, mixer and iterator rebuilt over the carried density]");
    MakeAccelerator(accType, params, /*restored*/false);   // a carried stage may have a new Hamiltonian / kT: not a restore

    auto prev = std::move(itsImp->scf);      // held ONLY until the new stage has copied its MOM reference
    itsImp->scf = carried
        ? std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
              carried.release(), itsImp->st.get(), itsImp->setup->opts.ortho, itsImp->setup->opts.orthoTol,
              qchem::SCFIterator::SeedRole::WorkingDensity, itsImp->setup->opts.numBands)
        : std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
              itsImp->setup->opts.seed, itsImp->st.get(), itsImp->setup->opts.ortho, itsImp->setup->opts.orthoTol, itsImp->setup->opts.numBands);
    // MOM continuation across TEMPERATURE: stage 0 self-adopts the seed's own freshly-filled occupied
    // subspace, every later stage adopts the stage before it -- so the CHARACTER the hot stage settled on
    // survives the fresh wavefunction, exactly as the density does.
    if (itsImp->setup->opts.momFromSeed)
        itsImp->scf->AdoptMOMReference(prev ? *prev->GetWaveFunction() : *itsImp->scf->GetWaveFunction());
    prev.reset();                            // the adoption copied what it needed
    // EVERY stage is observed, not just the first: this iterator is brand new, so the telemetry has to be
    // re-attached or an annealed run goes quiet after stage 0 -- silently, which is the worst kind.  The
    // DETECTORS ride the same hooks, so this is also what keeps the trajectory continuous across a stage
    // boundary (an order that died in stage 1 must not vanish because stage 2 got a fresh iterator).
    AttachProbes();
}

//---------------------------------------------------------------------------------------------------
//! ONE HOOK PAIR, serving the caller's telemetry AND the outcome detectors.
//!
//! The observer FILES what the iteration MEASURED: the integrated site moment rides the EnergyBreakdown
//! (R1.0h), so the detectors read it off the progress record -- no probe, no second sampling -- and the
//! caller's observer is composed behind the facade's own, so attaching telemetry late cannot disarm them.
void SolidCalculation::AttachProbes()
{
    auto userObs = itsImp->onIteration;
    itsImp->scf->SetObserver([this,userObs](const qchem::SCFIterator::IterationRow& row)
    {
        // The detectors read the observable off the breakdown the iteration produced (R1.0h) -- the same
        // number the trace reported, not a second pull.
        itsImp->diag.itsOrder.push_back(MaxSiteMoment(row.eb.charge.siteMoments, itsImp->diag.itsHasBasins));
        itsImp->diag.itsEee  .push_back(row.eb["Eee"]);
        // A Null accelerator computes NO [F,D]: its cell carries 0 (the trace hides it by tag), which
        // must not read as a MEASURED zero eigenvalue noise (found 2026-09-27, the NiO R0 gate) -- NaN instead.
        itsImp->lastCommutator = itsImp->stageAccel==SCFAccelerators::Type::Null
                               ? std::numeric_limits<double>::quiet_NaN() : row.ValueOf("[F,D]");
        itsImp->lastDensityChange = row.ValueOf("Δρ");
        if (const auto& r=row.eb.charge.samplingOverTau; r)   // the rho-tilde sampler's audit (the samp/τ column)
        {
            itsImp->samplingScfWorst=std::max(itsImp->samplingScfWorst, *r);
            if (*r>1.0 && itsImp->samplingScfFirstOver==0) itsImp->samplingScfFirstOver=row.iteration;
        }
        if (userObs) userObs(row);
    });
}

Outcome<SolidCalculation::Converged, SCFFailure>
SolidCalculation::Converge(const std::vector<SCFStage>& schedule)
{
    if (schedule.empty())
        throw std::runtime_error("SolidCalculation::Converge: an EMPTY anneal schedule has no meaning -- "
                                 "pass at least one stage, or use the single-SCFParams overload.");
    Outcome<Converged, SCFFailure> last = Outcome<Converged, SCFFailure>::Fail(SCFFailure{});
    for (size_t s=0; s<schedule.size(); ++s)
    {
        if (s>0)   // stage 0's graph is already standing (the ctor built it); later stages re-seed from it
            BuildStage("annealed schedule stage "+std::to_string(s), schedule[s].accelerator, std::move(itsImp->cd), schedule[s].params);
        last = Converge(schedule[s].params);
        EmitStageSummary(itsImp->setup->opts.label, s, schedule.size(), schedule[s].params.SmearingkT,
                         itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->scf->GetEnergy());
    }
    if (!itsImp->setup->rastersResolved && !itsImp->rastersSettled)
    {
        CalibrateRasters(schedule.back().params);   // C2 step 3: the pilot, after the FINAL stage
        return Outcome_();                          // the re-converged point's
    }
    return last;   // the FINAL stage's -- the earlier ones exist to feed it
}

namespace {
//! \brief THE PER-BAND +U SLOT WEIGHTS (QCHEM_DIAGNOSTICS=u_bands[=N], 2026-10-24): for the N bands either side of each spin
//! channel's frontier, each band's weight \f$w=\ell^\dagger P_{\rm slot}\ell\f$ on every U slot of every Hubbard manifold,
//! \f$\ell=\langle w_\mu|\psi\rangle\f$ through the +U term's OWN projectors (HubbardChannels) -- the comparison with QE's
//! projwfc.x on ortho-atomic projectors (doc/Records/HubbardUHistory.md, β-MnO₂ 2026-10-24: does our U'd band carry the
//! projector weight QE's does?).  Summed over the sites of one (element, l, slot): one column per such triple, so a band reads
//! "how much O-p-normal / Mn-d is it".  A manifold the term has no symmetry for is one column (its total).  Energies in eV
//! from the channel's LUMO (the lowest level less than half full).  Print only: no state, no return value.
void ReportBandSlotWeights(const std::string& label, const qchem::WaveFunction::cWaveFunction& wf,
                           const qchem::Hamiltonian::HubbardChannels& hub, const Structure& st, int nEach)
{
    const double eV=27.211386245988;
    std::vector<std::string> el; st.ForEachSite([&](int Z, const rvec3_t&, bool){ el.push_back(thePeriodicTable().GetSymbol(Z)); });
    const auto chans=hub.Channels();
    const char* lname="spdfgh";
    // Column per (element, l, slot label); remember which manifold slot feeds which column.
    std::vector<std::string> cols;
    std::vector<std::vector<std::pair<size_t,const rmat_t*>>> feeds(chans.size());   // per manifold: (column, slot projector | null = total)
    std::vector<qchem::Hamiltonian::HubbardSlots> slots(chans.size());
    auto col=[&](const std::string& k){ auto it=std::find(cols.begin(),cols.end(),k); if (it!=cols.end()) return size_t(it-cols.begin()); cols.push_back(k); return cols.size()-1; };
    for (size_t M=0;M<chans.size();M++)
    {
        slots[M]=hub.SlotProjectors(M);
        const std::string base=el.at(chans[M].site)+" "+std::string(1,lname[chans[M].l]);
        if (slots[M].slots.empty()) feeds[M].push_back({col(base),nullptr});
        else for (const auto& sl : slots[M].slots) feeds[M].push_back({col(base+" "+sl.label),&sl.P});
    }
    struct Band { double e, occ; std::vector<double> w; };
    std::map<int,std::vector<Band>> bySpin;            // Spin as int: None folds into one channel list
    std::map<int,double> cap;
    auto walk=[&]<class U>(const qchem::Orbitals::TOrbitals<U>& os, const qchem::Irrep& ir)
    {
        const auto* blk=dynamic_cast<const qchem::BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
        if (!blk) throw std::logic_error("ReportBandSlotWeights: an orbital block that is not an Orbital_DFT_IBS");
        std::vector<const qchem::Orbitals::TOrbital<U>*> all;
        for (const auto* o : os.template Iterate<qchem::Orbitals::TOrbital<U>>()) all.push_back(o);
        if (all.empty()) return;
        mat_t<U> C(blk->GetNumFunctions(), all.size());
        for (size_t i=0;i<all.size();i++) { const vec_t<U>& c=all[i]->GetCoeff(); for (size_t a=0;a<c.size();a++) C(a,i)=c[a]; }
        const auto amp=hub.ProjectorAmplitudes(*blk, C);
        const int s=int(ir.ms);
        cap[s] = ir.ms==qchem::Spin::None ? 2.0 : 1.0;
        for (size_t i=0;i<all.size();i++)
        {
            Band b{all[i]->GetEigenEnergy(), all[i]->GetOccupation(), std::vector<double>(cols.size(),0.0)};
            for (size_t M=0;M<amp.size();M++)
                for (const auto& [c,P] : feeds[M])
                {
                    double w=0.0;
                    const size_t m=amp[M].rows();
                    if (!P) for (size_t a=0;a<m;a++) w+=std::norm(amp[M](a,i));
                    else    for (size_t a=0;a<m;a++) for (size_t b2=0;b2<m;b2++) w+=std::real(blazem::conjs(amp[M](a,i))*(*P)(a,b2)*amp[M](b2,i));
                    b.w[c]+=w;
                }
            bySpin[s].push_back(std::move(b));
        }
    };
    for (const qchem::Irrep& ir : wf.GetQNs())
    {
        const qchem::Orbitals::Orbitals* os=wf.GetOrbitals(ir);
        if (const auto* c=dynamic_cast<const qchem::Orbitals::TOrbitals<dcmplx>*>(os)) walk(*c, ir);
        else if (const auto* r=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(os)) walk(*r, ir);
        else throw std::logic_error("ReportBandSlotWeights: an orbital set of unknown scalar");
    }
    const std::streamsize prec0=std::cout.precision();
    for (auto& [s,bands] : bySpin)
    {
        std::sort(bands.begin(), bands.end(), [](const Band& a, const Band& b){ return a.e<b.e; });
        size_t lumo=0; while (lumo<bands.size() && bands[lumo].occ>0.5*cap[s]) lumo++;
        if (lumo==bands.size()) continue;
        const double ref=bands[lumo].e;
        const char* sn = s==int(qchem::Spin::Up) ? "up" : s==int(qchem::Spin::Down) ? "down" : "both";
        std::cout<<"["<<label<<" u_bands] spin "<<sn<<": eV from the LUMO, occupation, weight on each +U slot (summed over sites)\n   e(eV)    occ";
        for (const auto& c : cols) std::cout<<" | "<<c;
        std::cout<<"\n";
        const size_t from = lumo>size_t(nEach) ? lumo-nEach : 0, to=std::min(bands.size(), lumo+size_t(nEach));
        for (size_t i=from;i<to;i++)
        {
            std::cout<<std::fixed<<std::setprecision(3)<<std::setw(8)<<(bands[i].e-ref)*eV<<std::setw(7)<<bands[i].occ;
            for (size_t c=0;c<cols.size();c++) std::cout<<" | "<<std::setw(std::max<int>(5,int(cols[c].size())))<<bands[i].w[c];
            std::cout<<(i+1==lumo ? "   <- HOMO" : i==lumo ? "   <- LUMO" : "")<<"\n";
        }
        std::cout<<std::defaultfloat<<std::setprecision(prec0)<<std::flush;
    }
}
} // namespace

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Converge(const SCFParams& params)
{
    assert(itsImp->scf);
    // THE RESTORED STAGE'S STEP (the Auto policy, MakeAccelerator): an Auto accelerator on a restored density pairs with a LINEAR
    // mixer, and a converged density wants a restart's step, not the ground state's (relax 1.0 diverged a restart, Si Γ 2026-09-29).
    SCFParams own=params;
    if (itsImp->stageRestored && itsImp->stageAccelAsked==SCFAccelerators::Type::Auto
        && itsImp->scf->ResolveMixing(params).Mixer==SCFParams::Mixing::Linear)
        own.StartingRelaxRo=std::min(own.StartingRelaxRo, 0.3);
    auto R=ConvergeStage(own);
    // THE PILOT'S START, CONDITIONING-AWARE (2026-11-01, beta-MnO2 on a dense-d basis, lambda_min(S) 6.75e-7).  An unresolved run starts on
    // the ball's own rasters (cheap; the pilot raises what the converged D shows it needs) -- but the orthonormal frame amplifies a
    // matrix-element error by up to 1/lambda_min, and there a cheap start's alias (x 8: e^-8 = 3e-4) wrecked the first diagonalization
    // (+4.19e4 Ha) before the pilot ever ran.  Measured: x 11 (alias 1.5e-5) and above gave the healthy first iterate.  So a pilot run
    // whose first iterate was REFUSED (UnphysicalStop) is rebuilt ONCE on rasters whose alias e^-x does not exceed lambda_min
    // (x = ln(1/lambda_min), never below what it had) and run again; the pilot then audits as always.  Only a refused run pays: every
    // healthy start is untouched.  A STATED raster (accuracy 0) is the user's and keeps its refusal.
    if (!itsImp->converged && !itsImp->scf->UnphysicalStop().empty() && !itsImp->setup->planeWave
        && !itsImp->setup->rastersResolved && !itsImp->rastersSettled && !itsImp->startRaisedForConditioning)
    {
        namespace L3 = BasisSet::Lattice;
        const double lmin=BasisSet::MinOverlapEigenvalue(*itsImp->bs);
        const BasisSet::Gaussian::RasterExponents had=itsImp->ladder->DeliveredRasterExponents();
        const double x=std::max({lmin>0.0 ? std::log(1.0/lmin) : 0.0, had.density, had.vlocLong, had.xcL0});
        BasisSet::Gaussian::GPWTolerances t=itsImp->setup->opts.tolerances;
        t.rasterX=std::max(t.rasterX,x); t.vlocRasterX=std::max(t.vlocRasterX,x); t.xcRasterX=std::max(t.xcRasterX,x);
        std::cout<<"["<<itsImp->setup->opts.label<<" raster pilot] the first iterate was REFUSED on the ball's own rasters (x density "
                 <<had.density<<"), and the overlap's smallest eigenvalue is "<<std::scientific<<std::setprecision(3)<<lmin<<std::defaultfloat
                 <<": restarting from the seed on rasters whose alias e^-x does not exceed it, x = "<<std::setprecision(4)<<x<<std::endl;
        itsImp->startRaisedForConditioning=true;
        auto sized=std::make_shared<SolidSetup>(*itsImp->setup);
        sized->opts.tolerances=t;
        const Lattice_3D lat=*itsImp->lat;
        const auto mol=itsImp->mol;
        const SCFAccelerators::Type asked=itsImp->stageAccelAsked;
        itsImp->spin.reset(); itsImp->scf.reset(); itsImp->accel.reset(); itsImp->cd.reset(); itsImp->ham.reset(); itsImp->ec.reset();
        itsImp->bs.reset();
        BuildPoint(std::move(sized), lat, mol);
        BuildRun(nullptr, own, asked, /*restored*/false, /*exactResume*/false);
        R=ConvergeStage(own);
    }
    if (itsImp->resumedFrom)   // an EXACT RESUME's first convergence: judged against the energy the state was saved at
    {
        const double E=itsImp->scf->GetEnergy().GetTotalEnergy();
        std::cout << "["<<itsImp->setup->opts.label<<" restart] resumed E="<<std::setprecision(10)<<E<<" vs saved "<<*itsImp->resumedFrom
                  << "  dE="<<std::setprecision(3)<<E-*itsImp->resumedFrom<<" after "<<itsImp->scf->GetIterationCount()<<" iterations"<<std::endl;
        itsImp->resumedFrom.reset();
    }
    return R;
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::ConvergeStage(const SCFParams& params)
{
    // ★ THE CONVERGE RESIDUE (doc/ParallelAndOraclePlan.md 1.1(a)).  The work AFTER Iterate returns is not
    // bookkeeping: GetChargeDensity() builds a fresh composite density and GetSpinDensity() rasters BOTH
    // spin channels, once per stage -- and an annealed run has a stage per schedule entry.
    qchem::report::Timed residue("scf: converge (residue -- final density + m(r) extraction)");
    // The banner and the recipe say what the iterator WILL build: SCFParams::Mixing::Auto resolved by the iterator itself.
    const SCFParams resolved=itsImp->scf->ResolveMixing(params);
    EmitSCFBanner(itsImp->setup->opts.label, resolved, itsImp->stageAccel, params.Mixer==SCFParams::Mixing::Auto,
                  itsImp->stageAccelAsked==SCFAccelerators::Type::Auto);
    { qchem::report::Section scf("scf"); ReportSCFRecipe(resolved, AccelName(itsImp->stageAccel)); }   // the run report's "what recipe am I running" header, as every facade files it
    itsImp->lastOccupation = {.useMOM=params.UseMOM, .momStartIter=(int)params.MOMStartIter,
                              .kT=params.SmearingkT, .momPenalty=params.MOMSmearPenalty};   // the iterator's own conversion
    itsImp->lastkT     = params.SmearingkT;
    itsImp->lastRecipe = std::string("accel=")+AccelName(itsImp->stageAccel)+" kT="+std::to_string(params.SmearingkT)
                       + " MOM="+(params.UseMOM?"on":"off")+" KerkerG0="+std::to_string(resolved.KerkerG0)
                       + " PulayDepth="+std::to_string(resolved.PulayDepth)+" NMaxIter="+std::to_string(params.NMaxIter);
    itsImp->samplingScfWorst=0.0; itsImp->samplingScfFirstOver=0;
    itsImp->scf->Iterate(params);
    if (itsImp->samplingScfWorst>1.0)   // never silent: reported here, acted on (once) at the next SCF boundary
    {
        std::cout<<"["<<itsImp->setup->opts.label<<" XC sampling] OVER TARGET in this SCF: worst delivered "<<itsImp->samplingScfWorst
                 <<" x tau, first at iteration "<<itsImp->samplingScfFirstOver
                 <<(itsImp->samplingTightened ? " (tau was already tightened once: not again)" : " (the next Hamiltonian build tightens tau once)")
                 <<std::endl;
        qchem::report::Scope here("xcSamplingOverTarget");
        qchem::report::Set("worstDeliveredOverTau", itsImp->samplingScfWorst);
        qchem::report::Set("firstIteration", itsImp->samplingScfFirstOver);
        qchem::report::Set("alreadyTightened", itsImp->samplingTightened);
        if (!itsImp->samplingTightened) itsImp->samplingPending=itsImp->samplingScfWorst;
    }
    itsImp->converged = itsImp->scf->Converged();
    itsImp->diag.itsConverged = itsImp->converged;
    // Take the density OUT of the wave function: it stays valid after the iterator's state moves on
    // because its basis block is `bs`, which this object owns and outlives it.
    auto cd = itsImp->scf->GetWaveFunction()->GetChargeDensity();   // BUILT for us; we take it
    itsImp->charge = cd->GetTotalCharge();
    itsImp->cd = std::move(cd);
    // m(r): only a run under the POLARIZED spin subgroup has one -- under imposed SU(2) it is identically
    // zero and is not built (V1.37: the subgroup is a property of the run, asked of the wave function, not a
    // type to cross-cast for).  RESET on the unpolarized branch: Converge runs once per anneal STAGE, so a
    // stale m(r) from an earlier stage must not survive into a run that no longer has one.
    const auto* wf = itsImp->scf->GetWaveFunction();
    itsImp->spin.reset();   // built on demand by Converged::SpinDensity()
    const bool polarizedRun = wf->GetSpinGroup()==SpinGroup::Polarized;
    ReportResult(itsImp->setup->opts.label, itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->charge, itsImp->scf->GetEnergy());
    {
        const qchem::EnergyBreakdown& E=itsImp->scf->GetEnergy();
        const std::streamsize prec0=std::cout.precision();
        // THE INTEGRATED SITE MOMENTS, a standard part of a polarized result (D-ENV step 3, user 2026-10-03: promoted from
        // the QCHEM_SITE_MOMENTS diagnostic).  An INTEGRATED observable -- the Becke-partitioned integral of w(rho_up -
        // rho_dn) per site -- never a point sample (doc/ feedback_integrated_observables).  A polarized run on an
        // atom-centred mesh that has NO site blocks is a DEFECT and is said so; on a uniform mesh the partition does not
        // exist, which is a configuration, said once.
        if (polarizedRun)
        {
            const rvec_t& mu=E.charge.siteMoments;
            double absSum=0.0, net=0.0; for (size_t a=0;a<mu.size();a++) { absSum+=std::fabs(mu[a]); net+=mu[a]; }
            if (mu.size()>0 && absSum>=1e-8)
            {
                std::cout<<"["<<itsImp->setup->opts.label<<"] site moments (Becke-partitioned Integral w (rho_up-rho_dn) d3r) [e]:";
                for (size_t a=0;a<mu.size();a++) std::cout<<"  "<<a<<":"<<std::showpos<<std::fixed<<std::setprecision(4)<<mu[a]<<std::noshowpos;
                std::cout<<"   net="<<std::showpos<<net<<std::noshowpos<<std::defaultfloat<<std::setprecision(prec0)<<std::endl;
            }
            else if (mu.size()==0)
                std::cout<<"["<<itsImp->setup->opts.label<<"] site moments: "
                         <<(itsImp->setup->xcMesh.cellKind==qcMesh::UnitCellKind::Becke
                            ? "UNAVAILABLE -- the atom-centred (Becke) XC mesh carries no site blocks: a DEFECT, not a configuration"
                            : "n/a (the uniform XC mesh has no site partition)")<<std::endl;
        }
    }
    if (auto v=qchem::Diagnostics::Value("u_bands"))
        if (const auto* hub=itsImp->ham->GetHubbardChannels())
        {
            const int n=std::atoi(v->c_str());
            ReportBandSlotWeights(itsImp->setup->opts.label, *itsImp->scf->GetWaveFunction(), *hub, *itsImp->st, n>1 ? n : 8);
        }
    // THE CHECKPOINT (CK-1), converged or not: the unattended run's safety net.  A failed WRITE warns and the run
    // goes on -- losing a finished SCF to a full disk would make the net the hazard (SaveState itself throws).
    if (!itsImp->setup->opts.saveStateTo.empty())
    {
        try { SaveState(itsImp->setup->opts.saveStateTo); }
        catch (const std::exception& e)
        {
            std::cout << "["<<itsImp->setup->opts.label<<"] WARNING: state NOT saved to "<<itsImp->setup->opts.saveStateTo
                      << ": "<<e.what()<<std::endl;
        }
    }
    return Outcome_();
}

//====================================================================================================
//  SAVED STATES (CK-1; the layout, the fingerprint classes and the reasons are in qchem.SolidState)
//====================================================================================================
StateFingerprint SolidCalculation::Fingerprint(double kT) const
{
    const auto& o=itsImp->setup->opts;
    std::ostringstream xc;
    xc<<MeshName(itsImp->setup->xcMesh.cellKind);
    if (itsImp->setup->xcMesh.cellKind==qcMesh::UnitCellKind::Becke) xc<<" nR="<<itsImp->setup->xcMesh.nRadial<<" L="<<itsImp->setup->xcMesh.angularDegree;
    else                                                      xc<<" eCut="<<itsImp->setup->xcMesh.eCut;
    // A free run on the PLAIN mesh (xcMeshSymmetric=false, what a SolidSession needs) quadratures on different points: name it, so a
    // restart across the choice says the quadrature changed.  Only then -- every state written before the option keeps its string.
    if (!itsImp->setup->imposed && !itsImp->setup->planeWave && !o.xcMeshSymmetric) xc<<" plain";
    RunIdentity id{.species=itsImp->setup->pp.Species(*itsImp->st), .Nelec=o.Nelec, .multiplicity=o.multiplicity, .functional=itsImp->setup->pp.Functional(), .ppSet=itsImp->setup->pp.Name(),
                         .hubbard=itsImp->hubbard, .densityEcut=o.densityEcut, .cutoffFactor=o.cutoffFactor,
                         .xcMesh=xc.str(), .ortho=int(o.ortho), .orthoTol=o.orthoTol, .kT=kT,
                         .orbitalBasis=itsImp->setup->planeWave ? "planeWave" : "gaussian", .planeWaveEcut=o.planeWaveEcut};
    // The numerics that move the ground state (OpenWork §3 "A restart's fingerprint is incomplete" gap 2).  The fit is the one the
    // Hamiltonian RESOLVES (Ham_PW_DFT: Auto -> Delta on a Becke mesh, else the plane-wave fit).
    const bool becke=itsImp->setup->xcMesh.cellKind==qcMesh::UnitCellKind::Becke;
    id.vxcFit= o.vxcFit==Hamiltonian::VxcFit::Delta || (o.vxcFit==Hamiltonian::VxcFit::Auto && becke) ? "Delta" : "PlaneWave";
    if (becke) id.beckeEps=itsImp->setup->xcMesh.beckeEps;
    if (!itsImp->setup->planeWave)
    {
        const auto& t=o.tolerances;
        {   // the raster demands (C2 step 3): what the run sampled on; every demand 0 is the ball's own rasters
            std::ostringstream r; r<<std::setprecision(6)<<"x="<<t.rasterX<<","<<t.vlocRasterX<<","<<t.xcRasterX;
            id.raster = t.rasterX==0.0 && t.vlocRasterX==0.0 && t.xcRasterX==0.0 ? std::string(BallRastersLabel) : r.str();
        }
        id.ladderFactor=o.ladderFactor;
        id.gpwTolerances={{"vlocEps",{t.vlocEps}}, {"localPPRelCutoff",{t.localPPRelCutoff}}, {"relFieldSharp",{t.relFieldSharp}},
                          {"screenEps",{t.screenEps}}, {"fieldSharp",{t.fieldSharp}}, {"relCutoff",{t.relCutoff}},
                          {"densityEps",{t.densityEps}}, {"mgridEcuts",t.mgridEcuts},
                          {"rasterX",{t.rasterX}}, {"vlocRasterX",{t.vlocRasterX}}, {"xcRasterX",{t.xcRasterX}}, {"accuracyEps",{o.accuracy.eps}}};
    }
    // The run's OWN policy table, resolved from its stated spec (not the process's installed one, which a later run may have replaced).
    for (const Deviation& d : RunPolicy(o.policy).Deviations())
        if (d.movesState) id.policy[d.knob]=d.value;
    return MakeStateFingerprint(itsImp->lat->GetUnitCell().GetCellMatrix(), *itsImp->st, *itsImp->bs,
                                o.multiplicity>=1 ? SpinGroup::Polarized : SpinGroup::UnPolarized, id);
}

void SolidCalculation::SaveState(const std::string& path) const
{
    const auto* wf=itsImp->scf ? itsImp->scf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("SolidCalculation::SaveState: no wave function yet");
    const StateSummary sum{.label=itsImp->setup->opts.label, .converged=itsImp->converged,
                           .iterations=itsImp->scf->GetIterationCount(),
                           .energy=itsImp->scf->GetEnergy().GetTotalEnergy(), .charge=itsImp->charge,
                           .commutator=itsImp->lastCommutator, .recipe=itsImp->lastRecipe,
                           .densityChange=itsImp->lastDensityChange};
    {
        qchem::report::Timed timed("state: save (HDF5)");
        WriteSolidState(path, Fingerprint(itsImp->lastkT), sum, *itsImp->bs, *wf);
    }
    std::cout << "["<<itsImp->setup->opts.label<<"] state saved to "<<path
              << (itsImp->converged ? "" : "  (NOT converged -- a restart will say so)") << std::endl;
}

Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>
SolidCalculation::Restart(const std::string& path, std::shared_ptr<const SolidSetup> setup, const SCFParams& params,
                          Observer onIteration)
{
    auto R=PrepareRestart(path, std::move(setup), params, std::move(onIteration));
    if (R)
    {
        (void)R.Value()->Converge(params);   // the restored stage, as PrepareRestart stood it up
        R.Value()->CalibrateRasters(params);  // C2 step 3: a restart on an unresolved setup is a pilot too
    }
    return R;
}

Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>
SolidCalculation::PrepareRestart(const std::string& path, std::shared_ptr<const SolidSetup> setup, const SCFParams& params,
                                 Observer onIteration)
{
    using O=Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>;
    const SCFAccelerators::Type accel=setup->opts.accelerator;
    const std::string tag="["+setup->opts.label+" restart] ";
    auto refuse=[&](RestartRefusal r)
    {
        std::cout << tag << "REFUSED ("<<path<<"): "<<r.details<<std::endl;
        return O::Fail(std::move(r));
    };
    auto saved=ReadSolidState(path);
    if (!saved) return refuse(saved.Error());

    // A STATE SIZED FOR THIS TARGET (C2 step 3): on an unresolved setup whose accuracy eps the saved state was made at, its rasters are the
    // ones a pilot settled for that eps -- adopt them, so an exact resume stays exact (the pilot would otherwise restart on the ball's
    // own rasters and re-size).  Any other state (another eps, a stated raster, an older file) leaves the setup to its own pilot.
    // The AUTO density ball is pilot-sized too (the ball ladder doubles cutoffFactor): adopting the rasters skips the pilot, so the
    // ball it settled must come with them -- else the resume runs the rasters on the starting ball and walks to that ball's fixed
    // point (MnO +U vendored twin, 2026-11-03: saved at cutoffFactor 4, resumed at 2, 3-4 iterations and 150 uHa off the cold gate).
    if (!setup->rastersResolved)
    {
        const auto& w=saved->fingerprint.warmNum;
        auto one=[&](const char* k)->std::optional<double>
        { auto it=w.find(k); return it!=w.end() && it->second.size()==1 ? std::optional<double>(it->second[0]) : std::nullopt; };
        const auto eps=one("gpw.accuracyEps"), rx=one("gpw.rasterX"), vx=one("gpw.vlocRasterX"), xx=one("gpw.xcRasterX");
        const auto cf=one("grids.cutoffFactor");
        if (eps && rx && vx && xx && *eps==setup->opts.accuracy.eps)
        {
            auto adopted=std::make_shared<SolidSetup>(*setup);
            adopted->opts.tolerances.rasterX=*rx; adopted->opts.tolerances.vlocRasterX=*vx; adopted->opts.tolerances.xcRasterX=*xx;
            const bool adoptBall = cf && !(setup->opts.densityEcut>0.0);  // a STATED densityEcut is the user's, never the pilot's
            if (adoptBall) adopted->opts.cutoffFactor=*cf;
            adopted->rastersResolved=true;
            std::cout<<tag<<"adopting the saved state's rasters x(density, V_long, XC L0)=("<<*rx<<", "<<*vx<<", "<<*xx<<")";
            if (adoptBall) std::cout<<" and its density ball (cutoffFactor "<<*cf<<")";
            std::cout<<": it was sized for this accuracy (eps "<<*eps<<")"<<std::endl;
            setup=std::move(adopted);
        }
    }
    std::unique_ptr<SolidCalculation> c(new SolidCalculation(BuildOnly{}));
    c->itsImp->onIteration=std::move(onIteration);
    // A restored density is at (or next to) the answer: a Ladder starts on its GDM rung (user 2026-10-21, OpenWork
    // §3 "Ladder restart starts on the DIIS+Kerker rung").  For THIS stage only -- reset after Converge below.
    // kT = 0 ONLY: a smeared stage restored from an integer state hands GDM an idempotent D', so GDM would ENGAGE and
    // minimise with integer occupations, never applying the smearing (measured 2026-10-21, Si k222 kT 0.02: 6 mHa off the
    // cold smeared run).  A smeared restart starts on DIIS, as it did before.
    const Lattice_3D& reference=setup->reference;
    auto mol=setup->mol;
    c->BuildPoint(std::move(setup), reference, std::move(mol));
    c->itsImp->accOpts.ladderStartOnGDM = params.SmearingkT==0.0;
    // Judged against the BUILT blocks (k, weight, size), not against a re-derivation of what they would be.
    auto warm=CompareFingerprints(saved->fingerprint, c->Fingerprint(params.SmearingkT));
    if (!warm) return refuse(warm.Error());
    auto cd=RestoreDensity(*saved, *c->itsImp->bs, c->itsImp->setup->polarized ? SpinGroup::Polarized : SpinGroup::UnPolarized);
    if (!cd) return refuse(cd.Error());

    const bool exact=warm->empty();
    {
        std::ostringstream os;
        os<<tag<<(exact ? "EXACT RESUME" : "WARM START")<<" from "<<path<<" (saved by '"<<saved->summary.label<<"', "
          <<(saved->summary.converged ? "converged" : "NOT converged")<<" after "<<saved->summary.iterations
          <<" iterations, E="<<std::setprecision(10)<<saved->summary.energy<<"; "<<saved->summary.recipe<<")";
        for (const std::string& d : *warm) os<<"\n"<<tag<<"  differs: "<<d;
        // THE SAVED STATE'S OWN RESIDUAL (gap 4, user 2026-10-29): even an exact resume starts where the state was SAVED, not at
        // the fixed point.  A density-MIXED run (Null accelerator: no [F,D] measured) stops with D a mixture, one undamped Fock
        // away from self-consistency, so its first iteration re-opens that residual (amplified on a soft system).
        const StateSummary& ss=saved->summary;
        if (std::isfinite(ss.densityChange))
        {
            os<<"\n"<<tag<<"  saved at Δρ "<<std::setprecision(3)<<ss.densityChange;
            if (std::isfinite(ss.commutator)) os<<", [F,D] "<<ss.commutator;
            else os<<" by a density-MIXED run (no [F,D]): D is a mixture, so the first Fock re-opens a residual of about that size -- "
                     "amplified on a soft system (β-MnO₂ r015: 40x), where even an exact resume iterates";
        }
        else os<<"\n"<<tag<<"  saved residual not recorded (a format-1.0 state, before 2026-10-30)";
        std::cout<<os.str()<<std::endl;
    }
    c->BuildRun(cd.TakeValue(), params, accel, /*restored*/true, exact);
    // The restored stage's accelerator is BUILT (BuildRun consumed the GDM-rung start); later stages on this calculation start as asked.
    c->itsImp->accOpts.ladderStartOnGDM=c->itsImp->setup->acc.ladderStartOnGDM;
    if (exact) c->itsImp->resumedFrom=saved->summary.energy;   // the first Converge reports the resumed energy beside it
    return O::Ok(std::move(c));
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Result() const {return Outcome_();}

qchem::EnergyBreakdown SolidCalculation::LastIterateTerms()  const {return itsImp->scf->GetEnergy();}
double                 SolidCalculation::LastIterateCharge() const {return itsImp->charge;}
// WHY A COMPLETE POOL (found 2026-09-21, when the projector was still contracted onto the site's shells -- the route
// HubbardProjectorPlan increment 5 retired).  A solid's valence basis is trimmed for the SOLID: the MnO VA span dropped
// the diffuse 0.18 d shell, and the free 4s2 3d5 atom run in those shells has eps(3d) 77 mHa too shallow and a too-compact
// "3d".  In a 16+16 even-tempered pool the same pseudo-atom reproduces CP2K's ATOM code to 3 / 0.3 / 0.1 mHa (E, eps 4s,
// eps 3d; UTCalculation ValenceBasisGen.MnQ7PseudoAtomInALargePool) -- so the projector is the pool's orbital, as is.
// THE PSEUDO-ATOMIC +U PROJECTOR (HubbardProjectorPlan D5/D6): the site species' pseudo-atom l orbital in its OWN
// pool, handed over as is -- no projection onto the site's shells, no knowledge of the orbital basis.
static std::vector<Math::Gaussian> PseudoAtomicFunction(const Structure& st, size_t site, int l,
                                                        const Pseudopotential::PseudopotentialSet& pp)
{
    int Z=0; { size_t i=0; st.ForEachSite([&](int z, const rvec3_t&, bool){ if (i==site) Z=z; i++; }); }
    const std::string element=thePeriodicTable().GetSymbol(Z);
    const int Zion=pp.Valence(element);
    const PseudoAtomOrbital& pa=[&]() -> const PseudoAtomOrbital&
    {
        qchem::SetupSCF setup;   // the pseudo-atom solve is setup, not the run's convergence (+122 on every beta-MnO2 panel case)
        return PseudoAtomOrbitalOf(element, Zion, l);
    }();
    std::ostringstream os;
    os<<"[+U projector] site "<<site<<" ("<<element<<", q"<<Zion<<") l="<<l<<": pseudo-atom "<<(l==0?"s":l==1?"p":l==2?"d":"f")
      <<" orbital eps="<<std::setprecision(4)<<pa.eps<<" Ha (16-exponent pool, E="<<std::setprecision(6)<<pa.Eatom
      <<"), in its own exponents";
    std::cout<<os.str()<<std::endl;
    return pa.Radial();
}

namespace {
//! \brief +U FROZEN for a scope (doc/LinearResponsePlan.md §3d Q6), and the previous state restored on EVERY exit
//! path -- the response kernel is applied lazily during the solve, so freezing only while it is BUILT would not do.
class FrozenHubbard
{
public:
    explicit FrozenHubbard(qchem::Hamiltonian::HubbardUTarget* t) : itsT(t), itsWas(t && t->OccupationsFrozen())
    {if (itsT) itsT->FreezeOccupations(true);}
    ~FrozenHubbard() {if (itsT) itsT->FreezeOccupations(itsWas);}
    FrozenHubbard(const FrozenHubbard&)=delete;
    FrozenHubbard& operator=(const FrozenHubbard&)=delete;
private:
    qchem::Hamiltonian::HubbardUTarget* itsT;
    bool                                itsWas;
};

//! The default perturbed set: the manifolds that CARRY a +U (hp.x's Hubbard sites), else all of them (§3d Q10).
std::vector<size_t> DefaultPerturbed(const qchem::Hamiltonian::HubbardChannels& hub)
{
    const auto ch=hub.Channels();
    std::vector<size_t> J;
    for (size_t i=0;i<ch.size();i++) if (ch[i].carriesU) J.push_back(i);
    if (J.empty()) for (size_t i=0;i<ch.size();i++) J.push_back(i);
    return J;
}
} // namespace

//---------------------------------------------------------------------------------------------------
//  THE GRID AUDIT (doc/GridEpsilonPlan.md §3)
//---------------------------------------------------------------------------------------------------
namespace
{
//! The raster demands of \a t, as the audit names a side.
std::string RasterName(const BasisSet::Gaussian::GPWTolerances& t)
{
    if (t.rasterX==0.0 && t.vlocRasterX==0.0 && t.xcRasterX==0.0) return "ball rasters";
    std::ostringstream os; os<<std::setprecision(4)<<"x(density, V_long, XC L0)=("<<t.rasterX<<", "<<t.vlocRasterX<<", "<<t.xcRasterX<<")";
    return os.str();
}
//! A Becke XC mesh as the pilot names a rung.
std::string BeckeName(const qcMesh::MeshParams& m)
{
    return "Becke (nR="+std::to_string(m.nRadial)+", L="+std::to_string(m.angularDegree)
           +(m.angular==qcMesh::AngularKind::Lebedev ? " Lebedev)" : " GL)");
}
//! One function of a block as a human reads it: its atom, its shell and l.  The shells' order is the block's function order (\c AoShellSource).
std::string FunctionName(const std::vector<Symmetry::Molecule::AoShell>& shells, const Structure& st, size_t f)
{
    for (size_t s=0;s<shells.size();s++)
    {
        const auto& sh=shells[s];
        if (f<sh.offset || f>=sh.offset+sh.nComponents()) continue;
        size_t atom=0, a=0; int z=0; double best=std::numeric_limits<double>::max();
        for (const Site& site : st.Sites()) { const double d=norm(site.R-sh.center); if (d<best) {best=d; atom=a; z=site.Z;} a++; }
        std::ostringstream os; os<<"atom "<<atom<<" Z="<<z<<" shell "<<s<<" l="<<sh.rep->L()<<" component "<<f-sh.offset;
        if (sh.exponents.size()>0)
        {
            double lo=sh.exponents[0], hi=lo;
            for (size_t p=1;p<sh.exponents.size();p++) { lo=std::min(lo,sh.exponents[p]); hi=std::max(hi,sh.exponents[p]); }
            os<<" alpha "<<lo<<".."<<hi;
        }
        return os.str();
    }
    throw std::logic_error("GridAudit: function "+std::to_string(f)+" is in no shell");
}
//! Two term lists of one (spin, block) channel, folded into \a terms (by position: one term order on both Hamiltonians): the largest
//! |a_ij - b_ij| and where it sits, and beside it the D-weighted Re Tr(D ΔH) and the S-metric norm.  \a frame = the run's own orthonormal
//! frame V for this block: an explicit S⁻¹ loses digits at small λ_min, and would count directions the run dropped.
template <class M> void FoldΔH(const std::vector<M>& a, const std::vector<M>& b, const SavedBlock& sb, const LASolver<dcmplx>& frame,
                               const std::string& spin, size_t block, std::vector<SolidCalculation::GridAuditTerm>& terms)
{
    if (a.size()!=terms.size() || b.size()!=terms.size()) throw std::logic_error("GridAudit: a block's term list differs from the energies'");
    for (size_t t=0;t<a.size();t++)
    {
        if (a[t].rows()!=b[t].rows() || a[t].columns()!=b[t].columns())
            throw std::logic_error("GridAudit: term "+std::to_string(t)+"'s block is "+std::to_string(a[t].rows())+" on the run, "
                                   +std::to_string(b[t].rows())+" on the reference -- the bases do not match function for function");
        const size_t n=a[t].rows();
        if (sb.D.rows()!=n) throw std::logic_error("GridAudit: block "+std::to_string(block)+"'s D does not match its term matrices");
        mat_t<dcmplx> ΔH(n,n);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) ΔH(i,j)=dcmplx(a[t](i,j))-dcmplx(b[t](i,j));
        auto& T=terms[t];
        double tr=0.0, largest=0.0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) tr+=std::real(sb.D(i,j)*ΔH(j,i));
        T.trDdH.push_back(sb.weight*tr);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
        {
            const double d=std::abs(ΔH(i,j));
            largest=std::max(largest,d);
            if (d>T.maxAbsΔH) { T.maxAbsΔH=d; T.maxElement.spin=spin; T.maxElement.block=block; T.maxElement.i=i; T.maxElement.j=j; T.maxElement.absD=std::abs(sb.D(i,j)); }
        }
        if (largest==0.0) continue;   // an analytic term: identical on both sides
        hmat_t<dcmplx> H(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) H(i,j)=ΔH(i,j);
        const hmat_t<dcmplx> X=frame.Transform(H);
        double f2=0.0;
        for (size_t i=0;i<X.rows();i++) for (size_t j=0;j<X.columns();j++) f2+=std::norm(X(i,j));
        T.sNormDH=std::max(T.sNormDH, std::sqrt(f2));
    }
}
} //anon

BasisSet::Gaussian::GPWTolerances SolidCalculation::ReferenceRasters(double margin) const
{
    namespace L3 = BasisSet::Lattice;
    if (!itsImp->ladder) throw std::logic_error("SolidCalculation: the basis has no GPW raster to refine");
    const std::optional<BasisSet::Gaussian::RasterExponents> x=itsImp->ladder->DeliveredRasterExponents();
    BasisSet::Gaussian::GPWTolerances t=itsImp->setup->opts.tolerances;
    t.rasterX    =x->density+margin;
    t.vlocRasterX=(x->vlocLong>0.0 ? x->vlocLong : x->density)+margin;   // no own V_long ladder (the kappa-sweep path): it rides the density's
    // The XC raster's error falls as a POWER of x, not e^{-x} (Si, C2 step 3), so x + margin would leave the reference's own error a
    // sizeable fraction of the run's at large x: at least twice the run's x there ((1/2)^k of it, k ~ 4 measured).
    t.xcRasterX  =std::max(x->xcL0+margin, 2.0*x->xcL0);
    return t;
}

SolidCalculation::GridAuditResult SolidCalculation::GridAudit(double referenceMargin) const
{
    const SolidSetup& S=*itsImp->setup;
    if (S.planeWave) throw std::invalid_argument("SolidCalculation::GridAudit: a plane-wave run has no raster to audit");
    if (!(referenceMargin>0.0)) throw std::invalid_argument("SolidCalculation::GridAudit: the reference must be FINER than the run (margin > 0)");
    qchem::report::Timed timed("postSCF: grid audit (reference basis + Hamiltonian + term matrices)");
    const BasisSet::Gaussian::GPWTolerances ref=ReferenceRasters(referenceMargin);
    return CompareAtFixedD(ref, rvec3_t(0,0,0), RasterName(ref));
}

ClientConsumes SolidCalculation::Consumes() const {return itsImp->consumes;}
SolidCalculation::RasterPilotCost  SolidCalculation::PilotCost() const {return itsImp->pilotCost;}

std::shared_ptr<const SolidSetup> SolidCalculation::SizedSetup() const
{
    if (itsImp->setup->rastersResolved) return itsImp->setup;
    if (!itsImp->rastersSettled) throw std::logic_error("SolidCalculation::SizedSetup: the pilot audit has not settled the rasters (the run did not converge?)");
    auto s=std::make_shared<SolidSetup>(*itsImp->setup);
    s->rastersResolved=true;
    return s;
}

void SolidCalculation::CalibrateRasters(const SCFParams& params)
{
    namespace G   = BasisSet::Gaussian;
    namespace L3  = BasisSet::Lattice;
    namespace rpt = qchem::report;
    if (itsImp->setup->rastersResolved || itsImp->rastersSettled) return;
    const std::string tag="["+itsImp->setup->opts.label+" raster pilot] ";
    if (!itsImp->converged)
    {
        std::cout<<tag<<"SKIPPED: the SCF did not converge, so there is no D to audit -- the rasters stay the ball's own"<<std::endl;
        return;
    }
    // THE TARGET, by what the client consumes (user ruling 2026-10-31): an energy run holds each raster term's |dE| to eps x the atom
    // count; a geometry run holds its max|dH_ij| to eps (the energy error cancels inside Tr D dH, M2c verdict 3).
    const G::ClientConsumes c=itsImp->consumes;
    const bool   onH   =c==G::ClientConsumes::MatrixElements;
    const size_t nAtoms=itsImp->lat->GetUnitCell().GetNumAtoms();
    const double target=itsImp->setup->opts.accuracy.Target(c, nAtoms);
    qchem::report::Timed timed("postSCF: raster pilot (reference basis + audit; re-converge when a raster rises)");
    const std::streamsize prec0=std::cout.precision();
    std::cout<<tag<<"sizing the rasters to "<<std::setprecision(3)<<target<<" Ha "
             <<(onH ? std::string("per matrix element (max|ΔH_ij|): a geometry run")
                    : "per raster term's energy (eps x "+std::to_string(nAtoms)+" atoms): an energy run")<<std::endl;
    rpt::json rounds=rpt::json::array();
    std::optional<std::pair<double,double>> xcBefore;   // the XC term's (x, delivered) one round back: the power-law fit's first point
    bool xcOver=false;                                  // the XC term was over its target in the last round audited
    std::optional<std::pair<std::string,double>> beckeBefore;   // (rung, its difference) one round back: the ladder's convergence ratio
    bool xcStatedOver=false;                            // a STATED Becke mesh was over its target (honoured, reported)
    std::optional<std::pair<std::string,double>> ballBefore;    // (ball, its difference) one round back
    bool ballStatedOver=false;                          // a STATED densityEcut was over its target (honoured, reported)
    bool metEveryTarget=false;                          // the SETTLED line, not the round cap, ended the ladder
    for (int round=0;; ++round)
    {
        const G::RasterExponents x=itsImp->ladder->DeliveredRasterExponents();
        const double xV = x.vlocLong>0.0 ? x.vlocLong : x.density;   // no own V_long ladder (the kappa-sweep path): the density's
        const GridAuditResult R=CompareAtFixedD(ReferenceRasters(kReferenceMargin), rvec3_t(0,0,0), "pilot reference (x + "+std::to_string(int(kReferenceMargin))+")");
        G::GPWTolerances t=itsImp->setup->opts.tolerances;
        double dDens=0.0, dV=0.0, dXC=0.0;
        rpt::json terms=rpt::json::array();
        std::cout<<std::scientific<<std::setprecision(3);
        for (const auto& T : R.terms)
        {
            double* demand=nullptr; double xRun=0.0; const char* which=nullptr;
            if      (T.label.find("LONG-range")!=std::string::npos) { demand = x.vlocLong>0.0 ? &dV : &dDens; xRun=xV; which="V_long"; }
            else if (T.entries=="Eee")                              { demand=&dDens; xRun=x.density; which="Hartree"; }
            else if (T.entries=="Exc" && T.maxAbsΔH>0.0)
            {
                // WHOSE error is the uniform XC's?  Below the density's rigorous x the collocated rho itself aliases, and v_xc[rho] carries
                // it: Na k222s at density x 6 read |dE_xc| 1.2e-4 with its XC raster at x 28 AND at x 60, and 5e-8 once the density raster
                // went to x 16 (C2 step 3).  So there the XC term raises the DENSITY raster, under the density's shape; above it only the
                // nonlinearity of v_xc is left, and the XC raster's own power law sizes L0 (Si Gamma, density x 27).
                if (x.density < itsImp->setup->opts.accuracy.Rigorous(G::kDensityTail)) { demand=&dDens; xRun=x.density; which="XC via rho (MEASURED)"; }
                else                                                                    { demand=&dXC;   xRun=x.xcL0;    which="XC (MEASURED: no bound)"; }
            }
            if (!demand) continue;   // an analytic term: nothing on a raster
            const double d   =onH ? T.maxAbsΔH : std::fabs(T.ΔE);
            double need=G::GridAccuracy::MeasuredDemand(d, xRun, target);   // 0 = meets the target here
            if (demand==&dXC)   // no bound and no shape: once two rounds are measured, the XC term's own power law (MeasuredPowerDemand)
            {
                if (xcBefore && need>0.0)
                    if (const double p=G::GridAccuracy::MeasuredPowerDemand(xcBefore->first, xcBefore->second, xRun, d, target); p>0.0) need=p;
                xcBefore=std::pair{xRun, d};
                // At most DOUBLE x per round: a slope from two nearby, noisy rounds extrapolates wildly (Si: k 1.4 from x 21 -> 30 would
                // ask x ~ 3000), and the next round measures again anyway.
                if (need>0.0) { need=std::min(need, 2.0*xRun); xcOver=true; } else xcOver=false;
            }
            const double C   =d/(std::pow(xRun,3.0)*std::exp(-xRun));             // the measured scale under the rigorous shape
            *demand=std::max(*demand, need);
            std::cout<<tag<<"round "<<round<<"  "<<std::setw(24)<<std::left<<which<<std::right<<" delivered "<<d<<" at x "<<std::defaultfloat
                     <<std::setprecision(4)<<xRun<<std::scientific<<std::setprecision(3)<<"  (C = "<<C<<" Ha)  ";
            if (need>0.0) std::cout<<"OVER: demands x "<<std::defaultfloat<<std::setprecision(4)<<need<<std::scientific<<std::setprecision(3)<<std::endl;
            else          std::cout<<"meets "<<target<<std::endl;
            terms.push_back({{"term",which},{"delivered",d},{"x",xRun},{"C",C},{"demand",need}});
        }
        // THE BECKE XC MESH on a run whose client consumes H (no raster: its own nested ladder).  BeckeXCParams() is calibrated on the
        // ENERGY; a displacement reads H (lane B, NaF/LiF 2026-11-01: (20, 2, 24) 4-10x wrong on a 0.1 bohr curvature), so the mesh is held
        // to max|dH_xc| <= target against the next rung at this D.  An atom-centred mesh moves with the atoms: the egg-box cannot see it.
        std::optional<qcMesh::MeshParams> xcMeshNext;
        const qcMesh::MeshParams& xcCur=itsImp->setup->xcMesh;
        if (onH && xcCur.cellKind==qcMesh::UnitCellKind::Becke)
        {
            const qcMesh::MeshParams ref=qcMesh::BeckeLadderAbove(xcCur);
            const GridAuditResult B=XCMeshAudit(ref);
            double d=0.0;
            for (const auto& T : B.terms) if (T.entries=="Exc") d=std::max(d, T.maxAbsΔH);
            // The rung difference UNDER-states the rung's own error when the ladder converges slowly: Si Gamma halves per rung (1.8e-5,
            // 8.9e-6), so its error is ~ the whole geometric tail d/(1 - r), r the measured ratio to the rung below -- twice d there.  One
            // rung measured: d alone (no ratio yet); a ratio >= 1 (not converging) also falls back to d, and then the next rung measures again.
            // The ratio is between two RUNGS: a round that re-audits the SAME rung (another ladder moved, this one did not) has no ratio --
            // d/d ~ 1 made the estimate blow up (NaF geometry, 2026-11-02: 4.7e-7 re-measured on the same mesh read as a 5.2e-2 tail).
            double est=d;
            if (beckeBefore && beckeBefore->first!=BeckeName(xcCur) && d<beckeBefore->second) est=d/(1.0-d/beckeBefore->second);
            beckeBefore=std::pair{BeckeName(xcCur), d};
            const bool stated=itsImp->setup->opts.xcMesh.cellKind!=qcMesh::UnitCellKind::Auto;
            std::cout<<tag<<"round "<<round<<"  "<<std::setw(24)<<std::left<<"XC on the Becke mesh"<<std::right<<" delivered "<<std::scientific
                     <<std::setprecision(3)<<d<<" (tail estimate "<<est<<") on "<<BeckeName(xcCur)<<" against "<<BeckeName(ref)<<"  ";
            if (est<=target) std::cout<<"meets "<<target<<std::endl;
            else if (stated) { std::cout<<"OVER, and the mesh is STATED: honoured (state xcMesh Auto to let the pilot size it)"<<std::endl; xcStatedOver=true; }
            else             { std::cout<<"OVER: climbs to "<<BeckeName(ref)<<std::endl; xcMeshNext=ref; }
            std::cout<<std::defaultfloat<<std::setprecision(prec0);
            terms.push_back({{"term","XC Becke mesh"},{"delivered",d},{"tailEstimate",est},{"mesh",BeckeName(xcCur)},{"reference",BeckeName(ref)},
                             {"over",est>target}});
        }
        // THE DENSITY BALL (2026-11-02, OpenWork F "density ball unjudged"): the G-ball the density is collocated into is NOT a raster --
        // the x + 8 reference above keeps the run's ball, so its truncation sat on BOTH sides and cancelled (lane B: H2O MOLOPT O, auto
        // densityEcut 24 Ha, 330 uHa from CP2K; 100 Ha closes it to 2.5 uHa).  A nested ladder like the Becke mesh's: the run against a
        // reference with TWICE the ball at the same rasters, at this D, judged on Hartree (below); an auto ball over the target doubles
        // (cutoffFactor), a STATED densityEcut is honoured and reported.  The target is the run's eps (pin 31), so the ladder's START
        // (cutoffFactor x alpha_max) is only where it begins: a better start -- a basis Resolution(eps) query -- changes no answer.
        std::optional<double> ballNext;
        std::string ballOtherOver;   // a term the ball audit finds OVER that a bigger ball does not cure (it blocks SETTLED, honestly)
        {
            const SolidCalcOptions& o=itsImp->setup->opts;
            const GridAuditResult Bl=CompareAtFixedD(itsImp->setup->opts.tolerances, rvec3_t(0,0,0), "density ball x 2", nullptr, 2.0);
            // THE BALL'S OWN TERM is Hartree (the ball is its quadrature): it alone decides whether the ball doubles.  EVERY OTHER TERM is judged
            // too (controller ruling 2026-11-02, the honest gate): a bigger ball also re-routes the local-PP LONG sweep over the multigrid levels,
            // and that row does NOT fall as the ball doubles -- it alternates (MnO geometry: +5.7e-5, -6.3e-5, +6.1e-5 at cutoffFactor 2, 4, 8,
            // against MnO's 4e-5; Si Gamma ~1e-5) -- so it is reported OVER and blocks SETTLED instead of driving a doubling it cannot answer.
            double d=0.0, dOther=0.0;
            std::string otherTerm;
            for (const auto& T : Bl.terms)
            {
                const double v = onH ? T.maxAbsΔH : std::fabs(T.ΔE);
                if (T.entries=="Eee") d=std::max(d, v);
                else if (v>dOther) { dOther=v; otherTerm=T.entries+": "+T.label.substr(0, T.label.find('('));}
            }
            const std::string worst="Eee";
            const bool   stated=o.densityEcut>0.0;
            const std::string ball = stated ? "densityEcut "+std::to_string(o.densityEcut)+" Ha" : "cutoffFactor "+std::to_string(o.cutoffFactor);
            double est=d;   // the geometric tail over the doubling rungs, as for the Becke mesh (and, as there, only between two DIFFERENT balls)
            if (ballBefore && ballBefore->first!=ball && d<ballBefore->second) est=d/(1.0-d/ballBefore->second);
            ballBefore=std::pair{ball, d};
            std::cout<<tag<<"round "<<round<<"  "<<std::setw(24)<<std::left<<"density ball"<<std::right<<" delivered "<<std::scientific
                     <<std::setprecision(3)<<d<<" ["<<worst<<"] (tail estimate "<<est<<") on "<<ball<<" against twice it  ";
            if (est<=target) std::cout<<"meets "<<target<<std::endl;
            else if (stated) { std::cout<<"OVER, and densityEcut is STATED: honoured (leave it auto to let the pilot size the ball)"<<std::endl; ballStatedOver=true; }
            else             { std::cout<<"OVER: the ball doubles (cutoffFactor "<<2.0*o.cutoffFactor<<")"<<std::endl; ballNext=2.0*o.cutoffFactor; }
            if (dOther>target)
            {
                ballOtherOver=otherTerm;
                std::cout<<tag<<"round "<<round<<"  "<<std::setw(24)<<std::left<<"  beside the ball"<<std::right<<" delivered "<<dOther<<" ["
                         <<otherTerm<<"] against twice the ball  OVER: not the ball's own term, a bigger ball does not cure it (OpenWork §3)"<<std::endl;
            }
            std::cout<<std::defaultfloat<<std::setprecision(prec0);
            terms.push_back({{"term","density ball"},{"delivered",d},{"worstTerm",worst},{"tailEstimate",est},{"ball",ball},{"over",est>target},
                             {"besideTheBall",dOther},{"besideTerm",otherTerm},{"besideOver",dOther>target}});
        }
        std::cout<<std::defaultfloat<<std::setprecision(prec0);
        rounds.push_back({{"round",round},{"rasters",RasterName(t)},{"terms",terms}});
        itsImp->pilotCost.rounds=round+1;
        const G::GPWTolerances was=t;
        t.rasterX    =std::max(t.rasterX,     dDens);
        t.vlocRasterX=std::max(t.vlocRasterX, dV);
        t.xcRasterX  =std::max(t.xcRasterX,   dXC);
        if (t==was && !xcMeshNext && !ballNext && !ballOtherOver.empty())
        {   // nothing the pilot sizes would change: another round would re-measure the same grids.  Stop and say so.
            std::cout<<tag<<"NOT SETTLED after "<<round<<" re-convergence(s): every grid the pilot sizes meets "<<target<<" Ha on "<<RasterName(t)
                     <<", but "<<ballOtherOver<<" is over it beside the ball (no knob of the pilot's answers it) -- the run keeps these grids and says so"
                     <<std::endl;
            break;
        }
        if (t==was && !xcMeshNext && !ballNext)
        {
            std::cout<<tag<<"SETTLED after "<<round<<" re-convergence(s): every raster term meets "<<target<<" Ha on "<<RasterName(t);
            if (onH && xcCur.cellKind==qcMesh::UnitCellKind::Becke) std::cout<<", the XC on "<<BeckeName(xcCur)<<(xcStatedOver ? " (STATED, OVER)" : "");
            std::cout<<", the density ball "<<(ballStatedOver ? "(STATED, OVER)" : "meets it");
            std::cout<<std::endl;
            metEveryTarget=true;
            break;
        }
        if (round+1>=kPilotRounds)
        {
            std::cout<<tag<<"NOT SETTLED after "<<kPilotRounds<<" rounds: a term still over its target on "<<RasterName(was)
                     <<(xcMeshNext ? ", the XC on "+BeckeName(xcCur) : std::string())<<" (the XC raster has no bound: its measured tail can fall slower than the shape the pilot fits) -- the run keeps "
                     <<"these rasters and says so"<<std::endl;
            if (xcOver && onH)
                std::cout<<tag<<"ADVICE: the uniform XC raster is the term over -- v_xc[rho] sampled on a fixed raster does not follow the atoms, "
                         <<"so a geometry run on it carries this as an egg-box; the atom-centred Becke mesh (xcMesh.cellKind Becke) follows "
                         <<"them"<<std::endl;
            break;
        }
        // REBUILD THIS POINT on the raised rasters and re-converge from the D just audited: a restored stage on the same Hamiltonian
        // terms (the Ladder enters on its GDM rung at kT = 0, as a restart does).  Torn down iterator-first: each borrows the next.
        std::cout<<tag<<"raising the rasters to "<<RasterName(t);
        if (xcMeshNext) std::cout<<" and the XC mesh to "<<BeckeName(*xcMeshNext);
        if (ballNext)   std::cout<<" and the density ball to cutoffFactor "<<*ballNext;
        std::cout<<" and re-converging from this D"<<std::endl;
        const SavedState st=Carry();
        auto sized=std::make_shared<SolidSetup>(*itsImp->setup);
        sized->opts.tolerances=t;
        if (xcMeshNext) sized->xcMesh=*xcMeshNext;
        if (ballNext)   sized->opts.cutoffFactor=*ballNext;
        const Lattice_3D lat=*itsImp->lat;
        const auto mol=itsImp->mol;
        // A restart's accelerator (CLAUDE.md: restarts run GDM or the DIIS->GDM Ladder at kT = 0): an Auto stage would resolve a restored
        // density to LINEAR mixing at relax 0.3 -- Si's energy pilot spent 39 iterations there.  A stated accelerator is obeyed.
        SCFAccelerators::Type asked=itsImp->stageAccelAsked;
        if (asked==SCFAccelerators::Type::Auto && params.SmearingkT==0.0) asked=SCFAccelerators::Type::Ladder;
        itsImp->spin.reset(); itsImp->scf.reset(); itsImp->accel.reset(); itsImp->cd.reset(); itsImp->ham.reset(); itsImp->ec.reset();
        itsImp->bs.reset();
        BuildPoint(std::move(sized), lat, mol);
        const SpinGroup g=itsImp->setup->polarized ? SpinGroup::Polarized : SpinGroup::UnPolarized;
        auto cd=RestoreDensity(st, *itsImp->bs, g);
        if (!cd) throw std::logic_error("SolidCalculation::CalibrateRasters: the audited D does not restore onto the raised rasters: "+cd.Error().details);
        itsImp->accOpts.ladderStartOnGDM = params.SmearingkT==0.0;
        BuildRun(cd.TakeValue(), params, asked, /*restored*/true, /*exactResume*/false);
        itsImp->accOpts.ladderStartOnGDM = itsImp->setup->acc.ladderStartOnGDM;
        (void)Converge(params);
        itsImp->pilotCost.reconvergeIterations+=itsImp->scf->GetIterationCount();
        if (!itsImp->converged)
        {
            std::cout<<tag<<"the re-convergence on the raised rasters FAILED: the run's result is that failure"<<std::endl;
            break;
        }
    }
    itsImp->rastersSettled=itsImp->converged;
    itsImp->pilotCost.settled=itsImp->rastersSettled && metEveryTarget;   // NOT SETTLED at the round cap is not settled (it read true)
    std::cout<<tag<<"cost: "<<itsImp->pilotCost.rounds<<" audit round(s), "<<itsImp->pilotCost.reconvergeIterations
             <<" re-convergence iteration(s) beyond the run's own"<<std::endl;
    rpt::EmitSection("rasterPilot", rpt::json{{"consumes", onH ? "matrixElements" : "energy"},{"target",target},{"rounds",rounds},
                                             {"rasters",RasterName(itsImp->setup->opts.tolerances)},
                                             {"xcMesh", itsImp->setup->xcMesh.cellKind==qcMesh::UnitCellKind::Becke ? BeckeName(itsImp->setup->xcMesh) : std::string("uniform")},{"settled",itsImp->rastersSettled},
                                             {"auditRounds",itsImp->pilotCost.rounds},{"reconvergeIterations",(long)itsImp->pilotCost.reconvergeIterations}});
}

SolidCalculation::GridAuditResult SolidCalculation::XCMeshAudit(const qcMesh::MeshParams& reference) const
{
    if (itsImp->setup->planeWave) throw std::invalid_argument("SolidCalculation::XCMeshAudit: a plane-wave run has no GPW basis to rebuild");
    qchem::report::Timed timed("postSCF: XC-mesh audit (reference Hamiltonian + term matrices)");
    return CompareAtFixedD(itsImp->setup->opts.tolerances, rvec3_t(0,0,0), reference.cellKind==qcMesh::UnitCellKind::Becke ? BeckeName(reference)
                                                                                                                         : std::string("uniform XC"), &reference);
}

SolidCalculation::GridAuditResult SolidCalculation::EggBox(const rvec3_t& shift) const
{
    const SolidSetup& S=*itsImp->setup;
    if (S.planeWave) throw std::invalid_argument("SolidCalculation::EggBox: a plane-wave basis is translation-invariant by construction -- nothing to gate");
    if (S.imposed)   throw std::invalid_argument("SolidCalculation::EggBox: an IMPOSED group's SALC phases follow the origin, so the shifted blocks "
                                                 "would not match function for function -- gate the free twin");
    if (!S.recipe)   throw std::invalid_argument("SolidCalculation::EggBox: the setup was resolved from a PRE-BUILT basis that cannot follow the "
                                                 "shifted atoms -- resolve it from a GaussianBasisRecipe");
    if (norm(shift)==0.0) throw std::invalid_argument("SolidCalculation::EggBox: a zero shift measures nothing");
    if (!itsImp->hubbard.empty()) throw std::invalid_argument("SolidCalculation::EggBox: the +U projectors are built at the run's positions and would "
                                                              "not follow the shift -- gate the U = 0 twin");
    qchem::report::Timed timed("postSCF: egg-box (shifted basis + Hamiltonian + term matrices)");
    std::ostringstream os; os<<"shift ("<<shift.x<<", "<<shift.y<<", "<<shift.z<<") bohr";
    return CompareAtFixedD(S.opts.tolerances, shift, os.str());
}

SolidCalculation::GridAuditResult SolidCalculation::CompareAtFixedD(const BasisSet::Gaussian::GPWTolerances& raster, const rvec3_t& shift,
                                                                    const std::string& referenceName, const qcMesh::MeshParams* xcMesh,
                                                                    double ballScale) const
{
    namespace rpt=qchem::report;
    const SolidSetup&       S    = *itsImp->setup;
    const SolidCalcOptions& opts = S.opts;
    const bool eggBox = norm(shift)>0.0;
    const auto* wf=itsImp->scf ? itsImp->scf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("SolidCalculation: no wave function yet (the audits run on a converged D)");
    const SpinGroup g = S.polarized ? SpinGroup::Polarized : SpinGroup::UnPolarized;
    // The grid-audit verdicts of each side's rasters, ASKED of the basis that owns them (V-REPORT-WRITEONLY: never read back out of the
    // run report): the run's own here, the reference's once its V_long has run below (both are printed, so the predicted alias of
    // whichever side is the coarse one can be read -- the run need not be the side on the ball's own rasters).
    auto verdictsOf=[](const BasisSet::Gaussian::GpwRasterLadder& ladder)
    {
        rpt::json v=rpt::json::object();
        for (const auto& r : ladder.DeliveredRasterVerdicts()) v[r.audit]={{"status",r.status},{"alias",r.alias}};
        return v;
    };
    const rpt::json verdicts=verdictsOf(*itsImp->ladder);
    rpt::json refVerdicts=rpt::json::object();

    const SavedState st=SnapshotState(*itsImp->bs, *wf);
    std::unique_ptr<BasisSet::Complex_BS> bsRef;
    std::shared_ptr<const BasisSet::Gaussian::GpwRasterLadder> ladderRef;
    std::unique_ptr<qchem::Hamiltonian::cHamiltonian> hamRef;
    // FENCED from here to the fold: the reference basis and Hamiltonian report their grids (and grid audits) as they BUILD AND FIRST
    // EVALUATE (the raster levels and the V_long sweep are lazy); a nested run keeps that out of this run's `grids` section, which must keep
    // describing the rasters the run USED.  (The run's own evaluator emitted its audits once, at its own first use, so nothing of the
    // run's lands in the fence.)  RAII: an exception mid-audit must not leave the report nested.
    struct Fence { explicit Fence(const std::string& n) {rpt::Begin(n);} ~Fence() {rpt::End();} };
    std::optional<Fence> fence(std::in_place, "gridAudit reference ("+referenceName+")");
    // The reference's GEOMETRY: the run's, or every atom translated rigidly by `shift` (the egg-box), its molecular basis remade from
    // the frozen recipe there (as SolidSession builds a moved point).  Held to the end: the reference basis and Hamiltonian borrow them.
    std::optional<Lattice_3D> latShifted;
    std::shared_ptr<const BasisSet::Real_BS> molRef=itsImp->mol;
    if (eggBox)
    {
        UnitCell cellR(S.reference.GetUnitCell());
        std::vector<rvec3_t> Rs; for (const Site& s : cellR.Sites()) Rs.push_back(s.R+shift);
        cellR.MoveSites(Rs);
        latShifted.emplace(cellR, S.reference.GetLimits());
        molRef=S.recipe->Make(*latShifted->GetStructure(), S.pp);
    }
    const Lattice_3D& latRef = eggBox ? *latShifted : *itsImp->lat;
    const auto        stRef  = eggBox ? latShifted->GetStructure() : itsImp->st;
    {
        namespace L3 = BasisSet::Lattice;
        // The density ball: the run's, or ballScale x it -- on the knob that sets it (an auto ball is cutoffFactor x alpha_max).
        const bool autoBall = opts.densityEcut<0.0;
        L3::GPWParams p{
            .densityEcut = autoBall ? opts.densityEcut : ballScale*opts.densityEcut, .cutoffFactor = autoBall ? ballScale*opts.cutoffFactor : opts.cutoffFactor,
            .images = opts.images, .kShift = opts.kShift, .ladderFactor = opts.ladderFactor,
            .imposeSymmetry = S.imposed, .siteSpins = S.siteSpins, .freeMeshInvariant = opts.xcMeshSymmetric, .tol = raster,
            .hamPreservesReal = !opts.forceComplex};
        ResolveGpwPoint(latRef, molRef, S.pp, p);   // the reference's own ladder (its rasters / ball differ), sized for the same PP
        ladderRef=p.ladder;
        bsRef.reset(L3::GPWFactory(latRef, molRef, p));
        if (rpt::Depth()>0)   // the reference's grids + verdicts, as setup files the run's
        {
            rpt::Section grids("grids");
            L3::EmitGpwGrids(*ladderRef, *molRef, *p.imposition);
        }
        rpt::Scope here(eggBox ? "eggBox" : "gridAudit");   // the REFERENCE Hamiltonian's meshes report one level into the audit, not over the run's grids
        hamRef.reset(qchem::Hamiltonian::Factory(g, stRef, bsRef.get(), S.pp, xcMesh ? *xcMesh : S.xcMesh, opts.vxcFit, itsImp->hubbard));
    }
    auto cdRun=RestoreDensity(st, *itsImp->bs, g);
    auto cdRef=RestoreDensity(st, *bsRef, g);
    if (!cdRun || !cdRef) throw std::logic_error("GridAudit: the converged D does not restore onto "+std::string(cdRun ? "the reference" : "its own")+" basis");
    const qchem::ChargeDensity::cDM_CD* dRun=cdRun->get();
    const qchem::ChargeDensity::cDM_CD* dRef=cdRef->get();
    itsImp->ham->RefreshForDensity(itsImp->bs.get(), dRun);
    hamRef     ->RefreshForDensity(bsRef.get(),       dRef);

    GridAuditResult R; R.reference=referenceName;
    auto toVerdicts=[](const rpt::json& v)
    {
        std::vector<GridAuditVerdict> out;
        for (auto it=v.begin(); it!=v.end(); ++it) out.push_back({it.key(), it.value().value("status",std::string()), it.value().value("alias",0.0)});
        return out;
    };
    const std::vector<qchem::Hamiltonian::TermEnergy> eRun=itsImp->ham->MakeTermEnergies(dRun), eRef=hamRef->MakeTermEnergies(dRef);
    if (eRun.size()!=eRef.size()) throw std::logic_error("GridAudit: the reference Hamiltonian has a different term list");
    for (size_t t=0;t<eRun.size();t++)
    {
        // Matched by the ENERGY ENTRIES, not the label: a raster term's description states its own point count, which is what differs.
        if (eRun[t].entries!=eRef[t].entries)
            throw std::logic_error("GridAudit: term "+std::to_string(t)+" reports '"+eRun[t].entries+"' here, '"+eRef[t].entries+"' on the reference");
        R.terms.push_back({eRun[t].label, eRun[t].entries, 0.0, eRun[t].E-eRef[t].E, eRun[t].E});
        R.ΔE+=eRun[t].E-eRef[t].E;
    }
    const std::vector<Spin> spins = g==SpinGroup::Polarized ? std::vector<Spin>{Spin::Up,Spin::Down} : std::vector<Spin>{Spin::None};
    size_t channel=0;   // the snapshot's canonical block order: spin irreps outer, basis blocks inner
    for (Spin s : spins)
        for (size_t i=0;i<itsImp->bs->GetNumIBS();++i,++channel)
        {
            const auto* rb=itsImp->bs->GetRealIBS(i);
            const auto* rr=bsRef->GetRealIBS(i);
            if (bool(rb)!=bool(rr)) throw std::logic_error("GridAudit: block "+std::to_string(i)+" is real on one basis and complex on the other");
            if (channel>=st.blocks.size() || st.blocks[channel].ms!=int(s)) throw std::logic_error("GridAudit: the snapshot's block order is not the audit's");
            const SavedBlock& sb=st.blocks[channel];
            hmat_t<dcmplx> S(sb.S.rows());
            for (size_t a=0;a<S.rows();a++) for (size_t b=a;b<S.rows();b++) S(a,b)=sb.S(a,b);
            const std::unique_ptr<LASolver<dcmplx>> frame(LASolver<dcmplx>::Factory(opts.ortho, opts.orthoTol));
            frame->SetBasisOverlap(S);
            const std::string spin = s==Spin::Up ? "up" : s==Spin::Down ? "down" : "none";
            if (rb) FoldΔH(itsImp->ham->MakeTermMatrices(rb,s,dRun), hamRef->MakeTermMatrices(rr,s,dRef), sb, *frame, spin, i, R.terms);
            else    FoldΔH(itsImp->ham->MakeTermMatrices((*itsImp->bs)[i],s,dRun), hamRef->MakeTermMatrices((*bsRef)[i],s,dRef), sb, *frame, spin, i, R.terms);
        }
    {   // the largest element of each term, named: block 0's shells speak for every block (AoShellSource's contract)
        const auto* src = itsImp->bs->GetRealIBS(0) ? dynamic_cast<const BasisSet::AoShellSource*>(itsImp->bs->GetRealIBS(0))
                                                    : dynamic_cast<const BasisSet::AoShellSource*>((*itsImp->bs)[0]);
        if (src)
        {
            const auto shells=src->GetAoShells();
            for (auto& t : R.terms) if (t.maxAbsΔH>0.0)
            {
                t.maxElement.functionI=FunctionName(shells, *itsImp->st, t.maxElement.i);
                t.maxElement.functionJ=FunctionName(shells, *itsImp->st, t.maxElement.j);
            }
        }
    }
    for (const auto& t : R.terms) R.maxAbsΔH=std::max(R.maxAbsΔH, t.maxAbsΔH);
    refVerdicts=verdictsOf(*ladderRef);
    fence.reset();
    R.runVerdicts=toVerdicts(verdicts); R.referenceVerdicts=toVerdicts(refVerdicts);
    // The rigorous bound per raster term, at the COARSE side's predicted alias (whichever side that is), and the measured eta
    // (doc/GridEpsilonPlan.md §1: the bound and the calibration isolated; §5 M2c verdict 3: forces' eta on max|dH|, never on E).
    auto coarseX=[&](const char* key)
    {
        double a=0.0;
        for (const rpt::json* v : {static_cast<const rpt::json*>(&verdicts), static_cast<const rpt::json*>(&refVerdicts)}) if (v->contains(key)) a=std::max(a, (*v)[key].value("alias",0.0));
        return a>0.0 ? -std::log(a) : 0.0;
    };
    namespace G = BasisSet::Gaussian;
    for (auto& t : R.terms)
    {
        if (t.label.find("LONG-range")!=std::string::npos)
        {
            t.x=coarseX("auditVlocLong"); if (t.x<=0.0) t.x=coarseX("auditDensity");   // the kappa-sweep path has no V_long audit
            t.bound=G::kVlocLongTail.Bound(t.x); t.eta=t.bound>0.0 ? t.maxAbsΔH/t.bound : 0.0; t.etaOn="H";
        }
        else if (t.entries=="Eee")
        {
            t.x=coarseX("auditDensity"); t.bound=G::kDensityTail.Bound(t.x);
            t.eta=(t.bound>0.0 && t.E!=0.0) ? std::fabs(t.ΔE/t.E)/t.bound : 0.0; t.etaOn="E";
        }
        else if (t.entries=="Exc" && t.maxAbsΔH>0.0) t.etaOn="measured";
    }

    // REPORT, contemporaneously (pin 17): console table + the run report's `gridAudit` section.
    const std::streamsize prec0=std::cout.precision();
    if (eggBox) std::cout<<"[egg box] rigid translation at the converged D (fixed D: the SCF egg-box differs at second order), "<<RasterName(opts.tolerances)
                         <<", "<<R.reference<<": analytic terms must read 0"<<std::endl;
    else        std::cout<<"[grid audit] delivered error at the converged D: "<<RasterName(opts.tolerances)<<" (this run) vs "<<R.reference<<" (reference)"<<std::endl;
    std::cout<<std::scientific<<std::setprecision(3);
    rpt::json rows=rpt::json::array();
    for (const auto& t : R.terms)
    {
        std::cout<<"[grid audit]   max|ΔH_ij| "<<std::setw(10)<<t.maxAbsΔH<<"  ΔE "<<std::setw(11)<<t.ΔE<<" Ha  ";
        if      (t.etaOn=="H" || t.etaOn=="E") std::cout<<"[x="<<std::defaultfloat<<std::setprecision(4)<<t.x<<std::scientific<<std::setprecision(3)
                                                       <<" bound "<<t.bound<<" eta_"<<t.etaOn<<" "<<t.eta<<"]  ";
        else if (t.etaOn=="measured")          std::cout<<"[MEASURED: no rigorous bound for the XC raster]  ";
        std::cout<<t.label<<" ["<<t.entries<<"]"<<std::endl;
        rows.push_back({{"term",t.label},{"entries",t.entries},{"maxAbsDH",t.maxAbsΔH},{"dE",t.ΔE},{"E",t.E},
                        {"x",t.x},{"bound",t.bound},{"eta",t.eta},{"etaOn",t.etaOn},{"trDdH",t.trDdH},{"sNormDH",t.sNormDH}});
        if (t.maxAbsΔH>0.0)
        {
            const auto& m=t.maxElement;
            double sum=0.0; for (double c : t.trDdH) sum+=c;
            std::cout<<"[grid audit]     Tr(D ΔH) "<<sum<<" = channels";
            for (double c : t.trDdH) std::cout<<" "<<c;
            std::cout<<"  S-norm "<<t.sNormDH<<"  max element: spin "<<m.spin<<" block "<<m.block<<" ("<<m.i<<","<<m.j<<") |D_ij| "<<m.absD
                     <<": ["<<m.functionI<<"] x ["<<m.functionJ<<"]"<<std::endl;
            rows.back()["maxElement"]={{"spin",m.spin},{"block",m.block},{"i",m.i},{"j",m.j},{"absD",m.absD},{"functionI",m.functionI},{"functionJ",m.functionJ}};
        }
    }
    std::cout<<"[grid audit]   max|ΔH_ij| "<<std::setw(10)<<R.maxAbsΔH<<"  ΔE "<<std::setw(11)<<R.ΔE<<" Ha  TOTAL"<<std::endl;
    for (const auto& [who,vs] : {std::pair<const char*,const rpt::json*>{"the run's own",&verdicts},{"the reference's",&refVerdicts}})
        for (auto it=vs->begin(); it!=vs->end(); ++it)
            std::cout<<"[grid audit]   "<<who<<" "<<it.key()<<": "<<it.value().value("status",std::string())
                     <<" (predicted alias "<<it.value().value("alias",0.0)<<")"<<std::endl;
    std::cout<<std::defaultfloat<<std::setprecision(prec0);
    rpt::EmitAt(eggBox ? "eggBox" : "gridAudit", eggBox ? R.reference : RasterName(opts.tolerances)+std::string("_vs_")+R.reference,
                    rpt::json{{"run",RasterName(opts.tolerances)},{"reference",R.reference},{"terms",rows},
                         {"maxAbsDH",R.maxAbsΔH},{"dE",R.ΔE},{"verdicts",verdicts},{"referenceVerdicts",refVerdicts}});
    return R;
}

std::vector<qchem::Hamiltonian::ACBN0Estimate> SolidCalculation::EstimateHubbardU() const
{
    std::unique_ptr<qchem::Hamiltonian::ACBN0Estimator> est=itsImp->ham->MakeACBN0Estimator();
    if (!est) throw std::logic_error("SolidCalculation::EstimateHubbardU: this run carries no Hubbard manifold (SolidCalcOptions::hubbard)");
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!wf) throw std::logic_error("SolidCalculation::EstimateHubbardU: no wave function yet");
    qchem::Response::FeedOccupiedOrbitals(*est, *wf);       // the orbital walk lives beside MakeReference (D4)
    est->Write(std::cout);                                   // the estimator reports at its own activity (pin 17)
    return est->Details();
}

//---------------------------------------------------------------------------------------------------
//  THE MANIFOLD WALK (hoisted from ConvergeHubbardU's [U table] snapshot, 2026-10-30, for panel #10): the frontier (homo, lumo)
//  and each +U manifold's spin-resolved occupation sum_i f_i |<w_mu|psi_i>|^2 (k-weighted), through the +U term's OWN projectors
//  -- so it works on any XC mesh (the Becke-basin site moment needs basins).  Integrated observables, never a point sample.
namespace {
struct ManifoldWalk
{
    double homo=-std::numeric_limits<double>::infinity(), lumo=std::numeric_limits<double>::infinity();
    std::vector<double> up, dn;   //!< per manifold
};
ManifoldWalk WalkManifolds(const qchem::WaveFunction::cWaveFunction& wf, const qchem::Hamiltonian::HubbardChannels& hub)
{
    const size_t nM=hub.Channels().size();
    ManifoldWalk W; W.up.assign(nM,0.0); W.dn.assign(nM,0.0);
    auto walk=[&]<class U>(const qchem::Orbitals::TOrbitals<U>& os, const qchem::Irrep& ir)
    {
        const auto* blk=dynamic_cast<const qchem::BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
        if (!blk) throw std::logic_error("SolidCalculation::WalkManifolds: an orbital block that is not an Orbital_DFT_IBS");
        const double cap = ir.ms==qchem::Spin::None ? 2.0 : 1.0;   // a level is OCCUPIED when more than half full
        std::vector<const qchem::Orbitals::TOrbital<U>*> occ;
        for (const auto* o : os.template Iterate<qchem::Orbitals::TOrbital<U>>())
        {
            if (o->GetOccupation()>0.5*cap) W.homo=std::max(W.homo, o->GetEigenEnergy());
            else                            W.lumo=std::min(W.lumo, o->GetEigenEnergy());
            if (o->GetOccupation()>0.0) occ.push_back(o);
        }
        if (occ.empty()) return;
        mat_t<U> C(blk->GetNumFunctions(), occ.size());
        for (size_t i=0;i<occ.size();i++) { const vec_t<U>& c=occ[i]->GetCoeff(); for (size_t a=0;a<c.size();a++) C(a,i)=c[a]; }
        const auto amp=hub.ProjectorAmplitudes(*blk, C);
        const double w=ir.sym->GetWeight();
        for (size_t M=0;M<nM;M++)
        {
            double n=0.0;
            for (size_t i=0;i<occ.size();i++) for (size_t mu=0;mu<amp[M].rows();mu++) n+=occ[i]->GetOccupation()*std::norm(amp[M](mu,i));
            n*=w;
            if      (ir.ms==qchem::Spin::Up)   W.up[M]+=n;
            else if (ir.ms==qchem::Spin::Down) W.dn[M]+=n;
            else { W.up[M]+=0.5*n; W.dn[M]+=0.5*n; }        // the folded doublet: half in each channel
        }
    };
    for (const qchem::Irrep& ir : wf.GetQNs())
    {
        const qchem::Orbitals::Orbitals* os=wf.GetOrbitals(ir);
        if (const auto* c=dynamic_cast<const qchem::Orbitals::TOrbitals<dcmplx>*>(os)) walk(*c, ir);
        else if (const auto* r=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(os)) walk(*r, ir);
        else throw std::logic_error("SolidCalculation::WalkManifolds: an orbital set of unknown scalar");
    }
    return W;
}
} // namespace

std::vector<double> SolidCalculation::MakeHubbardMoments() const
{
    const auto* hub=itsImp->ham ? itsImp->ham->GetHubbardChannels() : nullptr;
    if (!hub) throw std::logic_error("SolidCalculation::MakeHubbardMoments: this run carries no Hubbard manifold (SolidCalcOptions::hubbard)");
    const auto* wf=itsImp->scf ? itsImp->scf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("SolidCalculation::MakeHubbardMoments: no wave function yet");
    const ManifoldWalk W=WalkManifolds(*wf, *hub);
    std::vector<double> m(W.up.size());
    for (size_t M=0;M<m.size();M++) m[M]=W.up[M]-W.dn[M];
    return m;
}

SolidCalculation::HubbardLoopResult SolidCalculation::ConvergeHubbardU(const SCFParams& params, const HubbardLoop& loop)
{
    using Est=HubbardLoop::Estimator;
    const double eV=27.211386245988;
    const bool LR = loop.estimator==Est::LinearResponse;
    const char* who = LR ? "LR" : "ACBN0";
    HubbardLoopResult R;
    const auto* hub=itsImp->ham->GetHubbardChannels();
    auto* target=itsImp->ham->GetHubbardUTarget();
    if (!hub || !target) throw std::logic_error("SolidCalculation::ConvergeHubbardU: this run carries no Hubbard manifold");
    const size_t nM=hub->Channels().size();
    // Which manifolds the loop SETS: ACBN0 estimates every manifold; linear response the perturbed set (hp.x's Hubbard sites),
    // fixed ONCE here -- after step 1 the sites carry a U, and a default re-derived later would mean a different set.
    if (LR)
    {
        R.estimated = loop.perturbed.empty() ? DefaultPerturbed(*hub) : loop.perturbed;
        for (size_t J : R.estimated) if (J>=nM) throw std::out_of_range("SolidCalculation::ConvergeHubbardU: a perturbed manifold index past the manifolds");
    }
    else for (size_t M=0;M<nM;M++) R.estimated.push_back(M);
    std::vector<double> Ucur;                                 // the U each manifold currently runs with (eV)
    for (size_t M=0;M<nM;M++) Ucur.push_back(M<itsImp->hubbard.size() ? itsImp->hubbard[M].U*eV : 0.0);
    R.Uin_eV=Ucur;
    R.scfConverged=itsImp->converged;

    // THE STATE an estimate is made FROM, for the [U table]: the band gap and each manifold's integrated moment.  Both read the
    // converged orbitals (pin: integrated observables, never a point sample of m(r)); the moment goes through the +U term's OWN
    // projectors, so it works on any XC mesh (the Becke-basin site moment needs basins).
    auto snapshot=[&]()
    {
        const auto* wf=itsImp->scf->GetWaveFunction();
        if (!wf) throw std::logic_error("SolidCalculation::ConvergeHubbardU: no wave function");
        const ManifoldWalk W=WalkManifolds(*wf, *hub);
        R.gap_eV.push_back(std::isfinite(W.lumo) && std::isfinite(W.homo) ? (W.lumo-W.homo)*eV : std::numeric_limits<double>::quiet_NaN());
        std::vector<double> m(nM); for (size_t M=0;M<nM;M++) m[M]=W.up[M]-W.dn[M];
        R.moment.push_back(m);
    };

    for (size_t n=0; n<loop.maxOuter; n++)
    {
        snapshot();
        R.outer=n+1;
        // ESTIMATE from the current (converged) state, through the ONE strategy face (H5): each strategy's input is its own
        // construction business -- ACBN0 is fed the orbitals, linear response is a solved response (+U frozen for the solve).
        std::unique_ptr<qchem::Hamiltonian::HubbardUEstimator> est;
        if (LR)
        {
            auto r=HubbardLinearResponse(loop.nq, loop.krylov, R.estimated);   // prints its own tables
            if (!r) { R.failure=r.Error().detail; std::cout<<"[LR loop "<<n<<"] the response FAILED -- the loop stops: "<<R.failure<<std::endl; break; }
            est=std::make_unique<qchem::Response::LR_HubbardU>(std::move(*r), hub->Channels());
        }
        else
        {
            auto a=itsImp->ham->MakeACBN0Estimator();
            qchem::Response::FeedOccupiedOrbitals(*a, *itsImp->scf->GetWaveFunction());
            est=std::move(a);
        }
        est->Write(std::cout);
        R.last=est->Evaluate();
        std::vector<double> Unext=Ucur;
        for (const auto& e : R.last) Unext[e.manifold]=e.U*eV;
        double dmax=0.0;
        for (size_t M : R.estimated) dmax=std::max(dmax, std::abs(Unext[M]-Ucur[M]));
        R.U_eV.push_back(Unext);
        {
            std::ostringstream os;
            os<<"["<<who<<" loop "<<n<<"] U"<<(LR?"":"_eff")<<"(eV):"; for (size_t M : R.estimated) os<<" "<<std::fixed<<std::setprecision(4)<<Unext[M];
            os<<"  max|dU|="<<std::setprecision(5)<<dmax<<(itsImp->converged ? "" : "  (SCF NOT converged)");
            std::cout<<os.str()<<std::endl;
        }
        if (dmax<loop.tolU_eV) { R.converged=true; break; }
        if (n+1==loop.maxOuter) break;
        // APPLY on the same Hamiltonian (the write is the loop's, never the estimator's: H5), then re-converge from the current density.
        for (size_t M : R.estimated) target->SetU(M, Unext[M]/eV);
        for (size_t M : R.estimated)                          // the record a save reads (SetU = one U)
            if (M<itsImp->hubbard.size()) { itsImp->hubbard[M].U=Unext[M]/eV; itsImp->hubbard[M].Uirrep.clear(); }
        Ucur=Unext;
        // A NEW STAGE, not a continued Iterate: the converged mixer history (eight near-zero Pulay residuals)
        // would extrapolate the first post-U step straight back onto the old density and report "converged"
        // in one iteration (measured 2026-09-21: E moved by exactly E_U, the orbitals not at all).  Re-seeding
        // through BuildStage gives a fresh iterator, mixer and accelerator over the current density -- exactly
        // what an annealed schedule does between its stages.
        BuildStage("+U self-consistency: U changed, a fresh history over the current density", itsImp->stageAccelAsked, std::move(itsImp->cd), params);
        auto out=Converge(params);
        R.scfConverged=bool(out);
    }
    // THE [U table]: one greppable line per estimated manifold (NOTES, 2026-10-04: "make plucking trivial").  U_0 is the FIRST estimate
    // (the response of the starting state, at U_in), U_SC the last; the gap and moment are those of the states they were estimated from.
    if (!R.U_eV.empty())
    {
        const auto ch=hub->Channels();
        for (size_t M : R.estimated)
        {
            std::ostringstream os;
            os<<std::fixed<<std::setprecision(4)<<"[U table] "<<who<<" "<<itsImp->setup->opts.label<<(itsImp->setup->opts.label.empty()?"":" ")
              <<"site"<<ch[M].site<<" l="<<ch[M].l<<"  U_in="<<R.Uin_eV[M]<<"  U_0="<<R.U_eV.front()[M]<<"  U_SC="<<R.U_eV.back()[M]
              <<" eV  outer="<<R.outer<<(R.converged?" CONVERGED":" NOT-converged")
              <<"  gap "<<R.gap_eV.front()<<" -> "<<R.gap_eV.back()<<" eV  m "<<R.moment.front()[M]<<" -> "<<R.moment.back()[M];
            std::cout<<os.str()<<std::endl;
        }
    }
    return R;
}

const qchem::ChargeDensity::cDM_CD* SolidCalculation::LastIterateDensity() const {return itsImp->cd.get();}

Outcome<qchem::Response::ChannelResponse,qchem::Response::ResponseFailure>
SolidCalculation::IndependentResponse(ivec3_t Nq) const
{
    const auto* hub=itsImp->ham->GetHubbardChannels();
    if (!hub) throw std::logic_error("SolidCalculation::IndependentResponse: this run carries no Hubbard manifold -- the "
                                     "response channels ARE the +U manifolds (list them at U=0 to probe without +U)");
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!wf) throw std::logic_error("SolidCalculation::IndependentResponse: no wave function yet");
    // NaN (a recipe with no [F,D], e.g. the Null accelerator) is passed on as UNMEASURED, never as 0.
    const double noise=std::isfinite(itsImp->lastCommutator) ? std::fabs(itsImp->lastCommutator)
                                                             : std::numeric_limits<double>::quiet_NaN();
    qchem::Response::Reference ref=qchem::Response::MakeReference(*wf, itsImp->lastOccupation,
        {.acrossK=!itsImp->setup->opts.perKFermiLevels, .acrossSpin=itsImp->setup->opts.spinsShareFermi}, noise);
    qchem::Response::AmplitudeProbe probe=qchem::Response::MakeHubbardProbe(ref, *wf, *hub);
    auto r=qchem::Response::IndependentResponse(ref, probe, Nq);
    if (r) r->Write(std::cout, 27.211386245988, "1/eV");     // reported at its own activity (pin 17)
    else   std::cout << "[chi0] FAILED: " << r.Error().detail << std::endl;
    return r;
}


// THE SELF-CONSISTENT RESPONSE OVER A q-MESH (R2 at q = 0, R3 at q != 0): the same Reference and probe as
// IndependentResponse, plus the AO <-> MO frame and the Hamiltonian's analytic kernel, linearised about the converged
// ORBITALS' density (the state the Reference describes -- not the mixed iterate, which on a Kerker/Pulay recipe has
// no D).  +U is FROZEN for the whole solve (Timrov eq 20) and U is read from the PERTURBED set's own inverse (§3d Q6,
// Q10), over the q-mesh's supercell.
Outcome<qchem::Response::QMeshResponse,qchem::Response::ResponseFailure>
SolidCalculation::HubbardLinearResponse(ivec3_t Nq, const KrylovParams& krylov, std::vector<size_t> perturbed,
                                        qchem::Response::HubbardResolution res, bool qStars, bool littleGroupK) const
{
    const auto* hub=itsImp->ham->GetHubbardChannels();
    if (!hub) throw std::logic_error("SolidCalculation::HubbardLinearResponse: this run carries no Hubbard manifold -- the "
                                     "response channels ARE the +U manifolds (list them at U=0 to probe without +U)");
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!wf) throw std::logic_error("SolidCalculation::HubbardLinearResponse: no wave function yet");
    // (An imposed run is fine: the perturbation BREAKS the imposed group, and since R3 step 3 no density path of the
    //  response symmetrizes -- the transition density is (k+q, k) pairs contracted through the B2 faces, which have
    //  no fold or star-average in them.  The step-1 refusal of an imposed raster / Γ-only run is retired with the
    //  R2 composite route it guarded; gate GPW_Si_ImposedGamma_Raster_eqFreeChi.  §3d finding 5.)
    if (perturbed.empty()) perturbed=DefaultPerturbed(*hub);
    const FrozenHubbard frozen(itsImp->ham->GetHubbardUTarget());
    const double noise=std::isfinite(itsImp->lastCommutator) ? std::fabs(itsImp->lastCommutator)
                                                             : std::numeric_limits<double>::quiet_NaN();
    // THE REFERENCE MUST BE EIGENPAIRS OF ITS OWN FOCK (OpenWork §3, 2026-10-27): R0 is a sum over states.  GDM once handed out
    // virtuals labelled with other virtuals' energies (2.7 Ha off) and chi0 came out 400x too small with nothing complaining.
    // THE BOUND: chi0's relative error is ~ residual / response gap.  1e-4 Ha keeps it <~1e-3 on a 0.1 Ha gap and still refuses
    // a mislabel (O(0.1-1) Ha) by >1000x.  Measured: Si k211 1e-10 (diagonalizing) / 8e-9 (GDM-finished); NiO k222 GDM restart
    // 1.03e-5 (its [F,D] noise 2.4e-6) -- a 1e-5 bound refused that healthy state.
    {
        const double resid=qchem::Response::EigenpairResidual(*wf, *itsImp->ham, itsImp->bs.get());
        std::cout << "[response] reference eigenpair residual max|C^H F C - diag eps| = " << resid << " Ha" << std::endl;
        if (!(resid<=1e-4))
        {
            std::ostringstream os;
            os << "the wave function's orbitals are not eigenpairs of their own Fock (max|C^H F C - diag eps| = " << resid
               << " Ha > 1e-4): the sum-over-states response would be wrong.  Converge the stage to a diagonalized state first.";
            std::cout << "[response] REFUSED: " << os.str() << std::endl;
            return Outcome<qchem::Response::QMeshResponse,qchem::Response::ResponseFailure>::Fail(
                {qchem::Response::ResponseFailure::Why::NotConverged, os.str()});
        }
    }
    const qchem::Response::Reference ref=qchem::Response::MakeReference(*wf, itsImp->lastOccupation,
        {.acrossK=!itsImp->setup->opts.perKFermiLevels, .acrossSpin=itsImp->setup->opts.spinsShareFermi}, noise);
    const auto frame=qchem::Response::MakeOrbitalFrame(ref, *wf);
    const auto hp=qchem::Response::MakeHubbardProbe(ref, *wf, *hub, res);
    const bool resolved = res!=qchem::Response::HubbardResolution::Manifold;
    std::vector<size_t> channels;                                 // every channel of a perturbed manifold
    for (size_t c=0;c<hp.manifold.size();c++)
        if (std::find(perturbed.begin(), perturbed.end(), hp.manifold[c])!=perturbed.end()) channels.push_back(c);
    if (resolved)                                                 // say what each channel IS before the long solve
        for (size_t c=0;c<hp.manifold.size();c++)
            std::cout << "[response] channel " << c << ": " << hp.probe.Label(c) << "  [" << hp.composition[c] << "]"
                      << (std::find(channels.begin(), channels.end(), c)!=channels.end() ? "  perturbed" : "") << std::endl;
    const auto D0=wf->GetChargeDensity();
    const auto kernel=itsImp->ham->MakeResponseKernel(itsImp->bs.get(), D0.get());
    // Each q's table is printed as that q finishes; the supercell summary (real-space blocks, U) once all are in.
    // U is printed, never consumed here (R4 owns consumption).
    using O=Outcome<qchem::Response::QMeshResponse,qchem::Response::ResponseFailure>;
    // ONE q PER STAR (KAxisPlan P3, Q9), on request: it ASSERTS the ground state has the captured group, so the response
    // audits one member.  Manifold resolution only (a channel per manifold: a site permutation and a phase, no D^l).
    std::optional<qchem::Response::QStars> stars;
    if (qStars)
    {
        if (resolved) throw std::invalid_argument("SolidCalculation::HubbardLinearResponse: qStars needs HubbardResolution::Manifold "
                                                  "(a resolved channel rotates within its manifold -- D^l, KAxisPlan P4)");
        std::vector<qchem::Response::ChannelSite> cs;
        for (const auto& m : hub->Channels()) cs.push_back({m.site, m.l});
        // The ground state's group: the Shubnikov ops of the decoration the run was BUILT with (all-zero = the space group),
        // and the sites' fractional positions -- both asked of the held lattice.
        const rmat3d_t Ainv=Invert(itsImp->lat->GetUnitCell().GetCellMatrix());
        std::vector<rvec3_t> siteFrac;
        itsImp->st->ForEachSite([&](int, const rvec3_t& r, bool)
        { siteFrac.push_back(rvec3_t(Ainv(1,1)*r.x+Ainv(1,2)*r.y+Ainv(1,3)*r.z, Ainv(2,1)*r.x+Ainv(2,2)*r.y+Ainv(2,3)*r.z,
                                     Ainv(3,1)*r.x+Ainv(3,2)*r.y+Ainv(3,3)*r.z)); });
        const auto qmesh=ref.QMesh(Nq);
        if (qmesh) stars=qchem::Response::FoldQMesh(*qmesh, itsImp->lat->ShubnikovOps(itsImp->setup->decoration), siteFrac, cs, perturbed);
    }
    // THE LITTLE GROUP of each (q, J) (KAxisPlan P5): the ground state's group MEASURED on the orbitals (ruling 1), offered the
    // Shubnikov ops of the decoration the run was built with.  ⚠ COLLINEAR: a σ op relabels the up/down blocks (MakeOrbitalImages).
    std::optional<qchem::Response::OrbitalImages> images;
    if (littleGroupK)
    {
        images=qchem::Response::MakeOrbitalImages(ref, *wf, itsImp->lat->ShubnikovOps(itsImp->setup->decoration));
        images->Write(std::cout);
    }
    auto solved=qchem::Response::LinearResponse(ref, frame, *kernel, hp.probe, Nq, krylov, resolved ? channels : perturbed,
                                                stars ? &*stars : nullptr, 1e-6, images ? &*images : nullptr);
    if (!solved) { std::cout << "[response] FAILED: " << solved.Error().detail << std::endl; return solved; }
    qchem::Response::QMeshResponse r=solved.TakeValue();
    if (resolved) r.group=hp.manifold;
    r.Write(std::cout);
    const auto& lab=r.atQ[0].labels;
    for (const auto& v : InterSiteHubbardV(r))
        std::cout << "[V] " << lab[v.from] << " -> " << lab[v.to]
                  << " R=(" << v.R.x << "," << v.R.y << "," << v.R.z << ")  d=" << std::fixed << std::setprecision(4) << v.distance
                  << " bohr  V=" << v.V*27.211386245988 << " eV" << std::defaultfloat << std::endl;
    return O::Ok(std::move(r));
}

std::vector<SolidCalculation::InterSiteV> SolidCalculation::InterSiteHubbardV(const qchem::Response::QMeshResponse& r) const
{
    std::vector<InterSiteV> out;
    if (r.atQ.empty()) return out;
    const auto& J=r.atQ[0].perturbed;
    const size_t nJ=J.size();
    const ivec3_t Nq=r.Nq;
    const rmat_t W=r.Interaction();
    const auto ch=itsImp->ham->GetHubbardChannels()->Channels();
    auto manifold=[&](size_t c){ return r.group.empty() ? c : r.group[c]; };   // a resolved response's channel -> its manifold
    std::vector<rvec3_t> tau;                                     // the cell's site positions (Cartesian, a.u.)
    itsImp->st->ForEachSite([&](int, const rvec3_t& x, bool){ tau.push_back(x); });
    const auto A=itsImp->lat->GetUnitCell().GetCellMatrix();      // columns = a_1..a_3 (Matrix3D is 1-based)
    auto cart=[&](double f1, double f2, double f3)
    { return rvec3_t(A(1,1)*f1+A(1,2)*f2+A(1,3)*f3, A(2,1)*f1+A(2,2)*f2+A(2,3)*f3, A(3,1)*f1+A(3,2)*f2+A(3,3)*f3); };
    std::vector<ivec3_t> Rs;                                      // ToRealSpace's image order
    for (int x=0;x<Nq.x;x++) for (int y=0;y<Nq.y;y++) for (int z=0;z<Nq.z;z++) Rs.push_back(ivec3_t(x,y,z));
    for (size_t a=0;a<nJ;a++)                                     // a home-cell site
        for (size_t rr=0;rr<Rs.size();rr++)
            for (size_t b=0;b<nJ;b++)
            {
                if (rr==0 && b==a) continue;                      // the diagonal is U
                // The minimum image of (b, R) seen from (a, 0) in the supercell spanned by Nq.a_i.
                double best=std::numeric_limits<double>::infinity(); ivec3_t bestR;
                for (int n1=-1;n1<=1;n1++) for (int n2=-1;n2<=1;n2++) for (int n3=-1;n3<=1;n3++)
                {
                    const ivec3_t R(Rs[rr].x+n1*Nq.x, Rs[rr].y+n2*Nq.y, Rs[rr].z+n3*Nq.z);
                    const rvec3_t d=tau[ch[manifold(J[b])].site]+cart(R.x,R.y,R.z)-tau[ch[manifold(J[a])].site];
                    const double dist=std::sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
                    if (dist<best-1e-9) {best=dist; bestR=R;}
                }
                out.push_back({.from=J[a], .to=J[b], .R=bestR, .distance=best, .V=W(a, rr*nJ+b)});
            }
    std::stable_sort(out.begin(), out.end(), [](const InterSiteV& x, const InterSiteV& y){return x.distance<y.distance;});
    return out;
}

Outcome<qchem::Response::SelfConsistentResponse,qchem::Response::ResponseFailure>
SolidCalculation::HubbardLinearResponse(const KrylovParams& krylov, std::vector<size_t> perturbed) const
{
    using O=Outcome<qchem::Response::SelfConsistentResponse,qchem::Response::ResponseFailure>;
    auto r=HubbardLinearResponse(ivec3_t(1,1,1), krylov, std::move(perturbed));
    if (!r) return O::Fail(r.Error());
    return O::Ok(std::move(r->atQ[0]));
}

// THE FINITE-DIFFERENCE CROSS-CHECK (LR-cDFT): +-alpha on manifold J, +U frozen, three fresh stages.
Outcome<SolidCalculation::FiniteDifferenceChi,SCFFailure>
SolidCalculation::HubbardFiniteDifferenceChi(size_t J, double alpha, const SCFParams& params, bool reseed)
{
    using O=Outcome<FiniteDifferenceChi,SCFFailure>;
    const auto* hub=itsImp->ham->GetHubbardChannels();
    auto* target=itsImp->ham->GetHubbardUTarget();
    if (!hub || !target) throw std::logic_error("SolidCalculation::HubbardFiniteDifferenceChi: this run carries no Hubbard manifold");
    const auto channels=hub->Channels();
    if (J>=channels.size()) throw std::out_of_range("SolidCalculation::HubbardFiniteDifferenceChi: no such manifold");
    // ⛔ An IMPOSED run's SCF star-averages the XC density every iteration (the Becke (ρ, m) projector, the raster
    // average) -- so a +-alpha run on one site would respond with a SYMMETRIZED XC potential and an unsymmetrized
    // Hartree one: neither the physical response nor the linear response's.  Found 2026-09-29 on NiO AFM-II.
    if (itsImp->setup->imposed)
        throw std::invalid_argument("SolidCalculation::HubbardFiniteDifferenceChi: an IMPOSED run symmetrizes the XC density "
                                    "every SCF iteration, and the alpha*P perturbation breaks the group -- run the cross-check FREE");
    const double alpha0 = J<itsImp->hubbard.size() ? itsImp->hubbard[J].alpha : 0.0;   // restored at the end
    // The channel occupations of the CURRENT orbitals, through the same projector amplitudes the probe uses.
    auto occupations=[&]()
    {
        const auto* wf=itsImp->scf->GetWaveFunction();
        if (!wf) throw std::logic_error("SolidCalculation::HubbardFiniteDifferenceChi: no wave function");
        const qchem::Response::Reference ref=qchem::Response::MakeReference(*wf, itsImp->lastOccupation,
            {.acrossK=!itsImp->setup->opts.perKFermiLevels, .acrossSpin=itsImp->setup->opts.spinsShareFermi}, std::numeric_limits<double>::quiet_NaN());
        const auto probe=qchem::Response::MakeHubbardProbe(ref, *wf, *hub);
        qchem::Response::BlockPairs D;                      // the state in its own MO basis: diag(occupation)
        for (const Irrep& ir : wf->GetQNs())
        {
            const auto* os=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(wf->GetOrbitals(ir));
            if (!os) throw std::logic_error("SolidCalculation::HubbardFiniteDifferenceChi: a real TRIM block (run with forceComplex)");
            qchem::Response::cmat_t X(os->GetNumOrbitals(), os->GetNumOrbitals(), dcmplx(0.0));
            size_t i=0;
            for (const auto* o : os->Iterate<Orbitals::TOrbital<dcmplx>>()) {X(i,i)=o->GetOccupation(); i++;}
            D.m.push_back(X);
        }
        const cvec_t n=probe.Measure(Symmetry::Invariant(), D);
        rvec_t out(n.size());
        for (size_t I=0;I<n.size();I++) out[I]=n[I].real();
        return out;
    };
    FiniteDifferenceChi R;
    for (const auto& c : channels) R.labels.push_back("site"+std::to_string(c.site)+" l="+std::to_string(c.l));
    R.perturbed=J; R.alpha=alpha;
    std::optional<SCFFailure> failed;
    // RESEED keeps the converged ground-state density ASIDE (the +-alpha runs start from the seed, so they never
    // read it) and restores from it: the restore then starts AT the fixed point, rather than re-converging it from
    // the seed.
    std::unique_ptr<qchem::ChargeDensity::cDM_CD> ground;
    if (reseed) ground=std::move(itsImp->cd);
    // THE RESTART RECIPE (OpenWork §3 "A carried density is re-filled at kT = 0", 2026-10-21): every stage below starts from a
    // converged density (the restore always; the +-alpha runs unless reseeded), where the ground state's mixing step is wrong --
    // a relax of 1.0 diverged on a restart (Si Γ 2026-09-29) and a Pulay relax 1.0 amplified FeS2's first-iterate kick.  So:
    // the LINEAR step capped at 0.3 (CLAUDE.md: a restart's linear fallback runs at 0.2-0.3), and at kT = 0 a Ladder enters at its
    // GDM rung (a smeared stage starts on DIIS: see Restart).  The caller's params and accelerator options are restored after.
    // LINEAR D-mixing only: a Kerker/Pulay step is a different alpha (the k211 FD gate converges its +-alpha runs at 0.45 and
    // NOT at 0.3 -- measured 2026-10-21).
    assert(itsImp->scf);
    SCFParams restart=params;   // the caller's Mixer is kept (an Auto one is resolved per stage, against that stage's accelerator)
    if (itsImp->scf->ResolveMixing(params).Mixer==SCFParams::Mixing::Linear)
    {   // the START and the CEILING: capping only the start let the adaptive step grow back to 1.0 (SCFParams::MaxRelaxRo)
        restart.StartingRelaxRo=std::min(restart.StartingRelaxRo, 0.3);
        restart.MaxRelaxRo     =std::min(restart.MaxRelaxRo,      0.3);
    }
    const bool startedOnGDM=itsImp->accOpts.ladderStartOnGDM;
    itsImp->accOpts.ladderStartOnGDM = params.SmearingkT==0.0;   // kT = 0 only (see Restart)
    const SCFParams& stageParams = reseed ? params : restart;
    {
        const FrozenHubbard frozen(target);
        for (double a : {+alpha, -alpha})
        {
            target->SetPerturbation(J, alpha0+a);
            std::cout << "[LRT] manifold " << R.labels[J] << ": alpha = " << alpha0+a << " Ha, +U frozen -- re-converging" << std::endl;
            // A FRESH stage (see ConvergeHubbardU): from the SEED, or from the current state.  The +U term stays
            // frozen at the GROUND STATE's n either way -- the freeze lives in the term, not in the density.
            if (reseed) BuildStage("finite-difference LRT: alpha "+std::to_string(alpha0+a)+" Ha, from the seed", itsImp->stageAccelAsked, nullptr, stageParams);
            else        BuildStage("finite-difference LRT: alpha "+std::to_string(alpha0+a)+" Ha", itsImp->stageAccelAsked, std::move(itsImp->cd), stageParams);
            auto out=Converge(stageParams);
            if (!out) {failed=out.Error(); break;}
            (a>0 ? R.nPlus : R.nMinus)=occupations();
        }
        target->SetPerturbation(J, alpha0);
    }
    // RESTORE the unperturbed, unfrozen state -- the run is left as it was found, whatever happened above.
    BuildStage("finite-difference LRT: restore the unperturbed state", itsImp->stageAccelAsked, reseed ? std::move(ground) : std::move(itsImp->cd), restart);
    R.restored=bool(Converge(restart));
    itsImp->accOpts.ladderStartOnGDM=startedOnGDM;
    std::cout << "[LRT] unperturbed state restored: " << (R.restored ? "converged" : "NOT converged") << std::endl;
    if (failed) return O::Fail(*failed);
    R.chi=rvec_t(R.nPlus.size());
    std::cout << "[LRT] chi_IJ = dn_I/dalpha_J, J = " << R.labels[J] << ", alpha = " << alpha << " Ha (1/Ha):";
    for (size_t I=0;I<R.chi.size();I++)
    {
        R.chi[I]=(R.nPlus[I]-R.nMinus[I])/(2*alpha);
        std::cout << "  " << R.labels[I] << " " << std::setprecision(8) << R.chi[I];
    }
    std::cout << std::endl;
    return O::Ok(std::move(R));
}

// THE FD ORACLE'S DOOR (friend only, forward.H).
qchem::Hamiltonian::cHamiltonian& SolidCalculation::ResponseHamiltonian() const {return *itsImp->ham;}
const WaveFunction::cWaveFunction& SolidCalculation::ResponseWaveFunction() const
{
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!wf) throw std::logic_error("SolidCalculation: no wave function yet");
    return *wf;
}

// The caller's observer is SWAPPED IN behind the facade's own (AttachProbes composes the two), so
// attaching telemetry late cannot silently disarm the outcome detectors.
void   SolidCalculation::OnIteration(Observer obs)      { itsImp->onIteration=std::move(obs); AttachProbes(); }
bool   SolidCalculation::DidConverge()    const         { return itsImp->converged; }
size_t SolidCalculation::IterationCount() const         { return itsImp->scf->GetIterationCount(); }

const RunDiagnostics& SolidCalculation::Diagnostics() const { return itsImp->diag; }

const qcMesh::MeshParams& SolidCalculation::ResolvedXCMesh() const { return itsImp->setup->xcMesh; }
const BasisSet::Complex_BS& SolidCalculation::Basis() const { return *itsImp->bs; }


//====================================================================================================
//  THE GEOMETRY POINT and THE SOLID SESSION (doc/CalcNetworkReview.md §7, step 8)
//====================================================================================================
std::unique_ptr<SolidCalculation> SolidCalculation::AtGeometry(std::shared_ptr<const SolidSetup> setup, const Lattice_3D& lat,
                                                               std::shared_ptr<const BasisSet::Real_BS> mol, const SavedState* carry,
                                                               const std::vector<SCFStage>& schedule, bool lowdin, Observer onIteration)
{
    if (schedule.empty())
        throw std::runtime_error("SolidCalculation: an EMPTY schedule has no meaning -- pass at least one stage.");
    std::unique_ptr<SolidCalculation> c(new SolidCalculation(BuildOnly{}));
    c->itsImp->onIteration=std::move(onIteration);
    // A SURFACE POINT: its client consumes matrix elements (a geometry -- forces, a path), so a pilot sizes on max|dH_ij|.
    c->itsImp->consumes=BasisSet::Gaussian::ClientConsumes::MatrixElements;
    // ITS PILOT STARTS AT THE BOUND (the controller's condition 2): a geometry run settles near the rigorous L(eps) (eta = 1) anyway --
    // NaF's path settled at x 17.1 / 16.7 against 19.4 / 21.6 -- so starting there replaces a ball-raster SCF plus a re-convergence by one
    // SCF the audit usually just confirms (it can still raise; it never lowers).  An energy run keeps the ball start: it settles at x 0-7,
    // and the bound would cost it up to ~10x the raster on every run.
    if (!setup->rastersResolved && !setup->planeWave)
    {
        auto seeded=std::make_shared<SolidSetup>(*setup);
        BasisSet::Gaussian::GridAccuracy rig=setup->opts.accuracy; rig.eta=BasisSet::Gaussian::AccuracyCalibration{};
        seeded->opts.tolerances=rig.Stated(setup->opts.tolerances);
        setup=std::move(seeded);
    }
    // ITS XC ROUTE FOLLOWS THE CONSUMER (user ruling 2026-11-01, pin 31): v_xc[rho] on a fixed uniform raster does not follow the atoms and
    // cannot meet max|dH| <= eps at any affordable x (Si Gamma: 3e-5 per element at x 150, a power law), so a geometry run whose XC was the
    // facade's choice (Auto) takes the atom-centred Becke mesh, which pin 30 makes admissible on a solid (H_R = {e}); the pilot then sizes
    // it on the nested ladder.  A STATED mesh is obeyed, and so is the CP2K-parity veto (policy.beckeXC) -- a declared deviation.
    if (!setup->planeWave && setup->opts.xcMesh.cellKind==qcMesh::UnitCellKind::Auto && setup->xcMesh.cellKind!=qcMesh::UnitCellKind::Becke
        && theRunPolicy().BeckeXC())
    {
        auto becke=std::make_shared<SolidSetup>(*setup);
        becke->xcMesh=qcMesh::BeckeXCParams();
        std::cout<<"[XC grid choice] the client consumes H (a geometry run): a uniform XC raster does not follow the atoms -> BECKE "
                 <<BeckeName(becke->xcMesh)<<std::endl;
        setup=std::move(becke);
    }
    {
        qchem::report::Timed residue("setup: facade ctor (residue -- decisions between the named buckets)");
        c->BuildPoint(std::move(setup), lat, std::move(mol));
        if (!carry)
        {
            c->BuildRun(nullptr, schedule.front().params, schedule.front().accelerator, /*restored*/false, /*exactResume*/false);
            (void)c->Converge(schedule.front().params);
        }
        else
        {
            // THE CARRY: the previous point's D, block for block, onto this point's basis -- the same recipe at moved atoms, so the
            // block layout (k, size, realness) is the one RestoreDensity asserts.  A WorkingDensity, NOT a restored one: the
            // Hamiltonian moved with the atoms, so the Auto policy's carried-stage branch (Null + Pulay on D) applies.
            const SpinGroup g=c->itsImp->setup->polarized ? SpinGroup::Polarized : SpinGroup::UnPolarized;
            std::optional<SavedState> moved;
            if (lowdin) { moved=*carry; LowdinCarry(*moved, *c->itsImp->bs, g); carry=&*moved; }
            auto cd=RestoreDensity(*carry, *c->itsImp->bs, g);
            if (!cd) throw std::logic_error("SolidCalculation::AtGeometry: the carried density does not fit the new point: "+cd.Error().details);
            const SCFStage& last=schedule.back();
            c->BuildRun(cd.TakeValue(), last.params, last.accelerator, /*restored*/false, /*exactResume*/false);
            (void)c->Converge(last.params);
        }
    }
    if (!carry && schedule.size()>1) c->RunRestOf(schedule);
    c->CalibrateRasters(schedule.back().params);   // C2 step 3: on an unresolved setup (a session's first point) the pilot sizes the rasters
    return c;
}

SavedState SolidCalculation::Carry() const
{
    const auto* wf=itsImp->scf ? itsImp->scf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("SolidCalculation::Carry: no wave function yet");
    return SnapshotState(*itsImp->bs, *wf);
}

//! One geometry point: the structure its molecular basis is built over, that basis, and the calculation -- declared in that
//! order so the calculation (which holds the Bloch basis over \c mol) dies first and the atoms last.
struct SolidSession::Point
{
    std::shared_ptr<const Structure>         cell;
    std::shared_ptr<const BasisSet::Real_BS> mol;
    std::unique_ptr<SolidCalculation>        calc;
};

SolidSession::SolidSession(std::shared_ptr<const SolidSetup> setup, std::vector<SCFStage> schedule, SolidCalculation::Observer onIteration,
                           Carry carry)
    : itsSetup(std::move(setup)), itsSchedule(std::move(schedule)), itsObserver(std::move(onIteration)), itsCarry(carry)
{
    if (itsSchedule.empty()) throw std::invalid_argument("SolidSession: an EMPTY schedule has no meaning -- pass at least one stage.");
    if (!itsSetup->imposed && !itsSetup->planeWave && itsSetup->opts.xcMeshSymmetric)
        throw std::invalid_argument("SolidSession: a FREE setup with xcMeshSymmetric=true builds its XC mesh invariant under the group "
                                    "detected at EACH geometry, so the mesh switches as the atoms move and the surface jumps -- resolve "
                                    "the setup with SolidCalcOptions::xcMeshSymmetric=false (the plain atom-centred mesh)");
}
SolidSession::~SolidSession() = default;

PESGeometry SolidSession::Reference() const
{
    PESGeometry g;
    g.cell=itsSetup->reference.GetUnitCell().GetCellMatrix();
    itsSetup->cell->ForEachSite([&](int Z, const rvec3_t& r, bool){ g.R.push_back(r); g.Z.push_back(Z); });
    return g;
}
const SolidCalculation* SolidSession::Current() const {return itsPoint ? itsPoint->calc.get() : nullptr;}

Outcome<PESPoint, PESFailure> SolidSession::ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R,
                                                      const PESRequest& want)
{
    using O=Outcome<PESPoint, PESFailure>;
    using Why=PESFailure::Why;
    const SolidSetup& S=*itsSetup;
    const UnitCell&   ref=S.reference.GetUnitCell();
    const size_t      nAtoms=ref.GetNumAtoms();
    if (R.size()!=nAtoms)
        throw std::invalid_argument("SolidSession::ConvergeAt: "+std::to_string(R.size())+" positions for "+std::to_string(nAtoms)+" atoms");
    if (!cell)
        throw std::invalid_argument("SolidSession::ConvergeAt: a solid's geometry has a cell -- pass Reference().cell for a fixed-cell step");
    if (want.forces || want.stress)
        return O::Fail({Why::NoGradient, std::string(want.stress ? "stress" : "forces")+": no term implements its geometry gradient yet "
                        "(doc/CalcNetworkReview.md §3.3, step 9) -- difference ConvergeAt energies meanwhile"});
    // THE CELL IS THE SETUP'S, for now: the Bloch basis, the grids (sized from the cell) and the imposed group were all built for it.
    // The signature takes the cell so a variable-cell driver needs no new face when the basis/grid side can follow one.
    {
        const rmat3d_t& A=ref.GetCellMatrix();
        double dA=0.0, scale=0.0;
        for (size_t i=1;i<=3;i++) for (size_t j=1;j<=3;j++) { dA=std::max(dA, std::fabs((*cell)(i,j)-A(i,j))); scale=std::max(scale, std::fabs(A(i,j))); }
        if (dA>1e-12*scale)
            return O::Fail({Why::CellChange, "the cell moved (max |dA| = "+std::to_string(dA)+" bohr): a setup's basis and grids are built for "
                            "ONE cell, and a variable-cell step is not supported yet -- resolve a new setup at that cell"});
    }

    // THE POINT'S CELL: the setup's, atoms at R.
    UnitCell cellR(ref);
    double moved=0.0;
    for (size_t a=0;a<nAtoms;a++) moved=std::max(moved, norm(R[a]-cellR.Sites()[a].R));
    cellR.MoveSites(R);
    const Lattice_3D latR(cellR, S.reference.GetLimits());

    // THE IMPOSED GROUP must survive (Q4): a geometry that lowers it is a NEW setup, the user's choice -- never a silent re-fold.
    if (S.imposed && moved>0.0)
    {
        auto order=[&](const Lattice_3D& l){ return S.siteSpins.empty() ? l.GetSpaceGroup().Order() : l.ShubnikovOps(S.siteSpins).size(); };
        if (const size_t now=order(latR), was=order(S.reference); now<was)
            return O::Fail({Why::SymmetryBroken, "these positions keep "+std::to_string(now)+" of the imposed group's "+std::to_string(was)
                            +" operations -- symmetrize the step, or resolve a new setup with the lower group"});
    }

    auto point=std::make_unique<Point>();
    point->cell=latR.GetStructure();
    if (moved==0.0)          point->mol=S.mol;   // the reference itself: the setup's own basis, bit for bit
    else if (S.planeWave)    point->mol=S.mol;   // a plane-wave basis does not follow the atoms (the Gaussian one is unread)
    else if (!S.recipe)
        throw std::invalid_argument("SolidSession: the atoms moved, but the setup was resolved from a PRE-BUILT basis that cannot follow "
                                    "them -- resolve it from a GaussianBasisRecipe");
    else
    {
        // THE FROZEN RECIPE at R -- and the conditioning window it was trimmed in (Q3).  The vet loop, asked about this geometry
        // with the frozen trim as its base, must find NOTHING more to drop at the run's orthoTol; if it does, the point is refused
        // rather than silently run on a different basis (re-trimming per step would make the surface discontinuous).
        const GaussianBasisRecipe& recipe=*S.recipe;
        if (S.opts.orthoTol>0.0)
        {
            auto make=[&](const BasisSet::Gaussian::ShellTrim& t)
            {
                GaussianBasisRecipe r=recipe;
                r.trim.shells.insert(r.trim.shells.end(), t.shells.begin(), t.shells.end());
                return r.Make(*point->cell, S.pp);
            };
            auto vet=BasisSet::Lattice::VetStageTrim(latR, make, {.images=S.opts.images, .kShift=S.opts.kShift}, S.opts.orthoTol);
            if (!vet.trim.empty())
            {
                std::ostringstream os; vet.trim.Write(os);
                return O::Fail({Why::BasisConditioning, "at these positions the frozen basis is near-dependent at orthoTol "
                                +std::to_string(S.opts.orthoTol)+": the ortho step would drop "+os.str()+" -- the basis left the window "
                                "its trim was decided in (resolve a new setup at a geometry nearer this one)"});
            }
            point->mol=std::move(vet.mol);
        }
        else point->mol=recipe.Make(*point->cell, S.pp);
    }

    // THE CARRY, then release the previous point before building the next (one point's memory at a time).
    std::optional<SavedState> carry;
    if (itsPoint) carry=itsPoint->calc->Carry();
    itsPoint.reset();
    point->calc=SolidCalculation::AtGeometry(itsSetup, latR, point->mol, carry ? &*carry : nullptr, itsSchedule, itsCarry==Carry::Lowdin,
                                         itsObserver);
    itsPoint=std::move(point);

    auto out=itsPoint->calc->Result();
    if (!out) return O::Fail({Why::SCF, out.Error().details, out.Error()});
    // THE RASTERS ARE FROZEN at the first converged point (C2 step 3): every later point samples on the rasters this one's pilot settled,
    // as it runs on the basis the setup froze -- a raster re-sized per point would make the surface jump.
    if (!itsSetup->rastersResolved) itsSetup=itsPoint->calc->SizedSetup();
    return O::Ok(PESPoint{.energy=out->Energy(), .forces={}, .iterations=out->IterationCount()});
}

} //namespace qchem
