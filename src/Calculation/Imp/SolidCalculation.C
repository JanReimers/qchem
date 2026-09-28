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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
module qchem.SolidCalculation;

import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;   // GaussianSharpness::MaxExponent -- alpha_max (the NARROW face)
import qchem.BasisSet.Orbital_1E_IBS;          // Real_OIBS / Complex_OIBS (the per-irrep bases to iterate)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH local PP -> alpha_pp
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.Mesh.XCPolicy;                    // XCMeshSharpness / ResolveXCMesh (the grid decision)
import qchem.ElectronConfiguration.Crystal;    // Crystal_EC
import qchem.Symmetry.Irrep;                   // Irrep, Spin
import qchem.Blaze;                            // blazem::eigen (the atomic radial projection)
import qchem.Symmetry.Atom.Spherical;          // AtomicSymmetry::Getl -- the pseudo-atom's l-orbital (the atomic +U radial)
import qchem.AtomCalculation;                  // the pseudo-atom in the block's own primitives (HubbardManifold::atomicRadial)
import qchem.BasisSet.AoShellSource;           // the site's shells (their exponents) for that pseudo-atom
import qchem.Orbitals;                         // TOrbitals/TOrbital -- the ACBN0 feed (EstimateHubbardU)
import qchem.BasisSet.Orbital_DFT_IBS;         // the block an orbital set lives on (EstimateHubbardU)
import qchem.RunPolicy;                        // the declared CP2K deviations (doc/OpenWork.md N5/T5)
import qchem.Parallel;                         // WorkerThreads() -- half of the thread state a row must state
import qchem.Reporting;                       // report::Timed -- the facade's own setup buckets (1.1, 2026-09-06)

namespace qchem
{

//---------------------------------------------------------------------------------------------------
//  The run's SHARPNESS -- an above-SCFIterator decision input.
//
//  Both sources come off ABSTRACT capability faces via the sanctioned abstract->abstract cross-cast, so
//  nothing here touches a concrete basis or a concrete PP model:
//    alpha_max -- BasisSet::Gaussian::LatticeSum1E::MaxExponent(), documented there as "the GPW
//                 density-grid cutoff floor".
//    alpha_pp  -- BasisSet::SpeciesRadialField_Gaussian::AsGaussians(Z, Short), whose terms carry
//                 alpha = 1/(2 r_loc^2).  A model with no closed-Gaussian short part does not implement
//                 the face; that leaves alpha_pp at 0, which the selector reads as "not measurable" --
//                 NOT as "smooth".
//---------------------------------------------------------------------------------------------------
static qcMesh::XCMeshSharpness GatherSharpness(const Lattice_3D& lat, const BasisSet::Real_BS& mol,
                                               const SolidCalcOptions& o, bool imposed)
{
    qcMesh::XCMeshSharpness s;
    s.cellEdge = lat.GetUnitCell().GetMaximumCellEdge();
    s.nAtoms   = int(lat.GetUnitCell().GetNumAtoms());
    s.imposed  = imposed;   // the RESOLVED value: CP2K_COMPAT can veto it, and the grid cost follows
    for (auto ibs : const_cast<BasisSet::Real_BS&>(mol).Iterate<BasisSet::Real_OIBS>())
        // THE NARROW FACE, since the 2026-09-08 ISP split: this asks for ONE number, so it names the
        // three-method GaussianSharpness capability and not the seventeen-method periodic aggregate it
        // used to cross-cast to.  Nothing here wants a lattice sum, a collocation or a stream fold.
        if (const auto* sh=dynamic_cast<const BasisSet::Gaussian::GaussianSharpness*>(ibs))
            { s.alphaMax = sh->MaxExponent(); break; }
    for (const auto& [element, valence] : o.species)
    {
        const int Z = int(thePeriodicTable().GetZ(element));
        const Pseudopotential::HGH_LocalPotential loc = Pseudopotential::GetGTH(element,"LDA",valence).local;
        const auto* g = static_cast<const BasisSet::SpeciesRadialField_Gaussian*>(&loc);
        for (const auto& t : g->AsGaussians(Z, BasisSet::FieldRange::Short)) s.alphaPP = std::max(s.alphaPP, t.alpha);
    }
    return s;
}

//---------------------------------------------------------------------------------------------------
struct SolidCalculation::Imp
{
    SolidCalcOptions opts;
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
    qcMesh::MeshParams                          xcMesh;          // AFTER Auto resolution
    std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scf;
    std::unique_ptr<qchem::ChargeDensity::cDM_CD>         cd;    // the converged density (outlives the WF)
    //! m(r) of the converged state.  OWNED: WaveFunction::GetSpinDensity() BUILDS it and hands over the
    //! unique_ptr (V1.25).  EMPTY on an unpolarized run -- that WF does not implement the face at all (V1.17).
    std::unique_ptr<SolidCalculation::sf_t>               spin;
    SCFAccelerators::SolidAcceleratorOptions    accOpts;
    SCFAccelerators::Type                       stageAccel = SCFAccelerators::Type::DIIS;  //!< the CURRENT stage's, for the banner
    bool   converged = false;
    //! What the LAST converge stage filled with and how well it converged its eigenvalues -- the linear
    //! response's reference is built over the SAME occupancy rule, gated against the SAME measured noise
    //! (doc/LinearResponsePlan.md E1).  The commutator is the final iteration's [F,D].
    qchem::OccupationConfig lastOccupation;
    double lastCommutator = std::numeric_limits<double>::quiet_NaN();
    bool   imposed   = false;   //!< the run imposed a symmetry, so its solution must still carry it (T2)
    double charge    = 0.0;
    //! The OUTCOME DETECTORS' record (N1/T3-T4).  Accumulates across every attempt this object makes --
    //! an annealed schedule and a re-Converge are CONTINUATIONS of one SCF, and an order that died in
    //! stage 1 must not become invisible because stage 2 started from the corpse.
    RunDiagnostics diag;
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
             <<" C="<<o.cutoffFactor<<" raster="<<(o.raster==BasisSet::PlaneWave::RasterPolicy::BallOnly
                                                   ? "BallOnly" : "AliasFree")
             <<" xcMesh="<<MeshName(xc.cellKind);
    // The radial/angular pair describes an ATOM-CENTRED mesh; printing it beside "Uniform" would be
    // stating a number that had no effect on the run, which is the opposite of self-description.
    if (xc.cellKind==qcMesh::UnitCellKind::Becke) std::cout<<" (nR="<<xc.nRadial<<" L="<<xc.angularDegree<<")";
    else                                          std::cout<<" (eCut="<<xc.eCut<<" Ha)";
    std::cout<<std::endl;
    std::cout<<"["<<o.label<<" run] symmetry: "
             <<(imposed ? (o.greyImposition ? "IMPOSED (grey -- the erasure control)"
                                            : "IMPOSED (Shubnikov from the decoration)")
                        : "FREE")
             // A VETOED imposition must never be silent -- a caller who asked for one and did not get it
             // would otherwise read the row as though it had, which is the mirror of the hazard that made
             // imposeSymmetry opt-in in the first place.
             <<(o.imposeSymmetry && !imposed ? "  [asked for, VETOED by CP2K_COMPAT]" : "")
             <<";  threads: OMP_NUM_THREADS="<<(omp?omp:"unset")
             <<" GPW_OMP_THREADS="<<qchem::WorkerThreads()<<" (BLAS pinned to 1)"
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
            std::cout<<(M.radial.empty() ? ", radial: every shell (CP2K)" : ", radial: ONE contracted ("+std::to_string(M.radial.size())+" shells"+(M.atomicRadial ? ", pseudo-atom" : "")+(M.orthoAtomic ? ", ORTHO-atomic" : ", atomic")+")")
                     <<", site group "<<(M.siteOps.empty() ? 1 : M.siteOps.size())<<" ops / grey "
                     <<(M.greyOps.empty() ? (M.siteOps.empty() ? 1 : M.siteOps.size()) : M.greyOps.size())<<")";
        }
        std::cout<<"   [* = differs from CP2K]"<<std::endl;
    }
}
// WHAT THE SCF IS DOING: the mixer, the accelerator, the occupation machinery.  Emitted per Converge,
// because that is where these take effect -- and an anneal changes them stage by stage.
static void EmitSCFBanner(const std::string& label, const SCFParams& p, SCFAccelerators::Type acc)
{
    // ⚠ THE PRECONDITIONER AND THE HISTORY ARE INDEPENDENT, AND THE BANNER USED TO HIDE IT (2026-08-28).
    // It printed "Pulay(depth 8, start 5)" whenever PulayDepth>0 and never mentioned Kerker, so a reader
    // of a benchmark row would conclude the G-space preconditioner had been swapped OUT.  It has not:
    // PulayParams carries G0 into PulayMixer, so the two COMPOSE -- which matters
    // because that composition is exactly CP2K's own recipe (its &MIXING BETA is the Kerker damping
    // denominator, applied alongside BROYDEN_MIXING).  This line is the one doc/Benchmark.md tells people
    // to copy beside a row, so it has to say both.
    std::cout<<"["<<label<<" scf] mixer: "
             <<(p.KerkerG0>0.0 ? "Kerker(G0="+std::to_string(p.KerkerG0)+")"
                               : std::string("LINEAR D-mixing (no G-space preconditioner)"))
             <<(p.PulayDepth>0 ? " + Pulay history(depth "+std::to_string(p.PulayDepth)+", start "
                                 +std::to_string(p.PulayStart)+")"
                               : std::string())
             <<" alpha="<<p.StartingRelaxRo
             <<";  XC rho source: "<<(theRunPolicy().XCFromDM() ? "rho[D] WHOLESALE"
                                     : p.XCCuspDeficit         ? "rho_mix + cusp deficit"
                                     :                           "rho_mix")
             <<";  accel: "<<AccelName(acc)
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
static std::vector<double> AtomicRadial(const BasisSet::Real_BS& mol, const Structure& st, size_t site, int l,
                                        const std::vector<std::pair<std::string,int>>& species);

SolidCalculation::SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                   const SolidCalcOptions& opts, const SCFParams& params,
                                   const SCFAccelerators::SolidAcceleratorOptions& acc)
    : itsImp(std::make_unique<Imp>())
{
    // ★ THE CTOR'S RESIDUE BUCKET (doc/ParallelAndOraclePlan.md 1.1(a)).  Timed is EXCLUSIVE, so this
    // outer scope charges itself exactly what the named setup buckets below do NOT -- the decisions
    // between them (sharpness gather, magnetic decoration, IonicSAD targets, the irrep/EC build, the
    // banner).  With the same bracket on Converge and on Iterate, the ledger PARTITIONS the run and
    // "everything not in a bucket" stops being a place time can hide.
    qchem::report::Timed residue("setup: facade ctor (residue -- decisions between the named buckets)");
    itsImp->opts    = opts;
    itsImp->accOpts = acc;
    itsImp->st      = lat.GetStructure();

    namespace L3 = BasisSet::Lattice;
    // THE RESOLVED IMPOSITION (doc/OpenWork.md N5, user 2026-08-26).  The caller's flag AND the policy's
    // permission: CP2K parity forbids the capability outright, because CP2K does no symmetry work at all
    // (see RunPolicy::SymmetryImposition for the evidence), and every banked recipe asks for it -- so the
    // veto has to live here rather than in each recipe.  ONE name from here on; nothing below reads
    // opts.imposeSymmetry again, or the two would drift.
    const bool imposed = opts.imposeSymmetry && theRunPolicy().SymmetryImposition();
    // THE WORKING-TYPE DECISION (doc/RealComplexPlan.md §3, Step 3c-3): a block is real ⇔ its irrep is
    // (TRIM) ∧ every term preserves realness.  This composition root builds the LDA GPW stack --
    // kinetic, the PP trio, Hartree, XC -- every member of which PreservesReal(), so the term half is
    // TRUE here; it is asserted against the constructed Hamiltonian below, so a future term that
    // breaks realness must also flip this forecast.  forceComplex is the §6 ansatz-policy downgrade.
    const bool hamPreservesReal = !opts.forceComplex;
    // THE MAGNETIC DECORATION (S3).  An imposition is only SHUBNIKOV if the factory is told which sites
    // carry which spin; without it an "imposed AFM" run star-averages under the SPATIAL group and the
    // order is erased.  DERIVED here rather than passed in: every input is already an option this class
    // owns, so a caller cannot get it inconsistent with the seed it asked for (the IonicSAD targets are
    // the same resolution the seed itself uses -- see ChargeDensity::IonicSADTargets).
    //   unpolarized  -> {} (no channels to decorate)
    //   greyImposition -> {} DELIBERATELY: that arm is the erasure NEGATIVE CONTROL
    std::vector<int> siteSpins;
    if (!imposed) { /* nothing to decorate: siteSpins are meaningless without an imposition */ }
    else if (!opts.siteSpins.empty() && !opts.greyImposition)
        siteSpins = opts.siteSpins;                 // STATED: a specific ordering, not the seed's guess
    else if (opts.multiplicity>=1 && !opts.greyImposition)
    {
        const std::map<size_t,int> targets = (opts.seed==qchem::ChargeDensity::SeedStrategy::IonicSAD)
                                           ? qchem::ChargeDensity::IonicSADTargets(itsImp->st.get(), "LDA")
                                           : std::map<size_t,int>{};
        siteSpins = qchem::ChargeDensity::MagneticDecoration(itsImp->st.get(), "LDA", targets);
    }
    // ★ THE FACADE'S SETUP IS BUCKETED (doc/ParallelAndOraclePlan.md 1.1, 2026-09-06).  Every MnO row drives
    // SolidCalculation directly, and until now only the phases INSIDE these calls were timed -- so the basis
    // build, the Hamiltonian ctor's non-mesh half, the seed and the ortho landed in no bucket at all.  That
    // was ~48 s of a 119 s threaded MnO run, i.e. the largest non-scaling block in it, and unnamed.  Labels
    // match the GPW_SCF harness's so the two paths' ledgers read the same.
    {
        qchem::report::Timed timed("setup: GPW basis build");
        itsImp->bs.reset(L3::GPWFactory(lat, mol, L3::GPWParams{
            .densityEcut = opts.densityEcut, .cutoffFactor = opts.cutoffFactor, .raster = opts.raster,
            .images = opts.images, .kShift = opts.kShift, .ladderFactor = opts.ladderFactor,
            .imposeSymmetry = imposed, .siteSpins = siteSpins,
            .hamPreservesReal = hamPreservesReal}));
    }
    // THE RUN REPORTS ITS OWN BASIS AND GRIDS (TE phase 2, 2026-09-15; the reporting rule: each class emits
    // contemporaneously with its own activity).  These two sections used to be emitted by the integration
    // test driver around this same construction; the facade is the orchestrator now.  Order matters and is
    // the one the driver established: the conditioning pre-flight on the ANALYTIC overlap first (no grids
    // needed), then the grid ladder (EmitGpwGrids is what forces EnsureLevels).  Only when a run report is
    // open -- the emitters are inert otherwise -- and never an abort: a dependent basis is the caller's to
    // handle through a pivoted ortho, and an unhandled one fails loudly in the iterator.
    if (qchem::report::Depth()>0)
    {
        qchem::report::Timed timed("setup: vet basis (analytic S) + grids report");
        size_t nRemoved=0;
        {
            qchem::report::Log("vetting basis conditioning");
            qchem::report::Section basis("basis");
            nRemoved=L3::VetGpwConditioning(*itsImp->bs);
        }
        if (nRemoved>0)
            qchem::report::Log("basis is rank-deficient ("+std::to_string(nRemoved)+" redundant functions, see basis.removed)");
        qchem::report::Log("building grid ladder");
        L3::EmitGpwGrids(*itsImp->bs);
    }

    // DECISION 1 -- the XC quadrature.  Resolve Auto HERE, once, from facts about the run.  Downstream
    // consumers compare ==Becke, so an unresolved Auto would silently read as Uniform; resolving it at the
    // point the spec enters the Hamiltonian is what makes that impossible.
    // The XC-mesh route is a DECLARED CP2K deviation (RunPolicy::BeckeXC), resolved here beside the
    // imposition for the same reason: the policy belongs to the run, the sizing belongs to qcMesh.
    itsImp->xcMesh = qcMesh::ResolveXCMesh(opts.xcMesh, GatherSharpness(lat, *mol, opts, imposed),
                                           theRunPolicy().BeckeXC());

    // DECISION 2 -- the spin bookkeeping.  multiplicity -> (nUp,nDown), with the parity check that catches
    // a singlet asked of an odd electron count BEFORE integer division silently empties a channel.
    const int twoS = opts.multiplicity>1 ? opts.multiplicity-1 : opts.Nelec%2;
    const bool polarized = opts.multiplicity>=1;
    if ((opts.Nelec-twoS)%2!=0 || twoS>opts.Nelec)
        throw std::runtime_error("SolidCalcOptions: multiplicity "+std::to_string(opts.multiplicity)
                                 +" parity disagrees with Nelec "+std::to_string(opts.Nelec));
    auto irreps = itsImp->bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry Sum_k)
    itsImp->ec = std::make_unique<Crystal_EC>(irreps, (opts.Nelec+twoS)/2, (opts.Nelec-twoS)/2,
                                              opts.globalFermi, opts.spinsShareFermi);

    // THE HUBBARD MANIFOLDS' SITE GROUPS (step 5 increment 2): each manifold is labelled by the stabiliser of
    // its site in the DECLARED DECORATION's Shubnikov group (sigma=None ops), grey parentage beside it.  The
    // decoration is the same one the imposition would use, derived here whether or not the run imposes: the
    // site symmetry of the ordered state is the physical question either way.  A caller who filled siteOps
    // by hand keeps them.
    std::vector<qchem::Hamiltonian::HubbardManifold> hubbard=opts.hubbard;
    if (!hubbard.empty())
    {
        std::vector<int> decoration=siteSpins;
        if (decoration.empty() && polarized && !opts.greyImposition)
        {
            if (!opts.siteSpins.empty()) decoration=opts.siteSpins;
            else
            {
                const std::map<size_t,int> targets = (opts.seed==qchem::ChargeDensity::SeedStrategy::IonicSAD)
                                                   ? qchem::ChargeDensity::IonicSADTargets(itsImp->st.get(), "LDA")
                                                   : std::map<size_t,int>{};
                decoration=qchem::ChargeDensity::MagneticDecoration(itsImp->st.get(), "LDA", targets);
            }
        }
        for (auto& M : hubbard)
        {
            if (M.atomicRadial && M.radial.empty()) M.radial=AtomicRadial(*mol, *itsImp->st, M.site, M.l, opts.species);
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
        itsImp->ham.reset(qchem::Hamiltonian::Factory(
            polarized ? qchem::SpinGroup::Polarized : qchem::SpinGroup::UnPolarized,
            itsImp->st, itsImp->bs.get(), opts.species, "LDA", itsImp->xcMesh, opts.vxcFit, hubbard));
    }
    // The forecast crosscheck: the basis was built on the promise that every term preserves realness
    // (the AND's term half, above); the constructed Hamiltonian must agree, or real blocks were built
    // that its terms cannot serve.
    assert((!hamPreservesReal || itsImp->ham->PreservesReal()) &&
           "SolidCalculation: the term stack no longer preserves realness -- update the forecast above");

    // DECISION 3 -- the accelerator, by policy, through the public typed door.
    itsImp->stageAccel = opts.accelerator;
    itsImp->accel.reset(SCFAccelerators::Factory(opts.accelerator, acc));

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
    std::unique_ptr<qchem::ChargeDensity::cChargeDensity> seed;
    {
        qchem::report::Timed timed("setup: seed density (SAD/IonicSAD atomic solves)");
        seed.reset(qchem::ChargeDensity::MakeSeedDensity<dcmplx>(opts.seed, itsImp->bs.get(), itsImp->st.get(),
                                                                 itsImp->ec.get(), polarizedHam));
    }
    // MEASURED ONLY WHERE IT CAN MEAN SOMETHING.  SiteMoments rasters BOTH spin channels before it can
    // discover it has no basins to integrate over, and unlike the per-iteration probe -- which rides a
    // raster the Fock build has already made for this density serial -- nothing has been built yet at
    // seed time, so this one would be paid in full.  An unpolarized run has m == 0 identically and a
    // uniform XC mesh has no basins at all, and the facade knows both facts here without asking.
    if (seed && polarizedHam && itsImp->xcMesh.cellKind==qcMesh::UnitCellKind::Becke)
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
            seed.release(), itsImp->st.get(), opts.ortho, opts.orthoTol);   // the iterator consumes it in Init
    }

    // Observe from iteration ONE: the ctor converges, so an observer attached afterwards has already
    // missed stage 0 (see SolidCalcOptions::onIteration).
    AttachProbes();

    // MOM CONTINUATION FROM THE SEED (S0e): pin the reference to the seed's OWN freshly-filled occupied
    // subspace before iteration 1, so the CONFIGURATION the seed chose survives, not merely its density.
    // No-op unless SCFParams::UseMOM is also set.
    if (opts.momFromSeed) itsImp->scf->AdoptMOMReference(*itsImp->scf->GetWaveFunction());

    // T2 is gated on the run having IMPOSED something: an imposition is an ASSERTION about the answer,
    // and only then is losing the order a contradiction rather than physics (a FREE run that finds m=0
    // has found m=0, and must never be second-guessed for it).
    itsImp->imposed = imposed;

    EmitRunBanner(opts, itsImp->xcMesh, itsImp->st->GetNumAtoms(), imposed, hubbard);
    (void)Converge(params);   // the ctor ATTEMPTS; the caller faces the result via Result()
}

// THE ANNEALED CTOR.  Delegates to the single-stage one with the FIRST stage's parameters/accelerator --
// so the graph is built exactly once and stage 0 runs as the plain ctor's convergence -- then continues
// through the rest of the schedule.  (Delegating rather than duplicating the build is what keeps the two
// ctors from drifting; every DECISION in the build is made in one place.)
SolidCalculation::SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                   const SolidCalcOptions& opts, const std::vector<SCFStage>& schedule,
                                   const SCFAccelerators::SolidAcceleratorOptions& acc)
    : SolidCalculation(lat, mol,
                       schedule.empty() ? opts : [&]{ SolidCalcOptions o=opts;
                                                      o.accelerator=schedule.front().accelerator; return o; }(),
                       schedule.empty() ? SCFParams{} : schedule.front().params, acc)
{
    if (schedule.empty())
        throw std::runtime_error("SolidCalculation: an EMPTY anneal schedule has no meaning -- pass at "
                                 "least one stage, or use the single-SCFParams constructor.");
    // Stage 0 has already run (the delegated ctor converged it) -- report it, then continue.
    EmitStageSummary(itsImp->opts.label, 0, schedule.size(), schedule.front().params.SmearingkT,
                     itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->scf->GetEnergy());
    for (size_t s=1; s<schedule.size(); ++s)
    {
        BuildStage(schedule[s].accelerator, std::move(itsImp->cd));
        (void)Converge(schedule[s].params);
        EmitStageSummary(itsImp->opts.label, s, schedule.size(), schedule[s].params.SmearingkT,
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

const ScalarFunction<double>* SolidCalculation::Converged::SpinDensity() const {return itsImp->spin.get();}
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
    if (itsImp->imposed && d.OrderCollapsed())
        return fail(SCFFailure::Why::OrderLost,
                    "the run IMPOSED a magnetic symmetry and did not keep it: the integrated site moment "
                    "peaked at "+std::to_string(d.OrderPeak())+" e and was dead ( < "
                    +std::to_string(RunDiagnostics::kCollapseFraction)+" of that) from step "
                    +std::to_string(d.OrderDiedAt())+" onward, so this energy answers a different "
                    "question than the one asked");

    return O::Ok(Converged(itsImp.get()));
}

// ONE STAGE: a FRESH Hamiltonian + accelerator (a kT change must not carry stale DIIS history across the
// re-seed), an iterator seeded from \a carried when there is one, and MOM continuation from \a prev.
// Ownership: the iterator OWNS ham and accel and deletes them, so \a prev must outlive the adoption and is
// released immediately after it -- the same ordering the driver this replaces used.
void SolidCalculation::BuildStage(SCFAccelerators::Type accType,
                                  std::unique_ptr<qchem::ChargeDensity::cDM_CD> carried)
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
    itsImp->stageAccel = accType;
    itsImp->accel.reset(SCFAccelerators::Factory(accType, itsImp->accOpts));

    auto prev = std::move(itsImp->scf);      // held ONLY until the new stage has copied its MOM reference
    itsImp->scf = carried
        ? std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
              carried.release(), itsImp->st.get(), itsImp->opts.ortho, itsImp->opts.orthoTol)
        : std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
              itsImp->opts.seed, itsImp->st.get(), itsImp->opts.ortho, itsImp->opts.orthoTol);
    // MOM continuation across TEMPERATURE: stage 0 self-adopts the seed's own freshly-filled occupied
    // subspace, every later stage adopts the stage before it -- so the CHARACTER the hot stage settled on
    // survives the fresh wavefunction, exactly as the density does.
    if (itsImp->opts.momFromSeed)
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
    auto userObs = itsImp->opts.onIteration;
    itsImp->scf->SetObserver([this,userObs](const qchem::SCFIterator::SCFProgress& p)
    {
        // The detectors read the observable off the breakdown the iteration produced (R1.0h) -- the same
        // number the trace reported, not a second pull.
        itsImp->diag.itsOrder.push_back(MaxSiteMoment(p.eb.charge.siteMoments, itsImp->diag.itsHasBasins));
        itsImp->diag.itsEee  .push_back(p.eb["Eee"]);
        // A Null accelerator computes NO [F,D]: SCFProgress then carries 0 (the trace hides it by tag), which
        // must not read as a MEASURED zero eigenvalue noise (found 2026-09-27, the NiO R0 gate) -- NaN instead.
        itsImp->lastCommutator = itsImp->stageAccel==SCFAccelerators::Type::Null
                               ? std::numeric_limits<double>::quiet_NaN() : p.commutator;
        if (userObs) userObs(p);
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
            BuildStage(schedule[s].accelerator, std::move(itsImp->cd));
        last = Converge(schedule[s].params);
        EmitStageSummary(itsImp->opts.label, s, schedule.size(), schedule[s].params.SmearingkT,
                         itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->scf->GetEnergy());
    }
    return last;   // the FINAL stage's -- the earlier ones exist to feed it
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Converge(const SCFParams& params)
{
    assert(itsImp->scf);
    // ★ THE CONVERGE RESIDUE (doc/ParallelAndOraclePlan.md 1.1(a)).  The work AFTER Iterate returns is not
    // bookkeeping: GetChargeDensity() builds a fresh composite density and GetSpinDensity() rasters BOTH
    // spin channels, once per stage -- and an annealed run has a stage per schedule entry.
    qchem::report::Timed residue("scf: converge (residue -- final density + m(r) extraction)");
    EmitSCFBanner(itsImp->opts.label, params, itsImp->stageAccel);
    itsImp->lastOccupation = {.useMOM=params.UseMOM, .momStartIter=(int)params.MOMStartIter,
                              .kT=params.SmearingkT, .momPenalty=params.MOMSmearPenalty};   // the iterator's own conversion
    itsImp->scf->Iterate(params);
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
    if (wf->GetSpinGroup()==SpinGroup::Polarized)
        itsImp->spin = wf->GetSpinDensity();   // BUILT for us; we take it
    else
        itsImp->spin.reset();
    // THE RESULT LINE -- the run reports its own outcome, contemporaneously, at the stated precision
    // (doc/Benchmark.md compares codes at the 1e-5 Ha level; the terms stay at default width, they are read
    // for structure).  Was the integration-test driver's line until TE phase 2 (2026-09-15).
    {
        const qchem::EnergyBreakdown& E=itsImp->scf->GetEnergy();
        const std::streamsize prec0=std::cout.precision();
        std::cout << "["<<itsImp->opts.label<<"] "<<(itsImp->converged ? "CONVERGED" : "NOT converged")
                  << " iters="<<itsImp->scf->GetIterationCount()<<" charge="<<itsImp->charge
                  << " Eelec="<<E.GetElectronicEnergy()
                  << " Etot="<<std::setprecision(10)<<E.GetTotalEnergy()<<std::setprecision(prec0)
                  << "  (Ekin="<<E["Kinetic"]<<" Een="<<E["Een"]<<" Eee="<<E["Eee"]<<" Exc="<<E["Exc"]
                  << " Enn="<<E["Enn"]<<" E_alphaZ="<<E["E_alphaZ"]<<")" << std::endl;
    }
    return Outcome_();
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Result() const {return Outcome_();}

qchem::EnergyBreakdown SolidCalculation::LastIterateTerms()  const {return itsImp->scf->GetEnergy();}
double                 SolidCalculation::LastIterateCharge() const {return itsImp->charge;}
// THE ATOMIC RADIAL OF A HUBBARD MANIFOLD (DFT+U increment 3, 2026-09-21).  The physically meaningful +U
// manifold is ONE radial d (p) function, not every shell of that l on the site (CP2K's mechanism manifold,
// on which an ACBN0 U comes out near-bare because the KS states live entirely inside it).  The radial hp.x
// uses is the pseudo-atom's own l orbital, so that is what is built here -- the SAME GTH pseudo-atom the
// valence-basis generator validated, LDA, unpolarized -- and then PROJECTED onto the site's shells.
//
// ⚠ IN A COMPLETE POOL, NOT IN THE SITE'S OWN SHELLS (found 2026-09-21).  A solid's valence basis is trimmed
// for the SOLID: the MnO VA span dropped the diffuse 0.18 d shell, and the free 4s2 3d5 atom run in those
// shells has eps(3d) 77 mHa too shallow and a too-compact "3d" (the SR span, with two s exponents, has no 4s
// at all and puts the 3d at -1.04 Ha).  In a 16+16 even-tempered pool the same pseudo-atom reproduces CP2K's
// ATOM code to 3 / 0.3 / 0.1 mHa (E, eps 4s, eps 3d; UTCalculation ValenceBasisGen.MnQ7PseudoAtomInALargePool).
// So: the TRUE pseudo-atom orbital from the pool, projected in the overlap metric onto the site's l shells --
// the best the basis can hold of it -- with the captured norm printed (1 = the basis holds it exactly).
static std::vector<double> AtomicRadial(const BasisSet::Real_BS& mol, const Structure& st, size_t site, int l,
                                        const std::vector<std::pair<std::string,int>>& species)
{
    // The site: position, Z, element, valence.
    rvec3_t R; int Z=0; { size_t i=0; st.ForEachSite([&](int z, const rvec3_t& r, bool){ if (i==site) { R=r; Z=z; } i++; }); }
    const std::string element=thePeriodicTable().GetSymbol(Z);
    int Zion=0; for (const auto& [el,val] : species) if (el==element) Zion=val;
    if (Zion==0) Zion=Pseudopotential::GetGTH(element, "LDA").zion;
    // The site's shells by l, in the block's shell order (the order Hubbard_U::Select walks too).
    std::map<int,std::vector<double>> byL;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(mol).Iterate<BasisSet::Real_OIBS>())
    {
        const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(ibs);
        if (!src) continue;
        for (const auto& sh : src->GetAoShells())
        {
            if (norm(sh.center-R)>1e-8) continue;
            if (sh.exponents.size()!=1)
                throw std::runtime_error("SolidCalculation: the atomic +U radial needs UNCONTRACTED primitives on site "
                                         +std::to_string(site)+" (a contracted shell: not this increment)");
            byL[sh.rep->L()].push_back(sh.exponents[0]);
        }
        break;                                                  // every block carries the same shells
    }
    if (!byL.count(l)) throw std::runtime_error("SolidCalculation: no l="+std::to_string(l)+" shell on site "+std::to_string(site));
    // The pseudo-atom in a COMPLETE even-tempered pool, every l up to the highest the site carries (an
    // occupied l the site lacks would throw inside: a manifold on such a site makes no sense anyway).
    auto pool=[](int n, double emin, double emax){ std::vector<double> e; for (int i=0;i<n;i++) e.push_back(emin*std::pow(emax/emin, double(i)/(n-1))); return e; };
    AtomCalcOptions o;
    o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=Zion;
    int lmax=0; for (const auto& [ll,es] : byL) lmax=std::max(lmax,ll);
    for (int ll=0; ll<=lmax; ll++) o.exponentsByL.push_back({ll, pool(16, 0.05, 200.0)});
    SCFParams p; p.MinVirial=1e30;                              // no virial under a PP (the valence generator's rule)
    // 120, not the default 20: the 16-exponent pool is a big basis and a late-3d atom is a long descent.  Ni q10
    // (10 valence electrons against Mn's 7) was still at Δρ = 2e-4 on iteration 20 -- descending smoothly, one
    // gate away -- and the throw below turned that into "no atomic +U radial" for the whole NiO run (2026-09-22).
    // 60 was the first fix and is NOT enough margin: a d^8 cation needed 120 to clear a 2e-6 Ha degeneracy
    // between the last occupied and first empty minority orbital (2026-09-23).  The cost of an extra iteration
    // is milliseconds; the cost of being one short is a dead run, so buy margin.
    p.NMaxIter=120;
    AtomCalculation atom(Z, Z-Zion, o, p);
    if (!atom.IsConverged()) throw std::runtime_error("SolidCalculation: the "+element+" pseudo-atom did not converge -- no atomic +U radial");
    // Its lowest occupied l orbital, as coefficients over the pool.
    std::vector<double> cpool; double eBest=1e300;
    for (const Irrep& ir : atom.GetIrreps(Spin::None))
    {
        const auto* as=dynamic_cast<const Symmetry::Atom::AtomicSymmetry*>(ir.sym.get());
        if (!as || int(as->Getl())!=l) continue;
        const auto* os=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(atom.Orbitals(ir));
        if (!os) continue;
        for (const auto* orb : os->template Iterate<qchem::Orbitals::TOrbital<double>>())
            if (orb->IsOccupied() && orb->GetEigenEnergy()<eBest)
            {
                eBest=orb->GetEigenEnergy();
                const vec_t<double>& c=orb->GetCoeff();
                cpool.assign(c.size(), 0.0); for (size_t k=0;k<c.size();k++) cpool[k]=c[k];
            }
    }
    const std::vector<double> ep=pool(16, 0.05, 200.0);
    if (cpool.size()!=ep.size())
        throw std::runtime_error("SolidCalculation: the "+element+" pseudo-atom has no occupied l="+std::to_string(l)+" orbital -- an atomic +U radial needs one");
    // Overlap of unit-normalised r^l e^{-a r^2} radials: (2 sqrt(ab)/(a+b))^{l+3/2}.
    auto ov=[l](double a, double b){ return std::pow(2.0*std::sqrt(a*b)/(a+b), l+1.5); };
    // The pool orbital's own norm (the two codes' radial conventions must agree: unit-normalised primitives).
    double n2=0.0;
    for (size_t a=0;a<ep.size();a++) for (size_t b=0;b<ep.size();b++) n2+=cpool[a]*cpool[b]*ov(ep[a],ep[b]);
    if (std::abs(n2-1.0)>1e-6)
        throw std::logic_error("SolidCalculation: the pseudo-atom's l orbital is not unit-normalised over normalised primitives (|chi|^2="
                               +std::to_string(n2)+") -- the atomic and periodic radial conventions differ");
    // PROJECT onto the site's l shells (overlap metric): S_site r = b, b_s = <g_s|chi>.
    const std::vector<double>& es=byL[l];
    const size_t ns=es.size();
    rsmat_t Ss(ns); rvec_t b(ns);
    for (size_t s=0;s<ns;s++)
    {
        for (size_t t=s;t<ns;t++) Ss(s,t)=ov(es[s],es[t]);
        double bs=0.0; for (size_t k=0;k<ep.size();k++) bs+=cpool[k]*ov(es[s],ep[k]);
        b[s]=bs;
    }
    rvec_t w; rmat_t V; blazem::eigen(Ss, w, V);
    rvec_t r(ns, 0.0);
    for (size_t k=0;k<ns;k++)
    {
        if (w[k]<1e-10*w[ns-1]) continue;                       // a near-dependent direction of the site's shells: skipped
        double vb=0.0; for (size_t s=0;s<ns;s++) vb+=V(s,k)*b[s];
        for (size_t s=0;s<ns;s++) r[s]+=V(s,k)*vb/w[k];
    }
    double captured=0.0; for (size_t s=0;s<ns;s++) captured+=r[s]*b[s];   // |P chi|^2 = b^T S^-1 b
    std::ostringstream os;
    os<<"[+U radial] site "<<site<<" ("<<element<<", q"<<Zion<<") l="<<l<<": pseudo-atom "<<(l==0?"s":l==1?"p":l==2?"d":"f")
      <<" orbital eps="<<std::setprecision(4)<<eBest<<" Ha (16-exponent pool, E="<<std::setprecision(6)<<atom.Energy()
      <<") projected onto "<<ns<<" shells: captured "<<std::fixed<<std::setprecision(4)<<captured<<" of its norm; coefficients";
    for (size_t k=0;k<ns;k++) os<<" "<<std::setprecision(3)<<r[k]<<"@"<<es[k];
    std::cout<<os.str()<<std::endl;
    std::vector<double> out(ns); for (size_t s=0;s<ns;s++) out[s]=r[s];
    return out;
}

std::vector<qchem::Hamiltonian::HubbardEstimate> SolidCalculation::EstimateHubbardU() const
{
    std::unique_ptr<qchem::Hamiltonian::HubbardUEstimator> est=itsImp->ham->MakeHubbardUEstimator();
    if (!est) throw std::logic_error("SolidCalculation::EstimateHubbardU: this run carries no Hubbard manifold (SolidCalcOptions::hubbard)");
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!wf) throw std::logic_error("SolidCalculation::EstimateHubbardU: no wave function yet");
    // Block by block: the occupied orbitals' AO coefficients and physical occupations, the block's BZ
    // weight, and its spin (Spin::None = the folded doublet; the estimator halves it into both channels).
    auto feed=[&]<class U>(const qchem::Orbitals::TOrbitals<U>& os, const qchem::Irrep& ir)
    {
        const auto* blk=dynamic_cast<const qchem::BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
        if (!blk) throw std::logic_error("SolidCalculation::EstimateHubbardU: an orbital block that is not an Orbital_DFT_IBS");
        std::vector<const qchem::Orbitals::TOrbital<U>*> occ;
        for (const auto* o : os.template Iterate<qchem::Orbitals::TOrbital<U>>()) if (o->IsOccupied()) occ.push_back(o);
        mat_t<U> C(blk->GetNumFunctions(), occ.size());
        rvec_t f(occ.size());
        for (size_t i=0;i<occ.size();i++)
        {
            const vec_t<U>& c=occ[i]->GetCoeff();
            for (size_t a=0;a<c.size();a++) C(a,i)=c[a];
            f[i]=occ[i]->GetOccupation();
        }
        est->Accumulate(*blk, ir.ms, ir.sym->GetWeight(), C, f);
    };
    for (const qchem::Irrep& ir : wf->GetQNs())
    {
        const qchem::Orbitals::Orbitals* os=wf->GetOrbitals(ir);
        if (const auto* c=dynamic_cast<const qchem::Orbitals::TOrbitals<dcmplx>*>(os)) feed(*c, ir);
        else if (const auto* r=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(os)) feed(*r, ir);
        else throw std::logic_error("SolidCalculation::EstimateHubbardU: an orbital set of unknown scalar");
    }
    std::vector<qchem::Hamiltonian::HubbardEstimate> out=est->Evaluate();
    est->Write(std::cout);                                   // the estimator reports at its own activity (pin 17)
    return out;
}

SolidCalculation::HubbardLoopResult SolidCalculation::ConvergeHubbardU(const SCFParams& params, const HubbardLoop& loop)
{
    const double eV=27.211386245988;
    HubbardLoopResult R;
    std::unique_ptr<qchem::Hamiltonian::HubbardUEstimator> est=itsImp->ham->MakeHubbardUEstimator();
    if (!est) throw std::logic_error("SolidCalculation::ConvergeHubbardU: this run carries no Hubbard manifold");
    std::vector<double> Ucur;                                 // the U each manifold currently runs with (eV)
    for (const auto& M : itsImp->opts.hubbard) Ucur.push_back(M.U*eV);
    R.scfConverged=itsImp->converged;
    for (size_t n=0; n<loop.maxOuter; n++)
    {
        R.last=EstimateHubbardU();                            // feeds the estimator from the current orbitals, prints [ACBN0]
        R.outer=n+1;
        std::vector<double> Unext; double dmax=0.0;
        for (size_t M=0;M<R.last.size();M++) { Unext.push_back(R.last[M].Ueff()*eV); dmax=std::max(dmax, std::abs(Unext[M]-Ucur[M])); }
        R.U_eV.push_back(Unext);
        {
            std::ostringstream os;
            os<<"[ACBN0 loop "<<n<<"] U_eff(eV):"; for (double u : Unext) os<<" "<<std::fixed<<std::setprecision(4)<<u;
            os<<"  max|dU|="<<std::setprecision(5)<<dmax<<(itsImp->converged ? "" : "  (SCF NOT converged)");
            std::cout<<os.str()<<std::endl;
        }
        if (dmax<loop.tolU_eV) { R.converged=true; break; }
        if (n+1==loop.maxOuter) break;
        // Apply on the same Hamiltonian and re-converge from the current density: the estimator writes the
        // term's U; the facade re-runs the SCF with the caller's parameters.  (est was fed from these orbitals;
        // Apply resets it, and the next EstimateHubbardU builds a fresh one.)
        est->Apply(R.last);
        est=itsImp->ham->MakeHubbardUEstimator();
        Ucur=Unext;
        // A NEW STAGE, not a continued Iterate: the converged mixer history (eight near-zero Pulay residuals)
        // would extrapolate the first post-U step straight back onto the old density and report "converged"
        // in one iteration (measured 2026-09-21: E moved by exactly E_U, the orbitals not at all).  Re-seeding
        // through BuildStage gives a fresh iterator, mixer and accelerator over the current density -- exactly
        // what an annealed schedule does between its stages.
        BuildStage(itsImp->stageAccel, std::move(itsImp->cd));
        auto out=Converge(params);
        R.scfConverged=bool(out);
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
        {.acrossK=itsImp->opts.globalFermi, .acrossSpin=itsImp->opts.spinsShareFermi}, noise);
    qchem::Response::AmplitudeProbe probe=qchem::Response::MakeHubbardProbe(ref, *wf, *hub);
    auto r=qchem::Response::IndependentResponse(ref, probe, Nq);
    if (r) r->Write(std::cout, 27.211386245988, "1/eV");     // reported at its own activity (pin 17)
    else   std::cout << "[chi0] FAILED: " << r.Error().detail << std::endl;
    return r;
}

// THE SELF-CONSISTENT RESPONSE AT q = 0 (R2): the same Reference and probe as IndependentResponse, plus the
// AO <-> MO frame and the Hamiltonian's analytic kernel, linearised about the converged ORBITALS' density (the
// state the Reference describes -- not the mixed iterate, which on a Kerker/Pulay recipe has no D).
Outcome<qchem::Response::SelfConsistentResponse,qchem::Response::ResponseFailure>
SolidCalculation::HubbardLinearResponse(const KrylovParams& krylov) const
{
    using O=Outcome<qchem::Response::SelfConsistentResponse,qchem::Response::ResponseFailure>;
    const auto* hub=itsImp->ham->GetHubbardChannels();
    if (!hub) throw std::logic_error("SolidCalculation::HubbardLinearResponse: this run carries no Hubbard manifold -- the "
                                     "response channels ARE the +U manifolds (list them at U=0 to probe without +U)");
    const auto* wf=itsImp->scf->GetWaveFunction();
    if (!wf) throw std::logic_error("SolidCalculation::HubbardLinearResponse: no wave function yet");
    const double noise=std::isfinite(itsImp->lastCommutator) ? std::fabs(itsImp->lastCommutator)
                                                             : std::numeric_limits<double>::quiet_NaN();
    const qchem::Response::Reference ref=qchem::Response::MakeReference(*wf, itsImp->lastOccupation,
        {.acrossK=itsImp->opts.globalFermi, .acrossSpin=itsImp->opts.spinsShareFermi}, noise);
    const auto frame=qchem::Response::MakeOrbitalFrame(ref, *wf);
    const qchem::Response::AmplitudeProbe probe=qchem::Response::MakeHubbardProbe(ref, *wf, *hub);
    const auto D0=wf->GetChargeDensity();
    const auto kernel=itsImp->ham->MakeResponseKernel(itsImp->bs.get(), D0.get());
    auto qs=ref.QMesh(ivec3_t(1,1,1));
    if (!qs) return O::Fail(qs.Error());
    auto r=qchem::Response::LinearResponse(ref, frame, *kernel, probe,
                                           std::make_shared<qchem::Response::MeshShift>((*qs)[0]), krylov);
    if (!r)
    {
        std::cout << "[response] FAILED: " << r.Error().detail << std::endl;
        return r;
    }
    r->Write(std::cout);
    // U_I = (chi0^-1 - chi^-1)_II, in eV (hp.x's unit).  Printed, never consumed here (R4 owns consumption).
    const mat_t<dcmplx> X0i=blazem::inv(r->chi0), Xi=blazem::inv(r->chi);
    std::cout << "[response] U = (chi0^-1 - chi^-1)_II at q = 0 (this cell):";
    for (size_t I=0;I<r->labels.size();I++)
        std::cout << "  " << r->labels[I] << " " << std::setprecision(6) << (X0i(I,I)-Xi(I,I)).real()*27.211386245988 << " eV";
    std::cout << std::endl;
    return r;
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
void   SolidCalculation::OnIteration(Observer obs)      { itsImp->opts.onIteration=std::move(obs); AttachProbes(); }
bool   SolidCalculation::DidConverge()    const         { return itsImp->converged; }
size_t SolidCalculation::IterationCount() const         { return itsImp->scf->GetIterationCount(); }

const RunDiagnostics& SolidCalculation::Diagnostics() const { return itsImp->diag; }

const qcMesh::MeshParams& SolidCalculation::ResolvedXCMesh() const { return itsImp->xcMesh; }
const BasisSet::Complex_BS& SolidCalculation::Basis() const { return *itsImp->bs; }

} //namespace qchem
