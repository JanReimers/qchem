// File: Calculation/Imp/SolidCalculation.C  The periodic SCF front door.
//
// Imports ZERO .Internal. modules -- the same discipline Imp/Calculation.C and Imp/AtomCalculation.C
// already keep, and the reason Step 4 had to open two public doors before this file could exist.
module;
#include <algorithm>   // std::max (the T2 site-moment scan)
#include <cmath>       // std::fabs
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

import qchem.BasisSet.Molecule.LatticeSum1E;   // LatticeSum1E::MaxExponent -- alpha_max (a capability face)
import qchem.BasisSet.Orbital_1E_IBS;          // Real_OIBS / Complex_OIBS (the per-irrep bases to iterate)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH local PP -> alpha_pp
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.Mesh.XCPolicy;                    // XCMeshSharpness / ResolveXCMesh (the grid decision)
import qchem.ElectronConfiguration.Crystal;    // Crystal_EC
import qchem.Symmetry.Irrep;                   // Irrep, Spin
import qchem.RunPolicy;                        // the declared CP2K deviations (doc/OpenWork.md N5/T5)
import qchem.Parallel;                         // WorkerThreads() -- half of the thread state a row must state

namespace qchem
{

//---------------------------------------------------------------------------------------------------
//  The run's SHARPNESS -- an above-SCFIterator decision input.
//
//  Both sources come off ABSTRACT capability faces via the sanctioned abstract->abstract cross-cast, so
//  nothing here touches a concrete basis or a concrete PP model:
//    alpha_max -- BasisSet::Molecule::LatticeSum1E::MaxExponent(), documented there as "the GPW
//                 density-grid cutoff floor".
//    alpha_pp  -- Pseudopotential::LocalPotential_Gaussian::ShortRangeGaussian(Z), whose terms carry
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
        if (const auto* ls=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(ibs))
            { s.alphaMax = ls->MaxExponent(); break; }
    for (const auto& [element, valence] : o.species)
    {
        const int Z = int(thePeriodicTable().GetZ(element));
        const Pseudopotential::HGH_LocalPotential loc = Pseudopotential::GetGTH(element,"LDA",valence).local;
        const auto* g = static_cast<const Pseudopotential::LocalPotential_Gaussian*>(&loc);
        for (const auto& t : g->ShortRangeGaussian(Z)) s.alphaPP = std::max(s.alphaPP, t.alpha);
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
    qchem::Hamiltonian::cHamiltonian*           ham = nullptr;   // owned by the iterator once handed over
    qcMesh::MeshParams                          xcMesh;          // AFTER Auto resolution
    std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scf;
    std::unique_ptr<qchem::ChargeDensity::cDM_CD>         cd;    // the converged density (outlives the WF)
    //! m(r) of the converged state.  OWNED: WaveFunction::GetSpinDensity() hands back a raw `new`, so it
    //! goes straight into a unique_ptr here (CLAUDE.md ownership) -- null on an unpolarized run.
    std::unique_ptr<SolidCalculation::sf_t>               spin;
    SCFAccelerators::SolidAcceleratorOptions    accOpts;
    SCFAccelerators::Type                       stageAccel = SCFAccelerators::Type::DIIS;  //!< the CURRENT stage's, for the banner
    bool   converged = false;
    bool   imposed   = false;   //!< the run imposed a symmetry, so its solution must still carry it (T2)
    double charge    = 0.0;
    //! The OUTCOME DETECTORS' record (N1/T3-T4).  Accumulates across every attempt this object makes --
    //! an annealed schedule and a re-Converge are CONTINUATIONS of one SCF, and an order that died in
    //! stage 1 must not become invisible because stage 2 started from the corpse.
    RunDiagnostics diag;
    //! The integrated site moment of the CURRENT iterate, left here by the order probe for the observer
    //! to file.  Two hooks, one row: the probe is the only place the density is in hand, the observer is
    //! the only place that fires exactly once per iteration.
    double lastOrder = 0.0;
};

//---------------------------------------------------------------------------------------------------
//  The INTEGRATED order parameter, off the Hamiltonian's own atom-centred partition.  Empty basins (a
//  uniform XC mesh, or an unpolarized run) give an empty vector, which the caller reads as "not
//  measurable" -- never as "zero", because those are different facts and only one of them is a failure.
static double MaxSiteMoment(const qchem::Hamiltonian::cHamiltonian& ham,
                            const qchem::ChargeDensity::cChargeDensity& cd, bool& hasBasins)
{
    const rvec_t m = ham.SiteMoments(&cd);
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
                          bool imposed)
{
    const char* omp=std::getenv("OMP_NUM_THREADS");
    std::cout<<"["<<o.label<<" run] system: "<<nAtoms<<" atoms, "<<o.Nelec<<" valence e, multiplicity "
             <<o.multiplicity<<(o.multiplicity>=1 ? " (POLARIZED)" : " (unpolarized)")
             <<", seed="<<SeedName(o.seed)<<std::endl;
    std::cout<<"["<<o.label<<" run] grids: densityEcut="
             <<(o.densityEcut<0 ? std::string("auto") : std::to_string(o.densityEcut))
             <<" C="<<o.cutoffFactor<<" raster="<<(o.raster==BasisSet::Lattice_3D::RasterPolicy::BallOnly
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
             <<" GPW_OMP_THREADS="<<qchem::WorkerThreads()<<" (BLAS pinned to 1)"<<std::endl;
    std::cout<<"["<<o.label<<" run] "<<theRunPolicy().Banner()<<std::endl;
}
// WHAT THE SCF IS DOING: the mixer, the accelerator, the occupation machinery.  Emitted per Converge,
// because that is where these take effect -- and an anneal changes them stage by stage.
static void EmitSCFBanner(const std::string& label, const SCFParams& p, SCFAccelerators::Type acc)
{
    // ⚠ THE PRECONDITIONER AND THE HISTORY ARE INDEPENDENT, AND THE BANNER USED TO HIDE IT (2026-08-28).
    // It printed "Pulay(depth 8, start 5)" whenever PulayDepth>0 and never mentioned Kerker, so a reader
    // of a benchmark row would conclude the G-space preconditioner had been swapped OUT.  It has not:
    // MakeGSpaceMixer passes kerkerG0 straight into PulayMixer, so the two COMPOSE -- which matters
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
              << " A=E-TS="<<E.GetTotalEnergy()<<" -TS="<<E.MinusTS
              << " E(internal)="<<(E.GetTotalEnergy()-E.MinusTS)
              << std::setprecision(prec0) << std::endl;
}

//---------------------------------------------------------------------------------------------------
SolidCalculation::SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                   const SolidCalcOptions& opts, const SCFParams& params,
                                   const SCFAccelerators::SolidAcceleratorOptions& acc)
    : itsImp(std::make_unique<Imp>())
{
    itsImp->opts    = opts;
    itsImp->accOpts = acc;
    itsImp->st      = lat.GetStructure();

    namespace L3 = BasisSet::Lattice_3D;
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
    itsImp->bs.reset(L3::GPWFactory(lat, mol, L3::GPWParams{
        .densityEcut = opts.densityEcut, .cutoffFactor = opts.cutoffFactor, .raster = opts.raster,
        .images = opts.images, .kShift = opts.kShift, .ladderFactor = opts.ladderFactor,
        .imposeSymmetry = imposed, .siteSpins = siteSpins,
        .hamPreservesReal = hamPreservesReal}));

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

    itsImp->ham = qchem::Hamiltonian::Factory(
        polarized ? qchem::Hamiltonian::Pol::Polarized : qchem::Hamiltonian::Pol::UnPolarized,
        itsImp->st, itsImp->bs.get(), opts.species, "LDA", itsImp->xcMesh, opts.vxcFit);
    // The forecast crosscheck: the basis was built on the promise that every term preserves realness
    // (the AND's term half, above); the constructed Hamiltonian must agree, or real blocks were built
    // that its terms cannot serve.
    assert((!hamPreservesReal || itsImp->ham->PreservesReal()) &&
           "SolidCalculation: the term stack no longer preserves realness -- update the forecast above");

    // DECISION 3 -- the accelerator, by policy, through the public typed door.
    itsImp->stageAccel = opts.accelerator;
    auto* accel = SCFAccelerators::Factory(opts.accelerator, acc);

    // THE SEED IS BUILT HERE, not inside the iterator -- the same factory call the SeedStrategy ctor
    // would have made (ChargeDensity::MakeSeedDensity with the Hamiltonian's own polarization), handed
    // to the explicit-seed ctor instead.  WHY THE FACADE TAKES IT OVER (N1/T2, 2026-08-26): the
    // postcondition needs the order the run STARTED with, and the iterator consumes the seed inside
    // Init -- which builds a Fock from it, diagonalizes and FILLS, so the earliest density a caller can
    // reach is already one aufbau fill downstream.  For a run whose order dies in that very first fill
    // the difference is the whole measurement: a Na2 seed staggered at +/-1 e reads +/-0.07 e by the
    // time Init hands its density back, which is below any honest floor and made the postcondition
    // silently skip.  Measured before it is consumed, the baseline is the seed's own.
    const bool polarizedHam = itsImp->ham->IsPolarized();
    std::unique_ptr<qchem::ChargeDensity::cChargeDensity> seed(
        qchem::ChargeDensity::MakeSeedDensity<dcmplx>(opts.seed, itsImp->bs.get(), itsImp->st.get(),
                                                      itsImp->ec.get(), polarizedHam));
    // MEASURED ONLY WHERE IT CAN MEAN SOMETHING.  SiteMoments rasters BOTH spin channels before it can
    // discover it has no basins to integrate over, and unlike the per-iteration probe -- which rides a
    // raster the Fock build has already made for this density serial -- nothing has been built yet at
    // seed time, so this one would be paid in full.  An unpolarized run has m == 0 identically and a
    // uniform XC mesh has no basins at all, and the facade knows both facts here without asking.
    if (seed && polarizedHam && itsImp->xcMesh.cellKind==qcMesh::UnitCellKind::Becke)
        itsImp->diag.itsSeedOrder = MaxSiteMoment(*itsImp->ham, *seed, itsImp->diag.itsHasBasins);

    itsImp->scf = std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
        itsImp->bs.get(), itsImp->ec.get(), itsImp->ham, accel,
        seed.release(), itsImp->st.get(), opts.ortho, opts.orthoTol);   // the iterator consumes it in Init

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

    EmitRunBanner(opts, itsImp->xcMesh, itsImp->st->GetNumAtoms(), imposed);
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
    const bool polarized = itsImp->opts.multiplicity>=1;
    itsImp->ham = qchem::Hamiltonian::Factory(
        polarized ? qchem::Hamiltonian::Pol::Polarized : qchem::Hamiltonian::Pol::UnPolarized,
        itsImp->st, itsImp->bs.get(), itsImp->opts.species, "LDA", itsImp->xcMesh, itsImp->opts.vxcFit);
    itsImp->stageAccel = accType;
    auto* accel = SCFAccelerators::Factory(accType, itsImp->accOpts);

    auto prev = std::move(itsImp->scf);      // held ONLY until the new stage has copied its MOM reference
    itsImp->scf = carried
        ? std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham, accel,
              carried.release(), itsImp->st.get(), itsImp->opts.ortho, itsImp->opts.orthoTol)
        : std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham, accel,
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
//! The iterator has exactly one order-probe slot and one observer slot, and the facade needs both for
//! itself while still honouring whatever the caller asked for -- so it COMPOSES rather than competes.
//! The probe is the only place this iteration's density is in hand (that is where the integrated site
//! moment can be taken, for free, off a raster the XC term has already built for this density serial);
//! the observer is the only thing that fires exactly once per iteration (the probe can also be called
//! for the iteration-0 banner under Verbose).  So the probe MEASURES and the observer FILES.
//!
//! The caller's probe still owns the DISPLAY column: their order parameter is the campaign's question,
//! and the library has no business overwriting it.  When they set none and the run is polarized, the
//! integrated moment takes the column itself under the name "m_site" -- which is what a magnetic run
//! wants to see anyway, and it is an INTEGRATED observable rather than a point sample of m(r).
void SolidCalculation::AttachProbes()
{
    const bool polarized = itsImp->opts.multiplicity>=1;
    auto userProbe = itsImp->opts.orderProbe;
    if (userProbe || polarized)
        itsImp->scf->SetOrderParameter(
            userProbe ? itsImp->opts.orderName : std::string("m_site"),
            [this,userProbe](const qchem::ChargeDensity::cDM_CD& cd)->double
            {
                itsImp->lastOrder = MaxSiteMoment(*itsImp->ham, cd, itsImp->diag.itsHasBasins);
                return userProbe ? userProbe(cd) : itsImp->lastOrder;
            });
    auto userObs = itsImp->opts.onIteration;
    itsImp->scf->SetObserver([this,userObs](const qchem::SCFIterator::SCFProgress& p)
    {
        itsImp->diag.itsOrder.push_back(itsImp->lastOrder);
        itsImp->diag.itsEee  .push_back(p.eb.Eee);
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
    EmitSCFBanner(itsImp->opts.label, params, itsImp->stageAccel);
    itsImp->scf->Iterate(params);
    itsImp->converged = itsImp->scf->Converged();
    itsImp->diag.itsConverged = itsImp->converged;
    // Take the density OUT of the wave function: it stays valid after the iterator's state moves on
    // because its basis block is `bs`, which this object owns and outlives it.
    auto cd = itsImp->scf->GetWaveFunction()->GetChargeDensity();   // BUILT for us; we take it
    itsImp->charge = cd->GetTotalCharge();
    itsImp->cd = std::move(cd);
    // m(r) the same way: BUILT for us (a raw `new`), so it is adopted here rather than leaked.
    itsImp->spin.reset(itsImp->scf->GetWaveFunction()->GetSpinDensity());
    return Outcome_();
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Result() const {return Outcome_();}

qchem::EnergyBreakdown SolidCalculation::LastIterateTerms()  const {return itsImp->scf->GetEnergy();}
double                 SolidCalculation::LastIterateCharge() const {return itsImp->charge;}

// The caller's observer is SWAPPED IN behind the facade's own (AttachProbes composes the two), so
// attaching telemetry late cannot silently disarm the outcome detectors.
void   SolidCalculation::OnIteration(Observer obs)      { itsImp->opts.onIteration=std::move(obs); AttachProbes(); }
bool   SolidCalculation::DidConverge()    const         { return itsImp->converged; }
size_t SolidCalculation::IterationCount() const         { return itsImp->scf->GetIterationCount(); }

const RunDiagnostics& SolidCalculation::Diagnostics() const { return itsImp->diag; }

const qcMesh::MeshParams& SolidCalculation::ResolvedXCMesh() const { return itsImp->xcMesh; }
const BasisSet::Complex_BS& SolidCalculation::Basis() const { return *itsImp->bs; }

} //namespace qchem
