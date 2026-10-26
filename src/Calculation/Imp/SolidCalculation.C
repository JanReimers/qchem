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

import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;   // GaussianSharpness::MaxExponent -- alpha_max (the NARROW face)
import qchem.BasisSet.Orbital_1E_IBS;          // Real_OIBS / Complex_OIBS (the per-irrep bases to iterate)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH local PP -> alpha_pp
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.Mesh.XCPolicy;                    // XCMeshSharpness / ResolveXCMesh (the grid decision)
import qchem.ElectronConfiguration.Crystal;    // Crystal_EC
import qchem.Symmetry.Irrep;                   // Irrep, Spin
import qchem.Blaze;                            // blazem::conjs (the u_bands weights)
import qchem.PseudoAtomOrbital;                // the pseudo-atom l orbital: the +U projector (HubbardProjectorPlan D5)
import qchem.Orbitals;                         // TOrbitals/TOrbital -- the ACBN0 feed (EstimateHubbardU)
import qchem.Diagnostics;                      // u_bands: the per-band +U slot weights (QCHEM_DIAGNOSTICS=u_bands[=N])
import qchem.BasisSet.Orbital_DFT_IBS;         // the block an orbital set lives on (EstimateHubbardU)
import qchem.Environment;                   // WarnRetiredEnvironment (D-ENV step 6a)
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
    //! m(r), built on the FIRST \c Converged::SpinDensity() call (reset by each Converge): built eagerly it re-collocated
    //! both spin channels' D at the end of every stage, two n×n copies beside the iterator's own state -- the run's peak.
    mutable std::unique_ptr<SolidCalculation::sf_t>       spin;
    SCFAccelerators::SolidAcceleratorOptions    accOpts;
    SCFAccelerators::Type                       stageAccel = SCFAccelerators::Type::DIIS;  //!< the CURRENT stage's, RESOLVED (banner, [F,D])
    SCFAccelerators::Type                       stageAccelAsked = SCFAccelerators::Type::Auto;  //!< ...as ASKED: what an FD/U-loop stage re-resolves
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
    //! The +U manifolds AS THE HAMILTONIAN CARRIES THEM: resolved by BuildRun, and kept in step with every
    //! ACBN0 outer step's Apply -- so a saved state records the U its density was converged at, not the
    //! options' starting U.
    std::vector<qchem::Hamiltonian::HubbardManifold> hubbard;
    //! The last Converge's recipe (kT, and the banner facts, for the saved state's summary).
    double      lastkT = 0.0;
    std::string lastRecipe;
    rmat3d_t    cellMatrix;   //!< the lattice's A, for the fingerprint (Structure alone does not carry it)
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
    std::cout<<"["<<label<<" scf] mixer: "
             <<(p.Mixer==M::Linear ? std::string("LINEAR D-mixing (no G-space preconditioner)")
              : p.KerkerG0>0.0     ? "Kerker(G0="+std::to_string(p.KerkerG0)+")"
              :                      std::string("rho-mixing (no G-space preconditioner)"))
             <<(p.Mixer==M::Pulay ? " + Pulay history(depth "+std::to_string(p.PulayDepth)+", start "
                                    +std::to_string(p.PulayStart)+")"
                                  : std::string())
             <<(autoMixer ? " [mixer: Auto]" : "")
             <<" alpha="<<p.StartingRelaxRo
             <<";  XC rho source: "<<(theRunPolicy().XCFromDM() ? "rho[D] WHOLESALE"
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
//  THE MAGNETIC ORDERING IS STATED (OpenWork §3 "The magnetic decoration is DERIVED FROM THE SEED", 2026-10-26).
//  The message names the sites, not a derived guess: the ordering is the caller's chemistry.
void RequireStatedMagneticOrdering(const SolidCalcOptions& o, const Structure& st)
{
    std::string sites;
    st.ForEachSite([&](int Z, const rvec3_t&, bool){ sites+=(sites.empty()?"":" ")+thePeriodicTable().GetSymbol(Z); });
    const size_t nAtoms=st.GetNumAtoms();
    if (!o.siteSpins.empty())
    {
        if (o.siteSpins.size()!=nAtoms)
            throw std::invalid_argument("SolidCalcOptions: siteSpins has "+std::to_string(o.siteSpins.size())
                                        +" entries for "+std::to_string(nAtoms)+" sites ("+sites+")");
        for (int s : o.siteSpins)
            if (s<-1 || s>1) throw std::invalid_argument("SolidCalcOptions: siteSpins entries are -1, 0 or +1, got "+std::to_string(s));
        return;
    }
    const bool polarized = o.multiplicity>=1;
    if (polarized && !o.hubbard.empty() && !o.greyImposition)
        throw std::invalid_argument("SolidCalcOptions: a polarized run with a Hubbard manifold must STATE its magnetic ordering --"
                                    " siteSpins, one of +1/-1/0 per site ("+sites+"); all zeros for a non-magnetic cell."
                                    "  It is no longer derived from the seed (doc/OpenWork.md §3)");
}

//---------------------------------------------------------------------------------------------------
static std::vector<Math::Gaussian> PseudoAtomicFunction(const Structure& st, size_t site, int l,
                                                        const std::vector<std::pair<std::string,int>>& species);

SolidCalculation::SolidCalculation(BuildOnly) : itsImp(std::make_unique<Imp>()) {}

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
    const std::vector<int> siteSpins=BuildBasis(lat, mol, opts, acc);
    BuildRun(lat, *mol, siteSpins, nullptr, params);
    (void)Converge(params);   // the ctor ATTEMPTS; the caller faces the result via Result()
}

// THE FIRST HALF OF THE BUILD: everything the saved-state fingerprint is judged against (the Bloch blocks), and
// nothing that depends on the seed.  Restart stops here to compare before paying for the Hamiltonian.
std::vector<int> SolidCalculation::BuildBasis(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                              const SolidCalcOptions& opts,
                                              const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    qchem::SetRunPolicy(opts.policy);   // D-ENV 6b: the policy is part of the run's options, installed ONCE here, before any factory consults it
    itsImp->opts    = opts;
    itsImp->accOpts = acc;
    itsImp->st      = lat.GetStructure();
    RequireStatedMagneticOrdering(opts, *itsImp->st);
    itsImp->cellMatrix = lat.GetUnitCell().GetCellMatrix();

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
    // order is erased.  STATED (siteSpins) wherever a Hubbard manifold rides along -- RequireStatedMagneticOrdering
    // above refused the run otherwise.  The derivation from the seed below survives ONLY for an imposed polarized
    // run WITHOUT +U, and it is the defect OpenWork §3 names: with SeedStrategy::SAD it labels neutral O (2p4)
    // magnetic and MnO AFM-II imposes a ferrimagnet (12 Shubnikov ops, not 24).
    //   unpolarized  -> {} (no channels to decorate)
    //   greyImposition -> {} DELIBERATELY: that arm is the erasure NEGATIVE CONTROL
    std::vector<int> siteSpins;
    if (!imposed) { /* nothing to decorate: siteSpins are meaningless without an imposition */ }
    else if (!opts.siteSpins.empty() && !opts.greyImposition)
        siteSpins = opts.siteSpins;                 // STATED: a specific ordering, not the seed's guess
    else if (opts.multiplicity>=1 && !opts.greyImposition)
    {
        const std::vector<int> targets = (opts.seed==qchem::ChargeDensity::SeedStrategy::IonicSAD)
                                       ? qchem::ChargeDensity::IonicSADTargets(itsImp->st.get(), "LDA")
                                       : std::vector<int>{};
        siteSpins = qchem::ChargeDensity::MagneticDecoration(itsImp->st.get(), "LDA", targets);
    }
    // ★ THE FACADE'S SETUP IS BUCKETED (doc/ParallelAndOraclePlan.md 1.1, 2026-09-06).  Every MnO row drives
    // SolidCalculation directly, and until now only the phases INSIDE these calls were timed -- so the basis
    // build, the Hamiltonian ctor's non-mesh half, the seed and the ortho landed in no bucket at all.  That
    // was ~48 s of a 119 s threaded MnO run, i.e. the largest non-scaling block in it, and unnamed.  Labels
    // match the GPW_SCF harness's so the two paths' ledgers read the same.
    const bool planeWave = opts.orbitalBasis==SolidCalcOptions::OrbitalBasis::PlaneWave;
    if (planeWave)
    {
        // THE PLANE-WAVE FAMILY (OpenWork §2 PW-1).  The factory builds ONE Γ block and knows no symmetry, so
        // anything else is refused here rather than silently run as Γ / unimposed.
        if (!(opts.planeWaveEcut>0.0))
            throw std::invalid_argument("SolidCalcOptions: orbitalBasis PlaneWave needs planeWaveEcut > 0 (Ha)");
        if (lat.MakeKMesh(opts.kShift).size()!=1 || opts.kShift!=rvec3_t(0,0,0))
            throw std::invalid_argument("SolidCalcOptions: the plane-wave family is Γ-only so far (k-mesh 1x1x1, no shift)");
        if (imposed)
            throw std::invalid_argument("SolidCalcOptions: the plane-wave family has no symmetry imposition yet");
        // Plane waves ARE orthonormal: the default Cholesky would factor an identity S into two dense identity
        // factors per spin block.  A stated non-default ortho is obeyed (an A/B door).
        if (itsImp->opts.ortho==qchem::Cholesky) itsImp->opts.ortho=qchem::Orthonormal;
        qchem::report::Timed timed("setup: PW basis build");
        itsImp->bs.reset(L3::Factory(L3::Type::PW, lat, opts.planeWaveEcut));
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
        itsImp->bs.reset(L3::GPWFactory(lat, mol, L3::GPWParams{
            .densityEcut = opts.densityEcut, .cutoffFactor = opts.cutoffFactor, .raster = opts.raster,
            .images = opts.images, .kShift = opts.kShift, .ladderFactor = opts.ladderFactor,
            .imposeSymmetry = imposed, .siteSpins = siteSpins, .tol = tol,
            .hamPreservesReal = hamPreservesReal}));
    }
    // THE RUN REPORTS ITS OWN BASIS AND GRIDS (TE phase 2, 2026-09-15; the reporting rule: each class emits
    // contemporaneously with its own activity).  These two sections used to be emitted by the integration
    // test driver around this same construction; the facade is the orchestrator now.  Order matters and is
    // the one the driver established: the conditioning pre-flight on the ANALYTIC overlap first (no grids
    // needed), then the grid ladder (EmitGpwGrids is what forces EnsureLevels).  Only when a run report is
    // open -- the emitters are inert otherwise -- and never an abort: a dependent basis is the caller's to
    // handle through a pivoted ortho, and an unhandled one fails loudly in the iterator.
    // (Plane waves are orthonormal and grid-free: neither section has anything to say about them.)
    if (planeWave) {}
    else if (qchem::report::Depth()>0)
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
    else L3::PrintGpwConditioning(*itsImp->bs);   // no report open: still say how conditioned S(k) is

    // DECISION 1 -- the XC quadrature.  Resolve Auto HERE, once, from facts about the run.  Downstream
    // consumers compare ==Becke, so an unresolved Auto would silently read as Uniform; resolving it at the
    // point the spec enters the Hamiltonian is what makes that impossible.
    // The XC-mesh route is a DECLARED CP2K deviation (RunPolicy::BeckeXC), resolved here beside the
    // imposition for the same reason: the policy belongs to the run, the sizing belongs to qcMesh.
    // A plane-wave run's Auto is its own FFT raster (QE's route): the sharpness gather prices Gaussians, and a
    // delta fit would collocate every PW at every mesh point.  A STATED grid is obeyed, as for GPW.
    if (planeWave && opts.xcMesh.cellKind==qcMesh::UnitCellKind::Auto)
        itsImp->xcMesh = qcMesh::MeshParams{.cellKind=qcMesh::UnitCellKind::Uniform};
    else
        itsImp->xcMesh = qcMesh::ResolveXCMesh(opts.xcMesh, GatherSharpness(lat, *mol, opts, imposed),
                                               theRunPolicy().BeckeXC());

    // DECISION 2 -- the spin bookkeeping.  multiplicity -> (nUp,nDown), with the parity check that catches
    // a singlet asked of an odd electron count BEFORE integer division silently empties a channel.
    const int twoS = opts.multiplicity>1 ? opts.multiplicity-1 : opts.Nelec%2;
    const bool polarized = opts.multiplicity>=1;
    if ((opts.Nelec-twoS)%2!=0 || twoS>opts.Nelec)
        throw std::runtime_error("SolidCalcOptions: multiplicity "+std::to_string(opts.multiplicity)
                                 +" parity disagrees with Nelec "+std::to_string(opts.Nelec));
    // A partial eigensolve must leave EMPTY bands above the fuller channel (smearing and the gap need them).
    if (opts.numBands>0 && opts.numBands<=size_t((opts.Nelec+twoS)/2))
        throw std::invalid_argument("SolidCalcOptions: numBands "+std::to_string(opts.numBands)+" must exceed the larger spin channel's "
                                    +std::to_string((opts.Nelec+twoS)/2)+" electrons");
    auto irreps = itsImp->bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry Sum_k)
    itsImp->ec = std::make_unique<Crystal_EC>(irreps, (opts.Nelec+twoS)/2, (opts.Nelec-twoS)/2,
                                              !opts.perKFermiLevels, opts.spinsShareFermi);
    // T2 is gated on the run having IMPOSED something (see Outcome_); the resolved value, stored once.
    itsImp->imposed = imposed;
    itsImp->hubbard = opts.hubbard;   // the U the run will carry (BuildRun fills the groups and radials)
    (void)polarized;
    return siteSpins;
}

// THE SECOND HALF: the Hamiltonian, the seed (the caller's when given -- a restored state), the iterator.
void SolidCalculation::BuildRun(const Lattice_3D& lat, const BasisSet::Real_BS& mol, const std::vector<int>& siteSpins,
                                std::unique_ptr<qchem::ChargeDensity::cChargeDensity> seed, const SCFParams& params)
{
    const SolidCalcOptions& opts=itsImp->opts;
    const bool polarized        = opts.multiplicity>=1;
    const bool hamPreservesReal = !opts.forceComplex;   // the BuildBasis forecast, restated for the assert below
    const bool imposed          = itsImp->imposed;

    // THE HUBBARD MANIFOLDS' SITE GROUPS (step 5 increment 2): each manifold is labelled by the stabiliser of
    // its site in the STATED DECORATION's Shubnikov group (sigma=None ops), grey parentage beside it -- the
    // same decoration the imposition uses, whether or not the run imposes: the site symmetry of the ordered
    // state is the physical question either way.  No longer derived from the seed (OpenWork §3, 2026-10-26):
    // a polarized +U run was refused at BuildBasis without siteSpins.  A caller who filled siteOps by hand keeps them.
    std::vector<qchem::Hamiltonian::HubbardManifold> hubbard=opts.hubbard;
    if (!hubbard.empty())
    {
        std::vector<int> decoration=siteSpins;
        if (decoration.empty() && polarized && !opts.greyImposition) decoration=opts.siteSpins;
        assert(!polarized || opts.greyImposition || !decoration.empty());
        for (auto& M : hubbard)
        {
            if (M.projector==qchem::Hamiltonian::HubbardProjectorKind::PseudoAtomic && M.atomicFunction.empty())
                M.atomicFunction=PseudoAtomicFunction(*itsImp->st, M.site, M.l, opts.species);
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

    // DECISION 3 -- the accelerator, by policy (Auto pairs it with the mixer), through the public typed door.
    MakeAccelerator(opts.accelerator, params, /*carried*/bool(seed));

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
    const auto seedRole = seed ? qchem::SCFIterator::SeedRole::WorkingDensity : qchem::SCFIterator::SeedRole::Guess;   // restored = carried
    if (!seed)   // a RESTORED state arrives pre-built (Restart); otherwise the strategy the options name
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
            seed.release(), itsImp->st.get(), opts.ortho, opts.orthoTol, seedRole, opts.numBands);   // the iterator consumes it in Init
    }

    // Observe from iteration ONE: the ctor converges, so an observer attached afterwards has already
    // missed stage 0 (see SolidCalcOptions::onIteration).
    AttachProbes();

    // MOM CONTINUATION FROM THE SEED (S0e): pin the reference to the seed's OWN freshly-filled occupied
    // subspace before iteration 1, so the CONFIGURATION the seed chose survives, not merely its density.
    // No-op unless SCFParams::UseMOM is also set.
    if (opts.momFromSeed) itsImp->scf->AdoptMOMReference(*itsImp->scf->GetWaveFunction());

    // (T2's gate -- itsImp->imposed -- was stored by BuildBasis: an imposition is an ASSERTION about the answer,
    // and only then is losing the order a contradiction rather than physics.)

    EmitRunBanner(opts, itsImp->xcMesh, itsImp->st->GetNumAtoms(), imposed, hubbard);
    itsImp->hubbard = hubbard;        // RESOLVED (groups, radials): what the Hamiltonian carries and a save records
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
        BuildStage(schedule[s].accelerator, std::move(itsImp->cd), schedule[s].params);
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

// THE ACCELERATOR: the caller's, or the Auto policy's (SolidCalcOptions::accelerator).  The iterator then resolves an Auto mixer
// against the accelerator built here (tSCFIterator::ResolveMixing asks ExtrapolatesFock()), so the two halves cannot disagree.
void SolidCalculation::MakeAccelerator(SCFAccelerators::Type asked, [[maybe_unused]] const SCFParams& params,
                                       [[maybe_unused]] bool carried)
{
    using T=SCFAccelerators::Type;
    // THE Auto POLICY -- today NEUTRAL: DIIS, the default before Auto existed (and an Auto mixer stays linear beside it).  The
    // MEASURED PAIR (OpenWork §3 "Kerker/Pulay as the DEFAULT mixer", sweep P1) is: a carried stage at kT = 0 -> the Ladder on its
    // GDM rung (ao.ladderStartOnGDM) with a linear step capped at 0.3; a seeded stage -> Null, so ResolveMixing picks Kerker/Pulay.
    // It waits on N4 (pin 28: a Kerker/Pulay run's XC does not yet see rho[D]) and three more blockers listed there; \a params and
    // \a carried are the inputs it will read.
    const T t = asked==T::Auto ? T::DIIS : asked;
    // A direct minimiser rotates the occupied block into the WHOLE virtual space, which a partial eigensolve never built.
    if (itsImp->opts.numBands>0 && (t==T::GDM || t==T::Ladder))
        throw std::invalid_argument("SolidCalcOptions: numBands (a partial eigensolve) beside the GDM / Ladder accelerator -- "
                                    "they need every virtual; use Null or DIIS, or numBands 0");
    itsImp->stageAccelAsked = asked;
    itsImp->stageAccel      = t;
    itsImp->accel.reset(SCFAccelerators::Factory(t, itsImp->accOpts));
}

// ONE STAGE: a FRESH Hamiltonian + accelerator (a kT change must not carry stale DIIS history across the
// re-seed), an iterator seeded from \a carried when there is one, and MOM continuation from \a prev.
// Ownership: the iterator OWNS ham and accel and deletes them, so \a prev must outlive the adoption and is
// released immediately after it -- the same ordering the driver this replaces used.
void SolidCalculation::BuildStage(SCFAccelerators::Type accType,
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
    MakeAccelerator(accType, params, /*carried*/bool(carried));

    auto prev = std::move(itsImp->scf);      // held ONLY until the new stage has copied its MOM reference
    itsImp->scf = carried
        ? std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
              carried.release(), itsImp->st.get(), itsImp->opts.ortho, itsImp->opts.orthoTol,
              qchem::SCFIterator::SeedRole::WorkingDensity, itsImp->opts.numBands)
        : std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
              itsImp->bs.get(), itsImp->ec.get(), itsImp->ham.get(), itsImp->accel.get(),
              itsImp->opts.seed, itsImp->st.get(), itsImp->opts.ortho, itsImp->opts.orthoTol, itsImp->opts.numBands);
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
            BuildStage(schedule[s].accelerator, std::move(itsImp->cd), schedule[s].params);
        last = Converge(schedule[s].params);
        EmitStageSummary(itsImp->opts.label, s, schedule.size(), schedule[s].params.SmearingkT,
                         itsImp->converged, itsImp->scf->GetIterationCount(), itsImp->scf->GetEnergy());
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
    // ★ THE CONVERGE RESIDUE (doc/ParallelAndOraclePlan.md 1.1(a)).  The work AFTER Iterate returns is not
    // bookkeeping: GetChargeDensity() builds a fresh composite density and GetSpinDensity() rasters BOTH
    // spin channels, once per stage -- and an annealed run has a stage per schedule entry.
    qchem::report::Timed residue("scf: converge (residue -- final density + m(r) extraction)");
    // The banner and the recipe say what the iterator WILL build: SCFParams::Mixing::Auto resolved by the iterator itself.
    const SCFParams resolved=itsImp->scf->ResolveMixing(params);
    EmitSCFBanner(itsImp->opts.label, resolved, itsImp->stageAccel, params.Mixer==SCFParams::Mixing::Auto,
                  itsImp->stageAccelAsked==SCFAccelerators::Type::Auto);
    itsImp->lastOccupation = {.useMOM=params.UseMOM, .momStartIter=(int)params.MOMStartIter,
                              .kT=params.SmearingkT, .momPenalty=params.MOMSmearPenalty};   // the iterator's own conversion
    itsImp->lastkT     = params.SmearingkT;
    itsImp->lastRecipe = std::string("accel=")+AccelName(itsImp->stageAccel)+" kT="+std::to_string(params.SmearingkT)
                       + " MOM="+(params.UseMOM?"on":"off")+" KerkerG0="+std::to_string(resolved.KerkerG0)
                       + " PulayDepth="+std::to_string(resolved.PulayDepth)+" NMaxIter="+std::to_string(params.NMaxIter);
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
    itsImp->spin.reset();   // built on demand by Converged::SpinDensity()
    const bool polarizedRun = wf->GetSpinGroup()==SpinGroup::Polarized;
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
                std::cout<<"["<<itsImp->opts.label<<"] site moments (Becke-partitioned Integral w (rho_up-rho_dn) d3r) [e]:";
                for (size_t a=0;a<mu.size();a++) std::cout<<"  "<<a<<":"<<std::showpos<<std::fixed<<std::setprecision(4)<<mu[a]<<std::noshowpos;
                std::cout<<"   net="<<std::showpos<<net<<std::noshowpos<<std::defaultfloat<<std::setprecision(prec0)<<std::endl;
            }
            else if (mu.size()==0)
                std::cout<<"["<<itsImp->opts.label<<"] site moments: "
                         <<(itsImp->xcMesh.cellKind==qcMesh::UnitCellKind::Becke
                            ? "UNAVAILABLE -- the atom-centred (Becke) XC mesh carries no site blocks: a DEFECT, not a configuration"
                            : "n/a (the uniform XC mesh has no site partition)")<<std::endl;
        }
    }
    if (auto v=qchem::Diagnostics::Value("u_bands"))
        if (const auto* hub=itsImp->ham->GetHubbardChannels())
        {
            const int n=std::atoi(v->c_str());
            ReportBandSlotWeights(itsImp->opts.label, *itsImp->scf->GetWaveFunction(), *hub, *itsImp->st, n>1 ? n : 8);
        }
    // THE CHECKPOINT (CK-1), converged or not: the unattended run's safety net.  A failed WRITE warns and the run
    // goes on -- losing a finished SCF to a full disk would make the net the hazard (SaveState itself throws).
    if (!itsImp->opts.saveStateTo.empty())
    {
        try { SaveState(itsImp->opts.saveStateTo); }
        catch (const std::exception& e)
        {
            std::cout << "["<<itsImp->opts.label<<"] WARNING: state NOT saved to "<<itsImp->opts.saveStateTo
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
    const auto& o=itsImp->opts;
    std::ostringstream xc;
    xc<<MeshName(itsImp->xcMesh.cellKind);
    if (itsImp->xcMesh.cellKind==qcMesh::UnitCellKind::Becke) xc<<" nR="<<itsImp->xcMesh.nRadial<<" L="<<itsImp->xcMesh.angularDegree;
    else                                                      xc<<" eCut="<<itsImp->xcMesh.eCut;
    const RunIdentity id{.species=o.species, .Nelec=o.Nelec, .multiplicity=o.multiplicity, .functional="LDA",
                         .hubbard=itsImp->hubbard, .densityEcut=o.densityEcut, .cutoffFactor=o.cutoffFactor,
                         .xcMesh=xc.str(), .ortho=int(o.ortho), .orthoTol=o.orthoTol, .kT=kT};
    return MakeStateFingerprint(itsImp->cellMatrix, *itsImp->st, *itsImp->bs,
                                o.multiplicity>=1 ? SpinGroup::Polarized : SpinGroup::UnPolarized, id);
}

void SolidCalculation::SaveState(const std::string& path) const
{
    const auto* wf=itsImp->scf ? itsImp->scf->GetWaveFunction() : nullptr;
    if (!wf) throw std::logic_error("SolidCalculation::SaveState: no wave function yet");
    const StateSummary sum{.label=itsImp->opts.label, .converged=itsImp->converged,
                           .iterations=itsImp->scf->GetIterationCount(),
                           .energy=itsImp->scf->GetEnergy().GetTotalEnergy(), .charge=itsImp->charge,
                           .commutator=itsImp->lastCommutator, .recipe=itsImp->lastRecipe};
    {
        qchem::report::Timed timed("state: save (HDF5)");
        WriteSolidState(path, Fingerprint(itsImp->lastkT), sum, *itsImp->bs, *wf);
    }
    std::cout << "["<<itsImp->opts.label<<"] state saved to "<<path
              << (itsImp->converged ? "" : "  (NOT converged -- a restart will say so)") << std::endl;
}

Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>
SolidCalculation::Restart(const std::string& path, const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                          const SolidCalcOptions& opts, const SCFParams& params,
                          const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    using O=Outcome<std::unique_ptr<SolidCalculation>, RestartRefusal>;
    const std::string tag="["+opts.label+" restart] ";
    auto refuse=[&](RestartRefusal r)
    {
        std::cout << tag << "REFUSED ("<<path<<"): "<<r.details<<std::endl;
        return O::Fail(std::move(r));
    };
    auto saved=ReadSolidState(path);
    if (!saved) return refuse(saved.Error());

    std::unique_ptr<SolidCalculation> c(new SolidCalculation(BuildOnly{}));
    // A restored density is at (or next to) the answer: a Ladder starts on its GDM rung (user 2026-10-21, OpenWork
    // §3 "Ladder restart starts on the DIIS+Kerker rung").  For THIS stage only -- reset after Converge below.
    // kT = 0 ONLY: a smeared stage restored from an integer state hands GDM an idempotent D', so GDM would ENGAGE and
    // minimise with integer occupations, never applying the smearing (measured 2026-10-21, Si k222 kT 0.02: 6 mHa off the
    // cold smeared run).  A smeared restart starts on DIIS, as it did before.
    SCFAccelerators::SolidAcceleratorOptions restartAcc=acc;
    restartAcc.ladderStartOnGDM = params.SmearingkT==0.0;
    const std::vector<int> siteSpins=c->BuildBasis(lat, mol, opts, restartAcc);
    // Judged against the BUILT blocks (k, weight, size), not against a re-derivation of what they would be.
    auto warm=CompareFingerprints(saved->fingerprint, c->Fingerprint(params.SmearingkT));
    if (!warm) return refuse(warm.Error());
    auto cd=RestoreDensity(*saved, *c->itsImp->bs,
                           opts.multiplicity>=1 ? SpinGroup::Polarized : SpinGroup::UnPolarized);
    if (!cd) return refuse(cd.Error());

    const bool exact=warm->empty();
    {
        std::ostringstream os;
        os<<tag<<(exact ? "EXACT RESUME" : "WARM START")<<" from "<<path<<" (saved by '"<<saved->summary.label<<"', "
          <<(saved->summary.converged ? "converged" : "NOT converged")<<" after "<<saved->summary.iterations
          <<" iterations, E="<<std::setprecision(10)<<saved->summary.energy<<"; "<<saved->summary.recipe<<")";
        for (const std::string& d : *warm) os<<"\n"<<tag<<"  differs: "<<d;
        std::cout<<os.str()<<std::endl;
    }
    c->BuildRun(lat, *mol, siteSpins, cd.TakeValue(), params);
    (void)c->Converge(params);
    c->itsImp->accOpts.ladderStartOnGDM=acc.ladderStartOnGDM;   // later stages on this calculation start as asked
    if (exact)
        std::cout << tag << "resumed E="<<std::setprecision(10)<<c->itsImp->scf->GetEnergy().GetTotalEnergy()
                  << " vs saved "<<saved->summary.energy<<"  dE="<<std::setprecision(3)
                  << c->itsImp->scf->GetEnergy().GetTotalEnergy()-saved->summary.energy
                  << " after "<<c->itsImp->scf->GetIterationCount()<<" iterations"<<std::endl;
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
                                                        const std::vector<std::pair<std::string,int>>& species)
{
    int Z=0; { size_t i=0; st.ForEachSite([&](int z, const rvec3_t&, bool){ if (i==site) Z=z; i++; }); }
    const std::string element=thePeriodicTable().GetSymbol(Z);
    int Zion=0; for (const auto& [el,val] : species) if (el==element) Zion=val;
    if (Zion==0) Zion=Pseudopotential::GetGTH(element, "LDA").zion;
    const PseudoAtomOrbital& pa=PseudoAtomOrbitalOf(element, Zion, l);
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
        double homo=-std::numeric_limits<double>::infinity(), lumo=std::numeric_limits<double>::infinity();
        std::vector<double> up(nM,0.0), dn(nM,0.0);
        auto walk=[&]<class U>(const qchem::Orbitals::TOrbitals<U>& os, const qchem::Irrep& ir)
        {
            const auto* blk=dynamic_cast<const qchem::BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
            if (!blk) throw std::logic_error("SolidCalculation::ConvergeHubbardU: an orbital block that is not an Orbital_DFT_IBS");
            const double cap = ir.ms==qchem::Spin::None ? 2.0 : 1.0;   // a level is OCCUPIED when more than half full
            std::vector<const qchem::Orbitals::TOrbital<U>*> occ;
            for (const auto* o : os.template Iterate<qchem::Orbitals::TOrbital<U>>())
            {
                if (o->GetOccupation()>0.5*cap) homo=std::max(homo, o->GetEigenEnergy());
                else                            lumo=std::min(lumo, o->GetEigenEnergy());
                if (o->GetOccupation()>0.0) occ.push_back(o);
            }
            if (occ.empty()) return;
            mat_t<U> C(blk->GetNumFunctions(), occ.size());
            for (size_t i=0;i<occ.size();i++) { const vec_t<U>& c=occ[i]->GetCoeff(); for (size_t a=0;a<c.size();a++) C(a,i)=c[a]; }
            const auto amp=hub->ProjectorAmplitudes(*blk, C);
            const double w=ir.sym->GetWeight();
            for (size_t M=0;M<nM;M++)
            {
                double n=0.0;
                for (size_t i=0;i<occ.size();i++) for (size_t mu=0;mu<amp[M].rows();mu++) n+=occ[i]->GetOccupation()*std::norm(amp[M](mu,i));
                n*=w;
                if      (ir.ms==qchem::Spin::Up)   up[M]+=n;
                else if (ir.ms==qchem::Spin::Down) dn[M]+=n;
                else { up[M]+=0.5*n; dn[M]+=0.5*n; }        // the folded doublet: half in each channel
            }
        };
        for (const qchem::Irrep& ir : wf->GetQNs())
        {
            const qchem::Orbitals::Orbitals* os=wf->GetOrbitals(ir);
            if (const auto* c=dynamic_cast<const qchem::Orbitals::TOrbitals<dcmplx>*>(os)) walk(*c, ir);
            else if (const auto* r=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(os)) walk(*r, ir);
            else throw std::logic_error("SolidCalculation::ConvergeHubbardU: an orbital set of unknown scalar");
        }
        R.gap_eV.push_back(std::isfinite(lumo) && std::isfinite(homo) ? (lumo-homo)*eV : std::numeric_limits<double>::quiet_NaN());
        std::vector<double> m(nM); for (size_t M=0;M<nM;M++) m[M]=up[M]-dn[M];
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
        BuildStage(itsImp->stageAccelAsked, std::move(itsImp->cd), params);
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
            os<<std::fixed<<std::setprecision(4)<<"[U table] "<<who<<" "<<itsImp->opts.label<<(itsImp->opts.label.empty()?"":" ")
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
        {.acrossK=!itsImp->opts.perKFermiLevels, .acrossSpin=itsImp->opts.spinsShareFermi}, noise);
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
                                        qchem::Response::HubbardResolution res) const
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
    const qchem::Response::Reference ref=qchem::Response::MakeReference(*wf, itsImp->lastOccupation,
        {.acrossK=!itsImp->opts.perKFermiLevels, .acrossSpin=itsImp->opts.spinsShareFermi}, noise);
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
    auto solved=qchem::Response::LinearResponse(ref, frame, *kernel, hp.probe, Nq, krylov, resolved ? channels : perturbed);
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
    const auto& A=itsImp->cellMatrix;                             // columns = a_1..a_3 (Matrix3D is 1-based)
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
    if (itsImp->imposed)
        throw std::invalid_argument("SolidCalculation::HubbardFiniteDifferenceChi: an IMPOSED run symmetrizes the XC density "
                                    "every SCF iteration, and the alpha*P perturbation breaks the group -- run the cross-check FREE");
    const double alpha0 = J<itsImp->hubbard.size() ? itsImp->hubbard[J].alpha : 0.0;   // restored at the end
    // The channel occupations of the CURRENT orbitals, through the same projector amplitudes the probe uses.
    auto occupations=[&]()
    {
        const auto* wf=itsImp->scf->GetWaveFunction();
        if (!wf) throw std::logic_error("SolidCalculation::HubbardFiniteDifferenceChi: no wave function");
        const qchem::Response::Reference ref=qchem::Response::MakeReference(*wf, itsImp->lastOccupation,
            {.acrossK=!itsImp->opts.perKFermiLevels, .acrossSpin=itsImp->opts.spinsShareFermi}, std::numeric_limits<double>::quiet_NaN());
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
    if (itsImp->scf->ResolveMixing(params).Mixer==SCFParams::Mixing::Linear) restart.StartingRelaxRo=std::min(restart.StartingRelaxRo, 0.3);
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
            if (reseed) BuildStage(itsImp->stageAccelAsked, nullptr, stageParams);
            else        BuildStage(itsImp->stageAccelAsked, std::move(itsImp->cd), stageParams);
            auto out=Converge(stageParams);
            if (!out) {failed=out.Error(); break;}
            (a>0 ? R.nPlus : R.nMinus)=occupations();
        }
        target->SetPerturbation(J, alpha0);
    }
    // RESTORE the unperturbed, unfrozen state -- the run is left as it was found, whatever happened above.
    BuildStage(itsImp->stageAccelAsked, reseed ? std::move(ground) : std::move(itsImp->cd), restart);
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
void   SolidCalculation::OnIteration(Observer obs)      { itsImp->opts.onIteration=std::move(obs); AttachProbes(); }
bool   SolidCalculation::DidConverge()    const         { return itsImp->converged; }
size_t SolidCalculation::IterationCount() const         { return itsImp->scf->GetIterationCount(); }

const RunDiagnostics& SolidCalculation::Diagnostics() const { return itsImp->diag; }

const qcMesh::MeshParams& SolidCalculation::ResolvedXCMesh() const { return itsImp->xcMesh; }
const BasisSet::Complex_BS& SolidCalculation::Basis() const { return *itsImp->bs; }

} //namespace qchem
