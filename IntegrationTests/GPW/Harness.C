// File: IntegrationTests/GPW/Harness.C  THE GPW TEST/PROBE HARNESS (doc/TestSuitePlan.md phase 2, 2026-09-15).
//
// Everything the GPW integration tests AND the CLIapps/gpwprobe instruments share, and NOTHING that runs an
// SCF: every run is a qchem::SolidCalculation and every cell is a qchem::Materials entry.  What lives here is
// the RECIPE vocabulary -- the basis factories, the production-shaped SCF gates, the named per-material
// recipes (Al, NaF, the Mn box), the env A/B valves -- plus the instruments (run report bracket, trajectory
// fingerprint, the XC quadrature probes and the Becke ladder) that a gate and a hand-run probe both read.
// Built as its own library (qcGPW_Harness) so ITMain and gpwprobe compile one copy.  No gtest in here: a
// helper that finds something wrong THROWS (doc/CleanupCandidates.md: throw is a marker), and the caller
// decides whether that is a test failure or a probe abort.
module;
#include <memory>
#include <vector>
#include <cmath>
#include <cstdlib>   // std::getenv/std::atof (the env A/B valves)
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <string>
#include <iostream>
#include <iomanip>      // setprecision (the order-parameter trajectory line)
export module qchem.Tests.GPW_Harness;
export import qchem.SolidCalculation;                    // SolidCalculation, SolidCalcOptions, SCFStage
export import qchem.Materials;                           // Materials::Material / Get
export import qchem.Lattice_3D;                          // Lattice_3D, UnitCell
export import qchem.BasisSet;                            // Complex_BS, Real_BS
export import qchem.BasisSet.Gaussian.Point.Factory;     // Gaussian::Factory, BasisSetData/Engine/Angular
export import qchem.SCFParams;                           // SCFParams
export import qchem.SCFIterator;                         // SCFProgress, SolidSCFIterator::Observer
export import qchem.Reporting;                           // report:: (the run bracket)
export import qchem.Mesh;                                // qcMesh::MeshParams / UnitCellKind
export import qchem.Mesh.XCPolicy;                       // BeckeXCParams
export import qchem.ChargeDensity.DensitySampler;        // DensitySampler, MakeDensitySampler
export import qchem.Energy;                              // EnergyBreakdown
export import qchem.Outcome;                             // Outcome<Converged,SCFFailure>
export import qchem.Types;
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;  // MakeSphericalLatticeView (GPW_SPHERICAL=1)
import qchem.BasisSet.DeltaFit_IBS;              // DeltaFit_IBS -- the delta basis the singles strategy runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform -- the uniform probe's own point count
import qchem.BasisSet.Lattice.BasisSet;           // VetGpwConditioning / EmitGpwGrids (GpwReport's pre-flight faces)
import qchem.Hamiltonian.Internal.PWTerms;        // Vxc_Quadrature (the XC probes)
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian;                         // cDynamic_HT, Dynamic_HT_RealBlock
import qchem.ChargeDensity;                       // cDM_CD
import qchem.ChargeDensity.Seed;                  // SeedStrategy
import qchem.SCFAccelerator.Factory;              // SCFAccelerators::Type
import qchem.Symmetry.Factory;                    // BlochFactory
import qchem.LASolver;                            // qchem::Ortho
import qchem.Blaze;

export namespace qchem::tests::gpw
{
using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Complex_BS;
using qchem::BasisSet::Gaussian::BasisSetData;

//---------------------------------- the quadrature bundle ----------------------------------
// The delta fit basis the SINGLES quadrature runs on: it owns the mesh's points, weights and functions.
std::shared_ptr<const qchem::BasisSet::DeltaFit_IBS> DeltaFitOver(qchem::BasisSet::FitQuadrature q)
{
    return std::make_shared<const qchem::BasisSet::DeltaFit_IBS>(std::move(q),
               qchem::Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
}

// ...and the SINGLES engine over that same bundle -- which since 2026-08-24 needs it too, because the
// atomic partition, the orbit fold and the Shubnikov tags are INJECTED into the strategy rather than asked
// of the basis (BasisSet::FIT_SF_ABS lost Symmetrize/SymmetrizeSpin).  ONE bundle, handed to both
// collaborators, exactly as tBasisSet::CreateVxcFitBasisSet does it in production -- so a probe cannot
// accidentally give the basis one quadrature and the strategy another.
// ★ THROUGH THE FACTORY, NOT `make_shared` ON A STRATEGY (user, 2026-09-10: *"If there is a factory
// available the tests should exercise that interface instead of direct construction"*).  A delta fit basis
// carries points and nothing else, so `MakeDensitySampler` can only pick SINGLES -- which means naming the
// strategy here bought nothing except a dependency on an Internal module.  Exercising the factory also gates
// the capability decision itself, which direct construction silently skipped.
std::shared_ptr<const qchem::ChargeDensity::DensitySampler>
SinglesEngineOver(qchem::BasisSet::FitQuadrature q)
{
    auto fit=DeltaFitOver(q);   // NOT inline with the move below: argument evaluation order is unspecified
    return qchem::ChargeDensity::MakeDensitySampler(std::move(fit), std::move(q));
}

//---------------------------------- the basis factories ----------------------------------
// The valence Si Gaussian basis (SIPP, MnD-Cartesian) on ANY structure -- the L_PP / GPW_UT builder.
std::shared_ptr<const Real_BS> MakeBasis(const Structure& st)
{
    return std::shared_ptr<const Real_BS>(
        BasisSet::Gaussian::Factory(BasisSetData::SIPP, &st,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
}
// GPW_SPHERICAL=1 (doc/SphericalLatticePlan.md I1/I2): wrap the molecular basis in its spherical
// (contaminant-free) lattice view -- the span-matched A/B against CP2K's spherical-d convention.
// s/p-only bases are unchanged in SPAN (T = identity blocks), so Si/NaF runs double as null tests.
std::shared_ptr<const Real_BS> MaybeSpherical(std::shared_ptr<const Real_BS> bs)
{
    if (std::getenv("GPW_SPHERICAL"))
        return BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(std::move(bs));
    return bs;
}
// The SHORT-RANGE variant (most diffuse valence primitives dropped) -- well-conditioned Bloch overlap in a solid.
std::shared_ptr<const Real_BS> MakeBasisSR(const Structure& st)
{
    return MaybeSpherical(std::shared_ptr<const Real_BS>(
        BasisSet::Gaussian::Factory(BasisSetData::SIPP_SR, &st,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian)));
}
// The low-q GTH valence basis (valgen-generated; carries Al/Na/F) -- the Al block drives the FCC-Al metal test.
// GPW_BASIS_SPH=1 swaps in VALENCE_LOWQ_SPH (the Mn true-s window restored -- doc/SphericalLatticePlan.md I3);
// meaningful ONLY together with GPW_SPHERICAL=1 (under Cartesian d that file is contaminant-rank-deficient).
//
// GPW_BASIS_SPAN=sph|va|vb NAMES THE SPAN, and is the form doc/Benchmark.md's MnO rows use.  va/vb are the
// exact-span variants CP2K also holds function-for-function (VALENCE-LOWQ-V{A,B}; VA = 118 functions on the
// MnO magnetic cell, held FULL RANK by both codes -- VB = 128).  They exist as BASIS FILES because the span
// used to be produced by doc/scripts/bisect_valence_sph.py OVERWRITING the committed valence_lowq_sph.bsd in
// the working tree: the run could not say which span it ran, and the row could not be reproduced afterwards.
// GPW_BASIS_SPH=1 == GPW_BASIS_SPAN=sph, kept because the banked run recipes are written with it.
std::shared_ptr<const Real_BS> MakeBasisLowQ(const Structure& st, BasisSetData which=BasisSetData::VALENCE_LOWQ_SR)
{
    if (which==BasisSetData::VALENCE_LOWQ_SR)
    {
        if (std::getenv("GPW_BASIS_SPH")) which=BasisSetData::VALENCE_LOWQ_SPH;
        if (const char* s=std::getenv("GPW_BASIS_SPAN"))
        {
            const std::string span(s);
            if      (span=="sph") which=BasisSetData::VALENCE_LOWQ_SPH;
            else if (span=="va" ) which=BasisSetData::VALENCE_LOWQ_VA;
            else if (span=="vb" ) which=BasisSetData::VALENCE_LOWQ_VB;
            else if (span=="sr" ) which=BasisSetData::VALENCE_LOWQ_SR;
            else throw std::runtime_error("GPW_BASIS_SPAN: expected one of sr|sph|va|vb, got '"+span+"'");
        }
    }
    return MaybeSpherical(std::shared_ptr<const Real_BS>(
        BasisSet::Gaussian::Factory(which, &st,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian)));
}


//---------------------------------- run reporting + the trajectory instruments ----------------------------------
// Shared GPW run reporting.  ANY GPW driver -- RunGPW's fixed recipe OR a bespoke one (GPW_NaF.Γ_Imp_Anchor's
// multi-species PP + Ladder accelerator) -- gets automatic reporting by holding one of these: RAII brackets the
// run (Begin/End + optional console), and the driver calls VetBasis (the fail-fast conditioning pre-flight,
// BEFORE it builds the Hamiltonian/grids) then EmitGrids.  Keeps each driver's own Hamiltonian/accelerator/
// params -- only the reporting + setup-order discipline is shared.
struct GpwReport
{
    explicit GpwReport(const std::string& name, bool verbose)
    {
        qchem::report::Begin(name);
        // GPW_REPORT=1 forces the console on for ANY driver, whatever its hard-coded `verbose`.  The
        // report carries the timing ledger (the where-did-the-time-go table), and a cost measurement
        // must not require editing the test that happens to reproduce the cost.
        static const bool kEnvReport=[]{ const char* s=std::getenv("GPW_REPORT"); return s && std::atoi(s)!=0; }();
        if (verbose || kEnvReport) qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
    }
    ~GpwReport() { qchem::report::ClearConsole(); qchem::report::End(); }
    GpwReport(const GpwReport&) = delete;
    GpwReport& operator=(const GpwReport&) = delete;

    //! Fail-fast conditioning pre-flight: emit basis.perIrrep/removed, return the redundant-function count
    //! (0 == OK).  Call BEFORE building the Hamiltonian/grids -- a positive return is the cue to abort.
    size_t VetBasis(const Complex_BS& bs)
    {
        qchem::report::Log("vetting basis conditioning");
        qchem::report::Section basis("basis");
        return BasisSet::Lattice::VetGpwConditioning(bs);
    }
    //! Emit the grids section (this is where the ladder is actually built) -- only after VetBasis passed.
    void EmitGrids(const Complex_BS& bs)
    {
        qchem::report::Log("building grid ladder");
        BasisSet::Lattice::EmitGpwGrids(bs);
    }
};

// PROBE (dynamics fingerprint, doc/GPWPlan §0): classify an SCF trajectory captured via the Observer hook.
// Three pathologies have DISTINCT time-series signatures, so one line names which regime the run is in --
// separating "the iteration can't find the min" (dynamics: sloshing/divergence) from "the min is a fit floor"
// (functional).  Captured from qchem::SCFIterator::SCFProgress {iteration, energy, dE=|ΔE|, [F,D], Δρ}.
struct FpRow { size_t it; double E, dEabs, fd, drho, order=0, eee=0; };

// The ORDER-PARAMETER trajectory, one compact line (printed only when a probe was set).  The per-iteration
// SCF column already shows it live; this line is the POST-MORTEM -- the whole time series in one place, with
// the death iteration named, so a bounded diagnosis run answers "WHERE did the order die?" without anyone
// re-reading a 30-row table.
//
// The reference is the run's PEAK |order|, NOT iteration 1 (fixed 2026-08-07, first use).  Iteration 1 is a
// terrible baseline: its Fock comes from the SEED, so the order parameter there is whatever survived one
// crude step -- MnO reads 0.0046 at iteration 1, peaks at 0.1064 by iteration 7 as the self-consistent
// exchange splitting builds it back up, then decays to 7e-5.  Judged against iteration 1 that is "SURVIVED"
// (1.6%); judged against the peak it is a 1400x collapse, which is what actually happened.  A quantity that
// GROWS before it dies needs the high-water mark as its yardstick.
void OrderTrajectory(const std::vector<FpRow>& s, const std::string& name, const char* label)
{
    if (s.empty() || name.empty()) return;
    std::cout << "["<<label<<" "<<name<<"]";
    for (const auto& r : s) std::cout << " " << std::fixed << std::setprecision(4) << r.order;
    std::cout << std::defaultfloat << std::endl;
    size_t peak=0;                                          // the high-water mark and where it happened
    for (size_t i=1;i<s.size();++i) if (std::fabs(s[i].order)>std::fabs(s[peak].order)) peak=i;
    const double mMax=std::fabs(s[peak].order), dead=0.01*mMax;
    size_t died=s.size();                                   // first index from which |order| stays below dead
    while (died>0 && std::fabs(s[died-1].order)<=dead) --died;
    std::cout << "["<<label<<" "<<name<<"] iter1="<<s.front().order
              << " peak="<<s[peak].order<<"@iter"<<s[peak].it<<" final="<<s.back().order;
    if (mMax>0.0 && died<s.size() && died>peak)
        std::cout << "  ** DIED at iteration "<<s[died].it<<" (|"<<name<<"| < 1% of the peak from there on)";
    else if (mMax>0.0 && std::fabs(s.back().order) < 0.5*mMax)
        std::cout << "  ** DECAYING (final is "<<(100.0*std::fabs(s.back().order)/mMax)<<"% of the peak)";
    else
        std::cout << "  order SURVIVED the run";
    std::cout << std::endl;
}
void Fingerprint(const std::vector<FpRow>& s, const char* label)
{
    if (s.empty()) { std::cout << "["<<label<<" fp] (no iterations)"<<std::endl; return; }
    const size_t n=s.size(), w=std::min<size_t>(n,8);
    double emin=1e300, emax=-1e300;
    for (size_t k=n-w;k<n;++k){ emin=std::min(emin,s[k].E); emax=std::max(emax,s[k].E); }
    size_t flips=0;                                   // sign changes of successive SIGNED ΔE over the window
    for (size_t k=n-w+1;k+1<n;++k)
    {
        double d0=s[k].E-s[k-1].E, d1=s[k+1].E-s[k].E;
        if (d0*d1<0.0) ++flips;
    }
    const double Ef=s.back().E, drhoF=s.back().drho, amp=emax-emin;
    // IS Δρ ACTUALLY FLOORED, or merely SMALL and still falling?  The two look identical in the last row and
    // are opposite diagnoses -- a floor is grid/functional work, a slow descent is just a low iteration cap.
    // Measured by the geometric rate over the last few iterations: a floored Δρ has rate ~1, a descending one
    // is bounded away from it.  (MnO run 24 was labelled FIT-FLOOR STALL while Δρ was still falling 7% PER
    // ITERATION -- 2.30e-4 -> 1.44e-4 over its last seven, no plateau at all -- and that mislabel sent a
    // whole session's reading of the campaign to the grids instead of to the iteration cap.)
    double drhoRate=0.0;   // Δρ(last)/Δρ(5 back), ^(1/5): ~1 = floored, <1 = still descending
    if (s.size()>=6 && s[s.size()-6].drho>0.0)
        drhoRate = std::pow(drhoF/s[s.size()-6].drho, 0.2);
    const bool drhoFloored = (drhoRate > 0.98);   // <2% per iteration = not going anywhere
    const double relAmp=amp/std::max(std::fabs(Ef),1e-30);   // energy swing RELATIVE to the total
    // Verdict priority separates the three pathologies (+ the benign degenerate case) by their distinct
    // signatures.  KEY distinction: a degenerate open shell has the ENERGY settled (small relAmp) while Δρ
    // never falls (ρ rotates in the degenerate subspace) -- benign; charge-transfer SLOSHING swings BOTH.
    const char* verdict =
        (drhoF < 1e-5)                                   ? "CONVERGED" :
        (drhoF > 1e-3 && relAmp < 5e-3)                  ? "DENSITY-DEGENERATE (E settled, ρ rotates -- benign)" :
        (flips >= 3 && relAmp > 5e-3)                    ? "OSCILLATING (charge-transfer sloshing / mixing unstable)" :
        (drhoF > 1e-5 && s.back().dEabs < 1e-5 && drhoFloored)
                                                         ? "FIT-FLOOR STALL (Δρ floored, ΔE tiny -- functional/grid)" :
        (drhoF > 1e-5 && s.back().dEabs < 1e-5)          ? "UNSETTLED (Δρ still descending at the cap -- raise NMaxIter, NOT a floor)" :
        (std::fabs(Ef) > 3.0*std::fabs(s.front().E))     ? "DIVERGING" : "UNSETTLED (hit iter cap mid-descent)";
    // Efinal AT A STATED PRECISION.  This line is what the plan docs quote, and it used to inherit whatever
    // cout was left at: run 61 printed -61.4029762007 only because its verbose per-iteration table had set
    // fixed+10 upstream, while the same run without GPW_MNO_VERBOSE printed -61.4.  A number's precision
    // must not depend on which OTHER diagnostics were switched on, so it is set (and restored) here.
    const std::streamsize prec0=std::cout.precision();
    std::cout << "["<<label<<" fp] iters="<<n<<" Efinal="<<std::setprecision(10)<<Ef<<std::setprecision(prec0)
              << " lastΔρ="<<drhoF
              << " oscFlips(last"<<w<<")="<<flips<<" Eamp(last"<<w<<")="<<amp<<" relAmp="<<relAmp
              << "  => "<<verdict<<std::endl;
}

// Optional keep-alive handles for post-SCF term-level probes (the Becke XC gate): the basis and the
// converged density, which stays valid after the iterator tears down because its basis block is bs.
struct GpwHandles   // (kept as a name: the probes below take it) -- a NON-OWNING VIEW over the facade's objects
{
    const Complex_BS* bs=nullptr;
    const qchem::ChargeDensity::cDM_CD* cd=nullptr;
};

//---------------------------------- the facade harness ----------------------------------
using qchem::Materials::Material;
using qchem::SolidCalcOptions;

//! The lattice a material's cell is run on: the k-mesh is the TEST's axis, never the material's.
Lattice_3D LatticeOf(const Material& m, const ivec3_t& k=ivec3_t(1,1,1)) { return Lattice_3D(*m.cell, k); }

//! The run options a material dictates -- its electron count and its pseudopotential vocabulary -- and a
//! label.  Everything else (grid, k, symmetry, spin, kT, machinery) is stated at the call site: those are
//! the axes a test is a point in, and the facade's defaults are the elided ones.
SolidCalcOptions OptionsFor(const Material& m, const std::string& label)
{
    SolidCalcOptions o;
    o.label=label; o.Nelec=m.Nelec(); o.species=m.species;
    return o;
}

//! The SCF gates.  Two named recipes cover nearly every test: PRODUCTION (the anchors' gates -- energy AND
//! density, 60 iterations) and TIGHT (density-only to 1e-6, the old positional driver's default).  The
//! virial/FD gates are off on every pseudopotential run (the textbook -V/K=2 does not hold).
SCFParams Gates(size_t nmax, double minDrho, double minDE)
{
    SCFParams par;
    par.NMaxIter=nmax; par.MinΔρ=minDrho; par.MinΔE=minDE;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.StartingRelaxRo=0.3; par.MergeTol=1e-4;
    return par;
}
SCFParams ProductionGates() { return Gates(60, 1e-3, 1e-6); }
SCFParams TightGates(size_t nmax=120) { return Gates(nmax, 1e-6, 1e30); }

//! The A/B valves the old positional driver carried, applied to a stated recipe.  Diagnostics only: unset,
//! nothing changes.  GPW_IMPOSE=0/1 (imposition is the one part of a multi-k run that reconstructs the full
//! BZ from irreducible blocks, so it is the first thing to remove when TRIM and complex meshes disagree);
//! GPW_SMEAR=kT (an integer aufbau fill is ambiguous at a degenerate frontier); GPW_VERBOSE=1; GPW_REAL=0
//! (build every block complex: a defect that appears only with real-TRIM narrowing was a wrongly typed
//! block); GPW_SEED=coreguess|uniform|sad|ionicsad (CoreGuess separates the operators from the seed);
//! GPW_ORTHO=cholesky|eigen|svd (a defect under one ortho only is IN the ortho); GPW_KERKER_G0=g (the
//! density preconditioner the supercell ladder cannot run without).
void EnvOverrides(SolidCalcOptions& o, SCFParams& par)
{
    if (const char* im=std::getenv("GPW_IMPOSE")) o.imposeSymmetry=std::atoi(im)!=0;
    if (const char* kt=std::getenv("GPW_SMEAR"))  par.SmearingkT=std::atof(kt);
    if (std::getenv("GPW_VERBOSE"))               par.Verbose=true;
    if (const char* rl=std::getenv("GPW_REAL"))   o.forceComplex=std::atoi(rl)==0;
    if (const char* kg=std::getenv("GPW_KERKER_G0")) par.KerkerG0=std::atof(kg);
    if (const char* sd=std::getenv("GPW_SEED"))
    {
        const std::string v(sd);
        if      (v=="coreguess") o.seed=qchem::ChargeDensity::SeedStrategy::CoreGuess;
        else if (v=="uniform")   o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
        else if (v=="sad")       o.seed=qchem::ChargeDensity::SeedStrategy::SAD;
        else if (v=="ionicsad")  o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
        else throw std::runtime_error("GPW_SEED: expected coreguess|uniform|sad|ionicsad, got '"+v+"'");
    }
    if (const char* ot=std::getenv("GPW_ORTHO"))
    {
        const std::string v(ot);
        if      (v=="cholesky") o.ortho=qchem::Cholesky;
        else if (v=="eigen")    o.ortho=qchem::Eigen;
        else if (v=="svd")      o.ortho=qchem::SVD;
        else throw std::runtime_error("GPW_ORTHO: expected cholesky|eigen|svd, got '"+v+"'");
    }
}

//! The per-iteration trace a test can attach through \c SolidCalcOptions::onIteration -- live from the
//! constructor, so stage 0 is in it -- and read back as the fingerprint / order trajectory the campaign
//! instruments print.  (The assertions stay on the facade's answers; this is the human-readable side.)
struct Trace
{
    std::vector<FpRow> rows;
    qchem::SCFIterator::SolidSCFIterator::Observer Observer()
    {
        return [this](const qchem::SCFIterator::SCFProgress& p)
               { rows.push_back({p.iteration, p.energy, p.dE, p.commutator, p.drho, p.order}); };
    }
    void Print(const std::string& label, bool polarized=false) const
    {
        Fingerprint(rows, label.c_str());
        OrderTrajectory(rows, polarized ? std::string("m_site") : std::string(), label.c_str());
    }
};

//! The failure text of an outcome, for an ASSERT message ("" when it succeeded).
template <class R> std::string Why(const R& r) { return r ? std::string() : r.Error().details; }

//! The basis + converged density a term-level probe needs, from a CONVERGED facade run.
GpwHandles Handles(const qchem::SolidCalculation& calc)
{
    auto r=calc.Result();
    if (!r) throw std::runtime_error("Handles: the run did not converge -- "+r.Error().details);
    return {&calc.Basis(), &r->DensityMatrix()};
}
//! ...and from a BOUNDED run (a symmetry probe that stops at a fixed iteration): the last iterate's.
GpwHandles LastIterateHandles(const qchem::SolidCalculation& calc) { return {&calc.Basis(), calc.LastIterateDensity()}; }

//---------------------------------- the named per-material recipes ----------------------------------
SolidCalcOptions AlOptions(const Material& al, const std::string& label="Al FCC Gamma")
{
    SolidCalcOptions o=OptionsFor(al, label);
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    // Uniform STATED (V2.2): Al has no entry in atomic_valence_densities.json yet, so the IonicSAD default
    // would throw at seed time.  valgen can generate one; until then this is the explicit opt-in.
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    return o;
}
//! The Al block's gates: density-only to 1e-5 (a metal's energy settles long before its density does).
SCFParams AlGates(size_t nmax=60) { return Gates(nmax, 1e-5, 1e30); }

//! THE COMMITTED NaF PRODUCTION RECIPE (from the NaF rocksalt gate): the diffuse-trimmed SR2 basis, the
//! Fock DIIS->GDM Ladder on |ΔE/E|, IonicSAD, pivoted Cholesky, Kerker G0=1 against the low-G charge-transfer
//! slosh, delayed MOM.  One place, so the ladder probes and the Becke gate cannot drift from it.
std::shared_ptr<const Real_BS> MakeBasisNaFSR2(const Structure& st)
{
    return std::shared_ptr<const Real_BS>(BasisSet::Gaussian::Factory(
        BasisSetData::VALENCE_LOWQ_SR2, &st, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
}
SolidCalcOptions NaFOptions(const Material& naf, const std::string& label)
{
    SolidCalcOptions o=OptionsFor(naf, label);
    o.accelerator=qchem::SCFAccelerators::Type::Ladder;
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    return o;
}
SCFParams NaFGates()
{
    SCFParams par=Gates(200, 1e-4, 1e-8);
    par.StartingRelaxRo=0.45; par.KerkerG0=1.0;
    par.UseMOM=true; par.MOMStartIter=10;
    par.MergeTol=SCFParams{}.MergeTol;   // (the recipe never set it)
    return par;
}
//! The Mn sextet atom-in-box recipe (GPW_MnBox.Γ_M6_Smear_eqFinite): the finite-molecule mode, IonicSAD, pivoted ortho.
SolidCalcOptions MnBoxOptions(const Material& box, const std::string& label)
{
    SolidCalcOptions o=OptionsFor(box, label);
    o.multiplicity=6;                                  // S=5/2 Hund: nUp=6, nDown=1
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.imposeSymmetry=false;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    return o;
}

//---------------------------------- the XC quadrature probes ----------------------------------
// One XC quadrature's answers on a converged density: E_xc, the quadrature's rho integral vs Tr(DS),
// and the summed Dirac+VWN V_xc matrix per Bloch block.  Uniform and Becke probes share the shape, so
// gates can diff ANY pair (uniform-vs-Becke, or Becke-vs-Becke for internal convergence).
struct XCProbe
{
    std::string label;
    double Exc=0, rhoLost=0;
    //! The quadrature's OWN point count -- asked of the mesh/raster, never computed from a rule.  Every
    //! route is a FIT (user, 2026-09-06): Becke fits onto DELTA functions at its mesh points, the uniform
    //! route onto the plane-wave \f${G}\f$ raster -- both orthonormal metrics -- so "how many points" is
    //! the one cost axis they share, and the ladder can only compare rungs if each states its own.
    size_t nPts=0;
    std::vector<hmat_t<dcmplx>> M;
};

// Fill a probe from an exchange+correlation term pair (through the PUBLIC term faces).
// MIXED-AWARE (doc/RealComplexPlan.md 3c-3): since the harness flip a Γ block is REAL, and the
// same-scalar Iterate view THROWS on it -- so walk per index, take the cross-scalar view first, and
// drive the term's real-block face.  The probe's currency stays chmat_t: a real block's matrix WIDENS
// losslessly (its imaginary part is exactly zero, and the per-term gates pin the real assembly bitwise
// against the complex one), so uniform-vs-Becke diffs read exactly as they did before the flip.
hmat_t<dcmplx> ProbeBlockMatrix(Hamiltonian::cDynamic_HT& t, const BasisSet::Complex_BS& bs, size_t i,
                                const qchem::ChargeDensity::cDM_CD* cd)
{
    if (const auto* rb=bs.GetRealIBS(i))
    {
        const auto* face=dynamic_cast<const Hamiltonian::Dynamic_HT_RealBlock*>(&t);
        if (!face) throw std::runtime_error("ProbeBlockMatrix: a real block needs the term's real-block face (3c-1)");
        const rsmat_t R=face->GetMatrix(rb,Spin::None,cd);
        hmat_t<dcmplx> M(R.rows());
        for (size_t r=0;r<R.rows();r++)
            for (size_t s=r;s<R.columns();s++) M(r,s)=R(r,s);
        return M;
    }
    return t.GetMatrix(bs[i],Spin::None,cd);
}

XCProbe ProbeXC(const std::string& label, const GpwHandles& h,
                Hamiltonian::cDynamic_HT& x, Hamiltonian::cDynamic_HT& c,
                const qchem::EnergyBreakdown& e)
{
    XCProbe p; p.label=label; p.Exc=e["Exc"]; p.rhoLost=e.charge.lost;
    for (size_t i=0;i<h.bs->GetNumIBS();++i)
    {
        hmat_t<dcmplx> M=ProbeBlockMatrix(x,*h.bs,i,h.cd);
        M+=ProbeBlockMatrix(c,*h.bs,i,h.cd);
        p.M.push_back(std::move(M));
    }
    return p;
}

XCProbe UniformXCProbe(const GpwHandles& h, const std::shared_ptr<const Structure>& st)
{
    auto exch=std::make_shared<Hamiltonian::SlaterExchange>(2.0/3.0);
    auto corr=std::make_shared<Hamiltonian::VWN_Correlation>();
    qcMesh::MeshParams mp; mp.relCutoff=std::max(exch->GridCutoffFactor(), corr->GridCutoffFactor());
    ChargeDensity::fitbasis_t vfb(h.bs->CreateVxcFitBasisSet(st.get(), mp));
    auto pair=ChargeDensity::MakeDensitySampler(vfb);   // raster fit basis -> the pair/collocation strategy
    Hamiltonian::Vxc_Quadrature x(exch,pair,SpinGroup::UnPolarized), c(corr,pair,SpinGroup::UnPolarized);
    EnergyBreakdown e; x.GetEnergy(e,h.cd); c.GetEnergy(e,h.cd);
    XCProbe p=ProbeXC("uniform", h, x, c, e);
    // Its own raster's size -- the same face the Hartree energy asks for its quadrature rule.
    auto* rt=dynamic_cast<const BasisSet::G_RasterTransform*>(vfb.get());
    p.nPts = rt ? rt->RasterSize() : 0;
    return p;
}

XCProbe BeckeXCProbe(const GpwHandles& h, const std::shared_ptr<const Structure>& st,
                     const std::string& label, const qcMesh::MeshParams& mpB)
{
    auto exch=std::make_shared<Hamiltonian::SlaterExchange>(2.0/3.0);
    auto corr=std::make_shared<Hamiltonian::VWN_Correlation>();
    auto mesh=std::make_shared<const qcMesh::Mesh>(st->CreateIntegrationMesh(mpB));
    auto engine=SinglesEngineOver({mesh, {}});      // ONE quadrature, shared by the pair; free probe: no fold
    Hamiltonian::Vxc_Quadrature x(exch,engine,SpinGroup::UnPolarized), c(corr,engine,SpinGroup::UnPolarized);
    EnergyBreakdown e; x.GetEnergy(e,h.cd); c.GetEnergy(e,h.cd);
    XCProbe p=ProbeXC(label, h, x, c, e);
    p.nPts=mesh->size();                            // the mesh's own count, not a (nR x degree) rule
    return p;
}

// Print + return the elementwise V_xc gap between two probes; EXPECTs applied by the caller.
// (qcMesh::BeckeXCParams -- the calibrated Becke recipe these probes use -- is library policy, declared
//  beside MeshParams in src/Mesh/Mesh.C.)
double DiffXC(const XCProbe& A, const XCProbe& B)
{
    std::printf("[Becke gate] Exc %s=%.6f %s=%.6f  dExc=%+.3e  (rho-lost %+.3e / %+.3e)\n",
                A.label.c_str(), A.Exc, B.label.c_str(), B.Exc, B.Exc-A.Exc, A.rhoLost, B.rhoLost);
    double dworst=0;
    for (size_t blk=0; blk<A.M.size(); blk++)
    {
        const auto &MA=A.M[blk], &MB=B.M[blk];
        double dmax=0, amax=0;
        size_t im=0, jm=0;
        const size_t n=MA.rows();
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
            {
                if (std::abs(MA(i,j)-MB(i,j))>dmax) {dmax=std::abs(MA(i,j)-MB(i,j)); im=i; jm=j;}
                amax=std::max(amax, std::abs(MA(i,j)));
            }
        std::printf("[Becke gate] Vxc blk%zu n=%zu  max|%s|=%.3e  max|%s-%s|=%.3e at (%zu,%zu): %.6f vs %.6f\n",
                    blk, n, A.label.c_str(), amax, A.label.c_str(), B.label.c_str(), dmax, im, jm,
                    MA(im,jm).real(), MB(im,jm).real());
        dworst=std::max(dworst,dmax);
    }
    return dworst;
}

// Score a ladder of Becke meshes on one frozen density against a fine reference.  Prints one row per
// rung: the mesh size, dExc, and the worst V_xc matrix element deviation -- the same two numbers the
// Becke gate uses, so a rung that passes here would pass the gate.
void BeckeLadder(const GpwHandles& h, const std::shared_ptr<const Structure>& st, const char* system)
{
    auto mesh=[&](int nR, int deg){ qcMesh::MeshParams mp=qcMesh::BeckeXCParams(nR, 2.0, deg);
                                    mp.angular=qcMesh::AngularKind::GaussLegendre;   // arbitrary degree, no table gaps
                                    return mp; };
    // The reference must be a STRICT REFINEMENT of the rung family -- same mhl_alpha, more points on both
    // axes.  (The tree's existing REF80 probe uses alpha=1.0, which is fine for a one-off cross-check but
    // would put a fixed alpha-mismatch offset under every rung here and read as an error FLOOR the axes
    // never get below.  A ladder needs its reference inside its own family.)
    const XCProbe REF=BeckeXCProbe(h, st, "REF", mesh(100,41));
    std::printf("\n[V2.6 ladder %s] reference nR=100 GL-41 (same alpha=2.0 family): Exc=%.8f, %zu pts\n",
                system, REF.Exc, REF.nPts);
    // Worst |V_xc(i,j)| deviation from the reference, over all blocks -- the error in the matrix that
    // actually enters the Fock, which is what a quadrature/fit has to get right.
    auto dVxc=[&REF](const XCProbe& P)
    {
        double dV=0;
        for (size_t b=0; b<REF.M.size() && b<P.M.size(); b++)
            for (size_t i=0;i<REF.M[b].rows();++i)
                for (size_t j=0;j<REF.M[b].columns();++j)
                    dV=std::max(dV, std::abs(dcmplx(P.M[b](i,j))-dcmplx(REF.M[b](i,j))));
        return dV;
    };
    auto print=[&](const char* what, const XCProbe& P)
    {
        std::printf("[V2.6 ladder %s] %-14s %9zu pts  dExc=%+.3e  max|dVxc|=%.3e\n",
                    system, what, P.nPts, P.Exc-REF.Exc, dVxc(P));
    };
    // ★ THE YARDSTICK (user, 2026-09-06): the question is not "has this rung converged" but "does it match
    // what the UNIFORM route already delivers" -- because that is the accuracy the production runs get on
    // the cheap route, and the Becke mesh is only worth its setup if it is at least that good.  Both routes
    // are FITS of the same v_xc: the uniform one onto the plane-wave {G} raster, the Becke one onto DELTA
    // functions at its mesh points, both with orthonormal metrics.  So they are scored identically here,
    // against one reference, and the ANSWER TO READ OFF is the smallest (nR, degree) whose max|dVxc| is at
    // or below the uniform row's.
    const XCProbe UNI=UniformXCProbe(h, st);
    print("UNIFORM", UNI);
    std::printf("[V2.6 ladder %s] --- ANGULAR sweep (nRadial=40 fixed) ---\n", system);
    for (int deg : {5,7,9,11,15,17,21,23,29})
        { char b[32]; std::snprintf(b,sizeof b,"nR=40 GL-%d",deg);  print(b, BeckeXCProbe(h,st,"rung",mesh(40,deg))); }
    std::printf("[V2.6 ladder %s] --- RADIAL sweep (degree=29 fixed) ---\n", system);
    for (int nR : {10,15,20,25,30,40,60})
        { char b[32]; std::snprintf(b,sizeof b,"nR=%d GL-29",nR);   print(b, BeckeXCProbe(h,st,"rung",mesh(nR,29))); }
}

//---------------------------------- the MnO campaign arm (the probe's) ----------------------------------
// THE MnO ARM, ON THE FACADE (2026-08-25).  What used to be GpwResult+GpwHandles from a 268-line driver
// is now: the calculation (which owns the graph), its outcome, and the telemetry the campaign reports on.
// The whole recipe is stated in ONE block below instead of being assembled across MakeGpwAccelerator (449),
// RunGpw (608), RunGpwAnnealed (761) and this function -- four thousand lines for one run.
struct MnOArm
{
    std::shared_ptr<UnitCell>                cell;
    std::unique_ptr<qchem::SolidCalculation> calc;    //!< owns basis/Hamiltonian/iterator; outlives `result`
    qchem::Outcome<qchem::SolidCalculation::Converged, qchem::SCFFailure> result
        = qchem::Outcome<qchem::SolidCalculation::Converged, qchem::SCFFailure>::Fail(qchem::SCFFailure{});
    std::vector<FpRow>                       series;  //!< per-iteration telemetry, for Instrumentation()
    std::vector<double>                      stageKT; //!< each stage's kT, so the report can name its stages
};

//! THE CAMPAIGN INSTRUMENTATION, in ONE call (user's suggestion, 2026-08-25).  RunGpwAnnealed did two jobs
//! -- RUN and REPORT -- and only the first belongs in a library.  Splitting them is what let the facade take
//! the run while the fingerprint, the order trajectory and the magnetic diagnostics stayed in the test,
//! where a campaign's instruments belong.
void Instrumentation(const MnOArm& arm, const std::string& label)
{
    // PER STAGE, not over the concatenation.  The observer accumulates every stage into one series and the
    // iteration counter RESTARTS each stage, so a single fingerprint over the whole thing reads a
    // discontinuity as a trajectory -- and the campaign reads its stages separately anyway.  Split where
    // the iteration number goes backwards; that IS the stage boundary.
    size_t begin=0, nStage=0;
    for (size_t i=1; i<=arm.series.size(); ++i)
        if (i==arm.series.size() || arm.series[i].it<=arm.series[i-1].it)
        {
            const std::vector<FpRow> stage(arm.series.begin()+begin, arm.series.begin()+i);
            const std::string tag = label + " kT="
                                  + (nStage<arm.stageKT.size() ? std::to_string(arm.stageKT[nStage])
                                                               : std::string("?"));
            Fingerprint(stage, tag.c_str());
            OrderTrajectory(stage, "m_stag", tag.c_str());
            // The T3 instrument's raw time series (doc/OpenWork.md T3): Eee is the CHARGE channel's
            // health, and the threshold that convicts a runaway can only be calibrated against the
            // healthy arm's own trajectory -- so print it, per stage, beside the order.
            std::cout<<"["<<tag<<" Eee]";
            for (const auto& r : stage) std::cout<<" "<<std::fixed<<std::setprecision(3)<<r.eee;
            std::cout<<std::defaultfloat<<std::endl;
            begin=i; ++nStage;
        }
    // The LIBRARY's own post-mortem (doc/OpenWork.md T4): the detectors now live in SolidCalculation, so
    // this line is what any caller gets, not something this test computes for itself.
    std::cout<<"["<<label<<" diag] "<<arm.calc->Diagnostics().Summary()<<std::endl;
    if (!arm.result) { std::cout<<"["<<label<<"] NO ANSWER: "<<arm.result.Error().details<<std::endl; return; }
    // The INTEGRATED site moments (T2's instrument, the Becke basins) beside the historical point probe --
    // the honest observable, in electrons, rather than a spin DENSITY sampled at one guessed offset.
    std::cout<<"["<<label<<"] Etot="<<std::setprecision(10)<<arm.result->Energy()
             <<" charge="<<arm.result->TotalCharge()<<std::endl;
}

} // namespace qchem::tests::gpw
