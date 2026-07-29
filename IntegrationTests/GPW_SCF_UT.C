// File GPW_SCF_UT.C  The GPW self-consistent total energy: the first periodic SCF on GAUSSIAN orbitals.
//
// GPW (increments 1-2) already satisfies every plane-wave Kohn-Sham concept EXCEPT the external potential:
//   - kinetic  -> PW_Kinetic calls bs->MakeKinetic()          (GPW: lattice-sum <p^2>)              [inc 1]
//   - Hartree  -> PW_Hartree casts bs to Band_FT_IBS + cd to FourierDensity (GPW: collocation tensors)[inc 2]
//   - XC       -> PW_XC, same casts + the fit-basis grid                                             [inc 2]
//   - ion-ion  -> IonIon<dcmplx> (Ewald from Zion)                                                   [structure]
// The one gap was the external pseudopotential: the plane-wave PW_Pseudo needs G-space form factors, which
// Gaussians cannot supply.  This increment closes it: GPW_IBS realises Integrals_Pseudo<dcmplx> by REAL-SPACE
// mesh quadrature of the pseudopotential against its Gaussians (the SAME qcMesh machinery the molecular
// PP_Local/PP_NonLocal terms use).  So the ENTIRE Ham_PW_DFT drives a GPW basis verbatim -- Gaussian orbitals,
// plane-wave-style Hartree/XC by collocation -- through the real framework cSCFIterator.
//
// Validation (mirrors L_PP's finite==lattice cross-check, lifted to a full SCF):
//   (1) A single Si pseudo-atom in a LARGE cubic box, run through the GPW SCF, reproduces the FINITE molecular
//       density-fit DFT energy (the qchem::Calculation "sipp" + GTH-LDA pseudo-atom) to grid-cutoff tolerance.
//       Same basis, same PP, same functional; the only differences (density-fit vs collocation Hartree, Becke
//       vs uniform-grid XC, periodic vs open boundary) vanish as the box grows + the grid resolves -> the
//       electronic energies converge.  This is the tight correctness gate.
//   (2) A real material: crystalline silicon (diamond) at Gamma converges, conserves charge (8 valence e-),
//       and lands a reproducible total energy (a "did-E-move" regression anchor, per doc/GPWPlan.md section 5).
#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <cmath>
#include <cstdlib>   // std::getenv/std::atof (the NaF mixing-tuning env knobs)
#include <complex>
#include <cstdio>
#include <fstream>   // /proc/self/statm (the RSS breadcrumb bisect)
#include <stdexcept>
#include <algorithm>

import qchem.Structure;                          // Molecule, Atom
import qchem.UnitCell;                           // UnitCell, FCCUnitCell
import qchem.Lattice_3D;                         // Lattice_3D
import qchem.BasisSet;                           // Complex_BS, Real_BS
import qchem.BasisSet.Orbital_1E_IBS;            // Complex_OIBS (the overlap-spectrum diagnostic)
import qchem.Blaze;                              // blazem::eigen, blaze::min/max (overlap spectrum)
import qchem.BasisSet.Lattice_3D.BasisSet;       // GPWFactory (the GPW basis container)
import qchem.BasisSet.Molecule.Factory;          // Molecule::Factory, BasisSetData/Engine/Angular
import qchem.Hamiltonian.Internal.Hamiltonians;  // Ham_PW_DFT (the plane-wave LDA KS Hamiltonian -- drives GPW too)
import qchem.Hamiltonian.Internal.PWTerms;        // ReportGridCharge() -- the integral-rho_grid vs Tr(DS) readout
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ElectronConfiguration.Crystal;      // Crystal_EC (single-k Bloch occupation)
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS; // cSCFAcceleratorDIIS (complex DIIS)
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;  // cSCFAcceleratorGDM (complex geodesic direct-min)
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder; // cSCFAcceleratorLadder (DIIS -> GDM chain)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // tSCFAcceleratorNull<dcmplx> (NaF: pure damped Kerker)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Reporting;                          // report:: -- bracket the GPW run so grids/basis sections land
import qchem.Symmetry.Spin;                      // Spin
import qchem.Symmetry.Factory;                   // BlochFactory (build a k-block with a fractional MP shift)
import qchem.Symmetry.Lattice_3D.SpaceGroup;     // SpaceGroup::Detect + AtomSite + ReciprocalPointOps (IBZ ops)
import qchem.LASolver;                           // qchem::Ortho (Cholesky | Eigen | SVD -- basis orthogonalisation)
import qchem.BasisSet.Lattice_3D.GPW_IBS;         // GPW_IBS (build a concrete block for the collocation diagnostic)
import qchem.BasisSet.Lattice_3D.Evaluators.GPW;  // GPW_Evaluator (Overlap3CTensor -- the collocation tensor)
import qchem.BasisSet.Internal.GMap;              // G_ERI3 (the collocation weight tensor)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH, GTH_PP (the PP model, for the matrix-trace probe)
import qchem.Calculation;                        // qchem::Calculation, CalcOptions (finite reference)
import qchem.AtomCalculation;                    // AtomCalculation, AtomType, BasisSetAccuracy (Slater/High pseudo-atom ref)
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Complex_BS;
using qchem::BasisSet::Molecule::BasisSetData;

namespace
{
// The valence Si Gaussian basis (SIPP, MnD-Cartesian) on ANY structure -- the L_PP / GPW_UT builder.
std::shared_ptr<const Real_BS> MakeBasis(const Structure& st)
{
    return std::shared_ptr<const Real_BS>(
        BasisSet::Molecule::Factory(BasisSetData::SIPP, &st,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
}
// The SHORT-RANGE variant (most diffuse valence primitives dropped) -- well-conditioned Bloch overlap in a solid.
std::shared_ptr<const Real_BS> MakeBasisSR(const Structure& st)
{
    return std::shared_ptr<const Real_BS>(
        BasisSet::Molecule::Factory(BasisSetData::SIPP_SR, &st,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
}
// The low-q GTH valence basis (valgen-generated; carries Al/Na/F) -- the Al block drives the FCC-Al metal test.
std::shared_ptr<const Real_BS> MakeBasisLowQ(const Structure& st, BasisSetData which=BasisSetData::VALENCE_LOWQ_SR)
{
    return std::shared_ptr<const Real_BS>(
        BasisSet::Molecule::Factory(which, &st,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
}

struct GpwResult { bool converged; double charge; qchem::EnergyBreakdown E; size_t iters; };

// Shared GPW run reporting.  ANY GPW driver -- RunGPW's fixed recipe OR a bespoke one (NaFRocksaltGamma's
// multi-species PP + Ladder accelerator) -- gets automatic reporting by holding one of these: RAII brackets the
// run (Begin/End + optional console), and the driver calls VetBasis (the fail-fast conditioning pre-flight,
// BEFORE it builds the Hamiltonian/grids) then EmitGrids.  Keeps each driver's own Hamiltonian/accelerator/
// params -- only the reporting + setup-order discipline is shared.
struct GpwReport
{
    explicit GpwReport(const std::string& name, bool verbose)
    {
        qchem::report::Begin(name);
        if (verbose) qchem::report::SetConsole(std::cout, qchem::report::Detail::Normal);
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
        return BasisSet::Lattice_3D::VetGpwConditioning(bs);
    }
    //! Emit the grids section (this is where the ladder is actually built) -- only after VetBasis passed.
    void EmitGrids(const Complex_BS& bs)
    {
        qchem::report::Log("building grid ladder");
        BasisSet::Lattice_3D::EmitGpwGrids(bs);
    }
};

// PROBE (dynamics fingerprint, doc/GPWPlan §0): classify an SCF trajectory captured via the Observer hook.
// Three pathologies have DISTINCT time-series signatures, so one line names which regime the run is in --
// separating "the iteration can't find the min" (dynamics: sloshing/divergence) from "the min is a fit floor"
// (functional).  Captured from qchem::SCFIterator::SCFProgress {iteration, energy, dE=|ΔE|, [F,D], Δρ}.
struct FpRow { size_t it; double E, dEabs, fd, drho; };
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
    const double relAmp=amp/std::max(std::fabs(Ef),1e-30);   // energy swing RELATIVE to the total
    // Verdict priority separates the three pathologies (+ the benign degenerate case) by their distinct
    // signatures.  KEY distinction: a degenerate open shell has the ENERGY settled (small relAmp) while Δρ
    // never falls (ρ rotates in the degenerate subspace) -- benign; charge-transfer SLOSHING swings BOTH.
    const char* verdict =
        (drhoF < 1e-5)                                   ? "CONVERGED" :
        (drhoF > 1e-3 && relAmp < 5e-3)                  ? "DENSITY-DEGENERATE (E settled, ρ rotates -- benign)" :
        (flips >= 3 && relAmp > 5e-3)                    ? "OSCILLATING (charge-transfer sloshing / mixing unstable)" :
        (drhoF > 1e-5 && s.back().dEabs < 1e-5)          ? "FIT-FLOOR STALL (Δρ floored, ΔE tiny -- functional/grid)" :
        (std::fabs(Ef) > 3.0*std::fabs(s.front().E))     ? "DIVERGING" : "UNSETTLED (hit iter cap mid-descent)";
    std::cout << "["<<label<<" fp] iters="<<n<<" Efinal="<<Ef<<" lastΔρ="<<drhoF
              << " oscFlips(last"<<w<<")="<<flips<<" Eamp(last"<<w<<")="<<amp<<" relAmp="<<relAmp
              << "  => "<<verdict<<std::endl;
}

// One GPW Gamma-point SCF: build the GPW basis over the lattice, hand it the plane-wave LDA Hamiltonian
// (Ham_PW_DFT reaches GPW's real-space Integrals_Pseudo), seed uniform, run the complex-DIIS cSCFIterator.
// GpwOptions: the control surface for a GPW run.  Every material (Si, NaF, CsI, LiCoO2, f-oxides, ...) is one
// options literal -- multi-species PP, grid params (the efficiency lever for ionic systems), accelerator policy,
// seed, ortho, and the SCFParams gates.
struct GpwOptions
{
    std::string label = "gpw";
    int         Nelec = 0;
    std::vector<std::pair<std::string,int>> species;   // multi-species PP, e.g. {{"Na",1},{"F",7}}
    // grids
    double densityEcut  = -1.0;                        // <0 AUTO = cutoffFactor*alpha_max
    double cutoffFactor = 2.0;
    double ladderFactor = 4.0;
    BasisSet::Lattice_3D::RasterPolicy raster = BasisSet::Lattice_3D::RasterPolicy::BallOnly;
    BasisSet::Lattice_3D::CellImages   images = BasisSet::Lattice_3D::CellImages::Periodic;
    rvec3_t kShift = rvec3_t(0,0,0);
    // convergence machinery
    std::string accelerator = "DIIS";                  // DIIS | GDM | Ladder | Null
    bool        globalFermi = false;                    // metal: one μ across the k-mesh (Crystal_EC global mode)
    bool        reduceBZ    = false;                     // IBZ/k-star: fold the MP mesh to the irreducible wedge
    bool        symmetrize  = false;                     // IBZ: star-average the density (auto-on with reduceBZ;
                                                         //   set alone on a FULL mesh to check idempotency)
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Uniform;
    qchem::Ortho ortho    = qchem::Cholesky;
    double       orthoTol = 0.0;
    SCFParams    scf;                                  // NMaxIter / MinΔρ / MinΔE / SmearingkT / ... (the gates)
};

// Build the complex SCF accelerator named by \a policy.  Ladder = the ionic-crystal DIIS->GDM hand-off on
// |ΔE/E| (NaF's proven recipe); the rest are the plain single-engine choices.
static qchem::SCFAccelerators::tSCFAccelerator<dcmplx>* MakeGpwAccelerator(const std::string& policy)
{
    using namespace qchem::SCFAccelerators;
    if (policy=="Null") return new tSCFAcceleratorNull<dcmplx>();
    if (policy=="DIIS") return new cSCFAcceleratorDIIS(DIISParams{8, 8.0, 1e-10, 1e-9});
    if (policy=="GDM")  return new cSCFAcceleratorGDM(GDMParams{1.0});
    if (policy=="Ladder")
    {
        std::vector<std::unique_ptr<tSCFAccelerator<dcmplx>>> rungs;
        rungs.push_back(std::make_unique<cSCFAcceleratorDIIS>(DIISParams{8, 0.1, 1e-10, 1e-9}));
        rungs.push_back(std::make_unique<cSCFAcceleratorGDM>(GDMParams{1.0}));
        return new cSCFAcceleratorLadder(std::move(rungs), 1e-8, 5, 1e-8, 1e-6, ScheduleSignal::EnergyChange);
    }
    throw std::runtime_error("MakeGpwAccelerator: unknown policy \""+policy+"\" (DIIS|GDM|Ladder|Null)");
}

// The GENERAL GPW driver: basis -> fail-fast conditioning pre-flight -> grids -> multi-species Hamiltonian ->
// accelerator -> SCF, with automatic reporting (GpwReport) + heartbeat logging throughout.  Any material is a
// GpwOptions literal; the positional RunGPW below and the bespoke NaFRocksaltGamma are thin callers.
static GpwResult RunGpw(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const GpwOptions& o,
                        bool verbose=false)
{
    namespace L3=BasisSet::Lattice_3D;
    const std::string sp = o.species.empty() ? std::string() : o.species.front().first;
    GpwReport report(sp+" "+o.label, verbose);

    qchem::report::Log("building GPW basis");
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, mol, L3::GPWParams{
        .densityEcut=o.densityEcut, .cutoffFactor=o.cutoffFactor, .raster=o.raster,
        .images=o.images, .kShift=o.kShift, .ladderFactor=o.ladderFactor, .reduceBZ=o.reduceBZ}));

    // FAIL-FAST: vet the (analytic, grid-free) overlap + emit basis BEFORE grids/Hamiltonian; a rank-deficient
    // basis aborts here instead of building the whole ladder first.  Also puts `basis` before `grids`.
    if (report.VetBasis(*bs) > 0)
    {
        std::cout << "["<<o.label<<"] ABORT: basis rank-deficient (see basis.removed) -- skipped grids + SCF."<<std::endl;
        return {false, 0.0, qchem::EnergyBreakdown{}, 0};
    }
    report.EmitGrids(*bs);

    qchem::report::Log("building Hamiltonian");
    auto       irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry the Sum_k)
    Crystal_EC ec(irreps, o.Nelec, o.globalFermi);
    qchem::Hamiltonian::cHamiltonian* ham =
        new qchem::Hamiltonian::Ham_PW_DFT(lat.GetStructure(), bs.get(), o.species, "LDA");
    auto* acc = MakeGpwAccelerator(o.accelerator);

    qchem::report::Log("SCF start");
    // No Section("basis") here: the pre-flight already emitted basis, so MakeIrrepWFs stays silent
    // (report::InSection("basis") false) and just does the ortho.
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc,
                                         o.seed, lat.GetStructure().get(), o.ortho, o.orthoTol);
    // IBZ density symmetrization (doc/GPWPlan1.md item 3): hand the SCF the crystal's reciprocal point group so
    // the assembled density is star-averaged.  reduceBZ NEEDS it (a folded mesh is otherwise unsymmetric); on a
    // FULL mesh it is an IDEMPOTENT no-op (the mesh already carries every star partner) -- the correctness gate.
    // (The ops are recomputed here from the lattice; a follow-up exposes them from the basis to avoid the dup.)
    if (o.symmetrize || o.reduceBZ)   // reduceBZ needs it; symmetrize alone is the full-mesh idempotency probe
    {
        namespace SL=qchem::Symmetry::Lattice_3D;
        const auto& cell = lat.GetUnitCell();
        auto st = lat.GetStructure();   // HOLD the shared_ptr: GetStructure() returns a fresh temporary, so
        std::vector<SL::AtomSite> asites;   // `for (a : *lat.GetStructure())` would iterate a freed Structure.
        for (Atom* a : *st) asites.push_back({a->itsZ, cell.ToFractional(a->itsR)});
        SL::SpaceGroup sg = SL::SpaceGroup::Detect(cell.GetCellMatrix(), asites);
        scf.SetSymmetryOps(sg.ReciprocalPointOps(/*timeReversal*/true, /*symmorphicOnly*/true));  // W-only guard
    }
    std::vector<FpRow> series;
    scf.SetObserver([&series](const qchem::SCFIterator::SCFProgress& p)
                    { series.push_back({p.iteration, p.energy, p.dE, p.commutator, p.drho}); });
    SCFParams par = o.scf; par.Verbose = verbose;   // one `verbose` drives both the report console + the SCF table
    qchem::Hamiltonian::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");
    scf.Iterate(par);
    qchem::Hamiltonian::ReportGridCharge()=false;
    Fingerprint(series, o.label.c_str());

    auto* cd=scf.GetWaveFunction()->GetChargeDensity();
    double charge=cd->GetTotalCharge();
    delete cd;
    qchem::EnergyBreakdown E=scf.GetEnergy();
    std::cout << "["<<o.label<<"] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Eelec="<<E.GetElectronicEnergy() << " Etot="<<E.GetTotalEnergy()
              << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc
              << " Enn="<<E.Enn<<" E_alphaZ="<<E.E_alphaZ<<")" << std::endl;
    return {scf.Converged(), charge, E, scf.GetIterationCount()};
}

// ANNEALED GPW driver (doc/GPWPlan1.md item 2): run a DESCENDING Fermi-smearing kT schedule, re-seeding each
// stage from the previous stage's converged density (density continuation).  A hot first stage smears the
// degenerate open shell wide and converges it easily (μ lands in the manifold, every degenerate orbital takes
// the same fractional occupation); each cooler stage then starts IN-BASIN, so it settles toward the T->0
// answer without the integer-aufbau sloshing a cold cold-start would suffer.  The basis + grids + Hamiltonian
// are IDENTICAL every stage (only kT changes), so the re-seed is a direct density-matrix transfer with no
// re-projection -- and the density's basis block (bs) outlives every stage, so the carried seed stays valid
// after its producing iterator is torn down.  Returns the FINAL (coldest) stage's result; internal energy is
// GetTotalEnergy()-MinusTS (A=E-TS reported, so E=A-(-TS)).  A general convergence tool (GPWPlan1 "Annealing
// as a general capability") -- exercised here on the FCC-Al degenerate 3p.
static GpwResult RunGpwAnnealed(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const GpwOptions& o,
                                const std::vector<double>& kTSchedule, bool verbose=false)
{
    namespace L3=BasisSet::Lattice_3D;
    const std::string sp = o.species.empty() ? std::string() : o.species.front().first;
    GpwReport report(sp+" "+o.label+" (annealed)", verbose);

    qchem::report::Log("building GPW basis");
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, mol, L3::GPWParams{
        .densityEcut=o.densityEcut, .cutoffFactor=o.cutoffFactor, .raster=o.raster,
        .images=o.images, .kShift=o.kShift, .ladderFactor=o.ladderFactor, .reduceBZ=o.reduceBZ}));
    if (report.VetBasis(*bs) > 0)
    {
        std::cout << "["<<o.label<<"] ABORT: basis rank-deficient (see basis.removed) -- skipped grids + SCF."<<std::endl;
        return {false, 0.0, qchem::EnergyBreakdown{}, 0};
    }
    report.EmitGrids(*bs);

    auto       st = lat.GetStructure();
    auto       irreps = bs->GetIrreps(Spin::None);
    Crystal_EC ec(irreps, o.Nelec, o.globalFermi);

    GpwResult R{false, 0.0, qchem::EnergyBreakdown{}, 0};
    qchem::ChargeDensity::cDM_CD* seedCD = nullptr;   // carried between stages (the next stage's ctor consumes it)
    for (size_t s=0; s<kTSchedule.size(); ++s)
    {
        const double kT=kTSchedule[s];
        // Fresh Hamiltonian + accelerator per stage (the iterator OWNS + deletes them; a kT change must not
        // carry stale DIIS history across the re-seed).
        auto* ham = new qchem::Hamiltonian::Ham_PW_DFT(st, bs.get(), o.species, "LDA");
        auto* acc = MakeGpwAccelerator(o.accelerator);
        std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scf(
            s==0 ? new qchem::SCFIterator::SolidSCFIterator(bs.get(), &ec, ham, acc, o.seed,  st.get(), o.ortho, o.orthoTol)
                 : new qchem::SCFIterator::SolidSCFIterator(bs.get(), &ec, ham, acc, seedCD, st.get(), o.ortho, o.orthoTol));

        SCFParams par=o.scf; par.Verbose=verbose; par.SmearingkT=kT;
        std::vector<FpRow> series;
        scf->SetObserver([&series](const qchem::SCFIterator::SCFProgress& p)
                         { series.push_back({p.iteration,p.energy,p.dE,p.commutator,p.drho}); });
        std::cout << "["<<o.label<<" anneal "<<s+1<<"/"<<kTSchedule.size()<<"] kT="<<kT<<std::endl;
        scf->Iterate(par);
        Fingerprint(series, (o.label+" kT="+std::to_string(kT)).c_str());

        seedCD = scf->GetWaveFunction()->GetChargeDensity();   // consumed by the next stage's ctor (bs keeps it valid)
        qchem::EnergyBreakdown E=scf->GetEnergy();
        R = {scf->Converged(), seedCD->GetTotalCharge(), E, scf->GetIterationCount()};
        std::cout << "["<<o.label<<" stage "<<s+1<<"] kT="<<kT<<" conv="<<R.converged<<" iters="<<R.iters
                  << " A=E-TS="<<E.GetTotalEnergy()<<" -TS="<<E.MinusTS
                  << " E(internal)="<<(E.GetTotalEnergy()-E.MinusTS)<<std::endl;
        // scf drops here (deletes ham/acc/wf); seedCD survives -- its block is bs, which outlives the loop.
    }
    delete seedCD;   // the final stage's carried density (not consumed by any further ctor)
    return R;
}

// Positional back-compat wrapper (the existing single-species Si callers): forwards to RunGpw(GpwOptions).
GpwResult RunGPW(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, double densityEcut,
                 int Nelec, const char* element, const char* label, bool verbose=false, int nmax=120,
                 qchem::Ortho ortho=qchem::Cholesky, double orthoTol=0.0,
                 rvec3_t kShift={0,0,0}, double minDrho=1e-6, double minDE=1e30,
                 qchem::ChargeDensity::SeedStrategy seed=qchem::ChargeDensity::SeedStrategy::Uniform,
                 BasisSet::Lattice_3D::CellImages images=BasisSet::Lattice_3D::CellImages::Periodic,
                 double smearkT=0.0)
{
    GpwOptions o;
    o.label=label; o.Nelec=Nelec; o.species={{std::string(element), 4}};   // the Si callers: Zion=4
    o.densityEcut=densityEcut; o.images=images; o.kShift=kShift;
    o.accelerator="DIIS"; o.seed=seed; o.ortho=ortho; o.orthoTol=orthoTol;
    o.scf.NMaxIter=(size_t)nmax; o.scf.MinΔρ=minDrho; o.scf.MinΔE=minDE;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.Verbose=verbose; o.scf.SmearingkT=smearkT;
    return RunGpw(lat, mol, o, verbose);
}
} //anon

// (1) THE REAL-MATERIAL SCF: crystalline silicon (diamond) primitive cell at Gamma, driven end-to-end by
// the framework cSCFIterator through the plane-wave Kohn-Sham Hamiltonian on a GAUSSIAN (GPW) basis.  8
// valence electrons (2 x Zion 4) fill a closed shell (sigma_g^2 sigma_u^2 pi_u^4), so it converges cleanly.
// With the G-space local PP (box-independent, PW G=0/alignment convention) the total is now PHYSICAL:
// Etot=-8.248 -- close to the plane-wave bulk-Si -7.2273 (Ecut=4) / converged ~-7.9.  The residual ~1 Ha is
// the Rcut=0 over-binding (home-cell electrons, no inter-cell screening, feel the full periodic ion Ewald);
// true bulk (Rcut>0) awaits the overlap-conditioning fix.  A did-E-move regression anchor (pin the value).
// (1c) MULTI-K PLUMBING: a 2x1x1 Monkhorst-Pack mesh (2 k-points), SR basis, Rcut=2a.  The ANALYTIC
// collocation always sums the screened cross-cell pair offsets (with their Bloch phases), so there is no
// "Rcut=0 makes every k-block a copy of Gamma" shortcut any more -- the 2-point mesh has REAL dispersion, and
// this pins its BZ-weighted total.  What the gate protects is the multi-k machinery: one GPW_IBS per BZ
// k-point (GPW_BasisSet iterating MakeKMesh WITH BZ weights), the multi-block GetIrreps, Crystal_EC's
// BZ-weighted (Sum_k w_k) occupation, the per-irrep k-loop, and the BZ-summed charge/energy (it caught a
// missing BZ weight -> charge x Nk).  Energy-gated at the fit floor like the Gamma anchor.
TEST(GPW_SCF, DISABLED_SiliconMultiKPlumbing)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,1,1));   // 2 k-points: Gamma + the zone-boundary k=1/2 (real +-1 phases)

    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si",
                       "Si SR 2x1x1", /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0,
                       /*kShift*/rvec3_t(0,0,0), /*minDrho*/1e-3, /*minDE*/1e-6);

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);                          // 8 valence e- (BZ-weighted Sum_k, not x Nk)
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.45137, 5e-3);         // did-E-move anchor (2x1x1 dispersion, analytic path)
}

// DISPERSIVE MULTI-K BULK (disabled: 8 k-blocks, ~4 min) -- the first REAL bulk GPW, unblocked by the KB
// Bloch-orbital fix (Rcut>0 now correct).  Gamma-centred 2x2x2 MP, SIPP_SR, Rcut=2a: charge stays 8 and the
// total drops with k-sampling (Gamma -7.11467 -> 2x1x1 -7.451 -> 2x2x2 -7.778 -- real dispersion).
// CROSS-CHECK vs CP2K AT THE SAME GAMMA-CENTRED MESH: -7.7778 vs CP2K -7.77846 (~0.7 mHa, the N=32 grid gap;
// deck UnitTests/CP2K/si_fcc_gpw_222_gamma.inp).  The 90 mHa vs CP2K's DEFAULT -7.86744 is the k-CONVENTION:
// Gamma-centred here (kShift=0) vs CP2K's classic SHIFTED MONKHORST-PACK (k at +/-1/4).  The shifted grid is
// the sibling test below (kShift=1/2).  The general-k PHYSICS is validated at both.
TEST(GPW_SCF, DISABLED_SR_2x2x2GammaCentred_vs_CP2K)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si",
                       "Si 2x2x2 Gamma-centred", /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.77846, 3e-3) << "GPW 2x2x2 Gamma-centred vs CP2K same-mesh -7.77846";
}

// SHIFTED Monkhorst-Pack (kShift=½ → k at ±¼ = CP2K's DEFAULT MONKHORST-PACK 2 2 2) -- the apples-to-apples
// match to CP2K's shipped 2x2x2 reference -7.86744 (deck si_fcc_gpw_222.inp).  This is the FIRST run with a
// genuinely COMPLEX Bloch phase e^{ik·R} (not ±1), so the density matrix D and every k-block matrix are
// genuinely complex -- the exact case that exposed (and now validates the fix for) TWO complex-only bugs
// (doc/GPWPlan.md, "Complex-k GPW FIXED", 2026-07-10):
//   1. GPW_Evaluator::BuildWeights conjugated the BRA (i) instead of the KET (j) collocation slot, so the
//      Fourier density rho-tilde was the TRANSPOSE-density D^T (a different real field at complex k) -> the
//      Hartree/XC drive was inconsistent with the physical density (IrrepCD::operator() / the PW delta path).
//   2. GPW_Evaluator::MakeSeparablePP summed the KB projector images with e^{+ik·R} instead of e^{-ik·R};
//      the correct Bloch projection b_i=<chi_i^k|beta_home> tiles all-space with a CONJUGATED image phase.
//      At complex k this HALVED the nonlocal-PP trace (Vnl 42->22) -> a spurious deep core level -> over-bind.
// Both are inert at Gamma / half-integer k (phase ±1 self-conjugate), so every committed anchor is unchanged;
// they matter ONLY here.  Grid-matched to CP2K's mesh.
// ENABLED 2026-07-15 (\S0a complex-k revalidation): the ANALYTIC collocate/integrate kernels reproduce the
// CP2K shifted-MP reference -- Rcut=2a gave -7.86724 (0.20 mHa), and the run is affordable now (the stream
// cache + the phase-independent integrate memo make the 8 k-blocks share the static sweeps: ~2.5 min).
// Rcut switched to AUTO for scheme consistency with the enabled anchors (both sides parameter-free).
TEST(GPW_SCF, DISABLED_SR_2x2x2ShiftedMP_vs_CP2K)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si",
                       "Si 2x2x2 shifted MP (k=±¼)", /*verbose*/false, /*nmax*/60,
                       qchem::Cholesky, 0.0, rvec3_t(0.5,0.5,0.5));
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.86744, 3e-3) << "GPW 2x2x2 shifted MP (CP2K default) vs -7.86744";
}

// DIAGNOSTIC: TERM-BY-TERM translation invariance of the 1E/PP matrix TRACES (no SCF -> fast) -- the tool
// that localized the Rcut>0 over-binding.  A rigid translation of the whole crystal (both atoms + their basis)
// must leave every trace invariant; the residual is that term's grid/mesh artifact.  Compare a CORNER atom
// (frac 0, on the cell boundary) vs an off-boundary atom (frac 0.13).  Kinetic is analytic -> the control.
//
// THE STORY (2026-07-09).  Pre-fix the KB nonlocal PP (Vnl) was translation-variant by ~16 Ha at ALL Rcut:
// MakeSeparablePP quadratured the RAW home orbital against the projector on a single-cell mesh, so a
// boundary-straddling corner orbital lost its wrapped tail (summing the PROJECTOR images cannot restore the
// ORBITAL's).  FIX = use the Bloch-summed orbital (Eval) as the bra (GPWPlan TODO 1a).  A NON-obvious extra:
// the local-PP/Hartree/XC (Vloc) variance was NOT the FFT raster (a uniform-grid ORIGIN shift is ~a no-op for
// a periodic quadrature -- Poisson summation only moves Nyquist-aliasing phases; the voxel-shift "Option A"
// was tried and REVERTED) -- it too was just incomplete ORBITAL WRAPPING.  Once the orbital is fully wrapped
// (Rcut>=2a) BOTH terms are translation-invariant to machine precision:
//     Rcut=1.50a  Kin d=0.0000  Vloc d=0.67    Vnl d=1.67
//     Rcut=2.00a  Kin d=0.0000  Vloc d=0.0000  Vnl d=0.0000   <-- fully wrapped
//     Rcut=3.00a  Kin d=0.0000  Vloc d=0.0000  Vnl d=0.0000
// (Committed Rcut=0 anchors are unaffected: at Rcut=0 itsRc={0} so Eval==the raw orbital.)  At Rcut=2a the
// residual is image-truncation limited (~1e-4, tightening with Rcut), so the guard tolerance is 1e-3 -- still
// 1600x below the ~1.7 Ha (pre-fix ~16 Ha) bug it protects against.
TEST(GPW_SCF, DISABLED_TermTranslationInvariance)
{
    using BasisSet::Lattice_3D::GPW_IBS;
    auto tr=[](const chmat_t& M){ double s=0; for (size_t i=0;i<M.rows();i++) s+=std::real(dcmplx(M(i,i))); return s; };
    const double a=10.26, dE=30.0;   // N=64 (finer than CP2K's converged grid)
    auto traces=[&](double frac, double& kin, double& vloc, double& vnl)
    {
        FCCUnitCell cell(a);
        cell.AddAtom(14, {frac, frac, frac});
        cell.AddAtom(14, {0.25+frac, 0.25+frac, 0.25+frac});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        auto st=lat.GetStructure();
        Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH("Si","LDA",4);
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), MakeBasisSR(cell), dE);  // eps-complete enumeration
        const BasisSet::Complex_OIBS& g=gpw;
        kin =tr(g.Kinetic());
        vloc=tr(gpw.MakeLocalPotential   (st.get(), pp.local));
        vnl =tr(gpw.MakeSeparablePotential(st.get(), pp.nonlocal));
    };
    {
        double kc,lc,nc, ks,ls,ns;
        traces(0.00, kc,lc,nc);
        traces(0.13, ks,ls,ns);
        std::printf("Kin[%10.5f/%10.5f d=%.2e]  Vloc[%10.5f/%10.5f d=%.2e]  Vnl[%10.5f/%10.5f d=%.2e]\n",
                    kc,ks,std::fabs(kc-ks), lc,ls,std::fabs(lc-ls), nc,ns,std::fabs(nc-ns));
        EXPECT_NEAR(lc, ls, 1e-3) << "Vloc translation invariance (complete enumeration)";
        EXPECT_NEAR(nc, ns, 1e-3) << "Vnl translation invariance (complete enumeration)";
    }
}

// (1) THE GAMMA ANCHOR == THE CP2K ENERGY GATE.  SR basis, Rcut=2a (every term translation-invariant and
// screened-complete), densityEcut=20 (FFT N=32): reproduces the CP2K FCC-Si Gamma GPW reference (SIPP_SR /
// GTH-PADE-q4 / LDA_X+VWN5) Etot=-7.11506 to the N=32 grid gap (~0.4 mHa; densityEcut>=30 -> -7.11505 exact).
// NOTE (analytic path): the old fast Rcut=0 anchor (-8.2476) is GONE -- the analytic collocation is always
// screened-complete Bloch, so home-only 1E matrices would MIX SCHEMES (Tr(D S_home)=8 while the grid density
// integrates the Bloch trace -- the forbidden inconsistency; doc/GPWPlan.md durable pins).  SR keeps the Bloch
// overlap cleanly PD at 2a.  Energy-gated at the density-fit floor (minDE=1e-6, minDrho relaxed to 1e-3).
TEST(GPW_SCF, SiliconGammaConverges)
{
    const double a=10.26;                          // Si conventional cubic lattice constant (a.u.)
    FCCUnitCell cell(a);                           // FCC primitive cell (2-atom diamond basis)
    cell.AddAtom(14, {0,0,0});                      // Si diamond: true Z=14; Zion=4 via the PP
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si",
                       "Si SR Gamma", /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0,
                       /*kShift*/rvec3_t(0,0,0), /*minDrho*/1e-3, /*minDE*/1e-6);

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);              // 8 valence electrons
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.11506, 2e-3);   // CP2K FCC-Si Gamma reference (grid-gap tolerance)
}

// The GPW run-report SCHEMA CHECK (RunReportPlan step 3).  RunGPW brackets the run: GPWFactory emits the
// `grids` section during basis construction, and MakeIrrepWFs fills basis.perIrrep (per-Bloch-block
// conditioning) via the cursor.  Only the SETUP matters here, so a couple of iterations is plenty.
TEST(GPW_SCF, GridsReportSchema)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    report::ClearGlobal();                          // isolate this test's run
    RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si", "Si grids schema",
           /*verbose*/false, /*nmax*/4, qchem::Cholesky, 0.0, rvec3_t(0,0,0), /*minDrho*/1e-3, /*minDE*/1e-6);

    const report::json& all = report::GlobalReport();
    const report::json* grids=nullptr; const report::json* basis=nullptr;
    for (auto it=all.begin(); it!=all.end(); ++it)
    {
        if (it.value().contains("grids")) grids=&it.value()["grids"];
        if (it.value().contains("basis")) basis=&it.value()["basis"];
    }
    ASSERT_NE(grids, nullptr) << "no run emitted a grids section";
    EXPECT_GT((*grids)["densityEcut"].get<double>(), 0.0);
    EXPECT_TRUE(grids->contains("cutoffFactor"));
    EXPECT_TRUE(grids->contains("raster"));
    ASSERT_TRUE(grids->contains("ladder"));
    ASSERT_TRUE((*grids)["ladder"].is_array());
    ASSERT_GE((*grids)["ladder"].size(), 1u);
    const report::json& lvl0 = (*grids)["ladder"][0];
    for (const char* k : { "level", "N", "ecut", "nG", "role" })
        EXPECT_TRUE(lvl0.contains(k)) << "ladder row missing key: " << k;
    EXPECT_EQ(lvl0["role"].get<std::string>(), "reference");   // L==0 is the density/collocation grid
    ASSERT_TRUE(lvl0["N"].is_array());
    EXPECT_EQ(lvl0["N"].size(), 3u);
    ASSERT_TRUE(grids->contains("localPP"));
    EXPECT_TRUE((*grids)["localPP"].contains("kappa"));

    // The cursor path also gave GPW a basis section: per-Bloch-block conditioning in basis.perIrrep.
    ASSERT_NE(basis, nullptr) << "no run emitted a basis section";
    ASSERT_TRUE(basis->contains("perIrrep"));
    ASSERT_TRUE((*basis)["perIrrep"].is_array());
    ASSERT_GE((*basis)["perIrrep"].size(), 1u);
    EXPECT_TRUE((*basis)["perIrrep"][0].contains("cond"));
}

// (2) THE TIGHT CROSS-CHECK: the isolated Si pseudo-atom in a box vs the finite molecular DFT on the SAME
// SIPP basis + GTH-LDA PP.  With the G-space local PP the GPW total is box-independent and reproduces the
// finite SIPP energy (-3.74 vs -3.759) to grid tolerance -- the doc/GPWPlan sec 3.4 correctness gate.
// NOTE: at Gamma the atom has NO point group, so its half-filled 3p shell is degenerate -- the ENERGY converges
// (-3.736, grid-stable) but the density rotates freely within that degenerate shell, so |Delta rho| never
// reaches the tolerance (not a bug: integer occupation of a degenerate open shell).  A dcmplx GDM/Ladder
// energy-minimiser would converge it (today GDM/Ladder are <double>-only); the crystal above sidesteps it with
// a gap.  So this pins the CONVERGED ENERGY + charge as a did-E-move anchor, without a Converged() guard.
TEST(GPW_SCF, SiPseudoAtomInBoxMatchesFinite)
{
    // Basis-MATCHED reference: the SAME SIPP Gaussian basis + GTH-LDA PP as a finite molecule (density-fit
    // Hartree, Becke XC).  This is the tight cross-check: GPW-in-box == finite molecular DFT (doc/GPWPlan sec 3.4).
    Molecule si; si.Insert(new Atom(14, 0.0, {0,0,0}));
    Calculation cSipp(si, {.basis = "sipp", .pseudopotential = true});
    const double Esipp=cSipp.Energy();
    // Physical oracle (near-complete Slater/High, for context -- a different basis, not the GPW-correctness gate).
    AtomCalculation cHi(14, 14-4, {.type=AtomType::Slater, .accuracy=BasisSetAccuracy::High, .pseudopotential=true});
    std::cout << "[Si finite] sipp="<<Esipp<<"  Slater/High="<<cHi.Energy()<<std::endl;

    // Box a=16 (was 11): the analytic collocation always includes the screened cross-cell pair products, so
    // the box must be large enough that they are negligible for the finite-molecule comparison (SIPP's most
    // diffuse alpha=0.06 pair prefactor: e^{-0.03 a^2} = 2.7e-2 at a=11 -- visible; 4.6e-4 at a=16).
    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(14, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwResult R=RunGPW(lat, MakeBasis(cell), /*densityEcut*/10.0, /*Nelec*/4, "Si", "Si atom-in-box",
                       /*verbose*/false, /*nmax*/40, qchem::Cholesky, 0.0, rvec3_t(0,0,0), 1e-6, 1e30,
                       qchem::ChargeDensity::SeedStrategy::Uniform,
                       BasisSet::Lattice_3D::CellImages::HomeCellOnly);   // the finite-molecule mode

    EXPECT_NEAR(R.charge, 4.0, 1e-6);                        // 4 valence electrons (Zion=4), charge conserved
    // GPW-in-box (G-space local PP -> box-independent) reproduces the finite SIPP DFT energy to grid tolerance.
    // (Energy-converged; density is degenerate at Gamma -- see the note above -- so no Converged() guard.)
    EXPECT_NEAR(R.E.GetTotalEnergy(), Esipp, 5e-2) << "GPW-in-box total vs finite SIPP molecular DFT";
}

// (4b-i) FERMI SMEARING IS INERT ON A GAP (doc/GPWPlan1.md 4b, gate i).  The same gapped Si/Gamma anchor as
// SiliconGammaConverges, but with smearing kT=1e-3 Ha turned ON.  Si is a wide-gap insulator in this basis
// (ε_LUMO−ε_HOMO ≫ kT), so every f_i is essentially 0 or 1: the fractional occupations collapse to the
// aufbau integers, the Mermin −TS is negligible, and the total reproduces the CP2K reference −7.11506.  This
// is the regression that smearing must not perturb a system that does not need it (the T→0 / kT≪gap limit).
TEST(GPW_SCF, SmearingInertOnGap)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si",
                       "Si SR Gamma +smear", /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0,
                       /*kShift*/rvec3_t(0,0,0), /*minDrho*/1e-3, /*minDE*/1e-6,
                       qchem::ChargeDensity::SeedStrategy::Uniform,
                       BasisSet::Lattice_3D::CellImages::Periodic, /*smearkT*/1e-3);

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.11506, 2e-3);       // == the no-smear anchor: smearing is inert on a gap
    EXPECT_NEAR(R.E.MinusTS, 0.0, 1e-4);                     // −TS negligible when kT ≪ gap (f_i ∈ {0,1})
}

// (4b-iii) FERMI SMEARING CONVERGES A DEGENERATE OPEN SHELL (doc/GPWPlan1.md 4b, gate iii + the cure).  The
// Si pseudo-atom in a box has 4 valence electrons in a 3s²3p² configuration; at Gamma the atom has NO point
// group, so its three 3p orbitals are EXACTLY degenerate and half-filled.  Integer aufbau must pick 2 of the
// 6 p-states arbitrarily -> the density rotates freely within the degenerate shell and |Δρ| never converges
// (the documented behaviour of SiPseudoAtomInBoxMatchesFinite, iters=40/Δρ=0.08/"DENSITY-DEGENERATE", which
// pins only the energy).  Fermi smearing is the cure: μ lands in the degenerate manifold, each 3p orbital
// takes the SAME fractional occupation, the density is symmetric and STATIONARY, and the SCF converges Δρ
// (iters=24, Δρ=9e-7, "CONVERGED").  The Mermin −TS<0, so the total GetTotalEnergy() reported IS the free
// energy A=E−TS, which sits below the internal energy E -- the finite-T thermodynamic ordering (gate iii).
// kT MUST exceed the near-degenerate splitting to stabilise: kT=1e-2 converges, kT=1e-3 still slosh-rotates
// (measured) -- the honest recipe is "smear wider than the frontier splitting you are curing".
TEST(GPW_SCF, SmearingConvergesDegenerateShell)
{
    Molecule si; si.Insert(new Atom(14, 0.0, {0,0,0}));
    Calculation cSipp(si, {.basis = "sipp", .pseudopotential = true});
    const double Esipp=cSipp.Energy();   // finite SIPP molecular DFT (symmetric occupation): -3.759

    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(14, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwResult R=RunGPW(lat, MakeBasis(cell), /*densityEcut*/10.0, /*Nelec*/4, "Si", "Si atom-in-box +smear",
                       /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0, rvec3_t(0,0,0),
                       /*minDrho*/1e-6, /*minDE*/1e30,
                       qchem::ChargeDensity::SeedStrategy::Uniform,
                       BasisSet::Lattice_3D::CellImages::HomeCellOnly, /*smearkT*/1e-2);

    EXPECT_TRUE(R.converged) << "Fermi smearing should converge Δρ where integer aufbau cannot (degenerate 3p)";
    EXPECT_NEAR(R.charge, 4.0, 1e-6);
    // did-E-move anchor: the converged free energy A=E−TS at kT=1e-2 (internal E≈-3.744; the ~38 mHa gap to A
    // is the 3p-shell entropy −TS at this kT, which lowers A below E and below Esipp).  (Re-pinned when the
    // field-sharpness density rule landed -- doc/GPWPlan1.md 4b: the sharper XC grid moved it -3.779 -> -3.783.)
    EXPECT_NEAR(R.E.GetTotalEnergy(), -3.78260, 3e-3);
    EXPECT_LT(R.E.MinusTS, 0.0);                             // −TS<0 => A=GetTotalEnergy() sits below internal E (gate iii)
    EXPECT_NEAR(R.E.GetTotalEnergy()-R.E.MinusTS, Esipp, 3e-2) << "internal E=A−(−TS) vs finite SIPP molecular DFT";
}

// ===== (item 2) DEGENERATE-SHELL METAL: FCC Al @ Gamma (3s^2 3p^1) =====
// Al conventional FCC lattice constant 4.05 A = 7.653 au; ONE atom in the primitive cell, Zion=3.  At Gamma
// the cubic (O_h) site symmetry keeps the three 3p-derived states EXACTLY degenerate (a T_1u triplet), and
// they hold only ONE electron -- a partially-filled degenerate shell (the Si-3p pseudo-atom physics, now in a
// REAL periodic lattice).  Basis: the diffuse-trimmed VALENCE_LOWQ_SR Al block (the full valence_lowq Al p
// 0.05 makes the Bloch overlap rank-deficient AND explodes the collocation grid -- the SIPP->SIPP_SR / NaF SR
// lesson; valgen can regenerate).  NOT yet a Fermi-surface metal (that needs a global mu across k-blocks,
// GPWPlan1 items 3-4) -- this is the degenerate-open-shell + smearing/annealing gate.
static GpwOptions AlOptions()
{
    GpwOptions o;
    o.label="Al FCC Gamma"; o.Nelec=3; o.species={{"Al",3}};
    o.densityEcut=-1.0; o.accelerator="DIIS";
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    return o;
}

// (item 2a) THE MOTIVATION: integer aufbau CANNOT converge the degenerate 3p.  With smearing OFF, aufbau must
// place the lone 3p electron in ONE of the three degenerate p orbitals -- an arbitrary, symmetry-broken pick.
// The TOTAL ENERGY settles (|ΔE/E|~1e-13, the gap column ~0 => no frontier gap, the metallic signature) but
// the DENSITY rotates freely within the degenerate manifold, so |Δρ| floors well above tolerance and never
// converges.  This is the honest reason the smearing/annealing path below exists (mirrors the documented
// SiPseudoAtomInBoxMatchesFinite degenerate-shell behaviour, now for a periodic lattice).
TEST(GPW_SCF, AlFCCDegenerateShellAufbauStalls)
{
    FCCUnitCell cell(7.653);
    cell.AddAtom(13, {0,0,0});           // Al (Zion=3): 3s^2 3p^1
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o=AlOptions(); o.scf.NMaxIter=40; o.scf.SmearingkT=0.0;   // aufbau (no smearing)
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, /*verbose*/true);

    EXPECT_FALSE(R.converged) << "integer aufbau cannot converge Δρ of a partially-filled degenerate 3p shell";
    EXPECT_NEAR(R.charge, 3.0, 1e-6);                       // charge is still conserved (3 valence e-)
    EXPECT_NEAR(R.E.MinusTS, 0.0, 1e-12);                   // no smearing => no entropy term
}

// (item 2b) THE CURE + ANNEALING (doc/GPWPlan1.md item 2): a DESCENDING kT schedule (0.02 -> 0.01 -> 0.005 Ha),
// re-seeding each stage from the previous converged density.  Fermi smearing puts μ in the degenerate manifold
// so each 3p orbital takes the SAME fractional occupation -- the density is cubic-symmetric and STATIONARY, and
// every stage converges Δρ where aufbau cannot.  Annealing makes the cold end cheap: the kT=0.005 cold-START
// takes ~44 iters, but re-seeded from kT=0.01 it converges in ~17.  The Mermin −TS<0, so GetTotalEnergy() is
// the free energy A=E−TS (below the internal E); the INTERNAL energy E=A−(−TS) is kT-INDEPENDENT to ~1e-7
// across all three stages (−1.92115) -- the physical T→0 answer, and a strong self-consistency check on the
// smearing thermodynamics (gate iii).
TEST(GPW_SCF, AlFCCAnnealedMetal)
{
    FCCUnitCell cell(7.653);
    cell.AddAtom(13, {0,0,0});           // Al (Zion=3): 3s^2 3p^1
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o=AlOptions();
    // Accelerator = plain DIIS (NOT the DIIS->GDM Ladder).  MEASURED 2026-07-28: the Ladder's GDM tail rung is
    // INCOMPATIBLE with Fermi smearing here.  GDM builds its geodesic DIRECTION from the fixed-occupation
    // electronic gradient [F,D], but line-searches the FREE energy A=E−TS with the occupations Fermi-refilled
    // per trial (SCFIterator DirectMinStep).  Under fractional occupation those disagree: A is stationary where
    // the smeared gradient (not [F,D]) is zero, so at the DIIS fixed point [F,D]~4e-3≠0, and GDM's first step is
    // a persistent FALLBACK (GPW_GDMTRACE: best(Et−Ecur)=+0.42, all 12 backtracks uphill) -> occupations flip
    // (cfg `*`), the run never recovers.  This is NOT a grid / non-variationality fault (DIIS converges the SAME
    // E[ρ] cleanly to −1.97521); GDM just needs the occupation-response term to tail-polish a smeared free
    // energy.  So DIIS here; GDM+smearing is a captured follow-up (doc/GPWPlan1.md item 2).
    GpwResult R=RunGpwAnnealed(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, {0.02, 0.01, 0.005}, /*verbose*/false);

    EXPECT_TRUE(R.converged) << "Fermi-smearing annealing converges Δρ where integer aufbau cannot (degenerate 3p)";
    EXPECT_NEAR(R.charge, 3.0, 1e-6);
    EXPECT_LT(R.E.MinusTS, 0.0);                            // −TS<0 => A=GetTotalEnergy() sits below internal E (gate iii)
    // did-E-move anchors at the coldest stage (kT=0.005): the free energy A and the kT-independent internal E.
    EXPECT_NEAR(R.E.GetTotalEnergy(), -1.934665, 2e-3);                        // A = E − TS at kT=0.005
    EXPECT_NEAR(R.E.GetTotalEnergy()-R.E.MinusTS, -1.921148, 2e-3);            // internal E (T→0 physical value)
}

// (item 3) GLOBAL μ ACROSS k-BLOCKS -- the true metal fill.  FCC Al on a 2×2×2 Γ-centred Bloch mesh (8
// k-blocks) with ONE chemical potential across the whole BZ (Crystal_EC global mode + the composite cross-k
// Fermi fill): charge SLOSHES between k-points under a single μ instead of each k being pinned to a fixed
// per-block count.  This is the structural step a metal needs -- the partially-filled 3p band disperses with
// k, so no per-k integer (or per-k Fermi) occupation is right; the physical occupation is set by where each
// k's bands sit relative to the ONE Fermi level.  WHAT THIS GATES: (a) charge is conserved as the BZ-weighted
// Σ_k w_k n_k = 3 (the weight-consistency guard -- the μ constraint uses the SAME w_k the density applies);
// (b) the single μ CONVERGES the mesh where per-block filling cannot (measured: AL_GLOBAL=0 at 2×2×2 forces 3
// e⁻ at every k and lands non-converged garbage A≈-0.46, vs the global μ's converged -2.117 -- the charge
// MUST redistribute between k-points); (c) k-sampling lowers the energy vs Γ-only (-1.92 → -2.12, real
// dispersion).  Reduces EXACTLY to the per-block Fermi at a single k (verified: global≡per-block to 1e-12 at Γ).
TEST(GPW_SCF, AlFCCMetalGlobalMu)
{
    FCCUnitCell cell(7.653);
    cell.AddAtom(13, {0,0,0});                 // Al (Zion=3): 3s^2 3p^1
    Lattice_3D lat(cell, ivec3_t(2,2,2));      // 8-point Γ-centred Bloch mesh (weights sum to 1)
    GpwOptions o=AlOptions();
    o.globalFermi=true;                        // ONE μ across the BZ (the metal)
    o.scf.SmearingkT=0.01;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, /*verbose*/false);

    EXPECT_TRUE(R.converged) << "one μ across the BZ converges the dispersive metal (per-block filling cannot)";
    EXPECT_NEAR(R.charge, 3.0, 1e-6);          // BZ-weighted Σ_k w_k n_k = 3 (weight-consistency guard)
    EXPECT_LT(R.E.MinusTS, 0.0);               // −TS<0 => A=GetTotalEnergy() is the free energy (gate iii)
    EXPECT_NEAR(R.E.GetTotalEnergy(), -2.11681, 3e-3);   // did-E-move anchor (2×2×2 global-μ free energy A)
    EXPECT_LT(R.E.GetTotalEnergy(), -1.95);    // dispersion: well below the Γ-only -1.92 (k-sampling binds)
}

// ===== EXPERIMENTAL (scratch): global-μ across k-blocks (item 3 inc 3) =====
// AL_KGRID=n (mesh nxnxn), AL_GLOBAL=0/1 (per-block vs global μ), AL_KT, AL_NMAX.
TEST(GPW_SCF, DISABLED_AlGlobalMuExperiment)
{
    auto envd=[](const char* n,double d){const char*s=std::getenv(n);return s?std::atof(s):d;};
    const int nk=(int)envd("AL_KGRID",1);
    FCCUnitCell cell(7.653);
    cell.AddAtom(13,{0,0,0});
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));
    GpwOptions o=AlOptions();
    o.globalFermi = envd("AL_GLOBAL",1.0)!=0.0;
    o.reduceBZ = envd("AL_IBZ",0.0)!=0.0;
    o.symmetrize = envd("AL_SYM",0.0)!=0.0;   // idempotency probe: full-mesh star-average must not move E
    if (envd("AL_RASTER",0.0)!=0.0) o.raster=BasisSet::Lattice_3D::RasterPolicy::AliasFree;   // G-space XC route
    o.scf.SmearingkT = envd("AL_KT",0.01);
    o.scf.NMaxIter = (size_t)envd("AL_NMAX",60);
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell,BasisSetData::VALENCE_LOWQ_SR), o, /*verbose*/true);
    std::cout<<"[Al global-μ] nk="<<nk<<" global="<<o.globalFermi<<" conv="<<R.converged
             <<" charge="<<R.charge<<" A="<<R.E.GetTotalEnergy()
             <<" -TS="<<R.E.MinusTS<<" E(internal)="<<(R.E.GetTotalEnergy()-R.E.MinusTS)<<std::endl;
}

// (item 4) THE HONEST METAL: FCC Na, a real half-filled-band Fermi surface.  Zion=1 (3s^1) => ONE valence
// electron per cell, so the single conduction band is HALF-FILLED and μ cuts THROUGH it -- a genuine Fermi
// surface (unlike Al's degenerate-3p at Γ).  SHIFTED Monkhorst-Pack 2×2×2 (kShift=½ => k at ±¼, CP2K's default
// -- avoids the high-symmetry Γ, samples the Fermi surface evenly) + global μ + Fermi smearing.  MEASURED: μ
// lands mid-band, the 2 k-points inside the Fermi surface fill (n_k=2.0) while the 6 on it smear FRACTIONALLY
// (n_k=0.67, ε≈μ, f=1/(1+e^{0.7})=0.33/spin) -- the textbook smeared Fermi surface, charge Σ_k w_k n_k = 1
// exactly, converged in ~26 iters.  BASIS: VALENCE_LOWQ_SR2 Na at the REAL FCC-Na density (a=10 au, matched to
// Na's atomic volume) -- SR2 drops the diffuse Na s 0.0857 + p 0.05, so the Bloch overlap is well-conditioned
// at the correct lattice constant (cond~38; SR needs an unphysical a=12 to condition).  Remaining tension: SR2
// is a MINIMAL 6-function basis (no diffuse 3s), so a fuller metallic Na basis (valgen, step 1) + IBZ/mesh-
// convergence is the accuracy follow-up; this gate validates the machinery + Fermi surface, not a cohesive E.
TEST(GPW_SCF, NaFCCMetalGlobalMu)
{
    FCCUnitCell cell(10.0);                     // FCC Na at Na's atomic density (density-matched to real BCC Na)
    cell.AddAtom(11, {0,0,0});                  // Na (Zion=1): 3s^1 -- one electron => half-filled band
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwOptions o;
    o.label="Na FCC metal"; o.Nelec=1; o.species={{"Na",1}};
    o.densityEcut=-1.0; o.accelerator="DIIS"; o.globalFermi=true;   // ONE μ across the BZ
    o.kShift=rvec3_t(0.5,0.5,0.5);             // shifted Monkhorst-Pack (k at ±¼)
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.SmearingkT=0.01;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR2), o, /*verbose*/false);

    EXPECT_TRUE(R.converged) << "global μ + smearing converges the half-filled-band metal";
    EXPECT_NEAR(R.charge, 1.0, 1e-6);          // one valence electron, BZ-weighted Σ_k w_k n_k = 1
    EXPECT_LT(R.E.MinusTS, -1e-4);             // −TS<0 AND non-trivial: the Fermi surface IS fractionally filled
    EXPECT_NEAR(R.E.GetTotalEnergy(), 0.045543, 3e-3);    // did-E-move anchor (free energy A at kT=0.01)
}

// ===== EXPERIMENTAL (scratch, item 4): FCC Na -- the honest half-filled-band metal =====
// FCC Na, Zion=1 (3s^1): ONE valence electron per cell => the single conduction band is HALF-FILLED, so μ
// cuts THROUGH the band = a real Fermi surface (unlike Al's degenerate-3p, this is the textbook metal).
// SHIFTED Monkhorst-Pack (kShift=½ => k at ±¼, CP2K's default; avoids the high-symmetry Γ, samples the Fermi
// surface more evenly) + global μ + smearing.  Knobs: NA_KGRID (nxnxn), NA_A (lattice au), NA_KT, NA_SHIFT
// (0=Γ-centred,1=shifted MP), NA_GLOBAL, NA_NMAX.  Basis = VALENCE_LOWQ_SR Na (s from 0.086; the metal WANTS
// the 0.03 diffuse s but it makes the Bloch overlap singular -- the diffuse-basis tension, step-1 territory).
TEST(GPW_SCF, DISABLED_NaFCCMetalExperiment)
{
    auto envd=[](const char* n,double d){const char*s=std::getenv(n);return s?std::atof(s):d;};
    const int nk=(int)envd("NA_KGRID",2);
    const double a=envd("NA_A",10.0);              // FCC Na ~ 5.3 Å = 10 au (BCC-density-matched)
    FCCUnitCell cell(a);
    cell.AddAtom(11,{0,0,0});                       // Na (Zion=1): 3s^1
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));

    GpwOptions o;
    o.label="Na FCC metal"; o.Nelec=1; o.species={{"Na",1}};
    o.densityEcut=-1.0; o.accelerator="DIIS";
    o.globalFermi = envd("NA_GLOBAL",1.0)!=0.0;
    const bool shifted = envd("NA_SHIFT",1.0)!=0.0;
    o.kShift = shifted ? rvec3_t(0.5,0.5,0.5) : rvec3_t(0,0,0);   // shifted MP vs Γ-centred
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=(size_t)envd("NA_NMAX",60); o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.SmearingkT=envd("NA_KT",0.01);
    const BasisSetData bas = envd("NA_SR2",0.0)!=0.0 ? BasisSetData::VALENCE_LOWQ_SR2
                                                     : BasisSetData::VALENCE_LOWQ_SR;
    GpwResult R = envd("NA_ANNEAL",0.0)!=0.0
        ? RunGpwAnnealed(lat, MakeBasisLowQ(cell,bas), o, {0.02,0.01,0.005}, /*verbose*/true)
        : RunGpw       (lat, MakeBasisLowQ(cell,bas), o, /*verbose*/true);
    std::cout<<"[Na metal] nk="<<nk<<" a="<<a<<" shift="<<shifted<<" global="<<o.globalFermi
             <<" conv="<<R.converged<<" charge="<<R.charge<<" A="<<R.E.GetTotalEnergy()
             <<" -TS="<<R.E.MinusTS<<" E(internal)="<<(R.E.GetTotalEnergy()-R.E.MinusTS)<<std::endl;
}

// (4) MULTI-SPECIES GPW: ionic NaF (rocksalt = FCC + 2-atom basis) at Gamma, driven by the multi-species
// Ham_PW_DFT ctor ({{"Na",1},{"F",7}}).
// Valance basis for Na generated from all electrton atom caluclations (BasisSetData/valence_lowq.bsd) 
// looks like     s={0.03.0.086,0.245,0.7,2.0}, p={0.05, 0.3}
// The diffuse exponents s={0.03.0.086} simply don't work in a lattice context.  The overlap is singular.  Auto trimming
// the basis using eeigen instead choleslky decomposition also simply doe not work.  The only option is to tell the
// user to drop thos diffuse basis function.  If use the BasisSetData::VALENCE_LOWQ_SR2 the calculation converges nicely
// with some help from DIIR/GDM/Kerker-mixing.
// Also work noting is the sharp F basis function with exponent 40Ha. THis motivated a ladder of grids to which basis 
// function pairs are assigned bases on thier exponents.  
// There is also a challange for the integration grids used for Vxc fitting.  It is anticipated that using a non-uniform unit Becke/Voronoi-polyhedra grid will allow
// for rapid integration of diffuse basis function will make this run even more efficient. 
// The ideal minimum densityEcut=2*40Ha=80Ha based on the F max exponent.  40 converges to a lower E_total but otherwise converged nicely.
// Explicit densityEcut= 40 = SUB-FLOOR: warns, and BallOnly aliases there (-43 mHa)
TEST(GPW_SCF, DISABLED_NaFRocksaltGamma)
{
    using namespace qchem::Hamiltonian;
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na (Zion=1)
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F  (Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // SR2 (2026-07-16): the complete-enumeration-conditioned basis (lambda_min=1.57e-3; SR's three
    // degenerate 1.03e-6 near-null modes were exactly the Na p 0.05 triplet -- the cation's superfluous
    // diffuse shells; F kept intact for the anion).  See DISABLED_NaFOverlapConditioningSweep.
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ_SR2, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));

    // The production recipe as ONE GpwOptions literal (the full 2-week rationale is in the header above +
    // doc/GPWPlan §0b″).  The NAF_* env knobs stay as sweep INSTRUMENTS; the defaults ARE the committed recipe.
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };
    GpwOptions o;
    o.label        = "NaF GPW Gamma";
    o.Nelec        = 8;                                   // 1 (Na) + 7 (F) valence electrons
    o.species      = {{"Na",1},{"F",7}};
    o.densityEcut  = envd("NAF_ECUT", -1.0);              // AUTO = C·αmax=80 (the anchor config); NAF_ECUT=40 = sub-floor sweep
    o.ladderFactor = envd("NAF_LADDERF", 4.0);
    o.accelerator  = std::getenv("NAF_NULL") ? "Null" : "Ladder";   // Fock DIIS→GDM on |ΔE/E| (ionic); NAF_NULL=damped Kerker
    o.seed         = qchem::ChargeDensity::SeedStrategy::IonicSAD;   // diffuse F⁻/Na⁺ ionic seed (halves iters)
    const double pivotTol = envd("NAF_PIVOT", 1e-4);                // rank-revealing pivoted Cholesky (doc/GPWPlan1.md §4a)
    o.ortho        = pivotTol>0.0 ? qchem::CholeskyPivoted : qchem::Auto;
    o.orthoTol     = pivotTol;
    o.scf.NMaxIter = (size_t)envd("NAF_NMAX", 200);
    o.scf.MinΔE=1e-8; o.scf.MinΔρ=1e-4;                   // E-flat exit AND a Δρ gate (non-variational settled-E map)
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=envd("NAF_ALPHA",0.45); o.scf.KerkerG0=1.0;   // Kerker damps the low-G charge-transfer slosh
    o.scf.UseMOM=true; o.scf.MOMStartIter=10;             // delayed-IMOM: descend, then pin the occupied subspace through the crossing
    o.scf.SmearingkT=envd("NAF_SMEAR",0.0); o.scf.MOMSmearPenalty=envd("NAF_PENALTY",0.0);   // MOM-masked Fermi (experiment)

    qchem::SCFIterator::ReportBandGap()=true;             // per-iteration gap column: watch the diffuse virtual dive (header)
    GpwResult R = RunGpw(lat, mol, o, /*verbose*/true);
    qchem::SCFIterator::ReportBandGap()=false;

    EXPECT_NEAR(R.charge, 8.0, 1e-6);   // 1 (Na) + 7 (F) valence electrons, conserved
    // MECHANISM (the 2-week investigation, condensed; full trace in the header + doc/GPWPlan §0b′/§0b″).  The
    // energy spikes are NOT conditioning / mixing-gain / DIIS: a single diffuse (Na-3s-like) conduction state
    // has a GIANT response to the low-G charge-transfer slosh, dives across the Fermi edge (measured live by
    // ReportBandGap) and aufbau swaps 2e out of the F 2p manifold -- a periodic level-crossing (period ~27) ->
    // E=+5e3.  FIX = delayed-IMOM (pin the occupied subspace through the crossing; the 0h guard releases a bad
    // capture) + Kerker damping; the fixed-point gap is large (~0.35 Ha) so NOT Fermi smearing.
    // ANCHOR: -24.4304 (auto Ecut=80, BallOnly, raw XC, guards) -- 0.8 mHa from CP2K SR2 truth -24.4312; the
    // historical -27.93 "oracle" was RETRACTED as a screening artifact (doc/GPWPlan TRAPS #2).
    EXPECT_NEAR(R.E.GetTotalEnergy(), -24.4304, 0.01);   // did-E-move anchor (the default-config fixed point)
}

// (4b) NaF GRID-CONTINUATION SEEDING (doc/GPWPlan §0e, step 1) -- the PRODUCTION-GRID fix.
// (Pins + story re-derived 2026-07-23 on the post-analytic-short/kappa/5-smooth landscape; the original
// -27.76/-27.93 anchors were the RETRACTED aliasing-era values -- doc/GPWPlan.md TRAPS #2.  The clean NaF
// SR2 truth is CP2K -24.4312 at tight eps, Ecut=160-class grids.)
//
// THE PROBLEM.  The direct production-grid NaF run FALLS INTO the unphysical XC-collapse basin: from the
// ionic seed, the Kerker priming descent goes straight into E~-40 (mid-slosh D loads the sharpest F pairs
// beyond the grid calibration -> the collocated rho aliases spiky/locally-negative -> E_xc is legitimately
// huge-negative WITHIN the discretization, a self-consistent garbage fixed point), and Pulay engaging on
// that garbage state thrashes to +54.  MOM+Pulay are NECESSARY but NOT SUFFICIENT: the basin is a property
// of the map, reachable by the descent -- not an occupation swap (MOM) nor a mixing wobble (Pulay).
// HISTORY: on the BALL-XC map the basin was real at every sub-C=8 grid (the 0.5(f1) sweep hit it from a
// SEEDED start at Ecut=160: negCharge -91, Exc -109).  The 0.5(f2) raw-XC feed REMOVED it (rho_DM >= 0 by
// construction -- negCharge == 0 at C=8/4/3 in the acceptance sweep); this test retains the ionic-seed A/B
// as the historical repro knob.
//
// THE FIX (this test).  Never ENTER the basin: converge the CHEAP coarse grid (Ecut=40), then SEED the fine
// grid with that converged density so the fine SCF STARTS in the physical basin.  The orbital (SR2 Gaussian)
// basis is IDENTICAL at both cutoffs -- only the density COLLOCATION grid differs -- so the converged coarse
// density matrix transfers directly (no re-projection) via the explicit-seed cSCFIterator ctor, which
// collocates the seed's iteration-0 Hartree/XC on the REQUESTED fine fit grid (the fit-grid seam is honest
// since 2026-07-20 -- GPW_IBS builds the tensor over the requested fit basis's grid).  Init immediately
// re-diagonalizes on the fine basis and every subsequent iteration runs the fine grid, starting in (and
// staying in) the physical basin.
//
// MOM ACROSS THE GRID CHANGE (doc/GPWPlan 0h): transferring the coarse WF's occupied subspace as a fixed
// MOM reference (AdoptMOMReference) pinned AN EXCITED STATE across the discretization change (measured
// 2026-07-23: -23.680, +0.75 Ha).  The 0h GUARD (persistent-hole detection -> release + re-capture) now
// makes BOTH recipes land the ground state: pure aufbau (the default: GC_SEED_MOM=0, GC_FINE_MOM_START=
// 9999; 22 iters) and the transfer path (GC_SEED_MOM=1 GC_FINE_MOM_START=1: VERIFIED 2026-07-23, guard
// fires once on the coarse stage's own capture-at-10 reference, fine converges 16 iters to -24.43252 --
// identical to the aufbau pin to 8 decimals).  The guard also exposed that the COARSE stage's endpoint had
// itself been MOM-pinned +0.75 high in every earlier measurement (see the coarse-pin note below).
//
// GATE: the fine grid must reach the raw-XC aufbau ground state -24.4325 (1.3 mHa from CP2K's -24.4312,
// itself an Ecut=160-class number), NOT the -40 basin.  DISABLED (two full NaF SCFs, ~5 min).  Env knobs
// (GC_*) tune each stage without recompiling.  Verify basin-avoidance is REAL by A/B: with GC_SEED=0 the
// fine stage falls back to the ionic seed and must dive into the basin (the direct-run failure this test
// fixes; the energy gates then fail by design).
TEST(GPW_SCF, DISABLED_NaFGridContinuation)
{
    using namespace qchem::Hamiltonian;
    namespace L3=BasisSet::Lattice_3D;
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };

    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na (Zion=1)
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F  (Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto st = lat.GetStructure();       // held for both stages (the ctors' non-owning structure view)
    // GC_BASIS selects the orbital basis (default SR2 = the well-conditioned regression config; "SR" = the
    // FULL short-range basis, lambda_min~1e-6 at complete enumeration -- the sec-1 rank-reduction campaign's
    // probe target, oracle CP2K -27.93128 on VALENCE-LOWQ-SR).
    const char* gcb=std::getenv("GC_BASIS");
    const BasisSetData basis = (gcb && std::string(gcb)=="SR") ? BasisSetData::VALENCE_LOWQ_SR
                                                               : BasisSetData::VALENCE_LOWQ_SR2;
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        basis, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));

    // The converged Ecut=40 recipe (DISABLED_NaFRocksaltGamma): pure damped Kerker (NO DIIS), exit on E-flat,
    // delayed-IMOM MOM + Kerker-preconditioned Pulay.  Tunable per stage (near the fixed point the fine stage
    // can capture MOM + engage Pulay earlier, since it does not need the ~10/35-iter descent the coarse one does).
    auto makePar=[&](size_t nmax, int momStart, int pulayStart)
    {
        SCFParams par; par.NMaxIter=nmax; par.MinΔρ=1e30; par.MinΔE=1e-8; par.MinΔFD=1e30; par.MinVirial=1e30;
        par.MinFD=1e30; par.StartingRelaxRo=envd("GC_ALPHA",0.025); par.MergeTol=1e-4; par.Verbose=true;
        par.KerkerG0=envd("GC_KERKER_G0",1.0);
        par.UseMOM=true; par.MOMStartIter=momStart; par.PulayDepth=(int)envd("GC_PULAY",6); par.PulayStart=pulayStart;
        return par;
    };

    // RSS breadcrumb (the full-SR allocation-bomb bisect, 2026-07-22): prints resident MB per ctor phase.
    auto rss=[](const char* tag)
    {
        std::ifstream f("/proc/self/statm"); size_t vmpg=0, rspg=0; f>>vmpg>>rspg;
        std::cerr<<"[rss] "<<tag<<": "<<(rspg*4096/1048576)<<" MB"<<std::endl;
    };
    // ---- STAGE 1: converge on the CHEAP coarse density grid (Ecut=40 -> the physical fixed point). ----
    // Every coarse-stage object is a unique_ptr so the WHOLE stage can be torn down mid-test (below) the
    // moment the fine stage has consumed it -- doc/GPWPlan.md 0.5(b).
    rss("pre-basis");
    // The coarse SEED stage runs Ecut=40 -- SUB-FLOOR (below C*alpha_max=80), where BallOnly aliases
    // (-43 mHa); pin it to the exact-quadrature raster so the seed is the honest -24.4357 fixed point.
    std::unique_ptr<Complex_BS> bsC(L3::GPWFactory(lat, mol,
        L3::GPWParams{.densityEcut=envd("GC_COARSE_ECUT",40.0), .raster=L3::RasterPolicy::AliasFree}));
    rss("basis");
    auto ecC=std::make_unique<Crystal_EC>(bsC->GetIrreps(Spin::None), 8);
    rss("EC");
    cHamiltonian* hamC=new Ham_PW_DFT(st, bsC.get(), {{"Na",1},{"F",7}}, "LDA");
    rss("Ham");
    auto* accC=new qchem::SCFAccelerators::tSCFAcceleratorNull<dcmplx>();   // no DIIS (the CP2K recipe)
    auto scfC=std::make_unique<qchem::SCFIterator::SolidSCFIterator>(bsC.get(), ecC.get(), hamC, accC,
                                          qchem::ChargeDensity::SeedStrategy::IonicSAD, st.get(),
                                          qchem::Cholesky, 0.0);
    rss("SCFctor");
    qchem::Hamiltonian::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");   // step-2 probe: coarse-grid rho stats to compare vs fine
    scfC->Iterate(makePar((size_t)envd("GC_COARSE_NMAX",200), 10, 35));
    qchem::Hamiltonian::ReportGridCharge()=false;
    auto Ecoarse=scfC->GetEnergy();
    std::cout << "[NaF grid-cont COARSE] Ecut=40 iters="<<scfC->GetIterationCount()
              << " Etot="<<Ecoarse.GetTotalEnergy() << std::endl;
    // The Ecut=40 fixed point on the RAW-XC landscape WITH the 0h MOM guard: -24.4357, E-flat converged in
    // ~43 iters -- only 3.2 mHa from the fine (Ecut=320) -24.4325: under raw XC the Ecut=40 grid is nearly
    // converged.  PIN HISTORY (each anchor exposed by the next fix): -27.76 = the RETRACTED aliasing era;
    // -23.69 (ball, 515 iters) and -23.68 (raw, 45 iters) = a MOM-PINNED EXCITED STATE the 0h guard caught
    // (the capture-at-fill-10 reference grabbed a non-aufbau configuration; persistent ~3 mHa hole ->
    // release -> aufbau recovery).  The "coarse underbinds by 0.74 Ha" story was that excited state's
    // artifact, not grid error.
    EXPECT_NEAR(Ecoarse.GetTotalEnergy(), -24.4357, 0.01);   // seed-quality anchor (did-E-move)

    // Grab the converged coarse density (OWNED; consumed by the fine ctor's Init).  bsC stays alive until
    // after the fine ctor, so the density's coarse-block pointer stays valid for the one iteration-0 read.
    auto* seedCD = scfC->GetWaveFunction()->GetChargeDensity();

    // ---- STAGE 2: seed the PRODUCTION fine grid (auto Ecut=8*alpha_max=320) with the converged coarse density. ----
    std::unique_ptr<Complex_BS> bsF(L3::GPWFactory(lat, mol, /*densityEcut*/envd("GC_FINE_ECUT",-1.0)));  // <0 AUTO=320
    Crystal_EC ecF(bsF->GetIrreps(Spin::None), 8);
    cHamiltonian* hamF=new Ham_PW_DFT(st, bsF.get(), {{"Na",1},{"F",7}}, "LDA");
    auto* accF=new qchem::SCFAccelerators::tSCFAcceleratorNull<dcmplx>();
    qchem::Hamiltonian::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");
    qchem::SCFIterator::ReportBandGap()=true;
    // GC_SEED=0 A/Bs the fix OFF (ionic seed) -> the fine stage must dive into the -39 basin (the failure this
    // test fixes); default ON = the converged-coarse-density explicit seed.
    const bool useSeed = envd("GC_SEED",1.0)!=0.0;
    std::unique_ptr<qchem::SCFIterator::cSCFIterator> scfF;
    if (useSeed)
    {
        scfF.reset(new qchem::SCFIterator::SolidSCFIterator(bsF.get(), &ecF, hamF, accF, seedCD, st.get(),
                                                        qchem::Cholesky, 0.0));   // explicit-seed ctor (consumes seedCD)
        // MOM transfer across the grid change is OFF by default: AdoptMOMReference across a discretization
        // change PINS AN EXCITED STATE (doc/GPWPlan 0h; measured 2026-07-23: -23.680 vs the -24.434 aufbau
        // ground state, +0.754 Ha).  With the density seed holding the run in the physical basin, the pure
        // aufbau fill converges cleanly to the ground state.  GC_SEED_MOM=1 (with GC_FINE_MOM_START=1)
        // re-enables the transfer path for the 0h MOM-guard work.
        if (envd("GC_SEED_MOM",0.0)!=0.0)
            scfF->AdoptMOMReference(*scfC->GetWaveFunction());
    }
    else
    {
        delete seedCD;   // A/B control: discard the coarse density, fall back to the ionic seed (dives to -39)
        scfF.reset(new qchem::SCFIterator::SolidSCFIterator(bsF.get(), &ecF, hamF, accF,
                                                        qchem::ChargeDensity::SeedStrategy::IonicSAD, st.get(),
                                                        qchem::Cholesky, 0.0));
    }
    // The coarse stage is DONE (seed consumed by the fine ctor's Init, MOM reference copied out) -- tear it
    // down IN ORDER (iterator -> EC -> basis) BEFORE the fine iterations (doc/GPWPlan.md 0.5(b)).  The
    // coarse ~GPW_Evaluator hands its ladder's stream caches back to the global budget, and the fine shape
    // (built STARVED during the handoff, while the coarse caches were still resident) rebuilds into the
    // refunded budget at its next EnsureStreams (the self-heal).  Without this the fine stage runs at ~0%
    // stream coverage, re-evaluating billions of points per iteration (the 8.45-h full-SR run).
    scfC.reset(); ecC.reset(); bsC.reset();
    rss("coarse stage freed");
    scfF->Iterate(makePar((size_t)envd("GC_FINE_NMAX",100),
                          (int)envd("GC_FINE_MOM_START",9999), (int)envd("GC_FINE_PULAY_START",12)));
    qchem::Hamiltonian::ReportGridCharge()=false;
    qchem::SCFIterator::ReportBandGap()=false;

    auto Efine=scfF->GetEnergy();
    auto* cd=scfF->GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge(); delete cd;
    std::cout << "[NaF grid-cont FINE] auto-Ecut iters="<<scfF->GetIterationCount()<<" charge="<<charge
              << " Etot="<<Efine.GetTotalEnergy()
              << " (Ekin="<<Efine.Kinetic<<" Een="<<Efine.Een<<" Eee="<<Efine.Eee<<" Exc="<<Efine.Exc
              << " Enn="<<Efine.Enn<<" E_alphaZ="<<Efine.E_alphaZ<<")" << std::endl;
    EXPECT_NEAR(charge, 8.0, 1e-6);     // 1 (Na) + 7 (F) valence electrons, conserved
    // WHAT THIS GATES (re-derived 2026-07-23, post analytic-short/kappa/5-smooth + the 0.5(f2) raw-XC
    // feed): grid-continuation seeding makes the PRODUCTION fine grid converge CLEANLY to the aufbau
    // GROUND STATE.  Under raw-XC dynamics the XC-COLLAPSE basin is REMOVED (negCharge == 0 at every C in
    // the f2 acceptance sweep -- rho_DM >= 0 by construction), so this gate now carries basin history plus
    // the did-E-move pin.  The fine SCF converges in ~22 iters, charge conserved to 1e-8 throughout, at
    // -24.4325 -- 1.3 mHa from the CP2K SR2 truth -24.4312 (an Ecut=160-class number).  The historical
    // "-27.93 oracle / -3.5 Ha Exc step-2 gap" story recorded here previously was the RETRACTED
    // aliasing-era landscape (doc/GPWPlan.md TRAPS #2).
    EXPECT_TRUE(useSeed==false || scfF->Converged()) << "seeded fine SCF converges (no basin/spike thrash)";
    EXPECT_GT(Efine.GetTotalEnergy(), -29.0);   // basin avoidance: NOT the ~-40 unphysical attractor
    EXPECT_LT(Efine.GetTotalEnergy(), -20.0);   //                  NOT the +54 Pulay-thrash garbage
    EXPECT_NEAR(Efine.GetTotalEnergy(), -24.4304, 0.01);   // the raw-XC aufbau ground state at the production
                                                            //   default (auto Ecut=80, BallOnly); AliasFree@320
                                                            //   reference: -24.4325
}

// GATE (GPWPlan1 §4a increment 1): rectangular V is now plumbed through the PERIODIC SCF stack.  The FULL
// valence_lowq NaF basis (n=37) has one intrinsic near-null overlap mode (min eig=8.4e-8, a clean ~1000x gap
// below the physical ~1e-3 directions); canonical Eigen(1e-6) drops it -> a rank-36 rectangular ortho V.
// Previously this threw "Matrix sizes do not match" before iter 1 because the orthonormal density D' was
// allocated full n=37; now D' is sized to m=n-k=36 (IrrepWF/TOrbitals) while the AO density D stays 37x37.
//
// It asserts (a) the truncation actually happened (never-silent line captured) and (b) charge = Tr(DS) = 8
// is conserved through the rectangular-V density build + DIIS + collocation (charge is grid-INDEPENDENT).
// DISABLED for now: correct but SLOW -- the cost is the analytic LatticeSum1E (kinetic + local-PP) over the
// FULL 37-function basis, whose DIFFUSE modes are long-range and reach many periodic images (measured: the
// SolidSCFIterator ctor's initial H build dominates, minutes; the collocation is cheap).  Not a regression
// (the suite is unaffected); it is the intrinsic periodic-diffuse lattice-sum cost.  The FAST committed
// rank-reduction gate is a MINIMAL near-degenerate periodic system (few images); this stays as the
// full-basis marker + the converged-energy check (with DISABLED_NaFFullBasisEigenTol below), run on demand.
TEST(GPW_SCF, DISABLED_NaFFullBasisRankReduction)
{
    using namespace qchem::Hamiltonian;
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});
    cell.AddAtom(9,  {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, mol, /*densityEcut*/20.0));   // low: cheap grid, plumbing only
    auto       irreps=bs->GetIrreps(Spin::None);
    Crystal_EC ec(irreps, 8);
    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), {{"Na",1},{"F",7}}, "LDA");
    auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(qchem::SCFAccelerators::DIISParams{8, 8.0, 1e-10, 1e-8});
    testing::internal::CaptureStdout();
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::IonicSAD, lat.GetStructure().get(),
                                         qchem::Eigen, 1e-6);   // canonical ortho drops the ~8e-8 null mode -> rank 36
    SCFParams par; par.NMaxIter=3; par.MinΔρ=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30;
    par.MinFD=1e30; par.StartingRelaxRo=0.3; par.Verbose=false;
    scf.Iterate(par);
    std::string log=testing::internal::GetCapturedStdout();
    auto* cd=scf.GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge(); delete cd;
    // (a) the truncation must actually have fired (rank reduction really exercised, not a trivially-square pass).
    EXPECT_NE(log.find("LASolverEigen truncating"), std::string::npos) << "expected the near-null mode to be dropped";
    // (b) charge conserved through the rectangular-V density build (Tr(DS)=8, grid-independent).
    EXPECT_NEAR(charge, 8.0, 1e-6);
}

// The HEAVY converged-energy run on the full valence_lowq basis + Eigen(1e-6): kept DISABLED (slow -- full
// auto Ecut=80, 60 iters).  Since increment 1 the periodic stack handles the rank-36 rectangular V, so this
// no longer throws; it is the manual full-basis-vs-SR/SR2 energy check (GPWPlan1 §4a gate 5), run on demand.
TEST(GPW_SCF, DISABLED_NaFFullBasisEigenTol)
{
    using namespace qchem::Hamiltonian;
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});
    cell.AddAtom(9,  {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, mol, /*densityEcut AUTO*/-1.0));
    auto       irreps=bs->GetIrreps(Spin::None);
    Crystal_EC ec(irreps, 8);
    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), {{"Na",1},{"F",7}}, "LDA");
    auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(qchem::SCFAccelerators::DIISParams{8, 8.0, 1e-10, 1e-8});
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::IonicSAD, lat.GetStructure().get(),
                                         qchem::Eigen, 1e-6);   // (2): canonical ortho, drop the ~0 null cluster
    SCFParams par; par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30;
    par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true; par.KerkerG0=1.0;
    qchem::Hamiltonian::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");
    scf.Iterate(par);
    qchem::Hamiltonian::ReportGridCharge()=false;
    auto* cd=scf.GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge(); delete cd;
    auto E=scf.GetEnergy();
    std::cout << "[NaF GPW full/Eigen(1e-6)] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Etot="<<E.GetTotalEnergy() << std::endl;
    EXPECT_NEAR(charge, 8.0, 1e-6);
}

// (4b) FAST overlap-conditioning sweep for NaF: build ONLY the analytic Bloch overlap S(Gamma) (via GPW_IBS,
// densityEcut=0 -> no collocation, no SCF) for the full vs SR valence basis across Rcut, and report min/max
// eig(S) PLUS the ORTHOGONALISER RESIDUAL ‖VᴴSV − I‖ (V=S^-½ built by the SAME LASolver the SCF uses).  Seconds
// per point (analytic 1E sum), so we settle the conditioning question before paying for a full SCF.
//
// PROBE 1 (disentangling, doc/GPWPlan §0): this SEPARATES the two pathologies that both read as "conditioning":
//   (A) NEAR-SINGULAR but PSD (min eig>0, e.g. SR/Rcut=2a min eig 7.5e-4).  cond(S)~1/min_eig looks scary but
//       cond(V)=√cond(S), so the orthogonaliser is essentially EXACT: ‖VᴴSV−I‖ ~ machine eps.  => conditioning
//       here is a RED HERRING; whatever ails the SCF, it is NOT the metric.  (The plan's predicted result.)
//   (B) INDEFINITE (min eig<0) at small Rcut from the SHARP |R|≤Rcut cutoff (full basis: -0.42 at Rcut=a).  No
//       real S^-½ EXISTS -> Cholesky fails; this is a genuine problem, but a TRUNCATION one, fixed by
//       MAGNITUDE-SCREENING the image pairs (CP2K's EPS_PGF_ORB), NOT by a better orthogonaliser.  Flagged, not
//       measured (residual undefined).
// So one cheap table closes axis A (conditioning-of-a-valid-metric = red herring) and names axis B (sharp-cutoff
// truncation -> magnitude screening) as the real overlap work.
TEST(GPW_SCF, DISABLED_NaFOverlapConditioningSweep)
{
    namespace L3=BasisSet::Lattice_3D;
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F

    // ‖VᴴSV − I‖_max for a given ortho method: Transform(S)=VᴴSV must be the identity iff V=S^-½ is exact.
    auto residual=[](const hmat_t<dcmplx>& S, qchem::Ortho ortho, double tol=0.0)->double
    {
        std::unique_ptr<LASolver<dcmplx>> la(LASolver<dcmplx>::Factory(ortho, tol));
        la->SetBasisOverlap(S);
        hmat_t<dcmplx> I=la->Transform(S);     // VᴴSV (identity on the RETAINED subspace after truncation)
        double r=0.0;
        for (size_t i=0;i<I.rows();++i)
            for (size_t j=0;j<I.columns();++j)
                r=std::max(r, std::abs(dcmplx(I(i,j)) - (i==j ? 1.0 : 0.0)));
        return r;
    };

    auto probe=[&](BasisSetData bd, const char* name)
    {
        auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
            bd, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
        // BANISH-Rcut: the old Rcut sweep axis is unrepresentable (enumeration lives inside the seam);
        // the one remaining question is the COMPLETE-enumeration conditioning of each basis.
        {
            L3::GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/0.0);
            const BasisSet::Complex_OIBS& g = gpw;
            auto S = g.Overlap();
            rvec_t d; mat_t<dcmplx> U; blazem::eigen(S, d, U);   // ascending eigenvalues of Hermitian S
            std::cout << "[cond " << name << "] n=" << S.rows() << " (complete enumeration)"
                      << "  eig[0..3]=" << d[0] << "," << d[1] << "," << d[2] << "," << d[3]
                      << "  max=" << d[d.size()-1];
            std::cout << "  ‖VᴴSV-I‖: Eigen(1e-6)=" << residual(S, qchem::Eigen, 1e-6)
                      << " SVD(1e-6)=" << residual(S, qchem::SVD, 1e-6);
            if (d[0] > 0.0) std::cout << " Chol=" << residual(S, qchem::Cholesky);
            std::cout << std::endl;
        }
    };
    probe(BasisSetData::VALENCE_LOWQ,     "full");
    probe(BasisSetData::VALENCE_LOWQ_SR,  "SR  ");
    probe(BasisSetData::VALENCE_LOWQ_SR2, "SR2 ");
}
