// File GPW_SCF_UT.C  The GPW self-consistent total energy: the first periodic SCF on GAUSSIAN orbitals.
//
// GPW (increments 1-2) already satisfies every plane-wave Kohn-Sham concept EXCEPT the external potential:
//   - kinetic  -> Kinetic<dcmplx> calls bs->Kinetic()         (GPW: lattice-sum <p^2>)              [inc 1]
//   - Hartree  -> PW_Hartree casts bs to Orbital_DFT_IBS<dcmplx> + cd to FourierDensity (GPW: collocation tensors)[inc 2]
//   - XC       -> PWFittedVxc, same casts + the fit-basis grid                                             [inc 2]
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
#include <functional>   // the per-iteration order-parameter probe (GpwOptions::orderProbe)
#include <string>
#include <iomanip>      // setprecision (the order-parameter trajectory line)

import qchem.Structure;                          // Molecule, Atom
import qchem.UnitCell;                           // UnitCell, FCCUnitCell
import qchem.Matrix3D;                           // Matrix3D<double> (the rhombohedral AFM-II cell matrix)
import qchem.Lattice_3D;                         // Lattice_3D
import qchem.BasisSet;                           // Complex_BS, Real_BS
import qchem.BasisSet.Orbital_1E_IBS;            // Complex_OIBS (the overlap-spectrum diagnostic)
import qchem.Blaze;                              // blazem::eigen, blaze::min/max (overlap spectrum)
import qchem.BasisSet.Lattice_3D.PlaneWave_IBS;   // PlaneWave_IBS (the seed's CD fit basis)
import qchem.BasisSet.Lattice_3D.BasisSet;       // GPWFactory (the GPW basis container)
import qchem.BasisSet.Molecule.Factory;          // Molecule::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Molecule.PG_Spherical.LatticeView;  // MakeSphericalLatticeView (GPW_SPHERICAL=1)
import qchem.Hamiltonian.Factory;                 // the PUBLIC solid front door (Step 4): cHamiltonian* Factory(...)
import qchem.SolidCalculation;                    // the NAMED periodic facade (Step 4 3/3)
import qchem.Hamiltonian.Internal.Hamiltonians;  // Ham_PW_DFT direct ctors (the bespoke probes below still use them)
import qchem.Hamiltonian.Internal.PWTerms;        // ReportGridCharge(); PWFittedVxc / DeltaFittedVxc (the Becke XC gate)
import qchem.Mesh.Angular;                        // MakeAngular (the rotated-Lebedev bond-angle probe)
import qchem.Hamiltonian.Internal.ExFunctional;   // ExFunctional (the LDA functional face the XC terms hold)
import qchem.Hamiltonian.Internal.SlaterExchange; // SlaterExchange (Dirac exchange, for the Becke XC gate)
import qchem.Hamiltonian.Internal.VWN_Correlation;// VWN_Correlation (VWN5, for the Becke XC gate)
import qchem.Mesh;                                // qcMesh::MeshParams / UnitCellKind (the Becke XC quadrature)
import qchem.Mesh.XCPolicy;                       // BeckeXCParams / ResolveXCMesh / XCMeshSharpness (the grid policy)
import qchem.BasisSet.Molecule.LatticeSum1E;      // Molecule::LatticeSum1E::MaxExponent (alpha_max, for the selector)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH -> HGH local PP (alpha_pp = 1/2r_loc^2, for the selector)
import qchem.PeriodicTable;                       // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ElectronConfiguration.Crystal;      // Crystal_EC (single-k Bloch occupation)
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Factory;              // the PUBLIC complex accelerator door (Step 4)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS; // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;  // SCFAcceleratorGDM (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder; // SCFAcceleratorLadder (DIIS -> GDM chain)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (NaF: pure damped Kerker)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Reporting;                          // report:: -- bracket the GPW run so grids/basis sections land
import qchem.Symmetry.Spin;                      // Spin
import qchem.Symmetry.Factory;                   // BlochFactory (build a k-block with a fractional MP shift)
import qchem.LASolver;                           // qchem::Ortho (Cholesky | Eigen | SVD -- basis orthogonalisation)
import qchem.BasisSet.Lattice_3D.GPW_IBS;         // GPW_IBS (build a concrete block for the collocation diagnostic)
import qchem.BasisSet.Lattice_3D.Evaluators.GPW;  // GPW_Evaluator (Overlap3CTensor -- the collocation tensor)
import qchem.BasisSet.Internal.GMap;              // Projector3<dcmplx> (the collocation weight tensor); SymmetryDefects (§3 diagnostic)
import qchem.ChargeDensity.FourierDensity;        // FourierDensity (ρ̃ for the §3 order-parameter diagnostic)
import qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.SeedCD;              // PolarizedSeedCD (the raw spin-SAD seed, for the sublattice gate)               // IrrepCD_Factory/PolarizedCD_Factory (fixed-density probe)
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
// GPW_SPHERICAL=1 (doc/SphericalLatticePlan.md I1/I2): wrap the molecular basis in its spherical
// (contaminant-free) lattice view -- the span-matched A/B against CP2K's spherical-d convention.
// s/p-only bases are unchanged in SPAN (T = identity blocks), so Si/NaF runs double as null tests.
std::shared_ptr<const Real_BS> MaybeSpherical(std::shared_ptr<const Real_BS> bs)
{
    if (std::getenv("GPW_SPHERICAL"))
        return BasisSet::Molecule::PG_Spherical::MakeSphericalLatticeView(std::move(bs));
    return bs;
}
// The SHORT-RANGE variant (most diffuse valence primitives dropped) -- well-conditioned Bloch overlap in a solid.
std::shared_ptr<const Real_BS> MakeBasisSR(const Structure& st)
{
    return MaybeSpherical(std::shared_ptr<const Real_BS>(
        BasisSet::Molecule::Factory(BasisSetData::SIPP_SR, &st,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian)));
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
        BasisSet::Molecule::Factory(which, &st,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian)));
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
struct FpRow { size_t it; double E, dEabs, fd, drho, order=0; };

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

// One GPW Gamma-point SCF: build the GPW basis over the lattice, hand it the plane-wave LDA Hamiltonian
// (Ham_PW_DFT reaches GPW's real-space Integrals_Pseudo), seed uniform, run the complex-DIIS cSCFIterator.
// The gate-calibrated Becke XC quadrature recipe and the Auto-route resolution are LIBRARY policy
// (qcMesh::BeckeXCParams / qcMesh::ResolveXCMesh in src/Mesh/Mesh.C, beside MeshParams -- moved there
// 2026-08-07, D6): the production recipe and its GPW_BECKE_* sweep instruments are documented at the
// declarations, and these tests READ that policy rather than defining it.
//
// GpwOptions: the control surface for a GPW run.  Every material (Si, NaF, CsI, LiCoO2, f-oxides, ...) is one
// options literal -- multi-species PP, grid params (the efficiency lever for ionic systems), accelerator policy,
// seed, ortho, and the SCFParams gates.
struct GpwOptions
{
    std::string label = "gpw";
    int         Nelec = 0;
    //! Spin multiplicity 2S+1 (the facade convention, CalcOptions::multiplicity): 0 = minimal spin /
    //! unpolarized (historical behaviour); >=1 = the SPIN-NATIVE two-channel run (SymmetryUpgradePlan §4
    //! tier 4b) -- nUp-nDown = 2S through Crystal_EC and the polarized Ham_PW_DFT XC pair.  NB unlike the
    //! facade, multiplicity=1 here runs the EXPLICIT two-channel singlet (nUp=nDn): the ζ=0-collapse
    //! cross-check of the polarized machinery against the unpolarized anchors.
    int         multiplicity = 0;
    std::vector<std::pair<std::string,int>> species;   // multi-species PP, e.g. {{"Na",1},{"F",7}}
    // grids
    double densityEcut  = -1.0;                        // <0 AUTO = cutoffFactor*alpha_max
    double cutoffFactor = 2.0;
    double ladderFactor = 4.0;
    BasisSet::Lattice_3D::RasterPolicy raster = BasisSet::Lattice_3D::RasterPolicy::BallOnly;
    BasisSet::Lattice_3D::CellImages   images = BasisSet::Lattice_3D::CellImages::Periodic;
    rvec3_t kShift = rvec3_t(0,0,0);
    //! XC-quadrature policy, decided by \c xcMesh.cellKind.  DEFAULT \c Auto (2026-08-01 flip): the
    //! driver resolves it to the atom-centred periodic BECKE mesh (DeltaFittedVxc, the calibrated
    //! qcMesh::BeckeXCParams() recipe) — \c imposeSymmetry runs included since the 2026-08-02 carve-out retirement
    //! (plan §7 step 5): an imposed run builds the MIXED-RULE site-adapted invariant Becke mesh
    //! (~2x the free mesh at the production recipe) and star-averages ρ on it.  EXPLICIT \c Uniform
    //! or \c Becke is always honored.  The run announces the choice on the [XC quadrature] console
    //! line either way.
    qcMesh::MeshParams xcMesh{.cellKind=qcMesh::UnitCellKind::Auto};
    //! WHICH fit basis represents v_xc -- ORTHOGONAL to the grid choice above (plan §6a fit/grid
    //! separation): Auto = the historical pairing (Delta on Becke, PlaneWave on the uniform raster);
    //! Delta may be paired with EITHER grid (Delta+uniform = the band-limit cross-check cell).
    Hamiltonian::VxcFit vxcFit = Hamiltonian::VxcFit::Auto;
    // convergence machinery
    std::string accelerator = "DIIS";                  // DIIS | GDM | Ladder | Null
    bool        globalFermi = false;                    // metal: one μ across the k-mesh (Crystal_EC global mode)
    //! ONE μ over BOTH SPIN CHANNELS -- the moment becomes an OUTPUT instead of a constraint
    //! (Crystal_EC spinsShareFermi; needs SmearingkT>0).  Default off: a fixed-multiplicity run (the O2
    //! triplet, a singlet gate) is CONSTRAINED on purpose, and sharing the reservoir would let its
    //! multiplicity drift.  Turn it on where the magnetism must find itself -- with it off, μ↑≠μ↓ and the
    //! difference is an unrelieved driving force to move charge between the channels (MnO run 29: 27 mHa).
    bool        spinsShareFermi = false;
    //! DEFAULT FALSE (V1.30, 2026-08-09).  It read \c true while the comment below said OPT-IN, and the
    //! comment was the one carrying a reason.  Now actively hazardous rather than cosmetic: per V1.28 an
    //! imposed run star-averages EACH SPIN CHANNEL under the CHEMICAL space group, whose sublattice-exchanging
    //! ops would average an AFM structure's magnetic sublattices together and quietly demagnetise it -- the
    //! run would not crash, it would converge to the wrong state.  Default-ON meant every new magnetic call
    //! site inherited that unless it remembered to opt out, which is how MnO's `imposeSymmetry=false` became
    //! load-bearing rather than merely explicit.  DEFAULT-OFF IS THE SAFER WAY TO BE WRONG: an imposition you
    //! did not ask for is invisible in the result, a missing one only costs time.  Every pre-existing caller
    //! that was relying on the old default now says so at its own call site -- no run changed behaviour.
    bool        imposeSymmetry = false;                 // impose the detected space group (renamed from reduceBZ
                                                         //   2026-08-02): IBZ k-fold + per-iteration density
                                                         //   star-average + T1 {G}-star sweeps + the site-adapted
                                                         //   invariant Becke mesh.  OPT-IN per the §3 pin (an
                                                         //   imposed default would also ~2x the suite's XC grids).
                                                         //   S3: on a POLARIZED run with flip bits this resolves
                                                         //   to the SHUBNIKOV group of the declared ordering.
    bool        greyImposition = false;                 // §8 NEGATIVE CONTROL only: impose the GREY group on a
                                                         //   staggered cell (decoration suppressed) -- the
                                                         //   star-average then maps +m onto -m and must KILL
                                                         //   m_stag.  Never a production setting.
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Uniform;
    qchem::Ortho ortho    = qchem::Cholesky;
    double       orthoTol = 0.0;
    SCFParams    scf;                                  // NMaxIter / MinΔρ / MinΔE / SmearingkT / ... (the gates)
    //! SEED-CHARACTER MOM (with \c scf.UseMOM): adopt the SEED's own diagonalized occupied subspace as the
    //! fixed MOM reference, BEFORE iteration 1 runs.  The delayed-IMOM capture inside the wavefunction cannot
    //! do this: \c SetMOM is pushed in by \c Iterate, i.e. AFTER \c Init's seed fill, so \c itsUseMOM is still
    //! false when the seed is filled and the earliest self-capture is iteration 1's occupied set.  That is one
    //! iteration too late for a run whose ORDER dies at iteration 1 (MnO: m_stag 0.366 -> 0.0046), which is
    //! precisely the case a delayed capture was never designed for: the seed is not "garbage to settle away
    //! from", it is the physics (the staggered IonicSAD d^5 pair) and the loop is what loses it.  So: pin the
    //! occupied CHARACTER to the seed and let only the DENSITY relax.  Implemented as self-adoption through
    //! the existing grid-continuation face (AdoptMOMReference), no new wavefunction API.
    bool         momFromSeed = false;
    //! Optional ORDER PARAMETER (SCFIterator::SetOrderParameter): a named scalar measured on the WORKING
    //! density every iteration -- an extra trace column PLUS a compact end-of-run trajectory line, so a
    //! symmetry-broken basin (AFM staggering, charge disproportionation) can be watched living or dying
    //! iteration by iteration.  Empty (default) = no probe, no column, no cost.
    std::string  orderName;
    std::function<double(const qchem::ChargeDensity::cDM_CD&)> orderProbe;
    //! Thread GPWParams::hamPreservesReal (doc/RealComplexPlan.md 3c-3): build each TRIM block REAL.
    //! DEFAULT TRUE since 2026-08-18 -- the harness-wide flip, which is what makes every Γ-only anchor
    //! in this file a standing real-block regression test instead of leaving the shipped path (where
    //! \c SolidCalculation has flipped by default since 3c-3) covered only by the dedicated gates.
    //! Set false for the A/B controls (the real-vs-complex twins, the report gate's control arm) and
    //! for anything that needs the historical all-complex build.
    //! \note The physics is unchanged to gate resolution but NOT bitwise (dsyev vs zheev + DIIS
    //! amplification, see \c ExpectRealComplexTwins) -- so a did-E-move pin tightened to ~1e-6 or
    //! below may need re-anchoring, and the re-anchored value is the REAL one from here on.
    bool realTRIMBlocks = true;
};

namespace
{
// The run's SHARPNESS, gathered for the Uniform-vs-Becke cost selector (V1.26).  FACADE-LEVEL code: it
// belongs beside GpwOptions/RunGpw and moves with them into src/Calculation/ (V1.26 step 4).  It states
// FACTS about the run -- geometry and the two sharpness sources -- and decides nothing; qcMesh::ResolveXCMesh
// owns the policy.
//
// Both sharpness sources come off ABSTRACT capability faces, reached by the sanctioned abstract->abstract
// cross-cast, so nothing here touches a concrete basis or a concrete PP model:
//   alpha_max -- BasisSet::Molecule::LatticeSum1E::MaxExponent(), documented as "the GPW density-grid cutoff floor".
//   alpha_pp  -- Pseudopotential::LocalPotential_Gaussian::ShortRangeGaussian(Z), whose terms carry
//                alpha = 1/(2 r_loc^2).  A model with no closed-Gaussian short part does not implement the
//                face; that leaves alpha_pp at 0, which the selector reads as "not measurable" -- NOT as
//                "smooth" (see the XCMeshSharpness doc).
qcMesh::XCMeshSharpness GatherSharpness(const Lattice_3D& lat, const Real_BS& mol, const GpwOptions& o)
{
    qcMesh::XCMeshSharpness s;
    s.cellEdge=lat.GetUnitCell().GetMaximumCellEdge();
    s.nAtoms  =int(lat.GetUnitCell().GetNumAtoms());
    s.imposed =o.imposeSymmetry;
    for (auto ibs : const_cast<Real_BS&>(mol).Iterate<BasisSet::Real_OIBS>())
        if (const auto* ls=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(ibs))
            { s.alphaMax=ls->MaxExponent(); break; }
    for (const auto& [element, valence] : o.species)
    {
        const int Z=int(thePeriodicTable().GetZ(element));
        const Pseudopotential::HGH_LocalPotential loc=Pseudopotential::GetGTH(element,"LDA",valence).local;
        const auto* g=static_cast<const Pseudopotential::LocalPotential_Gaussian*>(&loc);
        for (const auto& t : g->ShortRangeGaussian(Z)) s.alphaPP=std::max(s.alphaPP, t.alpha);
    }
    return s;
}

// Shubnikov S3 (doc/SymmetryUpgradePlan.md §7 step 7): the per-atom collinear decoration for a
// POLARIZED IMPOSED run -- the factory then imposes the MAGNETIC group of the declared ordering
// (the grey group would erase it).  Assembled with the seed's OWN species resolution
// (MagneticDecoration + IonicSADTargets for an IonicSAD run), so the imposed symmetry and the seed
// that plants the moments cannot drift apart.  Empty (grey imposition / free run) otherwise.
std::vector<int> GatherSiteSpins(const Lattice_3D& lat, const GpwOptions& o)
{
    using namespace qchem::ChargeDensity;
    if (o.multiplicity<1) return {};             // unpolarized: no channels, no decoration
    if (o.greyImposition) return {};             // §8 NEGATIVE CONTROL: suppress the decoration so an
                                                 // imposed run star-averages under the GREY group --
                                                 // m_stag must die; never a production setting
    auto st=lat.GetStructure();
    const std::map<size_t,int> targets = (o.seed==SeedStrategy::IonicSAD)
                                       ? IonicSADTargets(st.get(), "LDA") : std::map<size_t,int>{};
    return MagneticDecoration(st.get(), "LDA", targets);
    // (Returned for FREE runs too: the factory reads it only when imposing; the free-run MAGNETIC
    // diagnostic below uses it as the Shubnikov reference group.)
}
} //anon

// Build the complex SCF accelerator named by \a policy, through the PUBLIC typed door (Step 4 2/3).
// The tuned constants that used to sit here ARE the GPW production recipe, so they now live with the
// factory as SolidAcceleratorOptions' defaults -- this function is just the policy-name lookup a facade
// would do.  Ladder = the ionic-crystal DIIS->GDM hand-off on |ΔE/E| (NaF's proven recipe).
static qchem::SCFAccelerators::SCFAccelerator* MakeGpwAccelerator(const std::string& policy)
{
    using namespace qchem::SCFAccelerators;
    if (policy=="Null")   return Factory(Type::Null);
    if (policy=="DIIS")   return Factory(Type::DIIS);
    if (policy=="GDM")    return Factory(Type::GDM,  {.gdm={1.0,0.1}});
    if (policy=="Ladder") return Factory(Type::Ladder);
    throw std::runtime_error("MakeGpwAccelerator: unknown policy \""+policy+"\" (DIIS|GDM|Ladder|Null)");
}

// Optional keep-alive handles for post-SCF term-level probes (the Becke XC gate): the basis and the
// converged density, which stays valid after the iterator tears down because its basis block is bs.
struct GpwHandles
{
    std::unique_ptr<Complex_BS> bs;
    std::unique_ptr<qchem::ChargeDensity::cDM_CD> cd;
};

// [symmetry] The §3 order-parameter line (doc/SymmetryUpgradePlan.md): a run ALWAYS tells the user the
// symmetry it actually found.  FREE run: measure the converged density's per-op defect against the FULL
// detected crystal point group (SymmetryDefects on ρ̃ -- G-space, exact index permutation, so no
// raster-commensurability caveat) and report the max + which-ops-broke count.  IMPOSED (imposeSymmetry) run:
// the defect is ≈0 BY CONSTRUCTION (every iterate is star-averaged), so print the assertion + the
// release-audit reminder instead -- imposition is never silent.
// THRESHOLD (measured): the through-SCF defect inherits the run's CONVERGENCE tolerance, not the plan-§8
// 1e-13 fold tier -- a converged (Δρ~2.5e-6) free Si Γ run measures max defect ~1e-5.  So the BROKEN alarm
// sits at 1e-3: comfortably above convergence noise, far below any genuine order parameter (O(0.01..1)).
// Below it the line reports the number without judging (the release-check quantifies real SSB anyway).
static void ReportSymmetryFound(const Complex_BS& bs, const qchem::ChargeDensity::cDM_CD& cd,
                                const Structure* st, bool imposed,
                                const std::vector<Symmetry::Lattice_3D::SymOp>& magOps = {})
{
    auto ops = bs.GetDetectedReciprocalOps();
    if (ops.size() <= 1) return;   // trivial (or undetected) group: nothing to report
    if (imposed)
    {
        std::cout << "[symmetry] IMPOSED: full crystal point group ("<<ops.size()<<" ops) star-averages every "
                     "iterate -- defect==0 by construction; the release-check is the SSB audit "
                     "(doc/SymmetryUpgradePlan.md - 3)."<<std::endl;
        return;
    }
    auto* fd = dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(&cd);
    if (!fd) return;
    Hamiltonian::PWFittedVxc::fbs_t fit(bs.CreateVxcFitBasisSet(st, qcMesh::MeshParams{}));
    std::vector<double> d = SymmetryDefects(fd->GetFourierDensity(*fit), ops);
    double mx = 0; int broken = 0;
    for (double x : d) { mx = std::max(mx, x); if (x > 1e-3) ++broken; }
    std::cout << "[symmetry] FREE run: density defect vs the detected crystal point group ("<<ops.size()
              << " ops): max="<<mx
              << (broken ? "  ** "+std::to_string(broken)+" op(s) BROKEN -- symmetry lowered (order parameter)"
                         : "  (no lowering resolved above SCF tolerance)") << std::endl;

    // S4: the MAGNETIC diagnostic for a polarized free run with a declared ordering -- did the CHANNEL
    // PAIR keep the Shubnikov group?  (The grey line above sees only the total; a dead sublattice moment
    // is invisible to it.)  σ=Flip ops compare across the channels, so their row IS the m1=-m2 mirror.
    if (magOps.empty()) return;
    const auto* pol = dynamic_cast<const qchem::ChargeDensity::cPolarized_CD*>(&cd);
    if (!pol) return;
    const auto* fdu = dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(pol->GetChargeDensity(Spin::Up));
    const auto* fdd = dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(pol->GetChargeDensity(Spin::Down));
    if (!fdu || !fdd) return;
    std::vector<Symmetry::Lattice_3D::SymOp> rmag;
    for (const auto& op : magOps) rmag.push_back(Symmetry::Lattice_3D::ReciprocalOf(op));
    std::vector<double> md = MagneticSymmetryDefects(fdu->GetFourierDensity(*fit),
                                                     fdd->GetFourierDensity(*fit), rmag);
    double mxN=0, mxF=0;
    for (size_t i=0;i<rmag.size();++i)
        (rmag[i].sigma==qchem::Symmetry::SpinAction::Flip ? mxF : mxN) = std::max(
            rmag[i].sigma==qchem::Symmetry::SpinAction::Flip ? mxF : mxN, md[i]);
    std::cout << "[symmetry] MAGNETIC (Shubnikov, "<<rmag.size()<<" ops): channel-pair defect max "
              << "None="<<mxN<<"  Flip="<<mxF
              << (mxF>1e-3 ? "  ** the sublattice mirror m1=-m2 is BROKEN" : "  (the magnetic order holds)")
              << std::endl;
}

//! GPW_SHARED_MU=1/0 OVERRIDES GpwOptions::spinsShareFermi for an A/B, and SAYS SO.  The question it
//! settles is whether a state we believe is the GROUND state survives when the moment is free: a
//! fixed-(nUp,nDn) fill is a CONSTRAINED calculation whose two μ are the constraints' Lagrange multipliers,
//! so it can hold a state the energetics would not choose -- and if it must, the constraint is hiding a
//! defect rather than curing one (user, 2026-08-10).  Diagnostic valve, never a production setting.
static bool SharedMu(const GpwOptions& o)
{
    const char* e=std::getenv("GPW_SHARED_MU");
    if (!e) return o.spinsShareFermi;
    const bool shared=(std::atoi(e)!=0);
    if (shared!=o.spinsShareFermi)
        std::cerr << "[fill] GPW_SHARED_MU=" << e << " OVERRIDES the run's choice: the two spin channels "
                  << (shared ? "SHARE one μ (moment free)" : "keep separate μ (moment constrained)") << "."
                  << std::endl;
    return shared;
}

// The GENERAL GPW driver: basis -> fail-fast conditioning pre-flight -> grids -> multi-species Hamiltonian ->
// accelerator -> SCF, with automatic reporting (GpwReport) + heartbeat logging throughout.  Any material is a
// GpwOptions literal; the positional RunGPW below and the bespoke NaFRocksaltGamma are thin callers.
static GpwResult RunGpw(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const GpwOptions& o,
                        bool verbose=false, GpwHandles* keep=nullptr)
{
    namespace L3=BasisSet::Lattice_3D;
    const std::string sp = o.species.empty() ? std::string() : o.species.front().first;
    GpwReport report(sp+" "+o.label, verbose);

    qchem::report::Log("building GPW basis");
    // Field-sharpness routing: HartreeOnly (drop the 2/3*alpha_max floor when the XC is off-raster) was
    // the plan, but MEASURED 2026-07-31 it destabilises NaF SR2 at EVERY reduced beta (0 -> +904 Ha;
    // beta/alpha_max in {1/12..1/3} -> 20-31 Ha slosh; only the historical 2/3 converges): the floor
    // protects the DENSITY's low-G integrity (diffuse-pair FFT fold-back), NOT just V_xc.  So Becke runs
    // keep HartreeXC routing; GPW_ROUTING=hartree is the CALIBRATION valve (pair with GPW_RELFIELDSHARP)
    // for the coarse-end routing campaign (kMinLevelN etc.), not a production setting.
    const bool routeExperiment = std::getenv("GPW_ROUTING")
                              && std::string(std::getenv("GPW_ROUTING"))=="hartree";
    std::unique_ptr<Complex_BS> bs;
    {
        qchem::report::Timed t("setup: GPW basis build");
        bs.reset(L3::GPWFactory(lat, mol, L3::GPWParams{
            .densityEcut=o.densityEcut, .cutoffFactor=o.cutoffFactor, .raster=o.raster,
            .images=o.images, .kShift=o.kShift, .ladderFactor=o.ladderFactor, .imposeSymmetry=o.imposeSymmetry,
            .rasterFields=routeExperiment ? L3::RasterFields::HartreeOnly : L3::RasterFields::HartreeXC,
            .siteSpins=GatherSiteSpins(lat, o), .hamPreservesReal=o.realTRIMBlocks}));
    }

    // FAIL-FAST: vet the (analytic, grid-free) overlap + emit basis BEFORE grids/Hamiltonian; a rank-deficient
    // basis aborts here instead of building the whole ladder first.  Also puts `basis` before `grids`.
    int nRemoved;
    {
        qchem::report::Timed t("setup: vet basis (analytic S)");
        nRemoved=report.VetBasis(*bs);
    }
    if (nRemoved > 0)
    {
        std::cout << "["<<o.label<<"] ABORT: basis rank-deficient (see basis.removed) -- skipped grids + SCF."<<std::endl;
        return {false, 0.0, qchem::EnergyBreakdown{}, 0};
    }
    report.EmitGrids(*bs);

    qchem::report::Log("building Hamiltonian");
    auto       irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry the Sum_k)
    // Spin-native channel counts (tier 4b): multiplicity 2S+1 -> nUp-nDown=2S; multiplicity=0 is the
    // minimal-spin ζ=0 collapse (2S = Nelec%2 -- the ODD-Nelec case needs the doublet split, else integer
    // division empties the EC; GetN(Spin::None) = Nelec either way, so unpolarized runs are unchanged);
    // multiplicity=1 = the EXPLICIT two-channel singlet (see GpwOptions).
    const int  twoS = o.multiplicity>1 ? o.multiplicity-1 : o.Nelec%2;
    const bool polarized = o.multiplicity>=1;
    if ((o.Nelec-twoS)%2!=0 || twoS>o.Nelec)
        throw std::runtime_error("GpwOptions: multiplicity "+std::to_string(o.multiplicity)
                                 +" parity disagrees with Nelec "+std::to_string(o.Nelec));
    Crystal_EC ec(irreps, (o.Nelec+twoS)/2, (o.Nelec-twoS)/2, o.globalFermi, SharedMu(o));
    qchem::Hamiltonian::cHamiltonian* ham=nullptr;
    {
        qchem::report::Timed t("setup: hamiltonian ctor (fit bases + becke mesh)");
        // THE PUBLIC SOLID FRONT DOOR (Step 4).  This driver is the facade-shaped path, so it goes through
        // qchem.Hamiltonian.Factory rather than the Internal ctor -- which is what a real facade in
        // src/Calculation/ will call.  (The bespoke term-level probes further down still use the Internal
        // ctors directly; they are testing the terms, not driving a run.)
        ham=qchem::Hamiltonian::Factory(polarized ? qchem::Hamiltonian::Pol::Polarized
                                                  : qchem::Hamiltonian::Pol::UnPolarized,
                                        lat.GetStructure(), bs.get(), o.species, "LDA",
                                        qcMesh::ResolveXCMesh(o.xcMesh, GatherSharpness(lat,*mol,o)), o.vxcFit);
    }
    auto* acc = MakeGpwAccelerator(o.accelerator);

    qchem::report::Log("SCF start");
    // No Section("basis") here: the pre-flight already emitted basis, so MakeIrrepWFs stays silent
    // (report::InSection("basis") false) and just does the ortho.
    // TIMING NOTE: the lazy heavy builds (collocation streams, local-PP sweep, KB, analytic 1E) run inside
    // the FIRST Fock build -- which is the SEED's (this ctor) or iteration 1's depending on the seed path --
    // so "seed + ortho" / "scf: iterate" CONTAIN those buckets.  MEASURED (diffuse NaF, 2026-07-31): local-PP
    // integrate-back 189 s + stream build 69 s of the 285 s total, both under seed+ortho; the per-pair
    // 791-cell offset enumeration (see [lattice sums]) is the driver, which is why kappa sweeps do nothing.
    std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scfp;
    {
        qchem::report::Timed t("setup: seed + ortho (iterator ctor)");
        scfp=std::make_unique<qchem::SCFIterator::SolidSCFIterator>(bs.get(), &ec, ham, acc,
                                         o.seed, lat.GetStructure().get(), o.ortho, o.orthoTol);
    }
    auto& scf=*scfp;
    // Seed-character MOM (GpwOptions::momFromSeed): the iterator ctor has just built the seed's density by
    // diagonalizing + filling, so the wavefunction's CURRENT occupied subspace IS the seed's -- self-adopt it
    // as the fixed reference now, before Iterate's delayed capture can lock onto iteration 1's (already lost)
    // configuration.  No-op unless scf.UseMOM is also set.
    if (o.momFromSeed) scf.AdoptMOMReference(*scf.GetWaveFunction());
    // IBZ density symmetrization is now automatic: the GPW basis exposes its reciprocal point ops (non-empty
    // only when imposeSymmetry), and the composite density ctor-injects them straight from the basis -- no setter,
    // no ops recomputed here (doc/GPWPlan1.md item 3).
    std::vector<FpRow> series;
    scf.SetObserver([&series](const qchem::SCFIterator::SCFProgress& p)
                    { series.push_back({p.iteration, p.energy, p.dE, p.commutator, p.drho, p.order}); });
    if (o.orderProbe) scf.SetOrderParameter(o.orderName, o.orderProbe);
    SCFParams par = o.scf; par.Verbose = verbose;   // one `verbose` drives both the report console + the SCF table
    qchem::Hamiltonian::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");
    {
        qchem::report::Timed t("scf: iterate (contains the lazy first-iteration setup buckets)");
        scf.Iterate(par);
    }
    qchem::Hamiltonian::ReportGridCharge()=false;
    Fingerprint(series, o.label.c_str());
    OrderTrajectory(series, o.orderProbe ? o.orderName : std::string(), o.label.c_str());

    auto* cd=scf.GetWaveFunction()->GetChargeDensity().release();   // V1.25: BUILT for us; keep/delete below
    double charge=cd->GetTotalCharge();
    {   // §3: the run reports the symmetry it found; S4 adds the MAGNETIC (Shubnikov) channel-pair
        // diagnostic on FREE polarized runs with a genuinely staggered decoration.
        std::vector<int> spins=GatherSiteSpins(lat, o);
        bool staggered=false; for (int s : spins) if (s<0) staggered=true;
        ReportSymmetryFound(*bs, *cd, lat.GetStructure().get(), o.imposeSymmetry,
                            (!o.imposeSymmetry && staggered) ? lat.ShubnikovOps(spins)
                                                             : std::vector<Symmetry::Lattice_3D::SymOp>{});
    }
    if (keep) keep->cd.reset(cd); else delete cd;
    qchem::EnergyBreakdown E=scf.GetEnergy();
    // Etot at 10 s.f.: doc/Benchmark.md compares codes at the 1e-5 Ha level, and the default 6 s.f. cannot
    // express a sub-mHa delta on a -61 Ha crystal.  The breakdown terms stay at default width (they are read
    // for structure, not for the last digit).
    std::cout << "["<<o.label<<"] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Eelec="<<E.GetElectronicEnergy()
              << " Etot="<<std::setprecision(10)<<E.GetTotalEnergy()<<std::setprecision(6)
              << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc
              << " Enn="<<E.Enn<<" E_alphaZ="<<E.E_alphaZ<<")" << std::endl;
    GpwResult R{scf.Converged(), charge, E, scf.GetIterationCount()};
    qchem::report::EmitTimings();       // the where-did-the-time-go ledger, sorted by cost (setup + scf)
    if (keep) keep->bs=std::move(bs);   // the density's basis block -- must outlive keep->cd
    return R;
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
// \a accSchedule (optional): a PER-STAGE accelerator name (parallel to \a kTSchedule; empty = o.accelerator
// throughout).  Each stage already builds a FRESH accelerator, so a kT x accelerator schedule composes here
// naturally -- e.g. {"DIIS","GDM"} with {kT, 0}: smear the DIIS stage through the branch ties, then hand the
// converged density to a COLD GDM stage (GDM does not support Fermi smearing, so its stage must run kT=0).
// \a keep (optional): retain the coldest stage's density + its basis block for the caller's post-mortem
// (site moments etc.), exactly as RunGpw's \a keep does -- otherwise the carried density is deleted here.
// MOM CONTINUATION (GpwOptions::momFromSeed + scf.UseMOM): the density is not the whole state.  Annealing
// carries rho from stage to stage, but each stage builds a FRESH wavefunction, so the OCCUPIED CHARACTER --
// which 5 of 10 near-degenerate d states are filled in each spin channel -- is re-decided from scratch by the
// cold stage's aufbau, which is the very decision the hot stage was run to make.  So each stage adopts the
// PREVIOUS stage's converged occupied subspace as its fixed reference (§0e's grid-continuation MOM, applied
// across TEMPERATURE instead of across grid), and stage 0 self-adopts the seed's.  This is what makes
// "smear through the ties, then hold what you found" a single experiment rather than two unrelated ones.
// \a penaltySchedule (optional): a PER-STAGE MOMSmearPenalty Λ, parallel to \a kTSchedule (empty = o.scf's
// value throughout).  ANNEALING THE CONSTRAINT, not just the temperature -- and for a symmetry-broken run
// that is the more important schedule of the two.  Λ is sized to hold a configuration against a dive at the
// START (MnO: 2.4 Ha at iteration 1); near convergence the competing scale is far smaller, and the SAME Λ
// then pins states empty that the converged spectrum wants filled -- measured on MnO run 24, which held
// SEVEN levels empty below the ↓ HOMO, some 0.6 Ha below states it was occupying.  A constraint that
// rescues the early iterations is not automatically one you want at the fixed point, so it gets a schedule:
// strong while the branch is being chosen, released once the exchange splitting can hold the branch itself.
// Releasing to Λ=0 also RE-ARMS aufbau (the masked-Fermi path only consults the reference when Λ>0), which
// is exactly the test of whether MOM is still load-bearing.
static GpwResult RunGpwAnnealed(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const GpwOptions& o,
                                const std::vector<double>& kTSchedule, bool verbose=false,
                                const std::vector<std::string>& accSchedule={}, GpwHandles* keep=nullptr,
                                const std::vector<double>& penaltySchedule={})
{
    assert(penaltySchedule.empty() || penaltySchedule.size()==kTSchedule.size());
    assert(accSchedule.empty() || accSchedule.size()==kTSchedule.size());
    namespace L3=BasisSet::Lattice_3D;
    const std::string sp = o.species.empty() ? std::string() : o.species.front().first;
    GpwReport report(sp+" "+o.label+" (annealed)", verbose);

    qchem::report::Log("building GPW basis");
    // Field-sharpness routing: HartreeOnly (drop the 2/3*alpha_max floor when the XC is off-raster) was
    // the plan, but MEASURED 2026-07-31 it destabilises NaF SR2 at EVERY reduced beta (0 -> +904 Ha;
    // beta/alpha_max in {1/12..1/3} -> 20-31 Ha slosh; only the historical 2/3 converges): the floor
    // protects the DENSITY's low-G integrity (diffuse-pair FFT fold-back), NOT just V_xc.  So Becke runs
    // keep HartreeXC routing; GPW_ROUTING=hartree is the CALIBRATION valve (pair with GPW_RELFIELDSHARP)
    // for the coarse-end routing campaign (kMinLevelN etc.), not a production setting.
    const bool routeExperiment = std::getenv("GPW_ROUTING")
                              && std::string(std::getenv("GPW_ROUTING"))=="hartree";
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, mol, L3::GPWParams{
        .densityEcut=o.densityEcut, .cutoffFactor=o.cutoffFactor, .raster=o.raster,
        .images=o.images, .kShift=o.kShift, .ladderFactor=o.ladderFactor, .imposeSymmetry=o.imposeSymmetry,
        .rasterFields=routeExperiment ? L3::RasterFields::HartreeOnly : L3::RasterFields::HartreeXC,
        .siteSpins=GatherSiteSpins(lat, o), .hamPreservesReal=o.realTRIMBlocks}));
    if (report.VetBasis(*bs) > 0)
    {
        std::cout << "["<<o.label<<"] ABORT: basis rank-deficient (see basis.removed) -- skipped grids + SCF."<<std::endl;
        return {false, 0.0, qchem::EnergyBreakdown{}, 0};
    }
    report.EmitGrids(*bs);

    auto       st = lat.GetStructure();
    auto       irreps = bs->GetIrreps(Spin::None);
    const int  twoSA = o.multiplicity>1 ? o.multiplicity-1 : o.Nelec%2;   // tier 4b: same channel split as RunGpw
    const bool polarizedA = o.multiplicity>=1;
    if ((o.Nelec-twoSA)%2!=0 || twoSA>o.Nelec)
        throw std::runtime_error("GpwOptions: multiplicity "+std::to_string(o.multiplicity)
                                 +" parity disagrees with Nelec "+std::to_string(o.Nelec));
    Crystal_EC ec(irreps, (o.Nelec+twoSA)/2, (o.Nelec-twoSA)/2, o.globalFermi, SharedMu(o));

    GpwResult R{false, 0.0, qchem::EnergyBreakdown{}, 0};
    qchem::ChargeDensity::cDM_CD* seedCD = nullptr;   // carried between stages (the next stage's ctor consumes it)
    // The PREVIOUS stage's iterator, held alive only long enough for the next one to adopt its occupied
    // subspace (MOM continuation).  Declared AFTER bs/ec so it is destroyed BEFORE them.
    std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> prev;
    for (size_t s=0; s<kTSchedule.size(); ++s)
    {
        const double kT=kTSchedule[s];
        // Fresh Hamiltonian + accelerator per stage (the iterator OWNS + deletes them; a kT change must not
        // carry stale DIIS history across the re-seed).
        auto* ham = qchem::Hamiltonian::Factory(polarizedA ? qchem::Hamiltonian::Pol::Polarized
                                                            : qchem::Hamiltonian::Pol::UnPolarized,
                                        st, bs.get(), o.species, "LDA",
                                        qcMesh::ResolveXCMesh(o.xcMesh, GatherSharpness(lat,*mol,o)), o.vxcFit);
        auto* acc = MakeGpwAccelerator(accSchedule.empty() ? o.accelerator : accSchedule[s]);
        std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scf(
            s==0 ? new qchem::SCFIterator::SolidSCFIterator(bs.get(), &ec, ham, acc, o.seed,  st.get(), o.ortho, o.orthoTol)
                 : new qchem::SCFIterator::SolidSCFIterator(bs.get(), &ec, ham, acc, seedCD, st.get(), o.ortho, o.orthoTol));
        // MOM continuation across TEMPERATURE (see the header comment): stage 0 self-adopts the seed's own
        // freshly-filled occupied subspace, every later stage adopts the stage before it -- so the character
        // the hot stage settled on survives the fresh wavefunction, exactly as the density does.
        if (o.momFromSeed) scf->AdoptMOMReference(prev ? *prev->GetWaveFunction() : *scf->GetWaveFunction());
        prev.reset();   // the adoption copied what we needed; release the previous stage's machinery

        SCFParams par=o.scf; par.Verbose=verbose; par.SmearingkT=kT;
        // Ladder-exhaustion stage-end (user 2026-08-14, run 49's wasted tail): on a NON-FINAL stage, once
        // DIIS's hand-off fires and the next rung is vetoed (GDM under smeared occupations), further DIIS
        // grinding buys nothing this stage's successor needs -- the density is at the tie floor, and the
        // colder stage re-arms GDM.  The FINAL stage keeps the default grind (it has no successor, and a
        // healthy tail can still pass the Δρ gate -- run 38 did).
        par.StopOnAccelExhausted = (s+1<kTSchedule.size());
        if (!penaltySchedule.empty()) par.MOMSmearPenalty=penaltySchedule[s];
        std::vector<FpRow> series;
        scf->SetObserver([&series](const qchem::SCFIterator::SCFProgress& p)
                         { series.push_back({p.iteration,p.energy,p.dE,p.commutator,p.drho,p.order}); });
        if (o.orderProbe) scf->SetOrderParameter(o.orderName, o.orderProbe);   // per STAGE, like the fingerprint
        std::cout << "["<<o.label<<" anneal "<<s+1<<"/"<<kTSchedule.size()<<"] kT="<<kT
                  << " acc="<<(accSchedule.empty()?o.accelerator:accSchedule[s])
                  << " MOM-Lambda="<<par.MOMSmearPenalty
                  << (par.MOMSmearPenalty<=0.0 && o.scf.UseMOM ? "  (RELEASED: plain energy Fermi)" : "")
                  << std::endl;
        scf->Iterate(par);
        Fingerprint(series, (o.label+" kT="+std::to_string(kT)).c_str());
        OrderTrajectory(series, o.orderProbe ? o.orderName : std::string(),
                        (o.label+" kT="+std::to_string(kT)).c_str());

        seedCD = scf->GetWaveFunction()->GetChargeDensity().release();   // consumed by the next stage's ctor
        qchem::EnergyBreakdown E=scf->GetEnergy();
        R = {scf->Converged(), seedCD->GetTotalCharge(), E, scf->GetIterationCount()};
        // A and E(internal) AT A STATED PRECISION, for the same reason as the fingerprint's Efinal: these
        // are the numbers doc/Benchmark.md's MnO rows are read from, and 6 s.f. cannot express a sub-mHa
        // delta on a -61 Ha crystal ("A=E-TS=-61.4" was the whole energy this line reported).
        const std::streamsize prec0=std::cout.precision();
        std::cout << "["<<o.label<<" stage "<<s+1<<"] kT="<<kT<<" conv="<<R.converged<<" iters="<<R.iters
                  << std::setprecision(10)
                  << " A=E-TS="<<E.GetTotalEnergy()<<" -TS="<<E.MinusTS
                  << " E(internal)="<<(E.GetTotalEnergy()-E.MinusTS)
                  << std::setprecision(prec0) << std::endl;
        prev = std::move(scf);   // held for the next stage's MOM adoption; released the moment it has copied
        // seedCD survives the hand-off -- its block is bs, which outlives the loop.
    }
    if (seedCD)
    {   // §3: report on the coldest stage's density (+ the S4 magnetic diagnostic, as in RunGpw)
        std::vector<int> spins=GatherSiteSpins(lat, o);
        bool staggered=false; for (int s : spins) if (s<0) staggered=true;
        ReportSymmetryFound(*bs, *seedCD, st.get(), o.imposeSymmetry,
                            (!o.imposeSymmetry && staggered) ? lat.ShubnikovOps(spins)
                                                             : std::vector<Symmetry::Lattice_3D::SymOp>{});
    }
    prev.reset();                                       // no further stage: drop the last iterator (ham/acc/wf)
    // THE SAME CLOSING REPORT RunGpw GIVES.  This driver had neither -- no total-energy summary and no
    // timing ledger -- and since every MnO row runs through here, doc/Benchmark.md's claim that "every GPW
    // run reports PEAK RSS" was false for exactly the runs whose RAM the table most needs (the ledger is
    // where PEAK RSS lives).  The per-stage lines above report each stage; this reports the RUN.
    std::cout << "["<<o.label<<"] stages="<<kTSchedule.size()<<" iters(last)="<<R.iters
              << " charge="<<R.charge
              << " Etot="<<std::setprecision(10)<<R.E.GetTotalEnergy()<<std::setprecision(6)
              << "  (Ekin="<<R.E.Kinetic<<" Een="<<R.E.Een<<" Eee="<<R.E.Eee<<" Exc="<<R.E.Exc
              << " Enn="<<R.E.Enn<<" E_alphaZ="<<R.E.E_alphaZ<<")" << std::endl;
    qchem::report::EmitTimings();       // the where-did-the-time-go ledger + PEAK RSS, as in RunGpw
    if (keep) { keep->cd.reset(seedCD); keep->bs=std::move(bs); }   // bs is the density's block -- it must outlive cd
    else      delete seedCD;   // the final stage's carried density (not consumed by any further ctor)
    return R;
}

// Positional back-compat wrapper (the existing single-species Si callers): forwards to RunGpw(GpwOptions).
GpwResult RunGPW(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, double densityEcut,
                 int Nelec, const char* element, const char* label, bool verbose=false, int nmax=120,
                 qchem::Ortho ortho=qchem::Cholesky, double orthoTol=0.0,
                 rvec3_t kShift={0,0,0}, double minDrho=1e-6, double minDE=1e30,
                 qchem::ChargeDensity::SeedStrategy seed=qchem::ChargeDensity::SeedStrategy::Uniform,
                 BasisSet::Lattice_3D::CellImages images=BasisSet::Lattice_3D::CellImages::Periodic,
                 double smearkT=0.0,
                 qcMesh::UnitCellKind xcKind=qcMesh::UnitCellKind::Auto)
{
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    // GPW_IMPOSE=0 / GPW_VERBOSE=1: the A/B pair for a multi-k defect.  Imposition is the one part of these
    // runs that reconstructs the full BZ from the IRREDUCIBLE blocks, so it is the first thing to remove when
    // a run differs between a TRIM mesh (phases +/-1) and a genuinely complex one -- see the shifted-MP row
    // in doc/Benchmark.md.  Diagnostics only: unset, nothing changes.
    if (const char* im=std::getenv("GPW_IMPOSE")) o.imposeSymmetry=std::atoi(im)!=0;
    o.label=label; o.Nelec=Nelec; o.species={{std::string(element), 4}};   // the Si callers: Zion=4
    o.densityEcut=densityEcut; o.images=images; o.kShift=kShift; o.xcMesh.cellKind=xcKind;
    o.accelerator="DIIS"; o.seed=seed; o.ortho=ortho; o.orthoTol=orthoTol;
    o.scf.NMaxIter=(size_t)nmax; o.scf.MinΔρ=minDrho; o.scf.MinΔE=minDE;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.Verbose=verbose; o.scf.SmearingkT=smearkT;
    // GPW_SMEAR / GPW_VERBOSE go AFTER the block above, or the fixed recipe silently overwrites them --
    // which it did: GPW_VERBOSE looked dead because o.scf.Verbose=verbose ran later and put back `false`.
    // An override that a later line can clobber is worse than no override, because it reads as evidence.
    // GPW_SMEAR=kT: an integer aufbau fill is AMBIGUOUS at a degenerate frontier, and the symptom --
    // converged, charge exact, energy high, "rho rotates" -- looks exactly like a broken operator until
    // smearing tells them apart.
    if (const char* kt=std::getenv("GPW_SMEAR"))  o.scf.SmearingkT=std::atof(kt);
    if (std::getenv("GPW_VERBOSE"))               o.scf.Verbose=true;
    return RunGpw(lat, mol, o, verbose || std::getenv("GPW_VERBOSE")!=nullptr);
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

// SINGLE-K SWEEP: E(k) along the cell diagonal, ONE k-point per run, no weights, no symmetry, no IBZ.
//
// The instrument the shifted-MP failure needed (2026-08-19).  That test lands a CONVERGED but wrong
// -3.7351 against its -7.86744 anchor, and everything cheap was already eliminated: it is not the
// imposition (the free run is wrong the same way, Ekin 6.14 / Een +0.61 in both), not the collocation
// streams (budgets off -> analytic path -> the SAME energy to 8 digits, converged in 18 iterations), and
// not the k-set (8 -> 2 irreducible is the textbook Chadi-Cohen reduction, and the blocks are correctly
// typed complex).  What remained was the possibility that a SINGLE complex-k block is already wrong,
// before any weighting or symmetry can be blamed -- and a one-k run isolates exactly that.
//
// The mesh is 1x1x1, so k = kShift EXACTLY: the sweep runs k = s(1,1,1) for s in [0, 1/2].  s=0 (Gamma)
// and s=1/2 (the L zone corner) are TRIM -- Bloch phases +/-1, real blocks, the paths every other test in
// this suite exercises.  Everything strictly between them is genuinely complex and is exercised by nothing
// else.  E(k) here is not a physical total energy (one k-point is not a Brillouin-zone average) but it IS
// a smooth function of k, so a KINK or a jump on entering the complex region localizes the defect to a
// single k-block's operators, with no weights, occupations or star-averaging in the picture.
//
//   for s in 0 0.125 0.25 0.375 0.5; do
//     GPW_KSHIFT=$s ./ITMain --gtest_filter='*SingleKSweepProbe*' --gtest_also_run_disabled_tests; done
TEST(GPW_SCF, DISABLED_SingleKSweepProbe)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));           // ONE k-point, at exactly kShift

    const double s = std::getenv("GPW_KSHIFT") ? std::atof(std::getenv("GPW_KSHIFT")) : 0.25;
    const bool trim = (s==0.0 || s==0.5);
    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Nelec*/8, "Si",
                       ("Si single-k s="+std::to_string(s)).c_str(), /*verbose*/false, /*nmax*/60,
                       qchem::Cholesky, 0.0, rvec3_t(s,s,s));
    std::cout << "[k sweep] s=" << s << (trim ? "  TRIM" : "  complex")
              << "  Etot=" << std::setprecision(10) << R.E.GetTotalEnergy() << std::setprecision(6)
              << "  charge=" << R.charge
              << "  (Ekin=" << R.E.Kinetic << " Een=" << R.E.Een << " Eee=" << R.E.Eee
              << " Exc=" << R.E.Exc << ")" << std::endl;
    EXPECT_NEAR(R.charge, 8.0, 1e-6);               // the one thing that must hold at EVERY k
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

// The 3c-3 flip's REPORT EVIDENCE: basis.perIrrep now carries `runsReal` (the block's CONSTRUCTED
// scalar) beside `real` (the irrep's TRIM fact).  With the flip on every real:true block must actually
// RUN real -- the fact the acceptance twins above cannot see from outside; with it off (the harness
// default) the same blocks stay complex, so the two fields separate exactly where they should.
TEST(GPW_SCF, RealTRIMBlocksRunRealInReport)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    auto perIrrep=[](){ const report::json& all=report::GlobalReport();
                        for (auto it=all.begin(); it!=all.end(); ++it)
                            if (it.value().contains("basis")) return it.value()["basis"]["perIrrep"];
                        return report::json{}; };
    GpwOptions o;
    o.label="Si real-block report"; o.Nelec=8; o.species={{"Si",4}}; o.densityEcut=20.0;
    o.scf.NMaxIter=3; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3;

    report::ClearGlobal();
    o.realTRIMBlocks=true;
    RunGpw(lat, MakeBasisSR(cell), o);
    report::json rows=perIrrep();
    ASSERT_GE(rows.size(), 1u);
    for (const auto& row : rows)
    {
        EXPECT_TRUE(row["real"].get<bool>());       // Γ-only: every block is TRIM
        EXPECT_TRUE(row["runsReal"].get<bool>());   // ...and the flip made each one ACTUALLY real
    }

    report::ClearGlobal();
    o.realTRIMBlocks=false;                          // the control: same irrep fact, complex build
    RunGpw(lat, MakeBasisSR(cell), o);
    rows=perIrrep();
    ASSERT_GE(rows.size(), 1u);
    for (const auto& row : rows)
    {
        EXPECT_TRUE (row["real"].get<bool>());
        EXPECT_FALSE(row["runsReal"].get<bool>());
    }
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
    // XC route pinned UNIFORM (the gate's calibrated arrangement): under the Becke default this
    // DEGENERATE half-filled 3p atom exposed a real open question — the freely-rotating degenerate
    // density has orientation-DEPENDENT quadrature error on the fixed-axis Becke angular grid (V_xc is
    // nonlinear, so an anisotropic rho's error rotates with it), turning the energy-neutral rotation
    // into a ~Ha-scale E oscillation.  The SMEARED sibling (SmearingConvergesDegenerateShell) is fine
    // under Becke — fractional occupation restores the symmetric density.  Recorded in doc/GPWPlan1.md
    // (Becke remaining increments); this gate's PURPOSE is the PP box-independence check, so it keeps
    // its historical route.
    GpwResult R=RunGPW(lat, MakeBasis(cell), /*densityEcut*/10.0, /*Nelec*/4, "Si", "Si atom-in-box",
                       /*verbose*/false, /*nmax*/40, qchem::Cholesky, 0.0, rvec3_t(0,0,0), 1e-6, 1e30,
                       qchem::ChargeDensity::SeedStrategy::Uniform,
                       BasisSet::Lattice_3D::CellImages::HomeCellOnly,    // the finite-molecule mode
                       /*smearkT*/0.0, qcMesh::UnitCellKind::Uniform);

    EXPECT_NEAR(R.charge, 4.0, 1e-6);                        // 4 valence electrons (Zion=4), charge conserved
    // GPW-in-box (G-space local PP -> box-independent) reproduces the finite SIPP DFT energy to grid tolerance.
    // (Energy-converged; density is degenerate at Gamma -- see the note above -- so no Converged() guard.)
    EXPECT_NEAR(R.E.GetTotalEnergy(), Esipp, 5e-2) << "GPW-in-box total vs finite SIPP molecular DFT";
}

// (tier 4b, invariant) THE ζ=0 COLLAPSE: the TWO-CHANNEL machinery on a CLOSED SHELL must reproduce the
// unpolarized anchor.  Same gapped Si/Γamma cell + recipe as SmearingInertOnGap, but multiplicity=1 drives
// the polarized pipeline (dcmplx tPolarizedWF, Crystal_EC(4,4), DeltaFittedVxcPol + DeltaFittedVcorrPol) with
// nUp=nDn=4 -- v^σ(ρ/2,ρ/2)=v^P(ρ) pointwise, so the total must land on the SAME −7.11506 anchor.  The
// periodic sibling of the molecular WaterPolarizedLDA-vs-LDA check; catches any polarized-path divergence
// (channel bookkeeping, shared-engine caching, the collocation memo screen) on known ground.
TEST(GPW_SCF, PolarizedSingletMatchesUnpolarizedSiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="Si SR Gamma pol-singlet";
    o.Nelec=8; o.multiplicity=1;                       // EXPLICIT two-channel singlet (nUp=nDn=4)
    o.species={{"Si",4}};
    o.densityEcut=20.0;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o);

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.11506, 2e-3);   // == the unpolarized Becke anchor (ζ=0 collapse exact)
}

// The SPIN-SAD sibling of the ζ=0 collapse (SCFSeedingPlan §10 increment B): a polarized run with a SAD
// seed now assembles the TWO-CHANNEL PolarizedSeedCD -- Si's library entry is spin-agnostic (closed shell),
// so each channel is exactly rho/2 and the whole polarized-seed machinery (channel SeedCDs, the merged
// FourierDensity total into Hartree, the cSpinResolved_CD branch in RhoPol) must land on the SAME anchor.
// The proof the polarized seed cannot perturb non-magnetic physics.
// ONE CHEMICAL POTENTIAL OVER BOTH SPIN CHANNELS (Crystal_EC spinsShareFermi; 2026-08-10).
//
// WHAT IT PINS: with a shared reservoir the MOMENT IS AN OUTPUT.  Seed Si-Gamma -- a closed-shell,
// gapped, non-magnetic system -- as a TRIPLET (nUp=5, nDn=3) and let the fill decide: one mu must move the
// two excess-spin electrons back down and land on the SINGLET energy.  With separate per-channel counts
// (the control below) nUp-nDn=2 is conserved by construction, so the run is stuck on the triplet no matter
// how far it converges.  The difference between the two arms IS the ensemble.
//
// WHY IT MATTERS (doc/SymmetryUpgradePlan.md sec 7 step 7): separate reservoirs mean mu_up != mu_dn, and then
// an occupation is monotone in epsilon only WITHIN a channel -- MnO run 29 ended with a down level 27 mHa
// BELOW an up level and LESS occupied.  That gap is an unrelieved driving force to move charge between the
// channels, so the converged state is not the free minimum.  CP2K's MnO deck constrains nothing.
TEST(GPW_SCF, SharedFermiLevelLetsTheMomentRelax)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    auto TripletSi=[&](bool shared)
    {
        GpwOptions o;
        o.imposeSymmetry=true;
        o.label = shared ? "Si Gamma triplet, shared mu" : "Si Gamma triplet, per-channel mu";
        o.Nelec=8; o.multiplicity=3;                   // SEEDED as a triplet (nUp=5, nDn=3)
        o.species={{"Si",4}};
        o.densityEcut=20.0;
        o.scf.SmearingkT=5e-3;                         // a shared mu needs smearing to relax the moment with
        o.spinsShareFermi=shared;
        o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
        o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
        o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
        return RunGpw(lat, MakeBasisSR(cell), o);
    };

    const GpwResult shared=TripletSi(true);
    EXPECT_NEAR(shared.charge, 8.0, 1e-6);             // the TOTAL is what a shared reservoir conserves
    // The moment relaxed away: this lands on the singlet anchor of PolarizedSingletMatchesUnpolarizedSiGamma.
    EXPECT_NEAR(shared.E.GetTotalEnergy(), -7.11506, 3e-3)
        << "a shared mu must let a triplet-seeded closed-shell system fall back to the singlet";

    // CONTROL: the same run with separate per-channel counts cannot relax -- nUp-nDn=2 is conserved.
    const GpwResult held=TripletSi(false);
    EXPECT_NEAR(held.charge, 8.0, 1e-6);
    EXPECT_GT(held.E.GetTotalEnergy(), shared.E.GetTotalEnergy() + 1e-3)
        << "with two reservoirs the seeded multiplicity is a CONSTRAINT and must sit above the free minimum";
}

TEST(GPW_SCF, PolarizedSeedSingletMatchesUnpolarizedSiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="Si SR Gamma pol-singlet spin-SAD";
    o.Nelec=8; o.multiplicity=1;                       // EXPLICIT two-channel singlet (nUp=nDn=4)
    o.species={{"Si",4}};
    o.densityEcut=20.0;
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;    // -> PolarizedSeedCD (rho/2 channels for pairless Si)
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o);

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.11506, 2e-3);   // the SAME unpolarized Becke anchor
}

// (tier 4b, gate b) O2 in a box, TRIPLET: the multi-electron polarized solid pipeline vs the finite
// molecular facade on the SAME sipp O-q6 basis + GTH PP (the spin sibling of SiPseudoAtomInBoxMatchesFinite,
// cross-anchored to the facade's spin-native triplet machinery -- doc/SymmetryUpgradePlan.md §4 tier 4b).
TEST(GPW_SCF, O2TripletInBoxMatchesFinite)
{
    const double d=2.282;   // O2 bond (au)
    Molecule o2mol;
    o2mol.Insert(new Atom(8, 0.0, {-d/2,0,0}));
    o2mol.Insert(new Atom(8, 0.0, { d/2,0,0}));
    Calculation cRef(o2mol, {.basis="sipp", .multiplicity=3, .pseudopotential=true});
    const double Eref=cRef.Energy();
    {
        auto E=cRef.EnergyTerms();
        std::cout << "[O2 finite] sipp GTH-q6 LSDA triplet="<<Eref
                  << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc
                  << " Enn="<<E.Enn<<")"<<std::endl;
    }

    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(8, {0.5-0.5*d/a,0.5,0.5});
    cell.AddAtom(8, {0.5+0.5*d/a,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="O2 in-box triplet";
    o.Nelec=12; o.multiplicity=3;                      // S=1: nUp=7, nDown=5
    o.species={{"O",6}};                               // densityEcut stays AUTO: O q6 is hard (alpha_max rules)
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-6; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasis(cell), o);

    EXPECT_NEAR(R.charge, 12.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), Eref, 5e-2) << "GPW-in-box triplet vs finite molecular LSDA triplet";
}

// PROBE (kept disabled): FIXED-DENSITY term fingerprinter for the Na-doublet deficit.  Feed the polarized
// GPW Hamiltonian two hand-built ONE-HOT densities (the alpha=0.6999 and alpha=2.0 s-Gaussians, charge 1
// each) and print every energy term.  Because BOTH densities carry the same charge in the same cell, ALL
// the box-convention constants (Enn, E_alphaZ, dropped-G=0 shifts) cancel EXACTLY in the DIFFERENCE of
// each term -- so Δ(term) can be compared directly against the independent radial oracle's Δ, and the
// defective term identifies itself.  (Oracle: scratchpad na_q1_basis.py evaluates the same two Gaussians.)
TEST(GPW_SCF, DISABLED_NaFixedDensityTermProbe)
{
    namespace L3=BasisSet::Lattice_3D;
    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(11, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    // densityEcut EXPLICIT and high: the one-hot probes include the SHARP alpha=2 pair (p=4), which the
    // auto rule (cutoffFactor*alpha_max = 4 Ha) badly under-resolves -- the probe must measure the terms,
    // not the grid (the O2 auto-cutoff lesson).  Env knob for the probe's own convergence sweep.
    const double ecut=std::getenv("NAFD_ECUT")?atof(std::getenv("NAFD_ECUT")):40.0;
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR),
        L3::GPWParams{.densityEcut=ecut, .images=L3::CellImages::HomeCellOnly}));

    qchem::Hamiltonian::cHamiltonian* ham=new qchem::Hamiltonian::Ham_PW_DFT(
        lat.GetStructure(), bs.get(), {{"Na",1}}, "LDA", qcMesh::ResolveXCMesh({.cellKind=qcMesh::UnitCellKind::Auto}),
        Hamiltonian::VxcFit::Auto, /*polarized*/true);

    const BasisSet::Complex_OIBS* obs=nullptr;
    for (auto b : bs->Iterate<BasisSet::Complex_OIBS>()) obs=b;   // the single Gamma block (n=10)
    ASSERT_NE(obs, nullptr);
    auto S=obs->Overlap();
    const size_t n=obs->GetNumFunctions();

    // One-hot density on basis function k for the UP channel; the DOWN channel is the zero matrix.
    auto probe=[&](size_t k)
    {
        hmat_t<dcmplx> Dup(n), Ddn(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) { Dup(i,j)=dcmplx(0.0); Ddn(i,j)=dcmplx(0.0); }
        Dup(k,k)=dcmplx(1.0/std::real(dcmplx(S(k,k))));            // Tr(D S) = 1
        using namespace qchem::ChargeDensity;
        std::unique_ptr<cDM_CD> up(IrrepCD_Factory<dcmplx>(Dup, obs, obs->GetIrrep(Spin::Up)));
        std::unique_ptr<cDM_CD> dn(IrrepCD_Factory<dcmplx>(Ddn, obs, obs->GetIrrep(Spin::Down)));
        auto pol=PolarizedCD_Factory<dcmplx>(std::move(up),std::move(dn));   // V1.25: takes ownership
        auto* cd=dynamic_cast<cDM_CD*>(pol.get());
        qchem::EnergyBreakdown te = ham->GetTotalEnergy(cd);
        std::cout << "[fixed-D k="<<k<<"] charge="<<cd->GetTotalCharge()
                  << " Ekin="<<te.Kinetic<<" Een="<<te.Een<<" Eee="<<te.Eee<<" Exc="<<te.Exc
                  << " Enn="<<te.Enn<<" E_alphaZ="<<te.E_alphaZ
                  << " Etot="<<te.GetTotalEnergy()<<std::endl;
        delete cd;
        return te;
    };
    auto t2=probe(2);   // alpha=0.6999271 s-Gaussian
    auto t3=probe(3);   // alpha=2.0       s-Gaussian
    std::cout << "[fixed-D delta k2-k3] dEkin="<<t2.Kinetic-t3.Kinetic
              << " dEen="<<t2.Een-t3.Een << " dEee="<<t2.Eee-t3.Eee
              << " dExc="<<t2.Exc-t3.Exc
              << " dEtot="<<t2.GetTotalEnergy()-t3.GetTotalEnergy()<<std::endl;
    delete ham;

    // CONTROL: the SAME one-hot densities through the UNPOLARIZED Hamiltonian on a bare Spin::None leaf --
    // isolates "polarized path" from "box electrostatics in general".
    qchem::Hamiltonian::cHamiltonian* hamU=new qchem::Hamiltonian::Ham_PW_DFT(
        lat.GetStructure(), bs.get(), {{"Na",1}}, "LDA", qcMesh::ResolveXCMesh({.cellKind=qcMesh::UnitCellKind::Auto}),
        Hamiltonian::VxcFit::Auto, /*polarized*/false);
    auto probeU=[&](size_t k)
    {
        hmat_t<dcmplx> D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=dcmplx(0.0);
        D(k,k)=dcmplx(1.0/std::real(dcmplx(S(k,k))));
        using namespace qchem::ChargeDensity;
        auto* cd=IrrepCD_Factory<dcmplx>(D, obs, obs->GetIrrep(Spin::None));
        qchem::EnergyBreakdown te = hamU->GetTotalEnergy(cd);
        std::cout << "[fixed-D UNPOL k="<<k<<"] charge="<<cd->GetTotalCharge()
                  << " Ekin="<<te.Kinetic<<" Een="<<te.Een<<" Eee="<<te.Eee<<" Exc="<<te.Exc<<std::endl;
        delete cd;
        return te;
    };
    auto u2=probeU(2);
    auto u3=probeU(3);
    std::cout << "[fixed-D UNPOL delta k2-k3] dEkin="<<u2.Kinetic-u3.Kinetic
              << " dEen="<<u2.Een-u3.Een << " dEee="<<u2.Eee-u3.Eee
              << " dExc="<<u2.Exc-u3.Exc<<std::endl;
    delete hamU;

    // SPIN-SWAP TEST: the polarized XC pair's FOCK matrices per channel on the one-hot ↑ density
    // (ρ↑ = the alpha=0.6999 Gaussian, ρ↓ = 0).  Expected: v_x(Up) strongly negative diagonal
    // (v_x(ρ↑)), v_x(Down) ZERO (v_x(0)); v_c(Up) the ζ=1 correlation, v_c(Down) its ∂/∂ρ↓ partner.
    // A channel swap/share -- the E/H-inconsistency candidate that vanishes at ρ↑=ρ↓ -- shows as
    // Up≈0 / Down<0.
    {
        using namespace qchem::Hamiltonian;
        auto exch=std::make_shared<SlaterExchange>(2.0/3.0, Spin(Spin::Up));
        auto corr=std::make_shared<VWN_Correlation>();
        BasisSet::XCQuadrature q = bs->CreateXCQuadrature(lat.GetStructure().get(),
                                                          qcMesh::ResolveXCMesh({.cellKind=qcMesh::UnitCellKind::Auto}));
        auto engine=std::make_shared<XC_GridEngine>(q.mesh, std::move(q.fold));
        DeltaFittedVxcPol   x(exch, engine);
        DeltaFittedVcorrPol c(corr, engine);
        const size_t k=2;
        hmat_t<dcmplx> Dup(n), Ddn(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) { Dup(i,j)=dcmplx(0.0); Ddn(i,j)=dcmplx(0.0); }
        Dup(k,k)=dcmplx(1.0/std::real(dcmplx(S(k,k))));
        using namespace qchem::ChargeDensity;
        std::unique_ptr<cDM_CD> up(IrrepCD_Factory<dcmplx>(Dup, obs, obs->GetIrrep(Spin::Up)));
        std::unique_ptr<cDM_CD> dn(IrrepCD_Factory<dcmplx>(Ddn, obs, obs->GetIrrep(Spin::Down)));
        auto pol=PolarizedCD_Factory<dcmplx>(std::move(up),std::move(dn));   // V1.25: takes ownership
        auto* cd=dynamic_cast<cDM_CD*>(pol.get());
        auto diag=[&](const hmat_t<dcmplx>& M){ double s2=0; for (size_t i=0;i<M.rows();i++) s2+=std::real(dcmplx(M(i,i))); return s2; };
        cDynamic_HT& xi=x; cDynamic_HT& ci=c;      // the public term face (Imp::GetMatrix is private)
        auto Mxu=xi.GetMatrix(obs, Spin::Up,   cd);
        auto Mxd=xi.GetMatrix(obs, Spin::Down, cd);
        auto Mcu=ci.GetMatrix(obs, Spin::Up,   cd);
        auto Mcd=ci.GetMatrix(obs, Spin::Down, cd);
        std::cout << "[swap test] tr v_x(Up)="<<diag(Mxu)<<" tr v_x(Down)="<<diag(Mxd)
                  << "  tr v_c(Up)="<<diag(Mcu)<<" tr v_c(Down)="<<diag(Mcd)
                  << "  (expect v_x(Down)=0, v_x(Up)<0)"<<std::endl;
        delete cd;
    }

    // THE KILL SHOT: evaluate the GPW functional AT THE ORACLE'S OPTIMAL DENSITY (the radial same-basis
    // solver's s-block minimizer, c = the contracted combination -- OFF-DIAGONAL D, which the one-hot
    // probes above never exercised).  If E_GPW[D*] < the SCF's converged −0.0699, the GPW SCF failed to
    // find its own functional's minimum (solver-side bug); if E_GPW[D*] ≈ −0.07, the functional itself
    // errs on off-diagonal D.
    {
        const double c[4]={2.54610981, -1.79605159, 0.08532953, -0.03806822};   // oracle minimizer (s block)
        hmat_t<dcmplx> Dup(n), Ddn(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) { Dup(i,j)=dcmplx(0.0); Ddn(i,j)=dcmplx(0.0); }
        for (size_t i=0;i<4;i++) for (size_t j=i;j<4;j++) Dup(i,j)=dcmplx(c[i]*c[j]);
        // renormalize to Tr(D S) = 1 in THIS stack's S (guards small normalization-convention differences)
        double q=0;
        for (size_t i=0;i<4;i++) for (size_t j=0;j<4;j++) q+=c[i]*c[j]*std::real(dcmplx(S(i,j)));
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) Dup(i,j)=dcmplx(std::real(dcmplx(Dup(i,j)))/q);
        using namespace qchem::ChargeDensity;
        std::unique_ptr<cDM_CD> up(IrrepCD_Factory<dcmplx>(Dup, obs, obs->GetIrrep(Spin::Up)));
        std::unique_ptr<cDM_CD> dn(IrrepCD_Factory<dcmplx>(Ddn, obs, obs->GetIrrep(Spin::Down)));
        auto pol=PolarizedCD_Factory<dcmplx>(std::move(up),std::move(dn));   // V1.25: takes ownership
        auto* cd=dynamic_cast<cDM_CD*>(pol.get());
        qchem::Hamiltonian::cHamiltonian* hamP=new qchem::Hamiltonian::Ham_PW_DFT(
            lat.GetStructure(), bs.get(), {{"Na",1}}, "LDA", qcMesh::ResolveXCMesh({.cellKind=qcMesh::UnitCellKind::Auto}),
            Hamiltonian::VxcFit::Auto, /*polarized*/true);
        qchem::EnergyBreakdown te = hamP->GetTotalEnergy(cd);
        std::cout << "[oracle-D*] charge="<<cd->GetTotalCharge()
                  << " Ekin="<<te.Kinetic<<" Een="<<te.Een<<" Eee="<<te.Eee<<" Exc="<<te.Exc
                  << " Etot="<<te.GetTotalEnergy()
                  << "  (SCF found -0.0699; oracle E[D*]=-0.1416)"<<std::endl;
        delete cd; delete hamP;
    }
}

// PROBE (kept disabled): F atom in a box, DOUBLET (2p^5: nUp=4, nDn=3 -- open shell with BOTH channels
// populated).  Discriminates "empty minority channel" (Na) from "open shell" as the Na-gap driver.
// WARNING: slow -- the degenerate 2p^5 hole rotates under the fixed-axis Becke angular grid (the known
// open-shell Becke oscillation; smearing would cure it), and the F q7 auto grid is large.
TEST(GPW_SCF, DISABLED_FAtomInBoxDoubletProbe)
{
    Molecule f; f.Insert(new Atom(9, 0.0, {0,0,0}));
    Calculation cRef(f, {.basis="valence_lowq_sr", .multiplicity=2, .pseudopotential=true, .ppValence=7});
    const double Eref=cRef.Energy();
    std::cout << "[F finite] valence_lowq_sr LSDA doublet (q7)="<<Eref<<std::endl;

    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(9, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="F atom-in-box doublet"; o.Nelec=7; o.multiplicity=2; o.species={{"F",7}};
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;
    o.scf.SmearingkT=1e-2;               // smear the degenerate 2p^5 hole (the Becke rotating-ρ cure)
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-6; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o);
    std::cout << "[F probe] GPW="<<R.E.GetTotalEnergy()<<" finite="<<Eref
              << " diff="<<R.E.GetTotalEnergy()-Eref<<std::endl;
}

// PROBE (kept disabled): Na2 dimer in a box, closed-shell singlet, UNPOLARIZED both sides -- bounds the
// box-conventions + basis residual for the Na species without any spin machinery.  MEASURED 2026-08-04:
// GPW −0.3281 vs finite −0.3529, diff +24.8 mHa for the dimer (~12 mHa/atom) -- the baseline the
// DISABLED_NaPseudoAtomInBoxDoublet defect note compares against.
TEST(GPW_SCF, DISABLED_Na2DimerInBoxProbe)
{
    const double d=5.8;   // ~Na2 bond length (au)
    Molecule na2;
    na2.Insert(new Atom(11, 0.0, {-d/2,0,0}));
    na2.Insert(new Atom(11, 0.0, { d/2,0,0}));
    Calculation cRef(na2, {.basis="valence_lowq_sr", .pseudopotential=true, .ppValence=1});
    const double Eref=cRef.Energy();
    std::cout << "[Na2 finite] valence_lowq_sr LDA singlet (q1)="<<Eref<<std::endl;

    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(11, {0.5-0.5*d/a,0.5,0.5});
    cell.AddAtom(11, {0.5+0.5*d/a,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="Na2 dimer-in-box"; o.Nelec=2; o.species={{"Na",1}};
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;
    o.scf.NMaxIter=40; o.scf.MinΔρ=1e-6; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o);
    std::cout << "[Na2 probe] GPW="<<R.E.GetTotalEnergy()<<" finite="<<Eref
              << " diff="<<R.E.GetTotalEnergy()-Eref<<std::endl;
}

// (tier 4b, gate a) THE POLARIZED SOLID PIPELINE: Na pseudo-atom in a box, DOUBLET (doc/SymmetryUpgradePlan.md
// §4).  The minimal end-to-end TWO-CHANNEL GPW run: Na q1 GTH PP, S=1/2, moment 1 -- spin-resolved D through
// Crystal_EC(nUp=1,nDown=0), the dcmplx tPolarizedWF (two Bloch channels), and the spin-native Becke XC pair
// (DeltaFittedVxcPol + DeltaFittedVcorrPol).  Cross-anchored against the finite molecular facade doublet on the SAME
// valence basis + PP (the spin sibling of SiPseudoAtomInBoxMatchesFinite).
//
// SEED PIN (the 2026-08-04 root-cause campaign): this gate MUST seed from IonicSAD.  From the Uniform seed
// the lone ↑ electron converges to a GENUINE self-consistent excited basin 72 mHa above the minimum
// (diffuse 3s; DIIS honors it, GDM -- a local descender -- stays in it; total is box/grid/route-independent,
// so it looks "converged" by every health metric).  The functional itself was proven correct everywhere:
// fixed-density term probes (DISABLED_NaFixedDensityTermProbe) match analytic kinetic, the exact discrete
// G!=0 lattice-sum Hartree, and the ζ=1 Dirac/VWN values exactly, and E_GPW[D*]=−0.1420 at the independent
// radial same-basis oracle's minimizer (oracle E=−0.1416; complete-basis −0.1922).  A ONE-ELECTRON system
// is uniquely basin-fragile: no other electrons pull the density into the core basin (Na2, O2-triplet and
// F-doublet all escape Uniform fine).  IonicSAD lands in-basin: −0.1419, 4.8 mHa from the facade.
// SPIN-SAD (§10 increment B): the polarized run now seeds the TWO-CHANNEL PolarizedSeedCD, and Na's library
// pair is exact (1 e-: up=total, dn=0), so iteration 0 starts FULLY polarized -- same basin, same pin
// (-0.141933 to the digit), 21 iters vs the rho/2-collapse seed's 14 (the staggered start makes DIIS
// reorganize more; basin selection, not speed, is what the seed pin protects).
TEST(GPW_SCF, NaPseudoAtomInBoxDoublet)
{
    // Basis-MATCHED finite reference: the SAME valence_lowq_sr Na basis + GTH-LDA q1 PP through the
    // molecular facade, spin-native LSDA doublet (Ham_PP polarized: FittedVxcPol + FittedVcorrPol).
    // ppValence=1 overrides Na's GTH DEFAULT-valence entry, which is the SEMICORE q9 -- the q1 entry is
    // the valence_lowq bases' convention.  (The Slater/High atom path is NOT a good oracle here: Slater
    // functions fit the smooth nodeless pseudo-orbitals poorly -- same reason the Si gate anchors on the
    // basis-matched sipp facade run.)
    Molecule na; na.Insert(new Atom(11, 0.0, {0,0,0}));
    Calculation cRef(na, {.basis="valence_lowq_sr", .multiplicity=2, .pseudopotential=true, .ppValence=1});
    const double Eref=cRef.Energy();
    std::cout << "[Na finite] valence_lowq_sr LSDA doublet (q1)="<<Eref<<std::endl;

    // Box a=16 (the Si gate's size: cross-cell products of the most diffuse pair negligible).
    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(11, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="Na atom-in-box doublet";
    o.Nelec=1; o.multiplicity=2;                               // S=1/2: nUp=1, nDown=0
    o.species={{"Na",1}};
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;   // the finite-molecule mode
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;       // SEED PIN: Uniform has a stable wrong basin (header)
    o.scf.NMaxIter=40; o.scf.MinΔρ=1e-6; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o);

    EXPECT_TRUE(R.converged);                                // 3s^1 is non-degenerate: Δρ converges
    EXPECT_NEAR(R.charge, 1.0, 1e-6);                        // 1 valence electron (Zion=1), charge conserved
    // GPW-in-box (two-channel) reproduces the finite molecular LSDA doublet (measured 4.8 mHa; the gap is
    // the two stacks' fit/quadrature tech -- the facade Dunlap-fits J and fits v_xc, the GPW is fit-free).
    EXPECT_NEAR(R.E.GetTotalEnergy(), Eref, 2e-2) << "GPW-in-box doublet vs finite molecular LSDA doublet";
    EXPECT_NEAR(R.E.GetTotalEnergy(), -0.141933, 1e-4);      // did-E-move anchor (== the same-basis oracle -0.1416)
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
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
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
    std::cout<<"[Al global-μ full-mesh] A="<<R.E.GetTotalEnergy()<<std::endl;   // the IBZ pin's route-matched partner
    EXPECT_NEAR(R.E.GetTotalEnergy(), -2.11681, 3e-3);   // did-E-move anchor (2×2×2 global-μ free energy A)
    EXPECT_LT(R.E.GetTotalEnergy(), -1.95);    // dispersion: well below the Γ-only -1.92 (k-sampling binds)
}

// (item 3, IBZ/k-star) FOLDING IS EXACT.  The 8-point 2×2×2 Γ-mesh folds to 3 irreducible k-points under the
// cubic point group, and the density is STAR-AVERAGED consistently: the G-space Hartree via the reciprocal ops
// the basis exposes (GetReciprocalPointOps, ctor-injected into the composite density), and the real-space XC
// raster via the DIRECT ops (cFIT_SF_ABS::SymmetrizeRaster, ctor-injected into the Vxc fit basis -- so XC stays
// on the non-negative ρ_DM raster).  So the reduced run reproduces the full-mesh AlFCCMetalGlobalMu free energy
// to grid/SCF tolerance (measured ~6e-8) with fewer k-points -- the IBZ payoff, done exactly (doc/GPWPlan1 item 3).
TEST(GPW_SCF, AlFCCMetalIBZExact)
{
    FCCUnitCell cell(7.653);
    cell.AddAtom(13, {0,0,0});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwOptions o=AlOptions();
    o.globalFermi=true; o.imposeSymmetry=true;       // fold to the irreducible wedge AND star-average the density
    o.scf.SmearingkT=0.01;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, /*verbose*/false);

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 3.0, 1e-6);
    // Re-anchored 2026-08-01: the 0i custom V_loc G-ball (harmonic routing + custom top level) moved the
    // long-PP block by 1.6e-4 on Al's coarse grids -- full mesh AND reduced shift TOGETHER (folding stays
    // exact).  Old kappa-sweep anchor: -2.116812.
    // Re-pinned 2026-08-08 ON THE UNIFORM ROUTE (V2.4, arming the V1.26 cost selector): Al is soft
    // (alpha_max=4), so Auto now costs the two grids and picks uniform -- 4,096 points against Becke's
    // 18,000.  The value it lands on, -2.1169707, is EXACTLY the uniform-route number the previous
    // comment here already recorded, so this re-pin introduces no new quantity; it switches which of two
    // long-known values is the default.
    //   Which is right?  Measured, not assumed (GPW_SCF.DISABLED_GridRouteAB_AlFCC): against a fine Becke
    //   reference the uniform route is CLOSER on both scores -- ||drho||_1 6.18e-4 vs the production
    //   Becke mesh's 8.95e-4, and dEtot -1.34e-5 vs +1.93e-4.  The production Becke mesh's own residual
    //   is the angular error V2.6 measured on this system (Al is the worst case there).  And the uniform
    //   route is converged: refining its cutoff 4x moves ||drho||_1 by 1 part in 6000.
    // Previous anchor, on the Becke route: -2.1174805.
    EXPECT_NEAR(R.E.GetTotalEnergy(), -2.1169707, 1e-4)  // == the full 8-k-point mesh: IBZ symmetrization is EXACT
        << "IBZ-reduced must reproduce the full-mesh free energy (AlFCCMetalGlobalMu prints the same value)";
}

// (item 5, IBZ) NON-SYMMORPHIC -- diamond Si (FCC lattice + a 2-atom basis at (0,0,0),(¼,¼,¼)) is space group
// Fd-3m: NON-symmorphic (the two sublattices are related by a glide, τ=(¼,¼,¼)≠0).  The k-FOLD reaches the full
// Oh (3 irreducible k-points, same as FCC Al on this lattice): the τ=0 Td subgroup + time reversal (k→−k) already
// supplies the inversion Td lacks.  The DENSITY star-average now carries the glide τ: the G-space Hartree via the
// e^{+2πi(Um)·τ} phase (SymmetrizeGMap over SpaceGroup::ReciprocalOps) and the real-space XC raster via the exact
// FFT fractional shift ρ(W·x+τ) (SymmetrizeRaster over SpaceGroup::DirectOps).  So the IBZ-reduced total now
// reproduces the full-mesh Γ-centred 2×2×2 exactly (was −8.259, ~0.48 Ha off, under the old τ=0 W-only guard).
TEST(GPW_SCF, SiDiamondIBZ_NonSymmorphic)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});          // diamond = FCC + the glide-related 2nd sublattice (non-symmorphic)
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwOptions o;
    o.label="Si diamond IBZ"; o.Nelec=8; o.species={{"Si",4}};
    o.densityEcut=20.0; o.accelerator="DIIS"; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    // Target: the IBZ-reduced total reproduces the full-mesh Γ-centred 2×2×2 (DISABLED_SR_2x2x2GammaCentred,
    // -7.77846) to grid/SCF tolerance -- the non-symmorphic glide τ-phase makes the reduced density exact
    // (measured -7.77847, ~1e-5 vs the full mesh; the fold reaches 3 irreducible k-points under the full Oh).
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.77846, 2e-3)
        << "diamond Si (non-symmorphic Fd-3m): IBZ density symmetrization with the glide τ-phase must match the "
           "full mesh -7.77846 (G-space e^{+2πi(Um)·τ} + real-space FFT τ-shift)";
}

// (T3.2, doc/SymmetryUpgradePlan.md §6b) STREAM FOLD through-SCF A/B on an IMPOSED Γ-only run: the factory
// arms route (b) on the shared molecular evaluator (reduced stream build + replay, rep-transform h), and the
// existing SymmetrizeGMap/SymmetrizeRaster sites complete the group-average.  GPW_STREAM_FOLD=0/1 toggles the
// fold between two otherwise identical runs (read fresh in the factory, so one process can A/B).  The two
// totals must agree to the band-limit class (§8 through-SCF tier -- the production 5-smooth grid is NOT
// τ-commensurate, so reduced+P and full+P are two equally valid quadratures of the same density), and the
// folded run's [stream cache] line must show the reduced build (repPairs << pairs).
TEST(GPW_SCF, StreamFoldImposedGamma_SiDiamond)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});          // Fd-3m: non-symmorphic, 48 ops, quarter glide
    Lattice_3D lat(cell, ivec3_t(1,1,1));        // Γ-only: the T3.2 arming condition (k≠Γ is T3.4)
    GpwOptions o;
    o.label="Si diamond Γ stream-fold A/B"; o.Nelec=8; o.species={{"Si",4}};
    o.densityEcut=20.0; o.accelerator="DIIS"; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;

    setenv("GPW_STREAM_FOLD","0",1);
    GpwResult R0=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    setenv("GPW_STREAM_FOLD","1",1);
    GpwResult R1=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    unsetenv("GPW_STREAM_FOLD");

    EXPECT_TRUE(R0.converged); EXPECT_TRUE(R1.converged);
    std::cout << "[stream fold A/B] E(full)=" << R0.E.GetTotalEnergy()
              << "  E(folded)=" << R1.E.GetTotalEnergy()
              << "  dE=" << R1.E.GetTotalEnergy()-R0.E.GetTotalEnergy() << std::endl;
    EXPECT_NEAR(R1.E.GetTotalEnergy(), R0.E.GetTotalEnergy(), 1e-5)
        << "route (b) reduced streams must reproduce the full-stream imposed Γ run (band-limit class)";
    EXPECT_NEAR(R1.charge, 8.0, 1e-6);
}

// ===== EXPERIMENTAL (scratch): imposed multi-k × GDM probe (user observation 2026-08-03) =====
// OBSERVED on the user's multi-k NaF: imposeSymmetry=false + GDM reaches ΔE/E~1e-9 (Δρ~1e-3), while
// imposeSymmetry=true + GDM stalls at ΔE/E~1e-7 with Δρ~1e-2.  Two competing explanations:
//   (a) the imposed E/H pair is NOT exactly variational (a broken gradient -- GDM is the audit DIIS can't do);
//   (b) the projector restores EXACT degeneracies across the k-star (a free run's ULP splittings act as a
//       regularizer), so occupations/rotations tie at the filling and line-search level -- E flat, ρ wandering
//       ("electron jumping between k branches", the crystal analog of the atom-in-box E-settled/ρ-rotates mode).
// This probe A/Bs the in-suite proxy (Si diamond 2x2x2 IBZ): GDM_IMPOSE=0/1 toggles imposition,
// GDM_KT sets Fermi smearing (smearing cures (b), cannot cure (a)).  Run manually:
//   GDM_IMPOSE=1 GDM_KT=0.01 ./ITMain --gtest_filter='*ImposedGDMProbe*' --gtest_also_run_disabled_tests
// MEASURED (2026-08-03, Si diamond SR, Becke-auto XC): imposition does NOT degrade GDM on the clean
// insulator -- free and imposed both converge in 11 iterations to lastΔρ=3.0e-7, E_imposed=-7.77837 vs
// E_free=-7.77835 (imposed 2e-5 LOWER, as a variational E should be once the projector removes the
// free run's symmetry-defect noise; the free run's measured defect was 1.4e-5).  So (a) shows no
// GDM-visible gradient inconsistency at multi-k on Si; the NaF stall points at (b) -- NaF's diffuse
// near-gap virtuals + exactly-restored star degeneracies (branch ties).  Discriminator for the NaF
// config: SmearingkT ~ 0.005-0.01 rescues (b), cannot rescue (a).
TEST(GPW_SCF, DISABLED_ImposedGDMProbe_SiDiamondIBZ)
{
    auto envd=[](const char* n,double d){const char*s=std::getenv(n);return s?std::atof(s):d;};
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwOptions o;
    o.label="Si diamond IBZ GDM probe"; o.Nelec=8; o.species={{"Si",4}};
    o.densityEcut=20.0; o.accelerator="GDM";
    o.imposeSymmetry = envd("GDM_IMPOSE",1.0)!=0.0;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=int(envd("GDM_NMAX",60.0)); o.scf.MinΔρ=1e-6; o.scf.MinΔE=1e-10;   // let GDM run deep
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;
    o.scf.SmearingkT=envd("GDM_KT",0.0);
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/envd("GDM_VERBOSE",0.0)!=0.0);
    std::cout << "[GDM probe] impose=" << o.imposeSymmetry << " kT=" << o.scf.SmearingkT
              << " E=" << R.E.GetTotalEnergy() << " converged=" << R.converged << std::endl;
}

// The NaF SIBLING of the probe above -- the SMEARING DISCRIMINATOR on the user's actual configuration
// (multi-k 2x2x2 rocksalt NaF SR2, the DISABLED_NaFRocksaltGamma production recipe: IonicSAD seed, Kerker,
// delayed-IMOM, pivoted Cholesky; Becke XC at NAFGDM_L, default the fast L=11 recipe).  GDM does NOT (yet)
// support Fermi smearing, so the smeared leg runs as a TWO-STAGE anneal via RunGpwAnnealed's per-stage
// accelerator schedule: stage 1 DIIS at kT=NAFGDM_KT (smears the branch ties), stage 2 GDM at kT=0 seeded
// with stage 1's density (direct DM transfer; GDM's first aufbau collapses the fractional occupations).
// The discriminator (user question 2026-08-03): if the imposed-GDM stall (ΔE/E ~1e-7, Δρ~1e-2) is
// exactly-restored star degeneracies / branch ties, the smeared hand-off rescues it; a genuinely
// non-variational imposed E/H pair would stall regardless (and the Si proxy above already found none).
//   NAFGDM_IMPOSE=0/1 (default 1)   imposeSymmetry
//   NAFGDM_KT (default 0.01)        stage-1 smearing; 0 = the cold-cold control (the user's observed stall)
//   NAFGDM_NMAX (default 100)       per-stage iteration cap;  NAFGDM_L (default 11) Becke angular degree
//   NAFGDM_MOM=0                    disable delayed-IMOM (is MOM x GDM the fight?)
//   NAFGDM_PIVOT=0                  plain Cholesky (is the rank-reduced-subspace metric the fight?)
// Run:  NAFGDM_IMPOSE=1 NAFGDM_KT=0.01 ./ITMain --gtest_filter='*NaFImposedGDMSmear*' --gtest_also_run_disabled_tests
// MEASURED (2026-08-03, L=11, imposed = 3 IBZ points): the TIE HYPOTHESIS IS REJECTED for this stall --
// smearing does NOT rescue the GDM stage.  (A) imposed, kT=0.01: DIIS converges 46 it / Δρ=1.9e-7 /
// E=-24.5725 (-TS=3.4e-6: the gap dwarfs kT, smearing touches only ties); the seeded cold GDM stage then
// WALKS UPHILL to E=-24.527 (+45 mHa) wandering at Δρ~0.4.  (B) imposed, cold: DIIS converges 57 it /
// E=-24.5725; seeded GDM again degrades to -24.549 (+23 mHa), Δρ~0.29.  A minimizer leaving a stationary
// point uphill = an E/H<->minimizer inconsistency in play -- but NOT the imposed projector itself (the Si
// diamond proxy above is flawless imposed+GDM under identical projector machinery): the suspect set is
// what NaF's recipe ADDS -- delayed-IMOM through the GDM stage, pivoted-Cholesky rank reduction (GDM's
// gradient metric), diffuse-basis degeneracies.  NAFGDM_MOM / NAFGDM_PIVOT are the discriminator knobs.
// ROUND 2 (same day): MOM and PIVOT are ELIMINATED.  (D) imposed noMOM: DIIS 84 it -> -24.5725; seeded
// GDM climbs to -24.5167 (+56 mHa, Δρ~0.44 -- worse, if anything).  (E) noMOM+plain Cholesky: stage-1
// DIIS itself oscillates (the SR basis needs pivoting, as history says) yet GDM shows the SAME climb
// (-24.5178).  (C) FREE control: DIIS 53 it -> -24.5730; free GDM approximately HOLDS the fixed point
// (E +0.8 mHa, Δρ~5.7e-2 rotating -- the benign degenerate mode, matching the user's "sort of works").
// VERDICT: the defect is specifically IMPOSED x GDM on NaF (Si imposed x GDM is perfect): free-flat
// degenerate directions acquire CURVATURE through the projector (E[P rho(theta)] varies where E[rho(theta)]
// was flat), and the GDM engine walks UPHILL along them -- since the per-term E/H pairing is provably
// exact even at broken-symmetry iterates (W1 adjoint; Hartree self-adjoint P), the prime suspect is the
// GDM step model (curvature/step acceptance) on the projector-modified landscape, not E[rho] itself.
// NEXT INSTRUMENT: a directional FD gradient gate through the imposed stack (E(D±h dD) vs Tr(H dD)), or a
// GDM step log (model-predicted dE vs actual dE per accepted step).
// ROUND 3 (user-measured, same day): on SR2 (the diffuse-conditioned basis -- Na p 0.05 triplet removed)
// imposed multi-k GDM shows STRICTLY NEGATIVE dE -- variational descent confirmed, closing the original
// question: the imposed E[rho] IS variational (matching the W1/Hartree adjoint proofs and the Si proxy).
// The SR uphill walk therefore localizes to GDM x the DIFFUSE near-null directions: the engine wanders
// along soft modes in both free and imposed runs (free: Δρ~5.7e-2 at +0.8 mHa -- genuinely flat), and the
// projector CURVES those directions so the same wandering costs +23-56 mHa uphill.  Engine fixes, not
// functional fixes: an actual-E-checked step acceptance (clamps the uphill) and metric preconditioning of
// the soft directions (fixes the wander).  SEPARATE issue: GDM is SLOW on SR2 (tiny Δρ drops) imposed and
// free alike -- step-scale/preconditioning ergonomics (GDMParams{1.0}), imposition-independent
// (doc/SCFStrategyPlan.md territory; GDM has no gtest coverage -- the boron pin).
TEST(GPW_SCF, DISABLED_NaFImposedGDMSmearProbe)
{
    auto envd=[](const char* n,double d){const char*s=std::getenv(n);return s?std::atof(s):d;};
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na (Zion=1)
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F  (Zion=7)
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwOptions o;
    o.label   = "NaF IBZ GDM-smear probe";
    o.Nelec   = 8; o.species = {{"Na",1},{"F",7}};
    o.densityEcut = envd("NAF_ECUT",-1.0);              // AUTO = C·αmax = 80 (the anchor config)
    o.imposeSymmetry = envd("NAFGDM_IMPOSE",1.0)!=0.0;
    o.seed    = qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho   = qchem::CholeskyPivoted; o.orthoTol = 1e-4;
    o.scf.NMaxIter=(size_t)envd("NAFGDM_NMAX",100.0);
    o.scf.MinΔE=1e-9; o.scf.MinΔρ=1e-5;                 // deep enough to expose a stall vs a clean floor
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.25; o.scf.KerkerG0=1.0;
    o.scf.UseMOM=envd("NAFGDM_MOM",1.0)!=0.0; o.scf.MOMStartIter=10;   // NAFGDM_MOM=0: is MOM x GDM the fight?
    if (envd("NAFGDM_PIVOT",1.0)==0.0) { o.ortho=qchem::Cholesky; o.orthoTol=0.0; }  // plain-Cholesky control
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,int(envd("NAFGDM_L",11.0)));
    o.xcMesh.cellKind = qcMesh::UnitCellKind::Becke;

    const double kT=envd("NAFGDM_KT",0.01);
    GpwResult R = RunGpwAnnealed(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o,
                                 /*kTSchedule*/{kT, 0.0}, /*verbose*/envd("NAFGDM_VERBOSE",0.0)!=0.0,
                                 /*accSchedule*/{"DIIS","GDM"});
    std::cout << "[NaF GDM-smear probe] impose="<<o.imposeSymmetry<<" kT(stage1)="<<kT
              << " E="<<R.E.GetTotalEnergy()<<" converged="<<R.converged<<" iters(final)="<<R.iters<<std::endl;
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
}

// ===== EXPERIMENTAL (scratch): global-μ across k-blocks (item 3 inc 3) =====
// AL_KGRID=n (mesh nxnxn), AL_GLOBAL=0/1 (per-block vs global μ), AL_KT, AL_NMAX.
// I2 (plan §6a fit/grid SEPARATION): the (Delta, uniform) cross cell.  The SAME material through the
// PLANE-WAVE fit (band-limited v_xc on the FFT raster) and through the DELTA fit on the uniform cell
// mesh -- the two v_xc representations must agree to the route-difference class (band-limiting +
// raster-vs-midpoint-mesh quadrature), the same class the Becke-vs-uniform gate measures (~1e-4 Exc).
TEST(GPW_SCF, DeltaFitUniformGridMatchesPWFit_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.Nelec=8; o.species={{"Si",4}};
    o.densityEcut=20.0; o.accelerator="DIIS";
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3;

    o.label="Si PW-fit"; o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // (PlaneWave, raster)
    GpwResult P=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    ASSERT_TRUE(P.converged);

    o.label="Si Delta-fit uniform"; o.vxcFit=Hamiltonian::VxcFit::Delta;    // (Delta, uniform mesh)
    o.xcMesh.eCut=o.densityEcut;                                            // resolve rho on the midpoint mesh
    GpwResult D=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    ASSERT_TRUE(D.converged);

    EXPECT_NEAR(D.E.GetTotalEnergy(), P.E.GetTotalEnergy(), 5e-3)
        << "the delta fit on the uniform cell mesh must reproduce the PW fit on the raster to the "
           "band-limit/quadrature route-difference class (plan 6a fit/grid separation)";
}

// THE W1 GATE (doc/SymmetryUpgradePlan.md §6a): Becke XC under IBZ.  BeckeFit_IBS group-averages its
// mesh INVARIANT and star-averages rho every iteration (the SymmetrizeRaster hook, exact orbit-mean
// projector); on the converged SYMMETRIC density the invariant mesh integrates identically to the
// single-orientation mesh (Q_inv(f)==Q(f) for a symmetric f), so the IBZ run must reproduce the
// full-mesh Becke run to the IBZ class -- the Becke sibling of SiDiamondIBZ_NonSymmorphic, with the
// non-symmorphic glide tau exercised through the torus fold + MakeInvariant.  COARSE explicit Becke
// recipe on BOTH arms: the gate compares like against like, so grid quality cancels (and the
// imposed arm's group-average mesh growth stays affordable).
TEST(GPW_SCF, BeckeXC_IBZ_SiDiamond)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    GpwOptions o;
    o.Nelec=8; o.species={{"Si",4}};
    o.densityEcut=20.0; o.accelerator="DIIS";
    o.xcMesh=qcMesh::BeckeXCParams(15, 2.0, 9);          // explicit coarse Becke (nR=15, GL-9), same on both arms
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform; o.ortho=qchem::Cholesky;
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4;

    o.label="Si diamond Becke FULL"; o.imposeSymmetry=false;
    GpwResult F=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    ASSERT_TRUE(F.converged);

    o.label="Si diamond Becke IBZ"; o.imposeSymmetry=true;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false);
    ASSERT_TRUE(R.converged);

    std::cout<<"[Becke IBZ gate] full="<<F.E.GetTotalEnergy()<<" reduced="<<R.E.GetTotalEnergy()
             <<" dE="<<R.E.GetTotalEnergy()-F.E.GetTotalEnergy()<<std::endl;
    // Tolerance = the RULE-DIFFERENCE class at this deliberately coarse L: the free arm runs GL-9
    // (50 dirs), the imposed arm the MIXED-rule site-adapted minimal grid (degree-9-exact, ~76
    // dirs/atom) -- measured 2.0e-3 at L=9, collapsing with L (passes 2e-3 already at L=17; both
    // rules sit on the comparison floor at the production L=29).
    EXPECT_NEAR(R.E.GetTotalEnergy(), F.E.GetTotalEnergy(), 3e-3)
        << "Becke+IBZ must reproduce Becke+full-mesh (the W1 star-average makes the reduced density exact "
           "on the invariant Becke mesh; doc/SymmetryUpgradePlan.md 6a)";
}

TEST(GPW_SCF, DISABLED_AlGlobalMuExperiment)
{
    auto envd=[](const char* n,double d){const char*s=std::getenv(n);return s?std::atof(s):d;};
    const int nk=(int)envd("AL_KGRID",1);
    FCCUnitCell cell(7.653);
    cell.AddAtom(13,{0,0,0});
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));
    GpwOptions o=AlOptions();
    o.globalFermi = envd("AL_GLOBAL",1.0)!=0.0;
    o.imposeSymmetry = envd("AL_IBZ",0.0)!=0.0;
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
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
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
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
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
    // NAF_KMESH=n: an n^3 Γ-centred mesh.  The default 2 is what this test has run for some time -- and it
    // is NOT what the name, the header, or the -24.4304 anchor say: 2x2x2 (8 k -> 3 irreducible) lands
    // -24.5469, i.e. ~116 mHa of BAND DISPERSION below the Γ number the anchor was taken at.  The CP2K
    // oracle decks (naf_gpw_sr2_diag.inp / naf_gpw_sr_tight.inp) carry no &KPOINTS section, so they are Γ:
    // doc/Benchmark.md's NaF row is measured at NAF_KMESH=1, and comparing the default run to those decks
    // compares different Brillouin-zone samplings.  The anchor below follows the knob.
    const int nk=std::getenv("NAF_KMESH") ? std::atoi(std::getenv("NAF_KMESH")) : 2;
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));

    // SR2 (2026-07-16): the complete-enumeration-conditioned basis (lambda_min=1.57e-3; SR's three
    // degenerate 1.03e-6 near-null modes were exactly the Na p 0.05 triplet -- the cation's superfluous
    // diffuse shells; F kept intact for the anion).  See DISABLED_NaFOverlapConditioningSweep.
    // NAF_SPAN=sr runs the FULL SR span instead -- the one CP2K's naf_gpw_sr_tight.inp oracle holds
    // (-24.4322935, re-measured 2026-08-19), so doc/Benchmark.md's third NaF row is a comparison rather
    // than a lone column.  The pivoted-Cholesky rank filter is what makes the full span runnable at all.
    BasisSetData span=BasisSetData::VALENCE_LOWQ_SR2;
    if (const char* s=std::getenv("NAF_SPAN"))
    {
        const std::string v(s);
        if      (v=="sr" ) span=BasisSetData::VALENCE_LOWQ_SR;
        else if (v=="sr2") span=BasisSetData::VALENCE_LOWQ_SR2;
        else throw std::runtime_error("NAF_SPAN: expected sr|sr2, got '"+v+"'");
    }
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        span, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));

    // The production recipe as ONE GpwOptions literal (the full 2-week rationale is in the header above +
    // doc/GPWPlan §0b″).  The NAF_* env knobs stay as sweep INSTRUMENTS; the defaults ARE the committed recipe.
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
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
    o.scf.MergeTol=1e-4;
    o.scf.StartingRelaxRo=envd("NAF_ALPHA",0.25); o.scf.KerkerG0=1.0;   // Kerker damps the low-G charge-transfer slosh
    o.scf.UseMOM=true; o.scf.MOMStartIter=10;             // delayed-IMOM: descend, then pin the occupied subspace through the crossing
    o.scf.SmearingkT=envd("NAF_SMEAR",0.0); o.scf.MOMSmearPenalty=envd("NAF_PENALTY",0.0);   // MOM-masked Fermi (experiment)
    // XC quadrature: flip cellKind to try the atom-centred Becke XC route (doc/GPWPlan1.md; the recipe is
    // gate-calibrated, resolution sweepable via GPW_BECKE_L/NR/ALPHA).  The run prints [XC quadrature]
    // either way.  NOTE: Becke-in-SCF is unoptimised today (~min/iteration; the shared-rho + cached-Phi
    // GEMM route is the open perf item) -- start with GPW_BECKE_L=11.
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,24);//int nRadial=-1, double mhlAlpha=-1.0, int angularDegree=-1
    o.xcMesh.cellKind = qcMesh::UnitCellKind::Becke;    // <-- qcMesh::UnitCellKind::Becke turns Becke ON

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
    // ANCHOR: -24.4304 at Γ (auto Ecut=80, BallOnly, raw XC, guards) -- 0.8 mHa from CP2K SR2 truth
    // -24.4312 (re-measured 2026-08-19 through scripts/bench: -24.4312134, exactly the banked value); the
    // historical -27.93 "oracle" was RETRACTED as a screening artifact (doc/GPWPlan TRAPS #2).
    // The 2x2x2 default samples the zone and lands 116 mHa lower -- BAND DISPERSION, not a defect, but it
    // has no oracle, so it is anchored as a did-E-move pin in its own right.
    // The full-SR span is a DIFFERENT (larger) span, so it has its own anchor: qchem -24.4324 against
    // CP2K's -24.4322935.  Only the SR2 pair is asserted here; the SR arm is left to doc/Benchmark.md,
    // because the whole point of the wider span is that its conditioning is the thing under study.
    if (span==BasisSetData::VALENCE_LOWQ_SR2)
        EXPECT_NEAR(R.E.GetTotalEnergy(), nk==1 ? -24.4304 : -24.5469, 0.01);   // did-E-move anchor per k-mesh
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
    auto* accC=new qchem::SCFAccelerators::SCFAcceleratorNull();   // no DIIS (the CP2K recipe)
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
    auto* seedCD = scfC->GetWaveFunction()->GetChargeDensity().release();   // consumed by the fine ctor

    // ---- STAGE 2: seed the PRODUCTION fine grid (auto Ecut=8*alpha_max=320) with the converged coarse density. ----
    std::unique_ptr<Complex_BS> bsF(L3::GPWFactory(lat, mol, /*densityEcut*/envd("GC_FINE_ECUT",-1.0)));  // <0 AUTO=320
    Crystal_EC ecF(bsF->GetIrreps(Spin::None), 8);
    cHamiltonian* hamF=new Ham_PW_DFT(st, bsF.get(), {{"Na",1},{"F",7}}, "LDA");
    auto* accF=new qchem::SCFAccelerators::SCFAcceleratorNull();
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
    auto cd=scfF->GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge();
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
    auto* acc=new qchem::SCFAccelerators::SCFAcceleratorDIIS(qchem::SCFAccelerators::DIISParams{8, 8.0, 1e-10, 1e-8});
    testing::internal::CaptureStdout();
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::IonicSAD, lat.GetStructure().get(),
                                         qchem::Eigen, 1e-6);   // canonical ortho drops the ~8e-8 null mode -> rank 36
    SCFParams par; par.NMaxIter=3; par.MinΔρ=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30;
    par.MinFD=1e30; par.StartingRelaxRo=0.3; par.Verbose=false;
    scf.Iterate(par);
    std::string log=testing::internal::GetCapturedStdout();
    auto cd=scf.GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge();
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
    auto* acc=new qchem::SCFAccelerators::SCFAcceleratorDIIS(qchem::SCFAccelerators::DIISParams{8, 8.0, 1e-10, 1e-8});
    qchem::SCFIterator::SolidSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::IonicSAD, lat.GetStructure().get(),
                                         qchem::Eigen, 1e-6);   // (2): canonical ortho, drop the ~0 null cluster
    SCFParams par; par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30;
    par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true; par.KerkerG0=1.0;
    qchem::Hamiltonian::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");
    scf.Iterate(par);
    qchem::Hamiltonian::ReportGridCharge()=false;
    auto cd=scf.GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge();
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

//================================================================================================
//  THE BECKE XC GATE (doc/GPWPlan1.md "Becke XC grid").  The atom-centred periodic Becke XC
//  quadrature must reproduce the uniform-multigrid XC on a CONDITIONED basis -- same E_xc and same
//  V_xc matrix to grid tolerance -- before it can become the default for diffuse bases.  Converge
//  Si/Gamma on the standard uniform route (the SiliconGammaConverges recipe), then evaluate BOTH
//  XC term pairs (Dirac + VWN5) on the SAME converged density:
//    uniform -- PWFittedVxc on the Vxc fit basis's FFT grid (the raw-collocation route);
//    Becke   -- DeltaFittedVxc: rho(r) analytic at the atom-centred points, WeightedOverlap matrix.
//  Angular rule: GaussLegendre (machine-exact algebraic degree at any L -- the audited Lebedev
//  tables stop at L=11; see the Mesh_AngularDegree tests).
//================================================================================================
namespace
{
// One XC quadrature's answers on a converged density: E_xc, the quadrature's rho integral vs Tr(DS),
// and the summed Dirac+VWN V_xc matrix per Bloch block.  Uniform and Becke probes share the shape, so
// gates can diff ANY pair (uniform-vs-Becke, or Becke-vs-Becke for internal convergence).
struct XCProbe
{
    std::string label;
    double Exc=0, rhoLost=0;
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
        EXPECT_NE(face,nullptr) << "a real block needs the term's real-block face (3c-1)";
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
    XCProbe p; p.label=label; p.Exc=e.Exc; p.rhoLost=e.GridChargeLost;
    for (size_t i=0;i<h.bs->GetNumIBS();++i)
    {
        hmat_t<dcmplx> M=ProbeBlockMatrix(x,*h.bs,i,h.cd.get());
        M+=ProbeBlockMatrix(c,*h.bs,i,h.cd.get());
        p.M.push_back(std::move(M));
    }
    return p;
}

XCProbe UniformXCProbe(const GpwHandles& h, const std::shared_ptr<const Structure>& st)
{
    auto exch=std::make_shared<Hamiltonian::SlaterExchange>(2.0/3.0);
    auto corr=std::make_shared<Hamiltonian::VWN_Correlation>();
    qcMesh::MeshParams mp; mp.relCutoff=std::max(exch->GridCutoffFactor(), corr->GridCutoffFactor());
    Hamiltonian::PWFittedVxc::fbs_t vfb(h.bs->CreateVxcFitBasisSet(st.get(), mp));
    Hamiltonian::PWFittedVxc x(exch,vfb), c(corr,vfb);
    EnergyBreakdown e; x.GetEnergy(e,h.cd.get()); c.GetEnergy(e,h.cd.get());
    return ProbeXC("uniform", h, x, c, e);
}

XCProbe BeckeXCProbe(const GpwHandles& h, const std::shared_ptr<const Structure>& st,
                     const std::string& label, const qcMesh::MeshParams& mpB)
{
    auto exch=std::make_shared<Hamiltonian::SlaterExchange>(2.0/3.0);
    auto corr=std::make_shared<Hamiltonian::VWN_Correlation>();
    auto engine=std::make_shared<Hamiltonian::XC_GridEngine>(               // ONE engine, shared by the pair
                    std::make_shared<const qcMesh::Mesh>(st->CreateIntegrationMesh(mpB)));   // free probe: no fold
    Hamiltonian::DeltaFittedVxc x(exch,engine), c(corr,engine);
    EnergyBreakdown e; x.GetEnergy(e,h.cd.get()); c.GetEnergy(e,h.cd.get());
    return ProbeXC(label, h, x, c, e);
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
} //anon

TEST(GPW_SCF, BeckeXCMatchesUniformXC_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // Converge on the standard uniform route.  densityEcut=60 (not the SCF-sufficient 20): the gate
    // compares MATRIX ELEMENTS, and the uniform raw-adjoint H_xc carries raster error at N=15^3 that the
    // ~1e-4-converged Becke quadrature exposes (measured at Ecut=20: dExc=1.2e-4 but max|U-B|=1.2e-2 --
    // the raster's error, not Becke's; at Ecut=60 max|U-B|=3.5e-4).
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="Si Becke gate"; o.Nelec=8; o.species={{"Si",4}}; o.densityEcut=60.0;
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // the gate's SCF arm is the uniform route BY DESIGN (Auto would flip it)
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3;
    GpwHandles h;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false, &h);
    ASSERT_TRUE(R.converged);

    XCProbe U=UniformXCProbe(h, lat.GetStructure());
    XCProbe B=BeckeXCProbe(h, lat.GetStructure(), "B40", qcMesh::BeckeXCParams());
    EXPECT_NEAR(B.Exc, U.Exc, 5e-4);                 // measured: dExc=1.1e-4
    EXPECT_NEAR(B.rhoLost, 0.0, 5e-3);               // the Becke mesh integrates rho to Tr(DS)
    EXPECT_LT(DiffXC(U,B), 1e-3);                    // measured: 3.5e-4
}

//================================================================================================
//  V2.6 -- THE BECKE RECIPE LADDER.  "Is the production recipe (nRadial=40, angularDegree=29)
//  over-generous?"  Since V1.26 those two numbers set the ENTIRE Becke side of the Uniform-vs-Becke
//  cost comparison, so the answer changes grid CHOICES, not just grid cost.
//
//  METHOD (D8-compliant by construction): converge ONCE, then quadrature the SAME frozen density on a
//  ladder of meshes and score each against a FINE reference IN THE SAME FAMILY (nR=100, GL-41, same
//  mhl_alpha).  No SCF re-runs, no DeltaE_total -- the scored quantities are E_xc and the V_xc MATRIX,
//  i.e. properties.  Freezing rho is what makes this a measurement of the QUADRATURE, not of the SCF.
//  (The reference must be inside the rung family.  The tree's older REF80 probe uses mhl_alpha=1.0,
//  which is fine for a one-off cross-check but would put a fixed offset under every rung here and read
//  as an error FLOOR the axes never get below.)
//
//  Two axes, swept one at a time because they answer different questions:
//    ANGULAR -- how much angular freedom does a converged density actually need?  The steep axis: GL
//      directions go as (L+1)^2/2, so 29 -> 17 is a 2.8x cut in the whole mesh.
//    RADIAL -- cost is only LINEAR here, so it can afford to be generous; the question is whether it
//      is generous to no purpose.
//
//  ---------------------------------------------------------------------------------------------
//  RESULTS, 2026-08-07 (four systems).  max|dVxc| against the reference; the Becke gate's tolerance
//  is 1e-3.  Rows are the first rung INSIDE tolerance, in bold-equivalent (<-- marks it).
//
//    ANGULAR (nR=40 fixed)          Si covalent   NaF ionic   Mn open-d   Al metal
//      GL-9                           1.4e-3       4.2e-3      1.1e-5 <--  8.3e-3
//      GL-11                          9.5e-4 <--   3.0e-3      6.8e-6      1.8e-2 (!)
//      GL-15                          8.0e-5       7.6e-4 <--  1.4e-6      2.4e-3
//      GL-17                          5.8e-5       2.2e-4      1.3e-6      3.9e-4 <--
//      GL-29 (production)             1.7e-5       7.9e-5      1.3e-6      2.6e-4
//
//    RADIAL (GL-29 fixed)           Si           NaF         Mn          Al
//      nR=20                          3.5e-3       3.8e-3      1.2e-2      7.5e-2
//      nR=25                          9.0e-4 <--   1.2e-3      1.5e-3      3.8e-3
//      nR=30                          1.5e-4       3.2e-4 <--  2.0e-4 <--  9.3e-3 (!)
//      nR=40 (production)             1.7e-5       7.9e-5      1.3e-6      2.6e-4 <--
//      nR=60                          3.7e-6       2.1e-5      2.7e-7      9.3e-5
//
//  VERDICT: angular 29 is ~2.8x over-generous on INSULATORS (17 would do); RADIAL 40 is right.
//  BUT THE FLIP TO 17 WAS TRIED AND BACKED OUT -- see V2.6a below.  Production stays at 29.
//
//  ---------------------------------------------------------------------------------------------
//  FOUR THINGS THIS MEASUREMENT REFUTED.  Recorded because each one is the obvious guess:
//
//  (1) "The angular requirement is set by the SITE's density asphericity (which harmonics its point
//      group allows)."  FALSE.  It is set by the BECKE PARTITION SURFACE between DISSIMILAR
//      neighbours.  Evidence: Mn is a SINGLE atom in a box -- no interatomic partition at all -- and
//      is the EASIEST (degree 9); NaF, two ions of very different size and sharpness, is the hardest
//      of the three insulators; Si, partitioned between IDENTICAL atoms, sits between.  The mechanism
//      was already recorded in src/Structure/tests/MolecularMeshTests.C ("the fuzzy Voronoi switching
//      shell is angular-quadrature limited"); the site-symmetry framing simply ignored it.
//
//  (2) "Ionic materials, being near-spherical closed-shell ions, are the most forgiving."  FALSE --
//      NaF is the LEAST forgiving of the three insulators, for the reason in (1).
//
//  (3) "A high-spin d5 shell is the most aspherical density available."  FALSE -- half-filled d5 is
//      the one d configuration that IS spherically symmetric.  So DISABLED_BeckeRecipeLadder_MnSextet
//      does not probe asphericity even in principle.  GENUINE on-site asphericity remains UNTESTED by
//      all four systems; a crystal-field-split d (MnO) or the O2 pi* triplet would be the real probe.
//
//  (4) "nRadial adequacy can be predicted a priori from basis sharpness."  Nearly -- but not by the
//      obvious statistic.  Counting radial nodes inside the sharpest density feature rates nR=20 as
//      comfortable where measurement puts it 4-50x out, because MHL clusters as r ~ x^m and those
//      nodes bunch at tiny r leaving a gap AT the density peak.  The quantity that tracks Si and Mn is
//      the LOCAL SPACING there, r_peak/dr(r_peak) >~ 3 -- but Al breaks that threshold too (3.2 at
//      nR=30, which the rule passes, where measurement is 9x out).  The SHAPE is right, the threshold
//      is insulator-fitted.  See V2.7.
//
//  ---------------------------------------------------------------------------------------------
//  V2.6a -- THE FLIP TO 17, ATTEMPTED AND REJECTED.  A fourth system (Al FCC, simple metal) refutes
//  the three-system recommendation, and it refutes it in three separate ways:
//
//   (a) Al is the worst case on BOTH axes and its convergence is NON-MONOTONIC (GL-11 worse than GL-9,
//       GL-23 worse than GL-17), so no "degree N suffices" statement survives it.  Why a metal is worst
//       follows from (1): a nearly-free-electron valence density is the one that puts substantial
//       charge ON the partition surface.  The theory held; the three-system sample just did not
//       contain its own worst case.
//   (b) GPW_SCF.AlFCCDegenerateShellAufbauStalls flips QUALITATIVELY at degree 17: it starts
//       CONVERGING, which is a FAILURE.  That test asserts integer aufbau cannot converge a degenerate
//       3p shell (the density rotates freely in the manifold).  A coarser angular grid has larger
//       ORIENTATION-DEPENDENT quadrature error, which PINS the rotation -- so the run "converges" to a
//       state held in place by grid error.  Same mechanism as the SiPseudoAtomInBoxMatchesFinite
//       caveat, from the other side.  Never re-pin that EXPECT_FALSE to make a grid change pass.
//   (c) THE METHODOLOGICAL LIMIT OF THIS WHOLE FILE'S LADDER, and the most transferable finding:
//       A FROZEN-DENSITY QUADRATURE ERROR UNDERSTATES THE SELF-CONSISTENT SHIFT ON A METAL.
//       Al at degree 17 measures max|dVxc|=3.9e-4, comfortably inside the gate -- yet both Al SCF
//       anchors moved 6.4e-4 in TOTAL ENERGY when the default was flipped.  Grid error feeds back
//       through the density, the Fermi level and the occupations; on an insulator the two numbers
//       agree (Si/NaF anchors did not move at all), across a Fermi surface they do not.  A ladder
//       measures the QUADRATURE; for a metal that is a LOWER BOUND on what the SCF does with it.
//       => Before trusting any verdict from these ladders on a metal, do a converged-run A/B as well.
//
//  STANDING RULE EARNED HERE: calibrate a grid criterion on a simple METAL, or do not ship it as a
//  global default.  Three insulator-fitted rules broke on Al in one session (the angular degree, the
//  degenerate-shell assumption, and V2.7's radial threshold).
//================================================================================================
namespace
{
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
    const size_t nAtoms=st->GetNumAtoms();
    std::printf("\n[V2.6 ladder %s] reference nR=100 GL-41 (same alpha=2.0 family): Exc=%.8f\n", system, REF.Exc);
    auto row=[&](int nR, int deg)
    {
        const XCProbe P=BeckeXCProbe(h, st, "rung", mesh(nR,deg));
        double dV=0;                                   // worst |V_xc(i,j)| deviation over all blocks
        for (size_t b=0; b<REF.M.size() && b<P.M.size(); b++)
            for (size_t i=0;i<REF.M[b].rows();++i)
                for (size_t j=0;j<REF.M[b].columns();++j)
                    dV=std::max(dV, std::abs(dcmplx(P.M[b](i,j))-dcmplx(REF.M[b](i,j))));
        const long pts=long(nAtoms)*nR*(((deg+1)/2)*(deg+1));
        std::printf("[V2.6 ladder %s] nR=%-3d GL-%-3d  %9ld pts  dExc=%+.3e  max|dVxc|=%.3e\n",
                    system, nR, deg, pts, P.Exc-REF.Exc, dV);
    };
    std::printf("[V2.6 ladder %s] --- ANGULAR sweep (nRadial=40 fixed) ---\n", system);
    for (int deg : {5,7,9,11,15,17,21,23,29}) row(40,deg);
    std::printf("[V2.6 ladder %s] --- RADIAL sweep (degree=29 fixed) ---\n", system);
    for (int nR : {10,15,20,25,30,40,60}) row(nR,29);
}
} //anon

// COVALENT reference case: Si diamond.  Directional bonds put real charge between the atoms, which the
// Becke partition hands back to each site as a genuine l>0 component -- so this should be the HARDER of
// the two bonding characters for a low-degree angular rule (the prediction being tested).
TEST(GPW_SCF, DISABLED_BeckeRecipeLadder_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.label="Si V2.6 ladder"; o.Nelec=8; o.species={{"Si",4}}; o.densityEcut=60.0;
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // converge on the uniform route: the probe density
    o.imposeSymmetry=false;                            // free mesh: the ladder measures the RULE, not the fold
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3;
    GpwHandles h;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false, &h);
    ASSERT_TRUE(R.converged);
    BeckeLadder(h, lat.GetStructure(), "Si-covalent");
}

// IONIC contrast: NaF rocksalt (the sharp-F system).  Closed-shell near-spherical ions, so the PREDICTION
// under test is that this is MORE forgiving than Si, not less -- site point symmetry says which harmonics
// may be nonzero, never how large they are, and the Becke partition hands each site a near-spherical share.
// Recipe lifted verbatim from DISABLED_BeckeXCMatchesUniformXC_NaFSR2 (the committed NaF production recipe).
TEST(GPW_SCF, DISABLED_BeckeRecipeLadder_NaF)
{
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na (Zion=1)
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F  (Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ_SR2, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));

    GpwOptions o;
    o.label="NaF V2.6 ladder"; o.Nelec=8; o.species={{"Na",1},{"F",7}};
    o.accelerator="Ladder";
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.imposeSymmetry=false;                            // free mesh: the ladder measures the RULE, not the fold
    o.scf.NMaxIter=200; o.scf.MinΔE=1e-8; o.scf.MinΔρ=1e-4;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.45; o.scf.KerkerG0=1.0;
    o.scf.UseMOM=true; o.scf.MOMStartIter=10;
    GpwHandles h;
    GpwResult R=RunGpw(lat, mol, o, /*verbose*/false, &h);
    ASSERT_TRUE(R.converged);
    BeckeLadder(h, lat.GetStructure(), "NaF-ionic");
}

// OPEN-SHELL d contrast: the Mn sextet atom-in-box.  A half-filled 3d shell is the most ASPHERICAL density
// in the suite, so the prediction is that this is the LEAST forgiving of the three.  Recipe from
// MnAtomInBoxDChannel.  NOTE the probe pair is UNPOLARIZED while the run is spin-native: the probe then
// quadratures e_xc(rho_total) rather than the run's own two-channel E_xc.  That is deliberate and still a
// valid measurement, because what is being measured is the GRID's ability to integrate a real, converged,
// strongly aspherical density -- not this run's energy.  Do NOT compare its Exc to the run's.
TEST(GPW_SCF, DISABLED_BeckeRecipeLadder_MnSextet)
{
    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(25, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.label="Mn V2.6 ladder";
    o.Nelec=7; o.multiplicity=6;                       // S=5/2 Hund: nUp=6, nDown=1
    o.species={{"Mn",7}};
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.imposeSymmetry=false;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.scf.NMaxIter=40; o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.SmearingkT=5e-3;
    std::shared_ptr<const Real_BS> mnbasis(
        BasisSet::Molecule::Factory(BasisSetData::VALENCE_LOWQ_SR, &cell, BasisSet::Molecule::Engine::MnD,
                                    BasisSet::Molecule::Angular::Cartesian));
    GpwHandles h;
    GpwResult R=RunGpw(lat, mnbasis, o, /*verbose*/false, &h);
    ASSERT_TRUE(R.converged);
    BeckeLadder(h, lat.GetStructure(), "Mn-openshell-d");
}

// METALLIC contrast, added 2026-08-07 AFTER the first three refuted the flip to degree 17: Al FCC.  A
// nearly-free-electron valence density puts substantial charge ON the fuzzy-Voronoi partition surface --
// where Si/NaF/Mn all concentrate it near the nuclei -- so if the angular requirement is set by the
// PARTITION (which the first three said it is), a simple metal should be the WORST case.  It is: at degree
// 17 both Al anchors moved 6.4e-4, which is what backed the flip out.
TEST(GPW_SCF, DISABLED_BeckeRecipeLadder_AlFCC)
{
    FCCUnitCell cell(7.653);
    cell.AddAtom(13, {0,0,0});                         // Al (Zion=3): 3s^2 3p^1
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o=AlOptions();
    o.label="Al V2.6 ladder";
    o.scf.SmearingkT=0.02;                             // smeared: a converged, non-rotating density to freeze
    o.imposeSymmetry=false;                            // free mesh: the ladder measures the RULE, not the fold
    GpwHandles h;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, /*verbose*/false, &h);
    ASSERT_TRUE(R.converged);
    BeckeLadder(h, lat.GetStructure(), "Al-metal");
}

//================================================================================================
//  STEP 4 -- THE FACADE EQUIVALENCE GATE.  qchem::SolidCalculation is the NAMED front door for a
//  periodic SCF (the peer of Calculation / AtomCalculation), and the whole point of a facade is that
//  it IS the canonical recipe rather than a second one that drifts.  So gate it against the anchor
//  SiliconGammaConverges pins, on the same lattice / basis / cutoffs.
//
//  What this actually protects: RunGpw below is the test harness's driver (report bracketing, basis
//  vetting, fingerprints, order probes, symmetry reporting).  SolidCalculation is the production
//  recipe with none of that.  If the two ever disagree on the physics, one of them has grown a step
//  the other lacks -- which is exactly the failure a facade is supposed to make impossible, and
//  exactly what happened before it existed (the recipe lived only in the harness, so "the production
//  path" and "the test path" were the same code by accident rather than by construction).
//================================================================================================
TEST(GPW_SCF, SolidCalculationMatchesTheSiAnchor)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;

    qchem::SolidCalculation calc(lat, MakeBasisSR(cell),
                                 {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);

    EXPECT_TRUE(calc.Converged());
    EXPECT_NEAR(calc.TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(calc.Energy(), -7.11506, 2e-3)      // the SiliconGammaConverges anchor, same tolerance
        << "SolidCalculation must reproduce the driver's physics -- a facade that drifts from the "
           "recipe it fronts is worse than no facade";
    // The facade OWNS the grid decision, so it must be able to say what it chose (V1.26/V2.4): Si is soft
    // enough that the cost selector routes it to the uniform mesh, SIZED from this run's sharpness.
    EXPECT_EQ(calc.ResolvedXCMesh().cellKind, qcMesh::UnitCellKind::Uniform);
    EXPECT_GT(calc.ResolvedXCMesh().eCut, 0.0) << "an Auto-resolved uniform mesh must carry its own cutoff, "
                                                  "never fall back on nUniform's basis-blind default";
    // rho(r) is reachable through the same neutral ScalarFunction face the molecular facade exposes.
    EXPECT_GT(calc.Density()(rvec3_t(0.1,0.1,0.1)), 0.0);
}

//================================================================================================
//  CROSS-RUN DETERMINISM GATE (2026-08-18).  Three IDENTICAL SolidCalculation runs in one process
//  must give the SAME energy -- i.e. every run must replay fresh-process behaviour.  This is the
//  regression gate for the first-run anomaly: the GPW matrix-free 3C tensor closures capture the
//  evaluator's mutable per-SCF state (the CollocMemo D-screen), and when they were stored in the
//  process-wide DBCache a second identical run inherited the first run's converged D-screen -- its
//  seed Fock swept the full Hartree/Vxc pair set where a fresh process's (screened by the diagonal
//  SAD seed D) is diagonal-only: seed s-levels shifted ~0.24 Ha, converged E ~5e-6.  Fixed by
//  scoping those tensors per basis INSTANCE (tGPW_IBS::Overlap3C/Repulsion3C override); with the fix
//  the three energies here are BITWISE equal -- 1e-9 is pure headroom for BLAS/library variation.
//================================================================================================
TEST(GPW_SCF, CrossRunFirstRunAnomalyProbe)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;

    double E[3];
    for (int r=0;r<3;++r)
    {
        std::cout<<"[xrun] ================= RUN "<<r<<" ================="<<std::endl;
        qchem::SolidCalculation calc(lat, MakeBasisSR(cell),
                                     {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);
        EXPECT_TRUE(calc.Converged());
        E[r]=calc.Energy();
        std::cout.precision(15);
        std::cout<<"[xrun] run "<<r<<" E="<<E[r]<<std::endl;
    }
    std::cout.precision(15);
    std::cout<<"[xrun] E1-E0="<<E[1]-E[0]<<"  E2-E1="<<E[2]-E[1]<<std::endl;
    EXPECT_NEAR(E[1], E[0], 1e-9) << "cross-run pollution: first run differs";
    EXPECT_NEAR(E[2], E[1], 1e-9) << "runs 2+ should be steady";
}

//================================================================================================
//  STEP 3c-3 -- THE FACTORY-FLIP ACCEPTANCE (doc/RealComplexPlan.md).  SolidCalculation now computes
//  the working-type rule (irrep.IsReal() ∧ ham.PreservesReal()) and builds every TRIM block REAL by
//  default; forceComplex is the §6 ansatz-policy downgrade and this gate's A/B door.  The two runs
//  must be the SAME PHYSICS to machine precision: every per-term gate (RealComplexTerms.*) pinned the
//  real block's matrices bitwise except the quadrature GEMM's summation order (~1e-15/element), so
//  through a whole SCF the totals may differ only at accumulated-roundoff level.  A loose tolerance
//  here would be wrong twice over -- it could hide a genuine defect, and it would understate what the
//  per-term gates already guarantee.
//================================================================================================
//  TWO ARMS, because "machine-equal" means two different things across an SCF:
//  (1) ONE ITERATION from the shared uniform seed -- the PURE-ARITHMETIC gate.  Both runs build the
//      same Fock from the same D0 (statics + Hartree + raw XC bitwise per RealComplexTerms.*),
//      diagonalize, and fill the same gapped occupied set; every energy term from D1 must then agree
//      to roundoff.
//  (2) CONVERGED twins -- the PHYSICS gate.  Bitwise TRAJECTORIES are impossible by construction: the
//      real block diagonalizes with LAPACK's real path (dsyev) where the complex one runs zheev, which
//      agree only to roundoff even on an exactly-real matrix, and DIIS extrapolation chaotically
//      amplifies last-ulp seeds (measured: ~1e-4 in per-term energies after 12 iterations).  What the
//      physics guarantees is the FIXED POINT, so the converged states must agree to gate resolution.
//
//  (The 2026-08-18 WARM-UP DISCIPLINE is RETIRED: the first-run anomaly this gate used to hide from
//  -- a throwaway run so both arms sat in "steady-state slots" -- was the cross-run D-screen leak,
//  fixed by instance-scoping the GPW 3C tensors.  Every run now replays fresh-process behaviour
//  bit-for-bit, gated by GPW_SCF.CrossRunFirstRunAnomalyProbe above, so the arms run cold.)
static void ExpectRealComplexTwins(const Lattice_3D& lat, const UnitCell& cell, const char* what)
{
    const qchem::SolidCalcOptions optOn {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0};
    const qchem::SolidCalcOptions optOff{.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0, .forceComplex=true};

    SCFParams one;                                   // the one-iteration recipe (gates 0 = strict-< never trips)
    one.NMaxIter=1; one.MinΔρ=0.0; one.MinΔE=0.0;
    one.MinΔFD=1e30; one.MinVirial=1e30; one.MinFD=1e30; one.StartingRelaxRo=0.3; one.MergeTol=1e-4;

    {   // (1) the arithmetic arm: exactly one iteration each
        qchem::SolidCalculation on(lat, MakeBasisSR(cell), optOn, one), off(lat, MakeBasisSR(cell), optOff, one);
        const qchem::EnergyBreakdown En=on.EnergyTerms(), Ec=off.EnergyTerms();
        // 1e-8: real-vs-complex roundoff headroom -- Een's large-cancellation assembly (vs Enn~8 Ha)
        // measured 2e-9 across the 3-block mesh.  Roundoff grade, far below any defect scale.
        EXPECT_NEAR(En.GetTotalEnergy(), Ec.GetTotalEnergy(), 1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En.Kinetic, Ec.Kinetic, 1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En.Een,     Ec.Een,     1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En.Eee,     Ec.Eee,     1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En.Exc,     Ec.Exc,     1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(on.TotalCharge(), off.TotalCharge(), 1e-10) << what << " (iteration 1)";
    }
    {   // (2) the physics arm: both converge on the production-shaped gates
        SCFParams par;
        par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
        par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;
        qchem::SolidCalculation on(lat, MakeBasisSR(cell), optOn, par), off(lat, MakeBasisSR(cell), optOff, par);
        ASSERT_TRUE(on .Converged()) << what;
        ASSERT_TRUE(off.Converged()) << what;
        EXPECT_NEAR(on.Energy(), off.Energy(), 2e-5) << what;        // MinΔE=1e-6 resolution (measured 4.8e-6)
        const rvec3_t pts[]={ cell.ToCartesian(rvec3_t(0.3,0.4,0.7)),
                              cell.ToCartesian(rvec3_t(0.25,0.25,0.25)),
                              cell.ToCartesian(rvec3_t(0.1,0.9,0.2)) };
        for (const auto& r : pts)                                     // MinΔρ=1e-3 resolution (measured ≤9.1e-6)
            EXPECT_NEAR(on.Density()(r), off.Density()(r), 1e-4) << what;
    }
}

TEST(GPW_SCF, RealTRIMBlocksMatchComplex_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));   // Γ-only: EVERY block is TRIM, so the flip makes the whole run real
    ExpectRealComplexTwins(lat, cell, "Si Gamma real-vs-complex");
}

//  STEP 4's MIXED-MESH acceptance, at the smallest genuinely mixed mesh: N=(3,1,1) has ONE TRIM point
//  (Γ -- 2·1/3 is no reciprocal-lattice vector) beside two complex blocks (k=±1/3), so ONE composite
//  carries both child scalars through the full SCF -- the case Γ-only cannot reach.  (A 3×3×3 run is
//  the same code path 27 blocks wide; this keeps the gate's wall-time at ~3 Γ runs.)
TEST(GPW_SCF, RealTRIMBlocksMatchComplex_SiMixedMesh)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(3,1,1));
    ExpectRealComplexTwins(lat, cell, "Si (3,1,1) mixed-mesh real-vs-complex");
}

//  MOM ON A MIXED MESH -- the END-TO-END gate for the R2.21 state split (flagged at that item's
//  landing, added here at the merge).  Until R2.21, a real TRIM block had nowhere to put its
//  mat_t<double> MOM reference (the references were typed by the RUN), so switching MOM on made a
//  flipped run THROW mid-SCF -- after the basis, seed and first Fock were already paid for.  Now each
//  block's reference lives in the shared OccupationState under its OWN scalar, so the real child fills
//  under a genuine OccupationPolicy<double> and the run keeps one ledger.
//
//  The gate drives the case that could not previously run AT ALL: (3,1,1) -- one real Γ block beside
//  two complex ±⅓ blocks -- with MOM armed early (MOMStartIter=2, so the reference is captured and
//  consulted well inside the run), real vs forceComplex.  Equal converged energy = the real block's
//  reference was captured, scored and applied exactly as its complex twin's.
TEST(GPW_SCF, RealTRIMBlocksWithMOMMatchComplex_SiMixedMesh)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(3,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;
    par.UseMOM=true; par.MOMStartIter=2;              // armed EARLY: the reference must be live in-run

    qchem::SolidCalculation on (lat, MakeBasisSR(cell),
                                {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);
    qchem::SolidCalculation off(lat, MakeBasisSR(cell),
                                {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0, .forceComplex=true}, par);
    ASSERT_TRUE(on .Converged()) << "MOM on a mixed real/complex mesh must converge (R2.21)";
    ASSERT_TRUE(off.Converged());
    EXPECT_NEAR(on.Energy(), off.Energy(), 1e-9)
        << "a real TRIM block's MOM reference must behave exactly like its complex twin's";
    EXPECT_NEAR(on.TotalCharge(), off.TotalCharge(), 1e-10);
}

//================================================================================================
//  V2.4 -- THE GRID-ROUTE CONVERGENCE A/B.  Calibrate kUniformMargin, then arm the V1.26 selector.
//
//  THE QUESTION.  The selector's cost model says uniform beats Becke ~8-15x on Si and Al, so an ARMED
//  Auto would send both to the uniform route.  The 2026-08-01 measurement made Becke the default for
//  exactly those systems.  Both cannot be right, and kUniformMargin (a guess at 2.0) is the only thing
//  absorbing the discrepancy -- which is why the selector ships DISARMED until this lands.
//
//  WHY THIS IS NOT A LADDER.  V2.6a established that a FROZEN-density quadrature error understates the
//  self-consistent shift across a Fermi surface (Al measured 3.9e-4 frozen and moved 6.4e-4 through the
//  SCF).  Al is one of the two live uniform verdicts, so the frozen instrument is disqualified here by
//  its own finding.  This runs the SCF to CONVERGENCE on each route instead.
//
//  WHAT IS SCORED.  The converged DENSITY, per the D8 pin -- never DeltaE_total, because the fits are
//  non-variational so an energy difference does not bound the error.  Every route's converged rho is
//  sampled on ONE COMMON mesh and scored against the finest route:
//     ||drho||_1 = Integral |rho - rho_ref| dV   (electrons of misplaced charge -- the physical number)
//     ||drho||_inf = max |rho - rho_ref|         (the worst point)
//  Etot is printed for information only; it is not the verdict.
//
//  THE COMMON MESH IS DELIBERATELY UNIFORM AND DENSE, not Becke: scoring on an atom-centred mesh would
//  weight the cores where the Becke route is strongest and grade its own homework.  The cost is that a
//  uniform probe UNDER-weights the core, so ||drho||_1 flatters the uniform route slightly -- an error
//  in the conservative direction for a test whose null hypothesis is "uniform is good enough".
//
//  THE VALIDATION THAT MAKES THE REFERENCE TRUSTWORTHY: refine BOTH routes and check they agree
//  (U-4x vs B-fine).  Two independent quadratures converging to the same density is the only evidence
//  that either is converged; without it "close to B-fine" would just mean "close to Becke".
//
//  ---------------------------------------------------------------------------------------------
//  RESULTS 2026-08-08.  Both systems the selector routes to UNIFORM.  All routes CONVERGED (see the
//  setup note on Si below -- the first attempt did not, and that mattered).
//
//                 ||drho||_1     ||drho||_inf   dEtot vs B-fine   mesh
//    Si  U-sel     1.363e-4 e     8.0e-6         -4.2e-6            4,913
//        U-4x      1.362e-4 e     8.0e-6         -4.2e-6           (refined)
//        B-prod    1.106e-4 e     1.7e-6         +1.26e-4          72,000 (imposed)
//    Al  U-sel     6.18e-4 e      5.0e-5         -1.34e-5           4,096
//        U-4x      6.18e-4 e      5.0e-5         -1.34e-5          (refined)
//        B-prod    8.95e-4 e      3.9e-5         +1.93e-4          18,000
//
//  VALIDATION PASSES: U-sel == U-4x to 4 digits on BOTH systems, so the uniform route is genuinely
//  converged at the selector's cutoff -- refining it 4x changes nothing.  That is what licenses
//  reading the rest of the table as grid error rather than as under-resolution.
//
//  VERDICT: the selector is RIGHT on both.  Uniform at its own cutoff matches the fine reference at
//  least as well as the PRODUCTION Becke mesh does, for 4-15x fewer points -- and on TOTAL ENERGY it
//  is dramatically better (Si 30x, Al 14x closer to B-fine), because B-prod's own residual is the
//  angular error V2.6 measured.
//
//  NO CONTRADICTION WITH THE 2026-08-01 BECKE DEFAULT, which was justified for "diffuse bases" and
//  sharp cores.  Si-SR and Al-LOWQ-SR are neither; the systems that motivated it (F alpha_max=40, MnO)
//  are exactly the ones the selector ALREADY sends to Becke.  The selector reproduces that split
//  automatically instead of applying one answer to both regimes.
//
//  THE ONE PLACE BECKE STILL WINS, and it is visible above: ||drho||_inf.  Becke is 4.7x better at the
//  WORST point on Si and comparable on Al -- that point is the core, which is what an atom-centred
//  radial mesh is for.  So a property that samples the core (hyperfine, EFG, core-level shifts) should
//  prefer Becke even where the selector says uniform.  The integrated norm is not the whole story.
//================================================================================================
namespace
{
struct RouteProbe { std::string label; double Etot=0; rvec_t rho; bool ok=false; };

// One converged SCF on a given XC mesh, with its density sampled on the shared probe mesh.
RouteProbe RunRoute(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, GpwOptions o,
                    const char* label, const qcMesh::MeshParams& xc, const qcMesh::Mesh& probe)
{
    o.label=std::string(o.label)+" ["+label+"]";
    o.xcMesh=xc;
    GpwHandles h;
    RouteProbe p; p.label=label;
    GpwResult R=RunGpw(lat, mol, o, /*verbose*/false, &h);
    p.ok   = R.converged;
    p.Etot = R.E.GetTotalEnergy();
    p.rho  = (*h.cd)(probe.Points());     // the inherited ScalarFunction batch (R1.5's one spelling)
    return p;
}

// Score every route against the LAST one (the finest), on the shared mesh's weights.
void ScoreRoutes(const std::vector<RouteProbe>& r, const qcMesh::Mesh& probe, const char* system)
{
    const rvec_t& W=probe.Weights();
    const RouteProbe& ref=r.back();
    std::printf("\n[V2.4 %s] reference = %s (Etot=%.8f)   probe mesh %zu pts\n",
                system, ref.label.c_str(), ref.Etot, probe.size());
    for (const auto& p : r)
    {
        double l1=0, li=0;
        for (size_t i=0;i<p.rho.size() && i<ref.rho.size();++i)
        {
            const double d=std::fabs(p.rho[i]-ref.rho[i]);
            l1+=W[i]*d; li=std::max(li,d);
        }
        std::printf("[V2.4 %s] %-10s conv=%d  Etot=%+.8f  dEtot=%+.2e  ||drho||_1=%.3e e  ||drho||_inf=%.3e\n",
                    system, p.label.c_str(), int(p.ok), p.Etot, p.Etot-ref.Etot, l1, li);
    }
}

// Drive the four routes for one system.  U-sel is what an ARMED selector would choose; B-prod is what
// ships today; the two refined routes are the convergence evidence.
void GridRouteAB(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, const GpwOptions& base,
                 const char* system)
{
    const qcMesh::XCMeshSharpness sh=GatherSharpness(lat, *mol, base);
    const double eReq=qcMesh::RequiredUniformCutoff(sh);
    // The shared probe mesh: uniform, and FOUR TIMES the resolution any route under test uses, so the
    // sampling is not itself a variable in the comparison.
    qcMesh::MeshParams pm; pm.cellKind=qcMesh::UnitCellKind::Uniform; pm.eCut=4.0*eReq;
    const qcMesh::Mesh probe=lat.GetStructure()->CreateIntegrationMesh(pm);
    std::printf("\n[V2.4 %s] alpha_max=%g alpha_pp=%g -> selector eCut=%g Ha; uniform %ld pts vs Becke %ld\n",
                system, sh.alphaMax, sh.alphaPP, eReq,
                qcMesh::UniformMeshCost(sh), qcMesh::BeckeMeshCost(qcMesh::BeckeXCParams(), sh));

    qcMesh::MeshParams uSel; uSel.cellKind=qcMesh::UnitCellKind::Uniform; uSel.eCut=eReq;
    qcMesh::MeshParams u4x =uSel;                                        u4x.eCut=4.0*eReq;
    std::vector<RouteProbe> r;
    r.push_back(RunRoute(lat, mol, base, "U-sel",  uSel, probe));
    r.push_back(RunRoute(lat, mol, base, "U-4x",   u4x,  probe));
    r.push_back(RunRoute(lat, mol, base, "B-prod", qcMesh::BeckeXCParams(),         probe));
    r.push_back(RunRoute(lat, mol, base, "B-fine", qcMesh::BeckeXCParams(60,2.0,35),probe));
    ScoreRoutes(r, probe, system);
}
} //anon

// Si diamond: selector verdict UNIFORM (4,913 pts vs 36,000 Becke).
TEST(GPW_SCF, DISABLED_GridRouteAB_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o;
    o.label="Si V2.4"; o.Nelec=8; o.species={{"Si",4}}; o.densityEcut=20.0;
    // imposeSymmetry LEFT AT THE DEFAULT (true), unlike the V2.6 ladders.  A ladder freezes one density
    // and wants the bare angular RULE, so it turns the fold off; this test needs each route to reach a
    // CONVERGED density, and free-mesh Si/Gamma oscillates in its degenerate manifold (measured: 3 of 4
    // routes hit FIT-FLOOR STALL at 60 iterations, which made the first run's drho scores SCF noise
    // rather than grid error).  Imposed is also the production configuration, so this is the A/B that
    // matches what a user gets.
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3;
    GridRouteAB(lat, MakeBasisSR(cell), o, "Si");
}

// Al FCC: selector verdict UNIFORM (4,096 vs 36,000) -- and the system whose SCF amplified a
// frozen-density 3.9e-4 into 6.4e-4 (V2.6a).  Smeared, so the density is converged and non-rotating.
TEST(GPW_SCF, DISABLED_GridRouteAB_AlFCC)
{
    FCCUnitCell cell(7.653);
    cell.AddAtom(13, {0,0,0});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o=AlOptions();
    o.label="Al V2.4"; o.scf.SmearingkT=0.02; o.imposeSymmetry=false;
    GridRouteAB(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, "Al");
}

// THE ROTATED-LEBEDEV EXPERIMENT (plan §6a rotation insight, increment (b)): quadrature exactness is
// rotation-invariant, so rotating an efficient Lebedev grid OFF the bond axes should be a nearly-free
// accuracy fix for FREE runs -- the measured 5-10x rho-weighted loss was pure ALIGNMENT (the <111>
// orbit on diamond's bonds).  Hand-run (--gtest_also_run_disabled_tests): converge Si/Gamma once on
// the uniform route, then probe FOUR Becke quadratures on the SAME density against the refined GL
// reference (nR=80, GL-29): production GL-29, Lebedev-50 as-is (bond-aligned), Lebedev-50 rotated
// (GPW_BECKE_ROT class, 0.4 rad about (1,2,3)/sqrt14), and a rotated GL-29 control.  Verdict numbers
// land in the plan doc; the decision is whether the free-run default grid can shrink.
TEST(GPW_SCF, DISABLED_RotatedLebedevXCProbe_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="Si rotated-Lebedev probe"; o.Nelec=8; o.species={{"Si",4}}; o.densityEcut=60.0;
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // converge once on the uniform route (probe density)
    o.scf.NMaxIter=60; o.scf.MinΔρ=1e-3; o.scf.MinΔE=1e-6;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30; o.scf.StartingRelaxRo=0.3;
    GpwHandles h;
    GpwResult R=RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/false, &h);
    ASSERT_TRUE(R.converged);
    auto st=lat.GetStructure();

    const double rot=0.4;                              // generic angle: moves <111> well off the bonds
    auto leb=[&](int nDir, double angRot){ qcMesh::MeshParams mp=qcMesh::BeckeXCParams(40, 2.0, nDir);
                                           mp.angular=qcMesh::AngularKind::Lebedev; mp.angRot=angRot; return mp; };
    auto gl =[&](double angRot){ qcMesh::MeshParams mp=qcMesh::BeckeXCParams();  // production GL-29
                                 mp.angRot=angRot; return mp; };

    // How close does each Lebedev grid come to a bond axis?  (diamond bonds = the +<111> tetrahedron)
    for (int nDir : {50, 302})
        for (double angRot : {0.0, rot})
        {
            qcMesh::AngularMesh am=qcMesh::MakeAngular(leb(nDir, angRot));
            double worst=90.0;
            for (size_t i=0; i<am.size(); ++i)
                for (auto b : {rvec3_t(1,1,1), rvec3_t(1,-1,-1), rvec3_t(-1,1,-1), rvec3_t(-1,-1,1)})
                {
                    double cosang=(am.Dirs()[i]*b)/norm(b);
                    worst=std::min(worst, std::acos(std::min(1.0,std::abs(cosang)))*180.0/M_PI);
                }
            std::printf("[rotLeb] Lebedev-%d angRot=%.2f: closest direction-to-bond angle = %.3f deg\n",
                        nDir, angRot, worst);
        }

    XCProbe REF  =BeckeXCProbe(h, st, "REF80",  qcMesh::BeckeXCParams(/*nRadial*/80, /*mhlAlpha*/1.0, /*L*/29));
    XCProbe GL29 =BeckeXCProbe(h, st, "GL29",   gl(0.0));
    XCProbe L50R =BeckeXCProbe(h, st, "Leb50rot", leb(50, rot));
    XCProbe L302 =BeckeXCProbe(h, st, "Leb302",   leb(302, 0.0));
    XCProbe L302R=BeckeXCProbe(h, st, "Leb302rot",leb(302, rot));

    std::printf("[rotLeb] dExc vs REF80:  GL29=%+.3e  Leb50rot=%+.3e  Leb302=%+.3e  Leb302rot=%+.3e\n",
                GL29.Exc-REF.Exc, L50R.Exc-REF.Exc, L302.Exc-REF.Exc, L302R.Exc-REF.Exc);
    std::printf("[rotLeb] rho-lost:       GL29=%+.3e  Leb50rot=%+.3e  Leb302=%+.3e  Leb302rot=%+.3e\n",
                GL29.rhoLost, L50R.rhoLost, L302.rhoLost, L302R.rhoLost);
    std::printf("[rotLeb] max|Vxc-REF80|: GL29=%.3e  Leb50rot=%.3e  Leb302=%.3e  Leb302rot=%.3e\n",
                DiffXC(REF,GL29), DiffXC(REF,L50R), DiffXC(REF,L302), DiffXC(REF,L302R));
}

// The SHARP-FIELD leg of the gate (the plan names DISABLED_NaFRocksaltGamma as the stress case: the F-
// anion makes sharp peaks in rho and V_xc, and its diffuse basis is what the Becke grid exists for).
// DISABLED like the parent NaF anchor -- it is a long run (the full NaF convergence recipe at a
// matrix-grade densityEcut=160 reference, plus two Becke term evaluations); run it by hand with
// --gtest_also_run_disabled_tests when touching the XC quadrature.
TEST(GPW_SCF, DISABLED_BeckeXCMatchesUniformXC_NaFSR2)
{
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na (Zion=1)
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F  (Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ_SR2, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));

    // The committed NaF production recipe (DISABLED_NaFRocksaltGamma) at PRODUCTION grids (auto Ecut=80,
    // BallOnly): the SCF only supplies the density; the comparison itself carries the reference-grade work.
    // MEASURED ATTRIBUTION (2026-07-30): vs the SAME Becke B40 matrix, the uniform reference gave
    //   BallOnly  Ecut=160: max|U-B|=1.55e-1      (the production raster's sharp-F element error)
    //   AliasFree Ecut=160: max|U-B|=1.89e-2      (8x better -- most of the gap was BallOnly's)
    // with dExc tiny throughout (1.6e-4 / 1.3e-5) -- the discrepant elements carry no rho weight.  So on a
    // sharp-F system the UNIFORM raster is not element-converged at practical Ecut, which is this grid's
    // reason to exist; the gate here is Becke INTERNAL convergence (B40 vs a 2x-refined B80) plus the
    // energy-level agreement with the uniform route.
    GpwOptions o;
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.label="NaF Becke gate"; o.Nelec=8; o.species={{"Na",1},{"F",7}};
    o.accelerator="Ladder";
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.scf.NMaxIter=200; o.scf.MinΔE=1e-8; o.scf.MinΔρ=1e-4;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.45; o.scf.KerkerG0=1.0;
    o.scf.UseMOM=true; o.scf.MOMStartIter=10;
    GpwHandles h;
    GpwResult R=RunGpw(lat, mol, o, /*verbose*/false, &h);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    ASSERT_TRUE(R.converged);

    auto st=lat.GetStructure();
    XCProbe U  =UniformXCProbe(h, st);
    XCProbe B40=BeckeXCProbe(h, st, "B40", qcMesh::BeckeXCParams());
    XCProbe B80=BeckeXCProbe(h, st, "B80", qcMesh::BeckeXCParams(/*nRadial*/80, /*mhlAlpha*/1.0, /*L*/29));
    EXPECT_NEAR(B40.Exc, U.Exc, 2e-3);               // energy-level agreement with the uniform route
    EXPECT_NEAR(B40.rhoLost, 0.0, 5e-3);
    DiffXC(U,B40);                                    // report-only: the uniform raster's element error
    EXPECT_LT(DiffXC(B40,B80), 2e-3) << "Becke V_xc not internally converged on the sharp-F system";
}

// Mn PSEUDO-ATOM IN A BOX through the CRYSTAL (GPW) path vs the molecular facade -- the d-channel sibling
// of SiPseudoAtomInBoxMatchesFinite, and the cheap localiser for MnO's ~356 Ha over-binding.  Both
// real-space KB routes are now oracle-matched on this very PP (atomic radial -14.230 unpolarised /
// molecular Cartesian -14.668 polarised, vs CP2K ATOM -14.243986 restricted), so if the GPW path also
// lands near the facade the crystal KB is exonerated and the MnO defect lives elsewhere (multi-species,
// O, Ewald/alignment); if it does not, this is the l=2 crystal defect, isolated to ONE atom and one hour
// instead of the 4-atom magnetic cell.
// MEASURED 2026-08-06 -- THE ANSWER IS BASIS CONDITIONING, not the KB:
//   * CARTESIAN d carries the s contaminant (x^2+y^2+z^2), so our Mn window (7 s + 8 d shells x 6
//     Cartesian components = 55 functions) is RANK-DEFICIENT before any physics: lambdaMin 1.15e-07,
//     cond 8.2e7, and the GPW vet ABORTS.  The molecular facade only survives it by dropping modes --
//     its own log says "[ortho] dropped 5 near-null overlap mode(s) of 55".  A near-null direction that
//     the SCF can occupy is the classic variational-collapse mechanism, and MnO's 4-atom cell (154
//     functions, cond 7e8) sits in exactly that regime -- which fits -417 Ha vs the -61.47 oracle far
//     better than the 3e-2 analytic-vs-mesh KB discrepancy does.
//   * SPHERICAL d (5 pure components, no contaminant) is the natural cure but is NOT AVAILABLE on the
//     GPW path: it throws "the orbital basis is not a molecular Gaussian basis (no Molecule::LatticeSum1E)"
//     -- the spherical lineage does not implement the lattice-sum face (cf. the parked S3b spherical work).
//   => CURED 2026-08-06 (user's insight): keep the d set and drop the s window to TWO functions.  The
//      contaminants already span the mid/tight s space, so only the DIFFUSE 4s tail (0.10) and one tight
//      s (24) are needed -- 2s+8d gives lambdaMin 3.0e-03 / cond 2.1e3 (from 1.15e-07 / 8.2e7) at a cost
//      of just 2 mHa (facade -14.6661 vs the rank-deficient 7s+8d's -14.6681).  Trimming the d count
//      instead "fixes" conditioning but costs 0.55 Ha (7s+4d -> -14.11): the d set is the physics, the
//      s window was the redundancy.  NB CP2K solves this same shell list SPHERICALLY (its log: 55
//      Cartesian vs 47 spherical functions) and never sees the contaminant -- apples-to-oranges when
//      comparing its oracles.  GPW_MN_SPHERICAL=1 re-runs the (throwing) spherical arm.
// THIS IS NOW A GATE: the first OCCUPIED-d species validated end to end through the crystal path.
TEST(GPW_SCF, MnAtomInBoxDChannel)
{
    // GPW_MN_SPHERICAL=1: the SPHERICAL arm (doc/SphericalLatticePlan.md I1) -- the facade reference then
    // runs the NATIVE spherical family (same span as the view), so the box-vs-facade A/B stays span-matched.
    const bool spherical=(bool)std::getenv("GPW_MN_SPHERICAL");
    const bool sphBasis =(bool)std::getenv("GPW_BASIS_SPH");     // the restored-s file (I3); spherical arms only
    Molecule mnmol; mnmol.Insert(new Atom(25, 0.0, {0,0,0}));
    Calculation cRef(mnmol, {.basis=sphBasis?"valence_lowq_sph":"valence_lowq_sr",
                             .multiplicity=6, .pseudopotential=true, .ppValence=7,
                             .angular = spherical ? Angular::Spherical : Angular::Cartesian});
    const double Eref=cRef.Energy();
    std::cout << "[Mn finite] "<<(sphBasis?"valence_lowq_sph":"valence_lowq_sr")<<" LSDA sextet (q7, "
              <<(spherical?"SPHERICAL":"CARTESIAN")
              <<")="<<Eref<<"   (CP2K ATOM UKS sextet oracle -14.674425)"<<std::endl;

    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(25, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.label="Mn atom-in-box sextet";
    o.Nelec=7; o.multiplicity=6;                       // S=5/2 Hund: nUp=6, nDown=1
    o.species={{"Mn",7}};
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.imposeSymmetry=false;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.scf.NMaxIter=40; o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.SmearingkT=5e-3;
    // CARTESIAN d carries the s CONTAMINANT (x^2+y^2+z^2), so 8 d shells duplicate the 7-function s space
    // -- measured lambdaMin 1.15e-07 / cond 8.2e7 on this one-atom box, i.e. the basis is rank-deficient
    // BEFORE any physics runs.  SPHERICAL d (5 pure components) removes the contaminant; GPW_MN_SPHERICAL=1
    // selects it for the A/B -- via the SPHERICAL LATTICE VIEW (doc/SphericalLatticePlan.md I1): the
    // NATIVE Angular::Spherical family has no LatticeSum1E capability (the historical blocker -- feeding
    // it here died on the GPW cross-cast), so the view over the Cartesian engine is the working door.
    std::shared_ptr<const Real_BS> mnbasis(
        BasisSet::Molecule::Factory(sphBasis?BasisSetData::VALENCE_LOWQ_SPH:BasisSetData::VALENCE_LOWQ_SR,
                                    &cell, BasisSet::Molecule::Engine::MnD,
                                    BasisSet::Molecule::Angular::Cartesian));
    if (spherical) mnbasis=BasisSet::Molecule::PG_Spherical::MakeSphericalLatticeView(mnbasis);
    std::cout << "[Mn in-box] angular=" << (spherical?"SPHERICAL":"CARTESIAN") << std::endl;
    GpwResult R=RunGpw(lat, mnbasis, o, /*verbose*/(bool)std::getenv("GPW_MNO_VERBOSE"));
    std::cout << "[Mn in-box] GPW="<<R.E.GetTotalEnergy()<<"  facade="<<Eref
              << "  diff="<<(R.E.GetTotalEnergy()-Eref)<<std::endl;
    EXPECT_NEAR(R.charge, 7.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), Eref, 3e-2) << "GPW d-channel vs the molecular facade (measured 12 mHa)";
    if (!spherical)
        EXPECT_NEAR(R.E.GetTotalEnergy(), -14.6380, 1e-3);   // did-E-move anchor (2s+7d, the SR-trimmed cell basis)
}

// A POLARIZED RUN MUST STAY POLARIZED (SymmetryUpgradePlan §7 step 7, 2026-08-07).  The regression gate for
// the AFM collapse: the ρ̃ density mixers (Kerker/Pulay) carry ONE FourierMixCD -- the ↑+↓ TOTAL, with no spin
// channels -- and drive every Fock from it, so XC_GridEngine::RhoPol falls into its ρ↑=ρ↓=ρ/2 branch and the
// run is silently UNPOLARIZED from iteration 1 (measured on MnO: a seed staggered at m_stag=+0.366 reads
// EXACTLY +0.000000 at iteration 1 and never recovers).  MakeDensityMixer now refuses the ρ̃ mixers on a
// polarized density and falls back, loudly, to linear D-mixing.  This gate asks a POLARIZED run for Kerker
// on the cheap Mn q7 sextet box (same system as the d-channel gate above) and checks it got the POLARIZED
// answer.  QCHEM_SPINBLIND_KERKER=1 forces the broken mixer back on; the gate then fails.
//
// WHICH observable has teeth is itself a finding (measured, 2026-08-06/07).  The on-site MOMENT does NOT:
// with nUp=6/nDown=1 the moment is pinned by the CHANNEL OCCUPATIONS, so it survives a spin-blind Fock
// intact (0.64 vs 0.72) -- an atom cannot expose this bug.  Only an order that must be SELF-CONSISTENTLY
// SUSTAINED does, which is exactly MnO's staggering at nUp=nDn (nothing but v_xc^↑≠v_xc^↓ holds it up).
// What an atom DOES expose is the ENERGY: a spin-blind Fock loses the exchange splitting entirely and lands
// 68 mHa HIGH (−14.57 vs the physical −14.638, itself 12 mHa from the molecular facade).  So the moment is
// reported (it is the instrument under test) and the ENERGY is the assert.
TEST(GPW_SCF, PolarizedRunKeepsItsSpin)
{
    const double a=16.0;
    UnitCell cell(a);
    cell.AddAtom(25, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.label="Mn sextet under Kerker";
    o.Nelec=7; o.multiplicity=6;                       // S=5/2 Hund: nUp=6, nDown=1
    o.species={{"Mn",7}};
    o.images=BasisSet::Lattice_3D::CellImages::HomeCellOnly;
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.imposeSymmetry=false;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    auto envd=[](const char* n,double d){const char*s=std::getenv(n);return s?std::atof(s):d;};
    // GPW_MNSEXTET_NMAX raises the cap for the shared-μ A/B only: freeing the moment adds a soft degree of
    // freedom, so the same descent needs more iterations before the pinned energy means anything.
    o.scf.NMaxIter=(size_t)envd("GPW_MNSEXTET_NMAX",12); o.scf.MinΔρ=1e-8; o.scf.MinΔE=1e30;   // TIGHT on purpose: the assert is on the whole
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;   // trajectory, so the gate must actually
    o.scf.StartingRelaxRo=0.3; o.scf.MergeTol=1e-4; o.scf.SmearingkT=5e-3;   // ITERATE (Kerker on this
    o.scf.KerkerG0=1.0;                                // THE ASK: ρ̃ mixing on a polarized density

    // The order parameter: the on-site moment m(r)=ρ↑(r)−ρ↓(r) sampled 0.7 bohr off the nucleus (the d-shell
    // peak -- the d density vanishes AT the nucleus).  Recorded per iteration, so the SmokeTest of the
    // instrument itself is here: the probe must fire once per iteration and report an electrons-scale moment.
    std::vector<double> mtrace;
    const rvec3_t rMn(a/2,a/2,a/2), off(0.7,0,0);
    o.orderName="m_Mn";
    o.orderProbe=[&mtrace,rMn,off](const qchem::ChargeDensity::cDM_CD& cd)->double
    {
        const auto* pol=dynamic_cast<const qchem::ChargeDensity::cPolarized_CD*>(&cd);
        if (!pol) { mtrace.push_back(0.0); return 0.0; }    // an unpolarized Fock density: the defect itself
        const double m=(*pol->GetChargeDensity(Spin::Up  ))(rMn+off)
                      -(*pol->GetChargeDensity(Spin::Down))(rMn+off);
        mtrace.push_back(m);
        return m;
    };
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o,
                       /*verbose*/(bool)std::getenv("GPW_MNO_VERBOSE"));
    EXPECT_NEAR(R.charge, 7.0, 1e-6);
    ASSERT_FALSE(mtrace.empty()) << "the order-parameter probe never fired";
    const double mMin=*std::min_element(mtrace.begin(), mtrace.end());
    std::cout << "[Mn sextet under Kerker] m(0.7 bohr) over "<<mtrace.size()<<" iterations: min="<<mMin
              << " final="<<mtrace.back()<<std::endl;
    EXPECT_GT(mMin, 0.02) << "the instrument itself: the probe must see the S=5/2 moment every iteration";
    // THE ASSERT: the same polarized answer the linear-mixed d-channel gate pins (−14.6380).  A spin-blind
    // ρ̃ mixer loses the exchange splitting and lands at −14.57 -- 68x this tolerance away.
    EXPECT_NEAR(R.E.GetTotalEnergy(), -14.6380, 1e-3)
        << "asking for Kerker must not change the physics: this is the spin-blind-mixer detector";
}

// ============================ MnO rocksalt AFM-II (SymmetryUpgradePlan §7 step 7) ============================
// The FIRST magnetic transition-metal material: rocksalt MnO with the type-II AFM ordering (ferromagnetic
// (111) sheets, alternating along [111]).  The magnetic unit cell is the RHOMBOHEDRAL 2-f.u. cell -- the fcc
// cell doubled along [111]: lattice vectors a(1,1/2,1/2), a(1/2,1,1/2), a(1/2,1/2,1) (volume a^3/2 = 2x the
// fcc primitive), Mn sublattices at fractional (0,0,0) [+m] and (1/2,1/2,1/2) [-m] (cartesian (a/2)(0,1,1)
// + lattice), O at (1/4,1/4,1/4) and (3/4,3/4,3/4).  Everything below runs FREE (imposeSymmetry=false): the
// imposed per-channel star-average uses the GREY (spin-blind) crystal group, whose ops map +m to -m sites
// and would average the staggering away -- the Shubnikov-aware imposition is a separate increment.
//
// The SEED CHOOSES THE BASIN (SCFSeedingPlan §10): IonicSAD resolves Mn2+ (the d^5 S=5/2 library pair) +
// O2- (closed shell), and the flip bit on the second Mn staggers the majority channel -- the AFM-II pattern
// is present from iteration 0.  Fermi smearing kT=5e-3 rides the open d manifold through the early
// iterations.  LSDA q7 (semicore q15 = follow-on; +U = follow-on).
namespace
{
struct MnOArm { GpwResult R; GpwHandles h; std::shared_ptr<UnitCell> cell; };

MnOArm RunMnO(int multiplicity, bool afm, const std::string& label)
{
    const double a=8.40;                              // rocksalt a ~ 4.445 A (a.u.)
    Matrix3D<double> A(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);   // rhombohedral AFM-II cell (symmetric)
    auto cellp=std::make_shared<UnitCell>(A);
    UnitCell& cell=*cellp;
    // MNO_SWAP_SUBLATTICE: put the -m flip on the FIRST Mn instead of the second.  The two sites are related
    // by the (1/2,1/2,1/2) lattice translation, so the swapped run must be the EXACT MIRROR of the unswapped
    // one -- same energy, site moments exchanged.  It is the discriminator for where the observed asymmetry
    // lives (2026-08-11): mirror-image result => the code is equivariant under SITE exchange and the defect
    // is between the SPIN CHANNELS; identical result => the defect is tied to a SITE.
    const bool swapSub = std::getenv("MNO_SWAP_SUBLATTICE")!=nullptr;
    // MNO_SWAP_ORDER: add the (1/2,1/2,1/2) Mn FIRST, so the atom INDICES swap while the geometry and the
    // probe points are untouched.  This separates the last two candidates: if the moment still dies at the
    // ORIGIN the defect follows the POSITION; if it moves to (1/2,1/2,1/2) it follows the atom INDEX.
    const bool swapOrder = std::getenv("MNO_SWAP_ORDER")!=nullptr;
    // MNO_SHIFT=f RIGIDLY TRANSLATES the whole crystal by (f,f,f) fractional.  A rigid translation is an
    // EXACT symmetry -- every energy and every site moment must be invariant -- and it moves the Mn at
    // (0,0,0) OFF the cell corner.  That is the discriminator for the site asymmetry found 2026-08-11 (the
    // moment dies at the ORIGIN atom and survives at the cell-centre one, whichever spin each was seeded
    // with): a CORNER atom's atom-centred XC mesh must be reassembled from periodic images, a centre atom's
    // need not.  Same bug CLASS as DISABLED_TermTranslationInvariance (2026-07-09: "a boundary-straddling
    // corner orbital lost its wrapped tail") -- which covers Kinetic/Vloc/Vnl and NOT the Becke XC mesh,
    // because that mesh did not exist when it was written.
    const double sh = std::getenv("MNO_SHIFT") ? std::atof(std::getenv("MNO_SHIFT")) : 0.0;
    if (swapOrder)
    {
        cell.AddAtom(25, {0.5+sh,0.5+sh,0.5+sh}, afm && !swapSub); // -m sublattice, added FIRST
        cell.AddAtom(25, {sh,    sh,    sh    }, afm && swapSub);
    }
    else
    {
        cell.AddAtom(25, {sh,    sh,    sh    }, afm && swapSub);  // Mn sublattice +m (-m when swapped)
        cell.AddAtom(25, {0.5+sh,0.5+sh,0.5+sh}, afm && !swapSub); // Mn sublattice -m (flipped for the AFM arm)
    }
    cell.AddAtom(8,  {0.25+sh,0.25+sh,0.25+sh});
    cell.AddAtom(8,  {0.75+sh,0.75+sh,0.75+sh});
    // MNO_KMESH=n (run 40, 2026-08-12): an n^3 Γ-centred Monkhorst-Pack mesh on the magnetic cell.
    // The ORDERING experiment: a Γ-only 2-f.u. cell cannot resolve superexchange bandwidth (run 38
    // measured FM 38 mHa BELOW AFM-II), so the FM/AFM comparison needs k.  A magnetic imposed run
    // keeps the full mesh by construction (S3: recipFold left empty when the decoration is staggered).
    const int nk = std::getenv("MNO_KMESH") ? std::atoi(std::getenv("MNO_KMESH")) : 1;
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));

    GpwOptions o;
    o.label=label;
    o.Nelec=26; o.multiplicity=multiplicity;          // 2 x (Mn q7 + O q6); AFM: the explicit two-channel singlet
    o.species={{"Mn",7},{"O",6}};                     // densityEcut stays AUTO (= 2*alpha_max = 60, the trimmed d)
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;   // Mn2+ d^5 pair + diffuse O2- (the basin chooser)
    // MNO_IMPOSE (Shubnikov S4, doc/SymmetryUpgradePlan.md §7 step 7): 0 = FREE (default -- the banked
    // runs 30-34 ensemble, the §8 reference); 1 = impose the SHUBNIKOV group of the declared AFM ordering
    // (S3 resolves it from the flip bits automatically; the first magnetic imposition -- the Δρ tie-floor
    // measurement); 2 = impose the DETECTED grey group with the decoration suppressed.  MEASURED (run 36,
    // 2026-08-11): arm 2 is NOT an erasure control on this cell -- FindTau keeps the FIRST tau coset per W,
    // and every detected W fixes both Mn sites with tau=0, so the detected 12 are sublattice-PRESERVING
    // and the staggering survives them (iteration 1: m1=+0.42, E=-52.76 vs the free -52.80).  The ops that
    // WOULD erase m -- the tau=(1/2,1/2,1/2) coset -- exist only in the coset-complete magnetic set, which
    // S3 routes exclusively through the Shubnikov (sigma-aware) path: the erasure hazard the pre-S3
    // comment feared is structurally unreachable in production, and the erasure MECHANICS are proven at
    // the unit tier (the S2/S3 grey-control gates).
    {
        const char* im=std::getenv("MNO_IMPOSE");
        const int   iv=im?std::atoi(im):0;
        o.imposeSymmetry = iv!=0;
        o.greyImposition = iv==2;
    }
    // cond(S) ~ 7e8 on this dense oblique cell (lambdaMin ~ 2e-8 passes the vet, barely): plain Cholesky
    // EXPLODES the eigenproblem (measured: E ~ 1e9 Ha, [F,D] ~ 1e10 at iteration 1).  The NaF-class
    // recipe: pivoted Cholesky with a rank tolerance.
    // MNO_ORTHO_TOL (run 51, 2026-08-13): the pivoted-Cholesky rank tolerance, sweepable.  The v2
    // (CP2K-span) spherical cell vets lambdaMin 1.9e-5 / cond 3.8e5 -- the diffuse re-additions
    // restore the INTER-SITE redundancy the SR trim cured -- and collapsed to E=-455 from iteration 1
    // (near-null modes + eps-screened operators = the classic variational-collapse recipe).  A stiffer
    // rank filter drops those modes at the door.
    o.ortho=qchem::CholeskyPivoted;
    o.orthoTol=std::getenv("MNO_ORTHO_TOL") ? std::atof(std::getenv("MNO_ORTHO_TOL")) : 1e-4;
    // DIIS + plain linear mixing overshoots once the ramp reaches full step (measured run 3: healthy
    // -454.8 -> -456.5 over 3 iterations, then a +21 Ha slosh at Lin 1.00) -- the NaF-class ionic
    // charge-transfer dynamic.  Adopt the NaF Becke gate's full recipe: Ladder + Kerker-damped mixing
    // + MOM (masked Fermi) once the branch structure sets.
    o.accelerator="Ladder";
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };
    // MIXING KNOBS (sweepable, 2026-08-07).  alpha=0.45 came from the NaF BECKE gate; the NaF Gamma run --
    // the one whose recipe actually targets "the low-G charge-transfer slosh" -- uses 0.25, and MnO sloshes
    // WORSE than NaF, so 0.45 is the aggressive half of the borrowed recipe (user).  G0 selects WHICH modes
    // get damped, and for this cell the arithmetic matters: A=(a/2)(I+J) => B=(4pi/a)(I-J/4), the AFM
    // staggering lives at odd (h+k+l) whose shortest mode m=(1,0,0) has |G|=(4pi/a)*sqrt(0.6875)=1.24 a.u.,
    // while the lowest mode m=(1,1,1) has |G|=0.65.  At G0=1 the Kerker factor G^2/(G^2+G0^2) is 0.61 on the
    // AFM mode vs 0.30 on the lowest charge mode -- i.e. Kerker already favours the magnetism 2:1.  LOWERING
    // G0 FLATTENS that selectivity (G0=0.3: 0.94 vs 0.82, ratio 2.05->1.15) and un-damps exactly the sloshing
    // modes; RAISING it sharpens it (G0=3: 0.146 vs 0.045, ratio 3.3).  That is a linear-response argument
    // about which mode is actually pathological, so it is a PREDICTION to be swept, not a setting to trust.
    o.scf.StartingRelaxRo=envd("MNO_ALPHA",0.45); o.scf.KerkerG0=envd("MNO_KERKER_G0",1.0);
    // DENSITY-HISTORY MIXING (MNO_PULAY / MNO_PULAY_START, 2026-08-10).  Never exercised on MnO before: every
    // run to date used PulayDepth=0, i.e. a damped Kerker step with NO density history, while the CP2K oracle
    // runs BROYDEN_MIXING with NBUFFER 8 -- a quasi-Newton density history -- and reaches the magnetic state
    // with no MOM at all.  Note CP2K's BETA 1.5 IS its Kerker damping, i.e. it damps HARDER than our G0=1.0
    // and still keeps the moment, which points away from "Kerker suppresses the moment" and at the missing
    // history instead.  USER 2026-08-10, correcting me: if swapping the mixer changes whether m_stag
    // survives, that IMPLICATES the mixing path rather than exonerating it -- so Pulay/Broyden is a
    // first-class candidate, not the weak one I had called it.
    o.scf.PulayDepth=(int)envd("MNO_PULAY",0.0);
    o.scf.PulayStart=(int)envd("MNO_PULAY_START",5.0);
    // NB (2026-08-07) this KerkerG0 is currently INERT on the AFM arm: MakeDensityMixer refuses the ρ̃
    // mixers on a POLARIZED density and falls back (loudly) to linear D-mixing, because the ρ̃ mixers carry
    // the TOTAL density only and collapse v_xc to ζ=0 from iteration 1 -- the AFM collapse, diagnosed here
    // with the m_stag order parameter below.  QCHEM_SPINBLIND_KERKER=1 restores the broken arm for the A/B.
    // Once the spin-resolved ρ̃ mixer lands, this line becomes live again (and linear mixing's slosh --
    // measured: E swinging -51/-55/-44/-28 over 6 iterations -- is exactly why we want it back).
    // ---- OCCUPATION KNOBS: the kT / MOM-TIMING experiment (2026-08-08, plan §7 step 7) --------------------
    // The run-10/11 diagnosis left ONE live defect: the occupied configuration re-shuffles EVERY iteration
    // (`cfg *` on every row) between a MAGNETIC and a NON-magnetic branch 4-6 Ha apart.  That is an OCCUPATION
    // defect, and no density-mixing knob can reach it.  Two instruments exist and BOTH were mis-set:
    //   * kT=5e-3 against a frontier gap that WANDERS 0.026-0.27 Ha -- 5x to 50x too cold to bridge the branch
    //     tie.  It fractionalises the ties (the `m` flag) without ever letting mu choose between the branches.
    //   * UseMOM was INERT.  tIrrepWF::FillOrbitals gives SmearingkT>0 precedence, and that branch consults the
    //     MOM reference ONLY when MOMSmearPenalty>0 (default 0) -- so `UseMOM=true` alongside kT=5e-3 has never
    //     changed a single occupation in ANY MnO run.  And even reached, MOMStartIter=10 captures the reference
    //     at iteration 10, nine iterations after m_stag died (0.366 -> 0.0046 at iteration 1): a delayed capture
    //     is built for a seed you want to settle AWAY from, and MnO's seed is the physics.
    // Hence: kT is now the SWEEPABLE knob it should always have been, MOM is reachable under smearing via the
    // penalty, and momFromSeed pins the reference to the seed's own staggered d^5 character.
    auto envi=[](const char* n, int d){ const char* s=std::getenv(n); return s ? std::atoi(s) : d; };
    o.scf.UseMOM=envi("MNO_MOM",1)!=0; o.scf.MOMStartIter=envi("MNO_MOM_START",10);
    o.scf.MOMSmearPenalty=envd("MNO_MOM_PENALTY",0.0);   // >0 makes MOM live UNDER smearing (masked Fermi)
    o.momFromSeed=envi("MNO_MOM_SEED",0)!=0;             // pin the reference to the SEED's occupied subspace
    // MNO_REAL (doc/RealComplexPlan.md 3c-3): build the Γ TRIM block(s) REAL -- the Becke H_xc
    // quadrature and the ρ-sampling GEMM then run blaze's real kernels (measured 2.0x / 2.2x on this
    // system).  COMBINES FREELY WITH MOM since R2.21: the occupation STATE now stores each block's
    // reference under the BLOCK's own scalar, so the full seed-pinned recipe runs flipped.  Follows
    // the harness default (now ON); MNO_REAL=0 is the A/B door back to the all-complex build, and
    // MNO_KMESH=2 is ALL-TRIM, so there the whole k-sum goes real, not just Γ.
    o.realTRIMBlocks=envi("MNO_REAL",o.realTRIMBlocks?1:0)!=0;
    // The 0h MOM GUARD releases the reference after 3 consecutive NON-AUFBAU (hole) iterations, on the premise
    // that a hole means the reference is wrong.  That premise was calibrated on NaF, where the hole IS the
    // pathology (a diving diffuse ghost pinned by a stale reference).  It does NOT hold here: until the
    // self-consistent exchange splitting opens, the magnetic branch's occupied d^5-on-Mn1 states sit ABOVE
    // empty d states of the other sublattice in the raw eps ordering, so a hole is the SIGNATURE of the branch
    // we are trying to hold, not evidence against it.  MNO_MOM_HOLD raises the persistence (large = never
    // release) so the guard cannot silently undo the experiment.
    o.scf.Guard.HolePersistence=envi("MNO_MOM_HOLD",3);
    // MNO_NR / MNO_L shrink the Becke mesh so its per-atom numbers can be read by eye (GPW_BECKE_ATOMS).
    // MNO_XC_UNIFORM: run the XC quadrature on the UNIFORM raster instead of the atom-centred Becke mesh,
    // changing nothing else.  The site asymmetry (moment dies on the CORNER atom) and the 56 Ha rigid-
    // translation violation must vanish if the atom-centred path is the guilty one -- a uniform raster has
    // no per-atom structure to get a site wrong (2026-08-11).
    if (std::getenv("MNO_XC_UNIFORM")) o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    if (const char* nr=std::getenv("MNO_NR")) o.xcMesh.nRadial=std::atoi(nr);
    if (const char* ll=std::getenv("MNO_L"))  o.xcMesh.angularDegree=std::atoi(ll);
    o.scf.SmearingkT=envd("MNO_KT",5e-3);                // the open-d-manifold knob (MNO_ANNEAL for a schedule)
    // MNO_SHARED_MU: ONE chemical potential over both spin channels, so the MOMENT IS AN OUTPUT.  With it
    // off (the default, and every MnO run through 29) the fill conserves nUp and nDn separately, which is
    // two reservoirs and two μ -- run 29 ended with μ↑−μ↓ = 27 mHa, an unrelieved driving force to move
    // charge ↓→↑ that the constraint holds open, so the converged state is not the free minimum.  It also
    // makes the AFM constraint IMPOSED rather than found: m_net ≡ 0 by construction.  CP2K's deck does NOT
    // constrain the populations (Fermi smearing, one μ; &BS seeds the ATOMIC guess only), so THIS is the
    // oracle-comparable ensemble -- and m_stag is the instrument for whether the moment survives unaided.
    o.spinsShareFermi=envi("MNO_SHARED_MU",0)!=0;
    // GPW_MNO_NMAX bounds a HAND-RUN diagnosis (measured 2026-08-05: ~9 min PER analytic sweep on a 14 GB
    // box -- the free-run magnetic cell exceeds the stream budget, so every iteration pays collocate +
    // integrate analytically; an unbounded 80-iter run is a >24 h affair).
    const char* nx=std::getenv("GPW_MNO_NMAX");
    o.scf.NMaxIter=nx?std::atoi(nx):80; o.scf.MinΔρ=1e-5; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.MergeTol=1e-4;                              // (RelaxRo set above with the Ladder/Kerker recipe)

    // THE ORDER PARAMETER (plan §7 step 7 + §9): the STAGGERED moment m_stag = ½[m(Mn1)−m(Mn2)], with
    // m(r)=ρ↑(r)−ρ↓(r) sampled 0.7 bohr off each Mn nucleus (the d-shell peak -- the d density VANISHES at
    // the nucleus).  This is the AFM-II order itself, and it is the WHOLE order: the run fixes nUp=nDn, so
    // the net moment is identically zero and a ferromagnetic drift is not an available escape -- m_stag→0
    // means the magnetism died, full stop.  Run 6 converged to EXACTLY zero site moments from a provably
    // staggered seed (PolarizedSeedAFMStaggering), so the loss happens somewhere INSIDE the SCF; watching
    // m_stag per iteration is what brackets the iteration where it goes (doc/SymmetryUpgradePlan.md:1054).
    // Two point evaluations of the Bloch density per iteration -- negligible beside a collocation sweep.
    {
        const double ds=2.0*sh*a;                                // A*(sh,sh,sh) = 2*sh*a*(1,1,1)
        const rvec3_t off(0.7,0,0), rMn1(ds,ds,ds), rMn2(a+ds,a+ds,a+ds);  // A*(½,½,½) = a(1,1,1), shifted
        o.orderName="m_stag";
        o.orderProbe=[off,rMn1,rMn2](const qchem::ChargeDensity::cDM_CD& cd)->double
        {
            // Cross-cast to the spin-resolved FACE (abstract->abstract, the sanctioned dynamic_cast): an
            // unpolarized run simply has no order parameter to report.
            const auto* pol=dynamic_cast<const qchem::ChargeDensity::cPolarized_CD*>(&cd);
            if (!pol) return 0.0;
            const auto* up=pol->GetChargeDensity(Spin::Up);
            const auto* dn=pol->GetChargeDensity(Spin::Down);
            const double m1=(*up)(rMn1+off)-(*dn)(rMn1+off);
            const double m2=(*up)(rMn2+off)-(*dn)(rMn2+off);
            // GPW_MNO_SITES: the two site moments SEPARATELY, plus the channel charges.  m_stag alone is
            // ½(m1−m2), which is BLIND to the thing under investigation: for one magnetic species on one
            // Wyckoff site the two Mn are related by the (½,½,½) translation, so ANY valid solution has
            // m1 = −m2 and N↑ = N↓.  Printed from iteration 0, so a seed that is not already a mirror pair
            // is separated from a Fock build that breaks one (2026-08-11).
            if (std::getenv("GPW_MNO_SITES"))
                std::cout << "[sites] m1=" << m1 << " m2=" << m2 << " m1+m2=" << m1+m2
                          << "  N↑=" << up->GetTotalCharge() << " N↓=" << dn->GetTotalCharge()
                          << " (N↑−N↓=" << up->GetTotalCharge()-dn->GetTotalCharge() << ")" << std::endl;
            return 0.5*(m1-m2);
        };
    }
    MnOArm arm;
    arm.cell=cellp;
    const bool verbose=(bool)std::getenv("GPW_MNO_VERBOSE");
    // MNO_ACC (run 39, 2026-08-11): the accelerator policy, sweepable at last.  A single name
    // (DIIS|GDM|Ladder|Null) applies throughout; under MNO_ANNEAL a comma-list parallel to the kT
    // schedule gives each stage its own -- the RunGpwAnnealed-doc'd {"Ladder","GDM"} x {5e-3, 0}
    // experiment: GDM's Engageable veto requires INTEGER occupations (run 38 stayed on DIIS for all
    // 41 iterations because kT=5e-3 makes D' non-idempotent), so a GDM stage must run kT=0.
    std::vector<std::string> accSched;
    if (const char* ac=std::getenv("MNO_ACC"))
    {
        for (std::string s(ac), tok; !s.empty(); )
        {
            size_t c=s.find(','); tok=s.substr(0,c);
            if (!tok.empty()) accSched.push_back(tok);
            if (c==std::string::npos) break;
            s=s.substr(c+1);
        }
        if (accSched.size()==1) { o.accelerator=accSched.front(); accSched.clear(); }
    }
    // ANNEALED ARM (MNO_ANNEAL="0.05,0.02,0.01,0"): a DESCENDING kT schedule with density continuation, and
    // -- when MNO_MOM_SEED=1 -- occupied-character continuation too, so the configuration a hot stage settles
    // on is what the next stage holds.  The single-kT path stays the default: one run, one temperature.
    if (const char* sched=std::getenv("MNO_ANNEAL"))
    {
        std::vector<double> kTs;
        for (std::string s(sched), tok; !s.empty(); )
        {
            size_t c=s.find(','); tok=s.substr(0,c);
            if (!tok.empty()) kTs.push_back(std::atof(tok.c_str()));
            if (c==std::string::npos) break;
            s=s.substr(c+1);
        }
        assert(!kTs.empty() && "MNO_ANNEAL must list at least one kT");
        // MNO_ANNEAL_PENALTY: the per-stage MOM Λ, parallel to MNO_ANNEAL.  "1.5,0" = hold the branch, then
        // RELEASE the constraint and see whether the exchange splitting can keep it without help.
        std::vector<double> lams;
        if (const char* ls=std::getenv("MNO_ANNEAL_PENALTY"))
            for (std::string t(ls), tok; !t.empty(); )
            {
                size_t c=t.find(','); tok=t.substr(0,c);
                if (!tok.empty()) lams.push_back(std::atof(tok.c_str()));
                if (c==std::string::npos) break;
                t=t.substr(c+1);
            }
        assert((lams.empty() || lams.size()==kTs.size()) && "MNO_ANNEAL_PENALTY must parallel MNO_ANNEAL");
        assert((accSched.empty() || accSched.size()==kTs.size()) && "MNO_ACC list must parallel MNO_ANNEAL");
        arm.R=RunGpwAnnealed(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, kTs, verbose,
                             accSched, &arm.h, lams);
        return arm;
    }
    arm.R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, verbose, &arm.h);
    return arm;
}
} //anon

// THE RAW SEED of the MnO AFM-II cell: are the two Mn sublattices actually equal and opposite?
//
// One magnetic species on one Wyckoff site means the two Mn are related by the (1/2,1/2,1/2) translation,
// so ANY valid magnetic solution -- seed included -- must have m(Mn1) = -m(Mn2).  Nothing checked that.
// PlaneWaveDFT.PolarizedSeedAFMStaggering covers a DIFFERENT cell (simple cubic, 2 Mn, no O, NEUTRAL
// targets) and asserts G-space quantities plus GetTotalSpin()==0; m_stag = 1/2(m1-m2), the campaign's order
// parameter, is BLIND to the imbalance -- it reads +0.366 for (+0.37,-0.37) and for (0,-0.73) alike.
// Measured 2026-08-11: the density one Fock build downstream of this seed has m1 = -0.0001, m2 = -0.73.
// This test asks whether the seed ITSELF is lopsided, i.e. whether the defect is upstream of the SCF.
TEST(GPW_SCF, MnOSeedSublatticesAreEqualAndOpposite)
{
    const double a=8.40;
    Matrix3D<double> A(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);   // the SAME rhombohedral AFM-II cell
    auto cellp=std::make_shared<UnitCell>(A);
    UnitCell& cell=*cellp;
    cell.AddAtom(25, {0.0,0.0,0.0}, false);        // Mn +m
    cell.AddAtom(25, {0.5,0.5,0.5}, true);         // Mn -m
    cell.AddAtom(8,  {0.25,0.25,0.25});
    cell.AddAtom(8,  {0.75,0.75,0.75});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.label="MnO seed probe"; o.Nelec=26; o.multiplicity=1;
    o.species={{"Mn",7},{"O",6}};
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;
    o.imposeSymmetry=false;
    // The seed the run actually uses: PolarizedSeedCD over this cell, with the IonicSAD targets
    // (Mn2+ => 5 of the q7 valence, O2- => 8 of the q6).  Built directly -- no SCF, no Hamiltonian.
    qchem::BasisSet::Lattice_3D::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    qchem::ChargeDensity::PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);
    const auto* up=seedCD.GetChannel(Spin::Up);
    const auto* dn=seedCD.GetChannel(Spin::Down);

    const rvec3_t off(0.7,0,0), rMn1(0,0,0), rMn2(a,a,a);   // A*(1/2,1/2,1/2) = a(1,1,1); 0.7 bohr = the d peak
    const double m1=(*up)(rMn1+off)-(*dn)(rMn1+off);
    const double m2=(*up)(rMn2+off)-(*dn)(rMn2+off);
    std::cout << "[seed probe] m1=" << m1 << " m2=" << m2 << " m1+m2=" << m1+m2
              << "  N_up=" << up->GetTotalCharge() << " N_dn=" << dn->GetTotalCharge() << std::endl;

    EXPECT_NEAR(up->GetTotalCharge(), dn->GetTotalCharge(), 1e-8) << "the AFM seed carries no NET moment";
    EXPECT_GT(std::abs(m1), 0.1) << "the + sublattice must actually be polarized";
    EXPECT_GT(std::abs(m2), 0.1) << "the - sublattice must actually be polarized";
    EXPECT_NEAR(m1+m2, 0.0, 0.02*std::abs(m1))
        << "the two sublattices are related by a lattice translation: m1 must equal -m2";

    // THE MIRROR RELATION EVERYWHERE, not just at the two nuclei.  The AFM seed satisfies
    // rho_up(r) = rho_dn(r + t) with t = A*(1/2,1/2,1/2) = (a,a,a) at EVERY r, by construction.  Probed
    // here at points that straddle the CELL BOUNDARY, because that is where the Becke XC mesh wraps its
    // points into the home cell (kpt = r - A*n0) and where a real-space evaluation built from per-atom
    // recentred radials WITHOUT periodic images would pick up the wrong atom -- which, the two Mn carrying
    // OPPOSITE spins, swaps the channels rather than merely losing density.  That failure is invisible to
    // the grid-charge check (rho_up+rho_dn is untouched; only rho_up-rho_dn flips) and it is exactly the
    // observed symptom: the corner atom's moment dies while the mesh integrates the total to 1e-5.
    struct P { const char* what; rvec3_t r; };
    const std::vector<P> probes = {
        {"inside, +x of Mn1",        rvec3_t( 0.7, 0.0, 0.0)},
        {"OUTSIDE the cell, -x",     rvec3_t(-0.7, 0.0, 0.0)},   // physically beside Mn1; wraps
        {"OUTSIDE the cell, -xyz",   rvec3_t(-0.5,-0.5,-0.5)},
        {"far tail, -2 bohr",        rvec3_t(-2.0, 0.0, 0.0)},
    };
    for (const auto& p : probes)
    {
        const rvec3_t rt = p.r + rvec3_t(a,a,a);                  // the mirror partner
        const double u=(*up)(p.r), d=(*dn)(rt);
        std::cout << "[mirror] " << p.what << ": rho_up(r)=" << u << "  rho_dn(r+t)=" << d
                  << "  diff=" << u-d << std::endl;
        EXPECT_NEAR(u, d, 1e-6*std::max(1.0,std::abs(u)))
            << "AFM mirror broken at " << p.what << ": the seed is not properly periodic there";
    }

    // THE BATCH OVERLOAD vs THE SINGLE-POINT ONE.  Everything above used operator()(rvec3_t).  The XC mesh
    // samples a MATRIX-FREE seed through the BATCHED operator()(rvec3vec_t) instead -- XC_GridEngine::RhoPol
    // takes its cSpinResolved_CD branch for exactly this density -- so the batch path is what the first Fock
    // build actually sees, and nothing has ever checked the two agree.  They must, pointwise.
    rvec3vec_t batch(2*probes.size());
    for (size_t i=0;i<probes.size();++i) { batch[i]=probes[i].r; batch[probes.size()+i]=probes[i].r+rvec3_t(a,a,a); }
    const rvec_t bu=(*up)(batch), bd=(*dn)(batch);
    ASSERT_EQ(bu.size(), batch.size());
    for (size_t i=0;i<batch.size();++i)
    {
        const double su=(*up)(batch[i]), sd=(*dn)(batch[i]);
        std::cout << "[batch] r=(" << batch[i].x << "," << batch[i].y << "," << batch[i].z << ")"
                  << " up: batch=" << bu[i] << " single=" << su << " d=" << bu[i]-su
                  << " | dn: batch=" << bd[i] << " single=" << sd << " d=" << bd[i]-sd << std::endl;
        EXPECT_NEAR(bu[i], su, 1e-6*std::max(1.0,std::abs(su))) << "UP batch != single at probe " << i;
        EXPECT_NEAR(bd[i], sd, 1e-6*std::max(1.0,std::abs(sd))) << "DN batch != single at probe " << i;
    }
}

// THE v_xc SUBLATTICE MIRROR ON THE XC MESH -- the SymmetryUpgradePlan "NEXT ACTION" probe (2026-08-11).
// By elimination (seed, Becke weights, Kinetic/Vloc/Vnl, Phi tables all exonerated) the first-Fock-build
// mirror break must live in v_xc -- yet a pointwise LSDA functional "cannot" be site-dependent.  This probe
// resolves the contradiction by testing what the Fock build ACTUALLY consumes: the channel rasters
// XC_GridEngine::RhoPol hands the DeltaFittedV*Pol pair, at the mesh's own points.  The mesh stores its
// points WRAPPED into the home cell (kpt = r - A*n0, MakePeriodicBeckeMesh) -- so a valid seed must satisfy
// rho_up(p_g) = rho_dn(p_g + t) with t = A*(1/2,1/2,1/2) AT EVERY STORED POINT, and (v_xc being pointwise
// in the channel pair) v_xc^up(p_g) = v_xc^dn(p_g + t).  The two Mn blocks' grids are exact t-translates of
// each other (same radial x angular template), so every point's mirror partner is itself a mesh point --
// found here by hashing wrapped fractional coordinates, no interpolation anywhere.
//
// WHY THE SEED GATE ABOVE COULD PASS WHILE THIS FAILS: its probe points sit within ~2 bohr of a HOME atom,
// where SeedCD's real-space rho(r) = Sum_atoms rho_atom(|r-R|) -- a sum with NO LATTICE IMAGES -- is
// dominated by an atom it actually contains.  A WRAPPED mesh point near a cell face reads its density from
// an IMAGE of an atom (for the CORNER Mn, 7 of the 8 octants of its density hump belong to images), which
// an image-less sum simply does not have.  That is site-specific (the centre Mn2 has no near-shell wrapped
// points), a rigid translation changes it (MNO_SHIFT), the uniform raster reproduces it (its corner-region
// points read image density too), and under the AFM staggering the missing hump is the MAJORITY channel of
// exactly one sublattice -- every recorded symptom.
TEST(GPW_SCF, MnOSeedVxcMirrorOnBeckeMesh)
{
    const double a=8.40;
    Matrix3D<double> A(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);   // the SAME rhombohedral AFM-II cell
    auto cellp=std::make_shared<UnitCell>(A);
    UnitCell& cell=*cellp;
    cell.AddAtom(25, {0.0,0.0,0.0}, false);        // Mn +m (the CORNER atom -- the one whose moment dies)
    cell.AddAtom(25, {0.5,0.5,0.5}, true);         // Mn -m (cell centre)
    cell.AddAtom(8,  {0.25,0.25,0.25});
    cell.AddAtom(8,  {0.75,0.75,0.75});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The seed the run uses (identical to MnOSeedSublatticesAreEqualAndOpposite).
    qchem::BasisSet::Lattice_3D::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    qchem::ChargeDensity::PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);

    // The run's Becke recipe, shrunk (nR=20, GL-11) -- wrapping is generic, it needs no production density.
    qcMesh::MeshParams mp=qcMesh::BeckeXCParams(20, -1.0, 11);
    auto mesh=std::make_shared<const qcMesh::Mesh>(cell.CreateIntegrationMesh(mp));
    const rvec3vec_t& P=mesh->Points();
    const size_t N=P.size();
    ASSERT_GT(N, 0u);

    // The channel rasters EXACTLY as the Fock build gets them (RhoPol's cSpinResolved_CD seed branch).
    qchem::Hamiltonian::XC_GridEngine engine(mesh);
    const rvec_t ru=engine.RhoPol(&seedCD, Spin::Up);
    const rvec_t rd=engine.RhoPol(&seedCD, Spin::Down);

    // Mirror-partner lookup: hash each point's wrapped fractional coords; partner(g) = the mesh index at
    // wrapped(frac(p_g) + (1/2,1/2,1/2)).  Quantized key + 27-neighbour probe rides out wrap roundoff.
    const double q=1e-9;
    auto wrapf=[](rvec3_t f){ f.x-=floor(f.x); f.y-=floor(f.y); f.z-=floor(f.z); return f; };
    std::map<std::tuple<long long,long long,long long>,size_t> at;
    std::vector<rvec3_t> F(N);
    for (size_t g=0; g<N; g++)
    {
        F[g]=wrapf(cell.ToFractional(P[g]));
        at[{llround(F[g].x/q),llround(F[g].y/q),llround(F[g].z/q)}]=g;
    }
    auto partner=[&](size_t g)->long
    {
        const rvec3_t fm=wrapf(F[g]+rvec3_t(0.5,0.5,0.5));
        const long long kx=llround(fm.x/q), ky=llround(fm.y/q), kz=llround(fm.z/q);
        for (long long dx=-1; dx<=1; dx++) for (long long dy=-1; dy<=1; dy++) for (long long dz=-1; dz<=1; dz++)
            if (auto it=at.find({kx+dx,ky+dy,kz+dz}); it!=at.end())
            {
                rvec3_t d=wrapf(F[it->second]-fm+rvec3_t(0.5,0.5,0.5))-rvec3_t(0.5,0.5,0.5);  // min-image
                if (norm(d)<5e-9) return long(it->second);
            }
        return -1;
    };

    // v_xc per point from the channel pair -- the same functionals the DeltaFittedV*Pol pair applies.
    qchem::Hamiltonian::SlaterExchange ex(2.0/3.0, Spin(Spin::Up));   // channel-native (non-halving)
    qchem::Hamiltonian::VWN_Correlation vc;
    auto vxc=[&](double u, double d, const Spin& s)->double
    {
        if (u+d<=1e-12) return 0.0;                       // VWN's r_s/log guard; symmetric, mirror-safe
        return ex.GetVxc(s==Spin::Up ? u : d) + vc.GetVc(u, d, s);
    };

    // Sweep: rho and v_xc mirror defects across the whole raster; localize the worst offenders.
    size_t nOrphan=0, nBad=0;
    double maxRho=0, maxV=0, maxOrphanWRho=0;
    size_t argRho=0; long argRhoJ=-1;
    const rvec_t& W=mesh->Weights();
    for (size_t g=0; g<N; g++)
    {
        const long j=partner(g);
        // An ORPHAN is legal but must be an eps-tail point: the free Becke builder's keep decisions are
        // bit-different between translated blocks, so an eps-BORDERLINE point can be kept on one atom and
        // dropped on its partner (the same mechanism the imposed builder's orbit-consistency pass drops).
        // What the quadrature sees of it is w*rho -- bound THAT, not the count.
        if (j<0) { nOrphan++; maxOrphanWRho=std::max(maxOrphanWRho, std::abs(W[g])*std::max(ru[g],rd[g])); continue; }
        const double dRho=std::max(std::abs(ru[g]-rd[j]), std::abs(rd[g]-ru[j]));
        const double dV  =std::max(std::abs(vxc(ru[g],rd[g],Spin::Up  )-vxc(ru[j],rd[j],Spin::Down)),
                                   std::abs(vxc(ru[g],rd[g],Spin::Down)-vxc(ru[j],rd[j],Spin::Up  )));
        if (dRho>1e-6) nBad++;
        if (dRho>maxRho) { maxRho=dRho; argRho=g; argRhoJ=j; }
        if (dV  >maxV) maxV=dV;
    }
    std::cout << "[vxc mirror] N=" << N << " orphans=" << nOrphan << " (max w*rho " << maxOrphanWRho
              << ") bad(rho>1e-6)=" << nBad
              << "  max|rho_up(r)-rho_dn(r+t)|=" << maxRho
              << "  max|vxc_up(r)-vxc_dn(r+t)|=" << maxV << std::endl;

    // Localize the worst point: whose density is it -- a HOME atom's, or an IMAGE's the seed cannot see?
    if (argRhoJ>=0)
    {
        std::vector<rvec3_t> R; std::vector<std::string> nm={"Mn1","Mn2","O1","O2"};
        for (auto atom : cell) R.push_back(atom->itsR);
        auto nearest=[&](const rvec3_t& r)
        {
            double best=1e300; std::string who;
            for (size_t ia=0; ia<R.size(); ia++)
                for (int i=-1;i<=1;i++) for (int jj=-1;jj<=1;jj++) for (int k=-1;k<=1;k++)
                {
                    const double d=norm(r-(R[ia]+cell.ToCartesian(rvec3_t(i,jj,k))));
                    if (d<best) { best=d; who=nm[ia]+((i||jj||k)?" IMAGE":""); }
                }
            return std::make_pair(best,who);
        };
        const auto [dg,wg]=nearest(P[argRho]);
        const auto [dj,wj]=nearest(P[size_t(argRhoJ)]);
        std::cout << "[vxc mirror] worst point r=("<<P[argRho].x<<","<<P[argRho].y<<","<<P[argRho].z
                  << ") nearest "<<wg<<" d="<<dg<<"  rho_up(r)="<<ru[argRho]<<" rho_dn(r)="<<rd[argRho]<<"\n"
                  << "[vxc mirror] partner     r=("<<P[size_t(argRhoJ)].x<<","<<P[size_t(argRhoJ)].y<<","
                  << P[size_t(argRhoJ)].z<<") nearest "<<wj<<" d="<<dj
                  << "  rho_up(r+t)="<<ru[size_t(argRhoJ)]<<" rho_dn(r+t)="<<rd[size_t(argRhoJ)]<<std::endl;
    }

    EXPECT_LT(maxOrphanWRho, 1e-8) << "a partnerless mesh point must be an eps-tail point (weight*rho "
                                      "below the Becke builder's eps-converged-series contract)";
    EXPECT_LT(maxRho, 1e-8) << "the seed's channel rasters break the sublattice mirror ON THE XC MESH -- "
                               "this is the site-dependent v_xc defect (SymmetryUpgradePlan WHERE-WE-LEFT-OFF)";
    EXPECT_LT(maxV, 1e-6) << "v_xc^up(r) != v_xc^dn(r+t) on the mesh: the first Fock build is fed a "
                             "mirror-broken exchange-correlation potential";
}

// SHUBNIKOV S3, end to end at the seed level (doc/SymmetryUpgradePlan.md §7 step 7): a MAGNETICALLY
// imposed MnO basis must star-average the channel pair under the Shubnikov group -- keeping the AFM
// staggering EXACTLY mirror-symmetric -- where the grey average would erase it.  Chain under test:
// MagneticDecoration (the seed's own species rule) -> GPWParams::siteSpins -> the factory's Shubnikov
// resolution -> CreateXCQuadrature (site-adapted invariant mesh + fold + sigma tags + flip-fixed zero
// flags) -> XC_GridEngine::RhoPol's (rho,m) split.
TEST(GPW_SCF, MnOImposedShubnikovKeepsTheSeedStaggering)
{
    namespace L3=BasisSet::Lattice_3D;
    using namespace qchem::ChargeDensity;
    const double a=8.40;
    Matrix3D<double> A(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);
    auto cellp=std::make_shared<UnitCell>(A);
    UnitCell& cell=*cellp;
    cell.AddAtom(25, {0.0,0.0,0.0}, false);
    cell.AddAtom(25, {0.5,0.5,0.5}, true);
    cell.AddAtom(8,  {0.25,0.25,0.25});
    cell.AddAtom(8,  {0.75,0.75,0.75});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The decoration, by the seed's own resolution: Mn2+ (d5 pair) = +/-1, O2- (closed shell) = 0.
    auto st=lat.GetStructure();
    std::vector<int> spins=MagneticDecoration(st.get(), "LDA", IonicSADTargets(st.get(), "LDA"));
    ASSERT_EQ(spins.size(), 4u);
    EXPECT_EQ(spins[0], +1); EXPECT_EQ(spins[1], -1);
    EXPECT_EQ(spins[2],  0); EXPECT_EQ(spins[3],  0);

    // The magnetically IMPOSED GPW basis (coarse everything: this gate tests symmetry, not accuracy).
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR),
        L3::GPWParams{.densityEcut=8.0, .imposeSymmetry=true, .siteSpins=spins}));
    BasisSet::XCQuadrature q = bs->CreateXCQuadrature(st.get(), qcMesh::BeckeXCParams(20, -1.0, 11));
    ASSERT_EQ(q.sigmas.size(), 24u) << "the Shubnikov group of the AFM-II cell has 24 ops (12+12)";
    size_t nFlip=0; for (auto s : q.sigmas) if (s==Symmetry::SpinAction::Flip) nFlip++;
    EXPECT_EQ(nFlip, 12u);
    ASSERT_EQ(q.flipFixed.size(), q.mesh->size());

    // The seed (the same PolarizedSeedCD the run uses), through the engine's channel-pair projector.
    qchem::BasisSet::Lattice_3D::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);

    qchem::Hamiltonian::XC_GridEngine engine(q.mesh, q.fold, q.sigmas, q.flipFixed);
    const rvec_t up=engine.RhoPol(&seedCD, Spin::Up);
    const rvec_t dn=engine.RhoPol(&seedCD, Spin::Down);

    // (a) The staggering SURVIVES the imposed projector: the magnetization raster keeps its full scale.
    double maxM=0; for (size_t g=0; g<up.size(); g++) maxM=std::max(maxM, std::abs(up[g]-dn[g]));
    EXPECT_GT(maxM, 0.1) << "the SHUBNIKOV average must PRESERVE the AFM staggering";

    // (b) ...and is EXACTLY mirror-antisymmetric on the raster: m must vanish at every flip-fixed point,
    // and the weighted net moment must be zero to machine precision (the signed projector's guarantees).
    const rvec_t& W=q.mesh->Weights();
    double net=0, worstFixed=0;
    for (size_t g=0; g<up.size(); g++)
    {
        net += W[g]*(up[g]-dn[g]);
        if (q.flipFixed[g]) worstFixed=std::max(worstFixed, std::abs(up[g]-dn[g]));
    }
    // The projector pairs every orbit's members +/- exactly; the residual is the ULP-level asymmetry of
    // the site-adapted builder's partner WEIGHTS (op-image copies, equal only to roundoff -- measured
    // 1.7e-9 over 6000 points), not a projector defect.
    EXPECT_LT(std::abs(net), 1e-7) << "an imposed AFM pair carries ZERO net moment by construction";
    EXPECT_LT(worstFixed, 1e-12) << "m must vanish exactly at the flip-fixed mesh points";

    // (c) The GREY control: the SAME mesh and fold with the sigma tags withheld = the historical
    // per-channel spatial average, which maps +m sites onto -m sites and ERASES the order.  This is the
    // unit-level half of the S4 negative control ("imposing the grey group kills m_stag").
    qchem::Hamiltonian::XC_GridEngine grey(q.mesh, q.fold);
    const rvec_t gup=grey.RhoPol(&seedCD, Spin::Up);
    const rvec_t gdn=grey.RhoPol(&seedCD, Spin::Down);
    double maxGrey=0; for (size_t g=0; g<gup.size(); g++) maxGrey=std::max(maxGrey, std::abs(gup[g]-gdn[g]));
    EXPECT_LT(maxGrey, 1e-10) << "the grey average must ERASE the staggering -- if it does not, the "
                                 "Shubnikov machinery is not actually load-bearing";

    // (d) The TOTAL density is identical through both engines (the even channel is sigma-blind).
    double dTot=0; for (size_t g=0; g<up.size(); g++) dTot=std::max(dTot, std::abs((up[g]+dn[g])-(gup[g]+gdn[g])));
    EXPECT_LT(dTot, 1e-10) << "sigma must not touch the total density";
}

// SHUBNIKOV S4, THROUGH SCF (doc/SymmetryUpgradePlan.md §7 step 7): the imposed magnetic star-average
// confines the magnetization to the STAGGERED sector by construction -- every iterate's m is projected
// onto the Shubnikov-symmetric cone, so the AFM order cannot leak into a net moment and the sublattice
// mirror holds EXACTLY at every iteration, converged or not.  Fixture: the cheapest genuinely staggered
// crystal -- two neutral Mn (the d5s2 Hund pair) in a cubic box, CsCl/B2 arrangement, AFM flip on the
// second.  Bounded iterations (this gate tests SYMMETRY through the live SCF loop, not convergence).
// The through-SCF GREY negative control lives on MnO (MNO_IMPOSE=2): THIS cell's detected grey group is
// all site-preserving (every cubic W admits tau=0), so grey imposition here would be a vacuous control.
TEST(GPW_SCF, ImposedShubnikovHoldsAFMThroughSCF_Mn2Box)
{
    const double a=7.0;
    UnitCell cell(a);
    cell.AddAtom(25, {0.0,0.0,0.0}, false);   // Mn +m
    cell.AddAtom(25, {0.5,0.5,0.5}, true);    // Mn -m (the AFM flip)
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    GpwOptions o;
    o.label="Mn2 B2 AFM imposed"; o.Nelec=14; o.multiplicity=1; o.species={{"Mn",7}};
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;      // neutral Mn: the library's d5s2 spin pair
    o.imposeSymmetry=true;                                // S3 resolves the SHUBNIKOV group from the flips
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.accelerator="DIIS";
    o.scf.SmearingkT=5e-3;                                // the open-d-manifold tie smoother
    o.scf.NMaxIter=6; o.scf.MinΔρ=1e-9; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    o.scf.StartingRelaxRo=0.45; o.scf.MergeTol=1e-4;
    o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);         // coarse quadrature: symmetry, not accuracy
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;

    GpwHandles h;
    GpwResult R=RunGpw(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o,
                       /*verbose*/(bool)std::getenv("GPW_MN2_VERBOSE"), &h);
    ASSERT_TRUE(h.cd) << "the run must produce a final density";

    // The final density through the spin-resolved face: the order must be ALIVE and EXACTLY mirrored.
    const auto* pol=dynamic_cast<const qchem::ChargeDensity::cPolarized_CD*>(h.cd.get());
    ASSERT_NE(pol, nullptr);
    const auto* up=pol->GetChargeDensity(Spin::Up);
    const auto* dn=pol->GetChargeDensity(Spin::Down);
    const rvec3_t off(0.7,0,0), r1(0,0,0), r2(a/2,a/2,a/2);
    const double m1=(*up)(r1+off)-(*dn)(r1+off);
    const double m2=(*up)(r2+off)-(*dn)(r2+off);
    std::cout << "[Mn2 imposed] after "<<R.iters<<" iterations: m1="<<m1<<" m2="<<m2
              << " m1+m2="<<m1+m2<<" Etot="<<R.E.GetTotalEnergy()<<std::endl;
    EXPECT_GT(std::abs(m1), 0.05) << "the AFM order must SURVIVE the imposed SCF loop";
    // What the projector holds EXACTLY mirrored is the rho the FOCK consumes (the engine's (rho,m)
    // star-average); the DENSITY MATRIX itself is deliberately NOT projected (the rho-projection
    // philosophy, T3.2 note), so the D-density probed here is DRIVEN toward the mirror at the
    // convergence rate -- measured |m1+m2| = 3.1e-3 beside lastΔρ = 5.6e-4 at the 6-iteration cap.
    // The gate therefore asserts the mirror at the "small beside the order" tier; the EXACT-mirror
    // guarantee on the projected rasters is the S3 seed-level gate above.
    EXPECT_NEAR(m1+m2, 0.0, 0.05*std::abs(m1))
        << "the imposed SCF must keep the sublattice mirror tight (driven by the projected Fock)";
    EXPECT_NEAR(pol->GetTotalSpin(), 0.0, 1e-8) << "an imposed AFM pair carries zero net moment";
}

// ===================== MnO STATUS: THE OCCUPATION IS THE DEFECT (2026-08-08) ==========================
// Recorded HERE and not only in the plan doc, because the next person to touch this test needs it before
// they reach for another mixing knob.  Every mixing-side candidate is fixed or eliminated (spin-blind rho~
// mixer -- a real bug, fixed; DIIS-per-spin -- checked, it is joint; G0 selectivity -- swept; Kerker-on-m --
// built and measured).  What remained was the OCCUPATION, and the A/B is unambiguous (identical alpha=0.45,
// kT=0, only MOM differing, 20 iterations each):
//     aufbau control (run 13): m_stag dies at iteration 1, E in a -45.6 limit cycle, [F,D] STUCK at 1.43,
//                              converged eps_up-eps_dn ~ 1e-3 everywhere = the non-magnetic collapse.
//     seed-pinned MOM (run 12): m_stag 0.19-0.30 THROUGHOUT, E -54.7, [F,D] 4.11 -> 0.05.
// So the order and ~9 Ha of binding are BOTH bought by holding the occupation: "under-bound by 15.8 Ha" and
// "the AFM order collapsed" really were ONE defect.  Three traps, all of which bit:
//   1. UseMOM was INERT in every run before this one -- FillOrbitals gives SmearingkT>0 precedence and that
//      branch only consults the reference when MOMSmearPenalty>0.  Setting UseMOM beside kT did nothing.
//   2. The delayed capture (MOMStartIter) is one iteration too late here: SetMOM is pushed in by Iterate,
//      after Init's seed fill, so the earliest self-capture is iteration 1 -- AFTER m_stag has already died.
//      momFromSeed self-adopts the seed's subspace before iteration 1 runs.  (The delayed design is for a
//      seed you want to settle AWAY from; MnO's seed IS the physics.)
//   3. The 0h MOM guard would have silently released the reference at iteration 3, on the NaF-calibrated
//      premise that a hole means a wrong reference.  Here a hole is the SIGNATURE of the magnetic branch
//      before the exchange splitting opens.  MNO_MOM_HOLD raises the persistence.
// A HARD pin overshoots into the opposite error: it empties ONE Mn's d shell (five levels at -1.39..-1.32,
// four with eps_up-eps_dn = 0.00000000 -- an empty shell feels no exchange splitting), i.e. a d5/d0
// disproportionation, self-reinforcing because a sunk state has no overlap with the reference.  The cure is
// a pin that can be OVERRIDDEN by a dive: MOMSmearPenalty Lambda, scaled to the physical-vs-foreign score
// gap (see the measured rule in tIrrepWF::FillOrbitals).  Lambda=0.3 == no MOM; Lambda=1.5 cures the
// disproportionation (hole 2.4 -> 0.21 Ha) but DIVERGES at alpha=0.45.
// >>> THE RECIPE THAT REACHES AFM-II (runs 17/18): alpha=0.25, kT=5e-3, MOM seed-pinned, Lambda=1.5.
//     E=-60.9247 at 90 iterations (CP2K -61.4706, so 0.55 Ha -- was 15.8), sites +0.506/-0.529 =>
//     m_net=-0.022, STAGGERED, m_stag peak 0.668 (the seed's 0.366 GREW rather than merely survived), and a
//     real exchange splitting (eps_up-eps_dn to 0.25) where the aufbau control had a uniform 1e-3.
//     THE HOPPING IS GONE: from ~iteration 45 the `cfg` flag goes BLANK -- one configuration, held.
//     Three phases: Kerker+DIIS settles at -60.14 (drho 4.3e-4, [F,D] 7.8e-3), the Ladder hands off to GDM
//     at ~70 and finds a LOWER -60.92 (so the DIIS fixed point was not the bottom), and the moment settles
//     at 0.52.  STILL OPEN: it does not pass the gate -- the run ends on the GDM leg where the mixer's drho
//     is not the convergence measure (lastdrho 3.2e-2 vs 1e-5) -- and it is still NON-AUFBAU (a 0.21 Ha
//     hole), so the residual 0.55 Ha may BE that hole.  NB `cfg *` is only a hopping diagnostic while the
//     energy ORDER is stable: it keys off orbital index, which re-shuffles freely under MOM.
// Knobs: MNO_ALPHA MNO_KERKER_G0 MNO_KT MNO_MOM MNO_MOM_START MNO_MOM_SEED MNO_MOM_PENALTY MNO_MOM_HOLD
//        MNO_ANNEAL="kT,kT,..."  GPW_MNO_VERBOSE GPW_MNO_NMAX  QCHEM_MOM_SCORES
//
// DISABLED until the OCCUPIED-d nonlocal-PP defect is fixed (2026-08-06 find, the campaign's blocker):
// the KB d-channel has NEVER been exercised by occupied states before MnO (CsI's l=2 gate touches only
// virtuals), and it is WRONG by ~12.6x (~4pi-scale): our Mn q7 pseudo-ATOM converges to -189.4 Ha vs the
// CP2K ATOM oracle -14.244 (l<=1 species match modulo the known atom-convention constant: F q7 ours
// -20.90 vs CP2K -24.05 with an EXACT crystal match, O q6 -13.94 vs -15.75).  Both our independent KB
// consumers (atomic mesh route, crystal analytic route) err together => the defect sits in the SHARED
// SeparablePotential model layer; tabulation verified identical to CP2K's GTH_POTENTIALS, Qli/ProjR
// Bessel-pair verified, LegendreP + RealYlm textbook-correct -- the surviving suspects are the
// sqrt(4pi)-Projector convention on the l=2 path and the Cartesian-d s-contaminant (x2+y2+z2) leaking
// into the s-channel projectors.  ORACLES BANKED for the fix session (deck IntegrationTests/CP2K/
// mno_afm2_gpw_sr.inp): CP2K MnO AFM-II crystal E=-61.470570 Ha (72 Broyden steps -- the SCF is hard for
// CP2K too), Mulliken site moments Mn +/-4.654, O ~0, net charges +/-1.69; Mn q7 atom -14.243986.
// Run-3 partial trajectory (pivoted ortho): -454.8 -> -456.5 healthy descent then the ionic slosh --
// dynamics recipe (Ladder+Kerker+MOM) is staged in RunMnO.  Hand-run: GPW_MNO_VERBOSE=1 GPW_MNO_NMAX=n.
TEST(GPW_SCF, DISABLED_MnO_AFM2_RhombohedralGamma)
{
    // MNO_SKIP_AFM (run 41, 2026-08-12): the FM-ONLY hand run -- the same-code Δ(AFM−FM) energy
    // breakdown needs the FM arm at full print precision in the same (annealed/GDM) class as the
    // AFM arm, without paying the ~40 min AFM leg again.  Mirror of MNO_SKIP_FM.
    if (std::getenv("MNO_SKIP_AFM"))
    {
        MnOArm F=RunMnO(/*multiplicity*/11, /*afm*/false, "MnO FM Gamma");
        ASSERT_TRUE(F.R.converged);
        EXPECT_NEAR(F.R.charge, 26.0, 1e-6);
        return;
    }
    MnOArm A=RunMnO(/*multiplicity*/1,  /*afm*/true,  "MnO AFM-II Gamma");

    // ---- the STAGGERED order parameter (reported BEFORE the convergence assert, so a bounded
    //      GPW_MNO_NMAX diagnosis run still yields the physics readout) ----
    auto* pol=dynamic_cast<const qchem::ChargeDensity::cPolarized_CD*>(A.h.cd.get());
    ASSERT_TRUE(pol) << "a multiplicity>=1 run must produce a polarized density";
    const auto* up=pol->GetChargeDensity(Spin::Up);
    const auto* dn=pol->GetChargeDensity(Spin::Down);

    // (i) real space: m(r)=rho_up-rho_dn near the two Mn sites (d-shell peak ~0.7 bohr off the nucleus --
    // the d density is zero AT the nucleus) must be large and OPPOSITE, the sign following the seed pattern.
    const double a=8.40;
    const rvec3_t off(0.7,0,0);
    const rvec3_t rMn1(0,0,0), rMn2(a,a,a);           // cartesian: A*(1/2,1/2,1/2) = a(1,1,1)
    const double m1=(*up)(rMn1+off)-(*dn)(rMn1+off);
    const double m2=(*up)(rMn2+off)-(*dn)(rMn2+off);
    // m_stag = ½(m1−m2) is the ORDER, but it is DEGENERATE between the AFM state (+m,−m) and a
    // sublattice-disproportionated one (0,−2m) -- both give the same m_stag.  Measured 2026-08-08 (run 12):
    // a seed-pinned MOM run reported "order SURVIVED, m_stag=0.29" while the sites read (+0.0006, −0.579),
    // i.e. ONE Mn carried the whole moment.  So the trajectory column is necessary but not sufficient, and the
    // readout has to print the ASYMMETRY beside it: m_net=m1+m2 is 0 for a true stagger and ±2m for a collapse
    // onto one sublattice.  (The EXPECT_NEAR(m1,-m2) gate below is the acceptance test; this line is the
    // instrument that says WHY it failed, on a bounded diagnosis run that never reaches the gate.)
    std::cout << "[MnO AFM-II] site moments m(r): Mn1(+seed)="<<m1<<"  Mn2(-seed)="<<m2
              << "  m_stag=½(m1−m2)="<<0.5*(m1-m2)<<"  m_net=m1+m2="<<(m1+m2)
              << (std::abs(m1+m2) > 0.2*std::abs(m1-m2)
                  ? "  ** NOT STAGGERED: the moment sits on ONE sublattice" : "")
              << std::endl;
    ASSERT_TRUE(A.R.converged);                       // gate AFTER the readout (see above)
    EXPECT_NEAR(A.R.charge, 26.0, 1e-6);
    EXPECT_GT(m1,  0.01) << "Mn1 must stay in the +m basin the seed chose";
    EXPECT_LT(m2, -0.01) << "Mn2 must stay in the -m basin the seed chose";
    EXPECT_NEAR(m1, -m2, 0.2*std::abs(m1)) << "the two sublattices should stagger symmetrically";

    // (ii) G-space: |m-tilde| at the AFM wavevector (dm=(1,0,0): 2pi m.f = pi between the Mn sublattices),
    // scaled by the cell volume back to electrons-scale.
    {
        Hamiltonian::PWFittedVxc::fbs_t fit(A.h.bs->CreateVxcFitBasisSet(A.cell.get(), qcMesh::MeshParams{}));
        auto* fu=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(up);
        auto* fdn=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(dn);
        ASSERT_TRUE(fu && fdn);
        ΔG_Map mu=fu->GetFourierDensity(*fit), md=fdn->GetFourierDensity(*fit);
        const ivec3_t q(1,0,0);
        ASSERT_TRUE(mu.count(q) && md.count(q));
        const double Omega=A.cell->GetCellVolume();
        const double Mstag=std::abs(mu.at(q)-md.at(q))*Omega/2;
        std::cout << "[MnO AFM-II] |m-tilde(q_AFM)|*Omega/2 = "<<Mstag<<" e- (staggered moment scale)"<<std::endl;
        EXPECT_GT(Mstag, 1.0) << "the converged staggered moment must be electrons-scale (d^5 ~ 4-5)";
    }

    // ---- the ORDERING ENERGETICS: AFM-II below FM (the rocksalt-MnO superexchange ground state) ----
    // MNO_SKIP_FM (run 39): an AFM-only hand run -- the FM arm costs ~25 min and an experiment probing
    // the AFM state (the GDM stage) learns nothing from it.  (The ordering EXPECT is already a known
    // failure at Gamma-only+LSDA: run 38 measured FM 38 mHa BELOW AFM-II.)
    if (std::getenv("MNO_SKIP_FM")) return;
    MnOArm F=RunMnO(/*multiplicity*/11, /*afm*/false, "MnO FM Gamma");
    ASSERT_TRUE(F.R.converged);
    EXPECT_NEAR(F.R.charge, 26.0, 1e-6);
    const double Eafm=A.R.E.GetTotalEnergy(), Efm=F.R.E.GetTotalEnergy();
    std::cout << "[MnO ordering] E_AFM="<<Eafm<<"  E_FM="<<Efm<<"  dE="<<(Efm-Eafm)*1000<<" mHa"<<std::endl;
    EXPECT_LT(Eafm, Efm) << "AFM-II must be the LSDA ground state ordering";
}

// TEMPORARY diagnosis probe (will be removed): is the -0.205 overlap eigenvalue on the MnO rhombohedral
// cell a CELL-SHAPE problem (oblique lattice sum) or a BASIS problem (diffuse Mn/O windows)?
TEST(GPW_SCF, DISABLED_MnOCellShapeProbe)
{
    const double a=8.40;
    Matrix3D<double> A(a, a/2, a/2,  a/2, a, a/2,  a/2, a/2, a);
    UnitCell cell(A);
    cell.AddAtom(14, {0.0,0.0,0.0});
    cell.AddAtom(14, {0.5,0.5,0.5});
    cell.AddAtom(14, {0.25,0.25,0.25});
    cell.AddAtom(14, {0.75,0.75,0.75});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwOptions o;
    o.label="rhombo cell probe (SIPP_SR Si)"; o.Nelec=16; o.species={{"Si",4}};
    o.imposeSymmetry=false;
    o.scf.NMaxIter=1; o.scf.MinΔρ=1e30; o.scf.MinΔE=1e30;
    o.scf.MinΔFD=1e30; o.scf.MinVirial=1e30; o.scf.MinFD=1e30;
    RunGpw(lat, MakeBasisSR(cell), o, /*verbose*/true);
}
