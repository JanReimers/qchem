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
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // tSCFAcceleratorNull<dcmplx> (NaF: pure damped Kerker)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Symmetry.Spin;                      // Spin
import qchem.Symmetry.Factory;                   // BlochFactory (build a k-block with a fractional MP shift)
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

struct GpwResult { bool converged; double charge; qchem::EnergyBreakdown E; size_t iters; };

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
GpwResult RunGPW(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, double densityEcut,
                 int Nelec, const char* element, const char* label, bool verbose=false, int nmax=120,
                 qchem::Ortho ortho=qchem::Cholesky, double orthoTol=0.0,
                 rvec3_t kShift={0,0,0}, double minDrho=1e-6, double minDE=1e30,
                 qchem::ChargeDensity::SeedStrategy seed=qchem::ChargeDensity::SeedStrategy::Uniform,
                 BasisSet::Lattice_3D::CellImages images=BasisSet::Lattice_3D::CellImages::Periodic)
{
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, std::move(mol), densityEcut, kShift, images));
    auto       irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry the Sum_k)
    Crystal_EC ec(irreps, Nelec);                  // multi-k ready; a single Gamma block is the length-1 case
    // Ham_PW_DFT drives GPW verbatim (kinetic + external-PP + Hartree + Dirac X + VWN5 + ion-ion Ewald).
    qchem::Hamiltonian::cHamiltonian* ham=new qchem::Hamiltonian::Ham_PW_DFT(
        lat.GetStructure(), bs.get(), element, "LDA", 4);
    auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(qchem::SCFAccelerators::DIISParams{8, 8.0, 1e-10, 1e-9});
    // \a ortho / \a orthoTol: Cholesky (default) needs S positive-definite; for a periodic lattice basis the
    // diffuse Gaussians can go linearly dependent -> Eigen/SVD with a small-eigenvalue cutoff.
    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc,
                                         seed, lat.GetStructure().get(),
                                         ortho, orthoTol);
    std::vector<FpRow> series;   // capture the SCF trajectory for the dynamics fingerprint (Observer telemetry)
    scf.SetObserver([&series](const qchem::SCFIterator::SCFProgress& p)
                    { series.push_back({p.iteration, p.energy, p.dE, p.commutator, p.drho}); });
    SCFParams par;
    // \a minDrho / \a minDE: the SCF convergence gate (AND of the active criteria).  Default = the historical
    // Δρ<1e-6 only (energy/virial/[F,D] off) -- the committed Rcut=0 gapped anchors converge there.  A fitted
    // GPW SCF at a converged Rcut cannot drive Δρ below the density-fit floor (~2e-4, NON-variational -- see
    // doc/GPWPlan.md; NOT a bug, and NOT fixable by a direct minimiser, which needs a variational energy), so
    // the physical gate there is ENERGY: pass minDE~1e-6 + a relaxed minDrho~1e-3 to stop once the total energy
    // settles (~iter 15) instead of chasing fit noise to nmax.
    par.NMaxIter=nmax; par.MinΔρ=minDrho; par.MinΔE=minDE; par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=verbose;
    qchem::Hamiltonian::ReportGridCharge()=verbose;   // per-iteration integral-rho_grid vs Tr(DS) (grid-truncation charge loss)
    scf.Iterate(par);
    qchem::Hamiltonian::ReportGridCharge()=false;      // process-wide flag -- reset so it does not leak to other tests
    Fingerprint(series, label);                        // classify the trajectory: converged / oscillating / stalled / diverging

    auto* cd=scf.GetWaveFunction()->GetChargeDensity();
    double charge=cd->GetTotalCharge();
    delete cd;
    qchem::EnergyBreakdown E=scf.GetEnergy();
    std::cout << "["<<label<<"] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Eelec="<<E.GetElectronicEnergy() << " Etot="<<E.GetTotalEnergy()
              << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc
              << " Enn="<<E.Enn<<" Ealign="<<E.Ealign<<")" << std::endl;
    return {scf.Converged(), charge, E, scf.GetIterationCount()};
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
TEST(GPW_SCF, SiliconMultiKPlumbing)
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
TEST(GPW_SCF, SR_2x2x2ShiftedMP_vs_CP2K)
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

// (4) MULTI-SPECIES GPW: ionic NaF (rocksalt = FCC + 2-atom basis) at Gamma, driven by the multi-species
// Ham_PW_DFT ctor ({{"Na",1},{"F",7}}) on the GENERATED valence_lowq Gaussian basis (Na 5s2p + F 8s6p).  The
// plane-wave sibling (PlaneWaveDFT.FrameworkNaFThroughSCFIterator, Ecut=6) lands -20.3293; GPW here is the
// SAME PP + functional on a Gaussian basis.  FIRST-LIGHT machinery check: does multi-species GPW converge and
// conserve the 8 valence e- (1 Na + 7 F)?  F's tight 40-a.u. exponent wants a fine density grid (densityEcut
// high -- F is the hard atom), so the total is grid-underconverged at a modest cutoff; charge + convergence
// are the anchors here (a did-E-move total once the cutoff is dialled in).
// PROGRESSION (all charge=8, Gamma, densityEcut=40):
//   Rcut=0, full basis  -> Etot=-25.086 (first light; Rcut=0 over-binds, no inter-cell screening)
//   Rcut=2a, SR basis   -> Etot=-23.556 (this test; the over-binding partly removed by the periodic images)
//   PW Gamma reference  -> -20.3293
// The SR basis is REQUIRED at Rcut=2a: the full basis' truncated Bloch overlap is INDEFINITE there
// (DISABLED_NaFOverlapConditioningSweep), SR is PSD (min eig +7.5e-4, reported at startup).  The residual
// -23.56 vs -20.33 is grid (densityEcut=40 is coarse for F's 40-a.u. exponent) + SR basis incompleteness +
// the SCF hit the 60-iter cap (not fully gate-converged).  DISABLED: long.  Robust anchor
// = charge conservation; TODO to close the gap: converge densityEcut, more iters, then CP2K.
//
// 2026-07-15 (analytic path, auto Ecut=160): the SAME-BASIS CP2K oracle is **Etot = -27.93128**
// (UnitTests/CP2K/naf_gpw_sr_diag.inp, doc/CP2Kresults.md) -- CP2K's ENERGY settles to 1e-6 by ~130
// damped-Broyden/diagonalization iterations while its DENSITY limit-cycles forever (RMS 0.03-0.12), the
// exact charge-transfer cycle we see: the disease is the system+basis, NOT either implementation (CP2K's
// OT run never settled E at all, swinging -25.7..+253).  OUR Kerker(G0=1)+DIIS at relax 0.3 does NOT
// settle E in 60 iters (iteration 60 lands essentially randomly: -24.03 and +887.55 across runs; charge
// stays 8.0000000000 exactly throughout) -- the DIIS mid-cycle Fock extrapolations ARE the +900 spikes
// (they land exactly on the Nproj=8 iterations).  The recipe + tuning outcome (2026-07-16) is documented
// at the mixing knobs inside the test body: converging at Ecut=40, production grid blocked on
// quasi-Newton (Pulay/Broyden) density mixing.
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

    namespace L3=BasisSet::Lattice_3D;
    // densityEcut<0 => AUTOMATIC: the grid is floored to 4*alpha_max from the basis (F alpha_max=40 -> 160),
    // so F's tight 40-a.u. exponent is resolved without the caller specifying the Hartree value.
    // Mixing TUNING knobs (env overrides; defaults = the committed recipe).  The charge-transfer limit cycle
    // is a mixing-stability problem, so alpha/G0 sweeps are how this test gets tuned -- e.g.
    //   NAF_ALPHA=0.1 NAF_KERKER_G0=2 NAF_ECUT=40 NAF_NMAX=100 ./UTMain --gtest_filter=*NaFRocksalt* ...
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };
    // TUNED (2026-07-15/16, Ecut=40 alpha scan): alpha=0.2/0.1 limit-cycle (+-75 Ha, period ~48, passing
    // THROUGH the fixed point); alpha=0.05 contained (+-1 Ha) but not decaying; **alpha=0.025 converges**
    // (Ecut=40 -> ~-27.75, residual +-0.04 wobble).  G0: 1.0 good; 2.5 WORSE (over-damping low-G
    // destabilizes), 0.5 worse.
    //
    // THE PRODUCTION (AUTO Ecut=160) GRID IS BLOCKED ON QUASI-NEWTON DENSITY MIXING: the fine grid grows a
    // second, UNPHYSICAL attractor (E ~ -39, Exc ~ -143, integral rho_grid swinging 5.1<->7.7 vs Tr(DS)=8 --
    // the mid-slosh D loads the sharpest F pairs beyond the grid calibration, the XC of that spiky/negative
    // rho feeds back) which captures plain damped Kerker at EVERY alpha tried (0.2 down to 0.01: the
    // trajectory passes -26 and slides in).  CP2K's SAME map converges under BROYDEN (quasi-Newton, 8-step
    // history) -- the ingredient our complex path lacks.  => next increment: complex Pulay/Broyden
    // rho-tilde mixing on the FourierMixCD infrastructure; until then this test pins the CONVERGING
    // Ecut=40 regime (a real mixing-regression anchor: either bad attractor lands ~+65 or ~-39, far
    // outside the gate).  NAF_ECUT=-1 runs the production grid for Pulay development.
    const double tuneEcut =envd("NAF_ECUT", 40.0);           // <0 AUTO (=160 for F); 40 = the convergent regime
    const double tuneAlpha=envd("NAF_ALPHA", 0.025);
    const double tuneG0   =envd("NAF_KERKER_G0", 1.0);
    const size_t tuneNMax =size_t(envd("NAF_NMAX", 200));
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, mol, tuneEcut));   // PERIODIC: eps-complete enumeration
                                                                          // (the old Rcut=2a truncated S -- the
                                                                          // measured -2.25 e scheme mismatch)
    auto       irreps=bs->GetIrreps(Spin::None);
    Crystal_EC ec(irreps, 8);
    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), {{"Na",1},{"F",7}}, "LDA");
    // THE CP2K RECIPE (2026-07-15, doc/CP2Kresults.md): pure DAMPED KERKER MIXING, NO DIIS, exit on E-FLAT.
    //  - CP2K's damped-Broyden/diagonalization run settles E to 1e-6 by ~130 iterations at alpha=0.2 while its
    //    DENSITY limit-cycles forever (RMS 0.03-0.12) -- so the density is NOT a convergence criterion here
    //    (MinΔρ=1e30) and the physical exit is the relative-energy gate MinΔE (the SCFParams field added for
    //    exactly this non-variational settled-E case).
    //  - NO DIIS: our +431/+910 Ha energy spikes landed EXACTLY on the iterations where DIIS re-engaged
    //    (Nproj=8 in the 2026-07-15 traces) -- mid-cycle Fock extrapolation is the spike generator.  Earlier
    //    finding (2026-07-12) stands: pure relax mixing limit-cycles too; Kerker damping is what contains the
    //    low-G charge-transfer slosh, and heavier damping (alpha=0.2, CP2K's value) settles E.
    //    Kerker holds alpha=StartingRelaxRo fixed (no [F,D] auto-tune).
    // NAF_DIIS=1: re-admit DIIS (quasi-Newton in Fock space).  Its mid-cycle extrapolations were banned
    // for spiking -- but that was measured on the CORRUPTED (Rcut=2a scheme-mismatch) map; on the honest
    // map the alpha-independent giant-response instability is exactly DIIS's job.  Env knob for the sweep.
    auto envDIIS=[&]{ const char* d=std::getenv("NAF_DIIS"); return d && std::atof(d)!=0.0; };
    qchem::SCFAccelerators::tSCFAccelerator<dcmplx>* acc = envDIIS()
        ? static_cast<qchem::SCFAccelerators::tSCFAccelerator<dcmplx>*>(
              new qchem::SCFAccelerators::cSCFAcceleratorDIIS(qchem::SCFAccelerators::DIISParams{8, 8.0, 1e-10, 1e-9}))
        : new qchem::SCFAccelerators::tSCFAcceleratorNull<dcmplx>();
    qchem::ReportOverlapConditioning()=true;   // report min eig(S)/min sv(S) at SetBasisOverlap (the ctor below)
    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::IonicSAD, lat.GetStructure().get(),
                                         qchem::Cholesky, 0.0);   // diffuse F- / Na+ ionic seed (halved PW iters)
    qchem::ReportOverlapConditioning()=false;  // process-wide flag -- reset so it does not leak to other tests
    SCFParams par; par.NMaxIter=tuneNMax; par.MinΔρ=1e30; par.MinΔE=1e-8; par.MinΔFD=1e30; par.MinVirial=1e30;
    par.MinFD=1e30; par.StartingRelaxRo=tuneAlpha; par.MergeTol=1e-4; par.Verbose=true;
    par.KerkerG0=tuneG0;   // Kerker screening wavevector (a.u.^-1): damp the low-G charge-transfer slosh
    // MOM FIX (doc/GPWPlan 0b'' -- the diagnosed cure): occupy by MAX OVERLAP onto a FIXED reference {F 2s,
    // F 2p} subspace so aufbau cannot capture the diving diffuse virtual (the measured occupation swap
    // F 2p 6->4 e).  Delayed IMOM: plain aufbau for MOMStartIter fills (descend to the physical fixed
    // point), then capture the reference ONCE and hold it.  NAF_MOM=0 A/Bs it off; NAF_MOM_START tunes.
    auto envMOM=[&]{ const char* d=std::getenv("NAF_MOM"); return !d || std::atof(d)!=0.0; };  // default ON
    par.UseMOM=envMOM();
    par.MOMStartIter=(int)envd("NAF_MOM_START", 10);
    qchem::Hamiltonian::ReportGridCharge()=true;   // F's tight 40-a.u. exponent: watch integral rho_grid vs Tr(DS)
    qchem::SCFIterator::ReportBandGap()=true;       // BAND-GAP INSTRUMENT (doc/GPWPlan 0b''): watch eps_HOMO/eps_LUMO/gap
                                                    // per iteration -- the near-degenerate-frontier (giant-response)
                                                    // hypothesis predicts the gap collapsing near the fixed point,
                                                    // with the spurious level diving across the Fermi edge before each spike.
    scf.Iterate(par);
    qchem::Hamiltonian::ReportGridCharge()=false;  // process-wide flag -- reset so it does not leak to other tests
    qchem::SCFIterator::ReportBandGap()=false;     // idem (MOM is per-run via SCFParams -- nothing to reset)

    auto* cd=scf.GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge(); delete cd;
    auto E=scf.GetEnergy();
    std::cout << "[NaF GPW Gamma] iters="<<scf.GetIterationCount()<<" charge="<<charge<<" Etot="<<E.GetTotalEnergy()
              << " (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc
              << " Enn="<<E.Enn<<" Ealign="<<E.Ealign<<")" << std::endl;
    EXPECT_NEAR(charge, 8.0, 1e-6);     // 1 (Na) + 7 (F) valence electrons, conserved
    // ENERGY PIN SUSPENDED (banish-Rcut + SR2, 2026-07-16).  The old anchor (-27.73 +- 0.05) was
    // measured on the CORRUPTED map (Rcut=2a truncated S vs complete collocation) and is void.  On the
    // HONEST map: the scheme mismatch is DEAD (iteration-1 grid charge -4.9 e -> -2.4e-6 e) and the
    // recipe descends SMOOTHLY to a genuine fixed point (SR ~ -28.00, SR2 ~ -27.73), BUT departs and
    // blows up (+5e3 Ha) each time it nears it -- period ~27.  MEASURED CLASSIFICATION (doc/GPWPlan 0b'):
    //   - NOT conditioning: SR2 (lambda_min 1.6e-3, cond 2715) shows the SAME spikes as SR (1.03e-6);
    //   - NOT plain linear-mixing gain: alpha-INDEPENDENT (10/10/13 spikes at 0.025/0.0125/0.00625);
    //   - NOT fixed by DIIS (NAF_DIIS=1): 51 excursions, En>EMax flapping;
    //   - the departure is a SMOOTH climb over ~5 iters (a growing mode, not an occupation swap).
    // BAND-GAP INSTRUMENT VERDICT (2026-07-17, ReportBandGap; Ecut=40/alpha=0.025).  The static
    // near-degeneracy version of the hypothesis is FALSE and REFINED to a TRANSIENT giant-response mode:
    //   - the FIXED-POINT gap is HEALTHY, eps_LUMO-eps_HOMO ~ 0.33-0.37 Ha (~9-10 eV): NaF/Gamma in
    //     this ionic basis IS a wide-gap insulator at convergence (iters 30-37, 55-66 plateau at gap~0.35);
    //   - each spike is preceded ONE iteration earlier by eps_LUMO DIVING ~0.2-0.5 Ha (a diffuse virtual
    //     with a giant response to the low-G charge-transfer slosh): gap collapses (iter 12 -> 2.8e-2 with
    //     eps_LUMO crashing +0.167 -> -0.077; iter 68 -> 1.2e-4, eps_H/eps_L DEGENERATE) as the virtual
    //     crosses the occupied manifold, then AUFBAU fractionally occupies it ([partial-occ HOMO] fires
    //     exactly on the spike iters 14, 41) -> a ~1/sqrt(lambda) diffuse vector enters D -> E=+5e3..+7e3,
    //     [F,D] 0.09 -> 130.  Deterministic period ~27 (spikes 14/41/68 in one run).
    // => the mechanism is a giant-response DIFFUSE VIRTUAL causing a periodic aufbau LEVEL-CROSSING, not a
    //    small static gap.  Fix selection (doc/GPWPlan 0b''): NOT Fermi smearing (the fixed-point gap is
    //    large) -- the direct guard is MOM (occupied-subspace continuity through the crossing; coded,
    //    inactive: tIrrepWF::MOMScores/CaptureMOMReference), and the CAUSE-side cure is 0c Pulay/Broyden
    //    rho-mixing to damp the slosh that drives the dive (matches CP2K converging THIS map with Broyden).
    //
    // FRONTIER-WINDOW REFINEMENT (same instrument, 2-occ/4-virt window per iteration).  Two sharper facts:
    //   - it is ONE ISOLATED hyper-responsive virtual, NOT a wide-band cluster: at the dive (iter 11->12)
    //     the LUMO crashes +0.167 -> -0.077 (0.24 Ha in ONE step) while its virtual NEIGHBOURS (+0.42,
    //     +0.79) barely move, and the LUMO sits ~0.25 Ha clear of the next virtual at the plateau -- so
    //     the giant response is a single diffuse (Na-3s-like) conduction state overlapping the charge-
    //     transfer region, NOT an over-complete diffuse-band cluster (argues 3b physical, not 3a ghost);
    //   - the spike IS an OCCUPATION SWAP (this CORRECTS the 0b' "growing mode, not a swap" note): at each
    //     spike the F 2p level drops from (6.0) to (4.0) electrons -- the diving virtual captures 2 e out
    //     of the F 2p manifold (iters 14 and 41).  The smooth dive (the "growing mode") TERMINATES in the
    //     aufbau swap; they are two phases of ONE event.  MOM (pin the {F 2s, F 2p} occupied subspace) is
    //     therefore the direct fix; it should be clean since it is an isolated single-state swap.
    //
    // MOM FIX WIRED UP + VALIDATED (2026-07-17, SCFParams::UseMOM; default ON here).  DELAYED IMOM
    // (SCFParams::MOMStartIter=10): run plain aufbau for ~10 fills so the SCF descends to the
    // physical fixed point, then CAPTURE the {F 2s, F 2p} occupied subspace ONCE and hold it FIXED.  Two
    // wrong variants were measured + rejected first: RUNNING MOM (re-capture every iter) DRIFTS (a spike
    // corrupts the reference -> locks a +0.74 Ha level occupied while a -50 Ha level is empty -> wrong
    // -24.4); IMOM-from-iter-0 anchors the raw seed (mid-transient) -> catastrophe (+5 Ha occupied,
    // -112 Ha empty).  Delayed IMOM WORKS: the occupation swaps VANISH (partial-occ count 0), the diving
    // virtual is banished (-45 Ha, UNOCCUPIED), and the SCF converges SMOOTHLY+MONOTONICALLY to the
    // physical fixed point -27.76 (Ecut=40, dRho 6e-4 at 150 iters; gap 0.50 Ha).  vs CP2K oracle
    // -27.93128 at 320 Ry -- the ~0.17 Ha is the Ecut=40 grid (the production grid + 0c mixing are the
    // remaining levers).  ONE residual spike survives (iter 19) but it is NOT an occupation swap
    // (partial-occ 0) -- a density-MIXING excursion (the charge-transfer slosh) => the 0c Pulay/Broyden
    // item, not MOM's job.  NAF_MOM=0 A/Bs it off (restores the eternal spikes); NAF_MOM_START tunes.
    EXPECT_NEAR(E.GetTotalEnergy(), -27.76, 0.1);   // MOM now CONVERGES the map: a real did-E-move anchor
                                                    // (loose: slow linear-Kerker tail + the iter-19 mixing
                                                    // transient; 0c Pulay will tighten it)
}

// VALIDATION (2026-07-13): can we drop the SR-basis hand-tuning and use the FULL valence_lowq basis, relying on
// (a) magnitude screening + a generous overlap Rcut=4a to CONVERGE S (removing the Gibbs/truncation artifacts, so
// the residual negative eigenvalues shrink to ~0), then (b) canonical Eigen ortho with a tol=1e-6 that drops the
// intrinsic ~1e-6 over-complete null cluster (which sits in a clean ~1000x spectral gap below the physical ~1e-3
// directions -- see DISABLED_NaFOverlapConditioningSweep)?
//
// RESULT: NO, not yet -- this hits an INTEGRATION WALL, not a convergence problem.  The truncating ortho reduces
// the working dimension 37 -> 33 (drops the 4 null directions), but the GPW/periodic SCF stack (Crystal_EC /
// cDM_CD / the density collocation) assumes the FULL basis dimension n=37, so the rank-33 ortho subspace collides
// with it: "C++ exception: Matrix sizes do not match", thrown after the ortho truncation, before iter 1.  (The
// MOLECULAR SCF path handles rectangular V -- C=V*U' is 37x33, density full 37x37 -- but the PERIODIC path does
// not.)  So at the OVERLAP level (2) is clean (‖VᴴSV-I‖=6.6e-11, gap x1000) but at the SCF level it is BLOCKED:
// dropping SR needs rank-reduction plumbed through the periodic stack (density in the truncated subspace,
// back-transform, band-count reconciliation) -- its own increment.  Until then SR (dimension-preserving, cleanly
// PD, no truncation) stays the GPW conditioning answer.  Kept DISABLED as the marker for that future work.
// Collocation stays at collRcut=2a (the un-screened orbital-grid axis is a separate, later increment).
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
    qchem::ReportOverlapConditioning()=true;
    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::IonicSAD, lat.GetStructure().get(),
                                         qchem::Eigen, 1e-6);   // (2): canonical ortho, drop the ~0 null cluster
    qchem::ReportOverlapConditioning()=false;
    SCFParams par; par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30;
    par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=true; par.KerkerG0=1.0;
    qchem::Hamiltonian::ReportGridCharge()=true;
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
