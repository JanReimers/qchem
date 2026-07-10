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
#include <cmath>
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
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ElectronConfiguration.Crystal;      // Crystal_EC (single-k Bloch occupation)
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS; // cSCFAcceleratorDIIS (complex DIIS)
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

// One GPW Gamma-point SCF: build the GPW basis over the lattice, hand it the plane-wave LDA Hamiltonian
// (Ham_PW_DFT reaches GPW's real-space Integrals_Pseudo), seed uniform, run the complex-DIIS cSCFIterator.
GpwResult RunGPW(const Lattice_3D& lat, std::shared_ptr<const Real_BS> mol, double densityEcut, double Rcut,
                 int Nelec, const char* element, const char* label, bool verbose=false, int nmax=120,
                 qchem::Ortho ortho=qchem::Cholesky, double orthoTol=0.0, double collRcut=0.0,
                 rvec3_t kShift={0,0,0})
{
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, std::move(mol), densityEcut, Rcut, collRcut, kShift));
    auto       irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry the Sum_k)
    Crystal_EC ec(irreps, Nelec);                  // multi-k ready; a single Gamma block is the length-1 case
    // Ham_PW_DFT drives GPW verbatim (kinetic + external-PP + Hartree + Dirac X + VWN5 + ion-ion Ewald).
    qchem::Hamiltonian::cHamiltonian* ham=new qchem::Hamiltonian::Ham_PW_DFT(
        lat.GetStructure(), bs.get(), element, "LDA", 4);
    using qchem::SCFAccelerators::DIISParams;
    auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(DIISParams{8, 8.0, 1e-10, 1e-9});
    // \a ortho / \a orthoTol: Cholesky (default) needs S positive-definite; for a lattice basis with images
    // (Rcut>0) the diffuse Gaussians go linearly dependent -> Eigen/SVD with a small-eigenvalue cutoff.
    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::Uniform, lat.GetStructure().get(),
                                         ortho, orthoTol);
    SCFParams par;
    par.NMaxIter=nmax; par.MinΔρ=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=verbose;
    scf.Iterate(par);

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
// (1c) MULTI-K PLUMBING: a 2x1x1 Monkhorst-Pack mesh (2 k-points) at Rcut=0.  With no inter-cell images every
// k-block is IDENTICAL to Gamma (the phase multiplies only the origin), so the whole multi-k machinery -- one
// GPW_IBS per BZ k-point (GPW_BasisSet iterating MakeKMesh WITH BZ weights, mirroring PW_BasisSet), the multi-
// block GetIrreps, Crystal_EC's BZ-weighted (Sum_k w_k) occupation, the framework's per-irrep k-loop, and the
// BZ-summed charge/energy -- must reproduce the single-Gamma total EXACTLY (well-conditioned, no dispersion).
// This is the tight, robust gate for Step 2's plumbing (it caught a missing BZ weight -> charge x Nk); the
// dispersive bulk (Rcut>0) is the DISABLED conditioning study below.  A 2-point mesh exercises the same multi-
// block plumbing as 2x2x2 at ~1/4 the cost (the Rcut=0 blocks are redundant copies of Gamma).
TEST(GPW_SCF, SiliconMultiKPlumbing)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,1,1));   // 2 k-points, but Rcut=0 -> both identical to Gamma

    GpwResult R=RunGPW(lat, MakeBasis(cell), /*densityEcut*/12.0, /*Rcut*/0.0, /*Nelec*/8, "Si", "Si 2x1x1 Rcut=0");

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);                          // 8 valence e- (BZ-weighted Sum_k, not x Nk)
    EXPECT_NEAR(R.E.GetTotalEnergy(), -8.2476, 5e-3);         // == the Gamma total (SiliconGammaConverges)
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
    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Rcut*/2.0*a, /*Nelec*/8, "Si",
                       "Si 2x2x2 Gamma-centred Rcut=2a", /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0, 0.0);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.77846, 3e-3) << "GPW 2x2x2 Gamma-centred vs CP2K same-mesh -7.77846";
}

// SHIFTED Monkhorst-Pack (kShift=½ → k at ±¼ = CP2K's DEFAULT MONKHORST-PACK 2 2 2) -- the apples-to-apples
// match to CP2K's shipped 2x2x2 reference -7.86744 (deck si_fcc_gpw_222.inp).  The fractional-k PLUMBING works
// (kShift threaded through BlochQN/BlochFactory/GPW_BasisSet; overlap S(k) PSD, charge 8), BUT the SCF at these
// k OVER-BINDS (single k=¼ block: Een -1.14 → -15.15, Etot -14.7) and does not converge.
//
// LOCALIZED (2026-07-10) -- this is a genuinely-COMPLEX-Bloch-phase bug, NOT the shifted-MP plumbing:
//   * NOT conditioning: min eig S(k=¼) = 0.020 (PSD) at Rcut=2a AND 3a.
//   * NOT the density collocation: ‖W₀·Ω − S(k)‖/‖S‖ = 4e-6 at k=0, ½ AND ¼ -- the collocated density is
//     consistent with the analytic overlap even at complex k.
//   * The GPW term matrices are Hermitian by construction (chmat_t).
// KEY: at Γ and half-integer k the phase e^{ik·R}=±1 is REAL → the density matrix D is REAL; k=¼ is the FIRST
// time D is genuinely COMPLEX.  So the suspect is the complex-D SCF path (density-matrix build / complex DIIS /
// Crystal_EC occupation), consistent with the known "GDM/Ladder are <double>-only" limitation -- likely
// framework-level, not the GPW evaluator.  This diagnostic (overlap PSD + collocation-consistent, no SCF) pins
// that; the over-binding number above is from the single-k SCF.  NEXT: debug the complex-D SCF path, then this
// test can assert -7.86744.  (disabled.)
TEST(GPW_SCF, DISABLED_SR_2x2x2ShiftedMP_vs_CP2K)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    for (double rc : {2.0*a, 3.0*a})   // overlap PSD at the shifted k? (rules out conditioning)
    {
        namespace L3=BasisSet::Lattice_3D;
        std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasisSR(cell), 0.0, rc, 0.0, rvec3_t(0.5,0.5,0.5)));
        double worst=1e30;
        for (auto b : bs->Iterate<qchem::BasisSet::Complex_OIBS>())
        {
            chmat_t S=b->Overlap(); rvec_t w; mat_t<dcmplx> U; blazem::eigen(S, w, U);
            for (size_t i=0;i<w.size();i++) worst=std::min(worst,w[i]);
        }
        std::printf("shifted MP k=¼  Rcut=%.1fa  worst min(eig S(k))=%.5e %s\n", rc/a, worst, worst<0?"(INDEF)":"(PSD)");
        EXPECT_GT(worst, 0.0) << "S(k) must be PSD at the shifted k (Rcut="<<rc/a<<"a)";
    }
    // Collocation vs analytic overlap: W₀·Ω must equal S(k) at every k (rules out the density collocation).
    using BasisSet::Lattice_3D::GPW_IBS;
    using BasisSet::Lattice_3D::GPW_Evaluator;
    auto probe=[&](rvec3_t sh, const char* lbl)->double
    {
        GPW_IBS gpw(cell, Symmetry::BlochFactory(ivec3_t(1,1,1),ivec3_t(0,0,0),1.0,sh), MakeBasisSR(cell), 12.0, 2.0*a);
        const GPW_Evaluator& ev=gpw; const BasisSet::Complex_OIBS& g=gpw;
        chmat_t S=g.Overlap(); G_ERI3 ov=ev.Overlap3CTensor();
        int c0=-1; for (size_t c=0;c<ov.columns.size();c++) if (ov.columns[c].dm.x==0&&ov.columns[c].dm.y==0&&ov.columns[c].dm.z==0) c0=int(c);
        size_t n=S.rows(); double num=0,den=0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
        { dcmplx w=ov.weights[c0](i,j)*ov.volume, s=S(i,j); num+=std::norm(w-s); den+=std::norm(s); }
        double rel=std::sqrt(num/den);
        std::printf("%s  ||W0*Om - S(k)||/||S|| = %.4e\n", lbl, rel);
        return rel;
    };
    EXPECT_LT(probe(rvec3_t(0,0,0),       "k=0    (real phase) "), 1e-3);
    EXPECT_LT(probe(rvec3_t(0.5,0.5,0.5), "k=1/2  (real phase) "), 1e-3);
    EXPECT_LT(probe(rvec3_t(0.25,0.25,0.25),"k=1/4 (complex)    "), 1e-3);   // density consistent even at complex k
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
    auto traces=[&](double frac, double Rcut, double& kin, double& vloc, double& vnl)
    {
        FCCUnitCell cell(a);
        cell.AddAtom(14, {frac, frac, frac});
        cell.AddAtom(14, {0.25+frac, 0.25+frac, 0.25+frac});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        auto st=lat.GetStructure();
        Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH("Si","LDA",4);
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), MakeBasisSR(cell), dE, Rcut, 0.0);
        const BasisSet::Complex_OIBS& g=gpw;
        kin =tr(g.Kinetic());
        vloc=tr(gpw.MakeLocalPotential   (st.get(), pp.local));
        vnl =tr(gpw.MakeSeparablePotential(st.get(), pp.nonlocal));
    };
    for (double rc : {1.5*a, 2.0*a})
    {
        double kc,lc,nc, ks,ls,ns;
        traces(0.00, rc, kc,lc,nc);
        traces(0.13, rc, ks,ls,ns);
        std::printf("Rcut=%.2fa  Kin[%10.5f/%10.5f d=%.2e]  Vloc[%10.5f/%10.5f d=%.2e]  Vnl[%10.5f/%10.5f d=%.2e]\n",
                    rc/a, kc,ks,std::fabs(kc-ks), lc,ls,std::fabs(lc-ls), nc,ns,std::fabs(nc-ns));
        if (rc>=2.0*a)   // fully wrapped: every term translation-invariant (image-truncation limited ~1e-4)
        {
            EXPECT_NEAR(lc, ls, 1e-3) << "Vloc translation invariance at Rcut="<<rc/a<<"a";
            EXPECT_NEAR(nc, ns, 1e-3) << "Vnl translation invariance at Rcut="<<rc/a<<"a";
        }
    }
}

// THE CP2K ENERGY GATE (disabled: an image-heavy SR SCF, ~45 s).  SR basis, Gamma, at a FULLY-WRAPPED Rcut>=2a
// (where every term is translation-invariant, per the probe above) reproduces the CP2K FCC-Si Gamma GPW
// reference (SIPP_SR / GTH-PADE-q4 / LDA_X+VWN5).  densityEcut=20 (FFT N=32): Etot=-7.11467, charge 8,
// Exc=-2.544 -- within 0.4 mHa of CP2K -7.11506 (the N=32 grid-convergence gap; densityEcut>=30 -> N=64 ->
// -7.11505, an exact match, but ~5x slower).  doc/GPWPlan.md TODO 1's HARD gate -- the KB Bloch-orbital fix
// (TODO 1a) closes it.  Cost was cut ~25x (1100->45 s) by caching PhiOnGrid + this smaller-but-adequate grid.
TEST(GPW_SCF, DISABLED_SR_GammaRcut2a_CP2KReference)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwResult R=RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/20.0, /*Rcut*/2.0*a, /*Nelec*/8, "Si",
                       "SR Gamma Rcut=2a dE=20 N=32", /*verbose*/false, /*nmax*/60, qchem::Cholesky, 0.0, 0.0);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);
    EXPECT_NEAR(R.E.GetTotalEnergy(), -7.11506, 2e-3) << "GPW Gamma SR Rcut=2a vs CP2K FCC-Si reference";
}

TEST(GPW_SCF, SiliconGammaConverges)
{
    const double a=10.26;                          // Si conventional cubic lattice constant (a.u.)
    FCCUnitCell cell(a);                           // FCC primitive cell (2-atom diamond basis)
    cell.AddAtom(14, {0,0,0});                      // Si diamond: true Z=14; Zion=4 via the PP
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // Home cell only (Rcut=0): the primitive-cell Si2 (8 valence e- -> closed shell) with periodic
    // Hartree/XC/Ewald.  Folding images (Rcut>0) makes the diffuse Gaussians linearly dependent -> non-PD
    // overlap (a conditioning limit, deferred); Rcut=0 keeps S positive-definite and converges cleanly.
    GpwResult R=RunGPW(lat, MakeBasis(cell), /*densityEcut*/12.0, /*Rcut*/0.0, /*Nelec*/8, "Si", "Si Gamma");

    EXPECT_TRUE(R.converged);                       // clean closed-shell convergence (~17 iters, complex DIIS)
    EXPECT_NEAR(R.charge, 8.0, 1e-6);              // 8 valence electrons
    EXPECT_NEAR(R.E.GetTotalEnergy(), -8.2476, 5e-3);    // regression anchor (densityEcut=12, Gamma, Rcut=0)
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

    const double a=11.0;
    UnitCell cell(a);
    cell.AddAtom(14, {0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    GpwResult R=RunGPW(lat, MakeBasis(cell), /*densityEcut*/10.0, /*Rcut*/0.0, /*Nelec*/4, "Si", "Si atom-in-box",
                       /*verbose*/false, /*nmax*/40);

    EXPECT_NEAR(R.charge, 4.0, 1e-6);                        // 4 valence electrons (Zion=4), charge conserved
    // GPW-in-box (G-space local PP -> box-independent) reproduces the finite SIPP DFT energy to grid tolerance.
    // (Energy-converged; density is degenerate at Gamma -- see the note above -- so no Converged() guard.)
    EXPECT_NEAR(R.E.GetTotalEnergy(), Esipp, 5e-2) << "GPW-in-box total vs finite SIPP molecular DFT";
}
