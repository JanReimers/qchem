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
                 qchem::Ortho ortho=qchem::Cholesky, double orthoTol=0.0, double collRcut=0.0)
{
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, std::move(mol), densityEcut, Rcut, collRcut));
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

// (1d) BULK DISPERSION (Rcut>0) -- DISABLED: blocked by SIPP basis conditioning, not by the general-k code.
// Turning on inter-cell images to get real k-dispersion needs the analytic Bloch overlap to be positive-
// definite, which for the diffuse SIPP valence basis happens only at Rcut>=3a (min eig(S(k))=-0.0016 @2a ->
// +4.3e-6 @3a, measured across the 2x2x2 mesh).  But +4.3e-6 is NEAR-SINGULAR and does NOT improve with more
// images (+4.3e-6 @4a too) -- an INTRINSIC near-linear-dependence of the SIPP Gaussians folded periodically,
// not a truncation artifact.  Cholesky then amplifies noise by ~1/4e-6 and the SCF diverges (Etot blows up);
// SVD/Eigen truncation drops the near-null directions but a DIFFERENT count per k-block -> non-uniform per-k
// dimensions, which the framework's cross-k assembly cannot yet handle ("Matrix sizes do not match").
//
// This is the documented basis-conditioning wall (doc/GPWPlan.md sec 3b/sec 5), the deferred blocker for
// self-consistent bulk GPW.  The general-k MACHINERY itself is correct and validated at the matrix level by
// the Bloch invariants in GPW_UT (Hermiticity, the Bloch translation law, k->-k conjugation) and by the
// multi-k plumbing gate above.  The DECOUPLED image sets (overlap Rcut vs collocation collRcut) that make a
// dispersive run affordable are in place and exercised here; what remains for a converged bulk energy is the
// conditioning fix -- a less-diffuse/duals-cleaned valence basis, OR framework support for per-k orbital
// dimensions (canonical orthogonalisation dropping the near-null directions independently per k-block).
TEST(GPW_SCF, DISABLED_SiliconBulk2x2x2Dispersive)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    // overlap Rcut=3a (PSD but near-singular); collocation collRcut=1.5a (the local density reach -- the
    // decouple keeps the per-iteration collocation ~15x cheaper).  Diverges at present (see the header note).
    RunGPW(lat, MakeBasis(cell), /*densityEcut*/8.0, /*Rcut*/3.0*a, /*Nelec*/8, "Si", "2x2x2 rc=3a cc=1.5a",
           /*verbose*/true, /*nmax*/60, qchem::SVD, 1e-5, /*collRcut*/1.5*a);
}

// DIAGNOSTIC (disabled): the near-singular Bloch overlap is a BASIS problem (SIPP's diffuse 0.09s/0.06p
// primitives, RMS radius ~5 a.u., go near-linearly-dependent when Bloch-summed), NOT a GPW-code problem.
// Compare min eig(S(k)) over the 2x2x2 mesh at Rcut=3a for SIPP vs the short-range SIPP_SR (diffuse dropped).
TEST(GPW_SCF, DISABLED_Bulk2x2x2ConditioningBasisCompare)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    namespace L3=BasisSet::Lattice_3D;
    auto worstMinEig=[&](std::shared_ptr<const Real_BS> mol)->double
    {
        std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, std::move(mol), 0.0, 3.0*a));
        double worst=1e30;
        for (auto b : bs->Iterate<qchem::BasisSet::Complex_OIBS>())
        {
            chmat_t S=b->Overlap();
            rvec_t w; mat_t<dcmplx> U; blazem::eigen(S, w, U);
            for (size_t i=0;i<w.size();i++) worst=std::min(worst,w[i]);
        }
        return worst;
    };
    std::cout << "min eig(S(k)) @Rcut=3a  SIPP="<<worstMinEig(MakeBasis(cell))
              << "  SIPP_SR="<<worstMinEig(MakeBasisSR(cell)) << std::endl;
}

// Find the SMALLEST Rcut at which SIPP_SR's Bloch overlap is PSD -- so overlap AND collocation can use ONE
// consistent (and cheap) image reach (the decouple with collRcut != Rcut breaks the density normalisation:
// collocated charge = Tr(D S_collRcut) != Tr(D S_Rcut) = Nelec).
TEST(GPW_SCF, DISABLED_SR_OverlapSpectrumSweep)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    namespace L3=BasisSet::Lattice_3D;
    for (double rc : {1.0, 1.5, 2.0, 3.0})
    {
        std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasisSR(cell), 0.0, rc*a));
        double worst=1e30;
        for (auto b : bs->Iterate<qchem::BasisSet::Complex_OIBS>())
        {
            chmat_t S=b->Overlap(); rvec_t w; mat_t<dcmplx> U; blazem::eigen(S, w, U);
            for (size_t i=0;i<w.size();i++) worst=std::min(worst,w[i]);
        }
        std::cout << "SIPP_SR Rcut="<<rc<<"a  worst min(eig S(k))="<<worst<<(worst<0?" (INDEF)":" (PSD)")<<std::endl;
    }
}

// The well-conditioned dispersive bulk run: SIPP_SR (min eig(S(k))~0.016 -> plain Cholesky) at Rcut=1.5a,
// where the SR overlap is PSD AND converged (so overlap + collocation share ONE consistent reach, collRcut=0;
// the density then integrates to the right charge Tr(D S)=Nelec).  FINDINGS (2026-07-09):
//  - The SCF is now STABLE and CONVERGES (no divergence) -- the conditioning blocker is GONE with a proper
//    basis.  Etot settles at ~-15.15 (decoupled Rcut=3a/collRcut=1.5a gives the SAME ~-15.15, so it is NOT a
//    decouple/normalisation artifact -- the two configs agree).
//  - BUT ~-15.15 is ~2x the plane-wave bulk -7.7613 and MORE negative than the Rcut=0 GPW -8.25 -- a real
//    OVER-BINDING (a smaller basis would raise, not lower, the energy), pointing to a bug in the Rcut>0
//    EXTERNAL-PP / nuclear assembly: MakeLocalPP/MakeSeparablePP still use the cell's OWN atoms with no
//    periodic-IMAGE potential (GPWPlan sec 1 Increment 3 limit (c) / sec 4 "rigorous periodic external
//    potential").  At Rcut=0 the home-cell approximation is exact; with inter-cell hopping on it is not.
//  - NEXT: sum the PP over lattice images (image-atom potentials), then validate the bulk total against an
//    INDEPENDENT GPW code (CP2K, same GTH PP) -- our own PW -7.7613 is an Ecut=4 internal number, not an
//    absolute literature value.  (Runs slow: ~25 s/iter for 8 k-points; bump the timeout when re-enabling.)
TEST(GPW_SCF, DISABLED_SiliconBulk2x2x2SR)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    RunGPW(lat, MakeBasisSR(cell), /*densityEcut*/8.0, /*Rcut*/1.5*a, /*Nelec*/8, "Si", "2x2x2 SR rc=cc=1.5a",
           /*verbose*/true, /*nmax*/60, qchem::Cholesky, 0.0, /*collRcut*/0.0);
}

// DIAGNOSTIC (disabled): is the COLLOCATED density consistent with the analytic Bloch overlap once images are
// on?  The G=0 collocation weight W_0(i,j)*Omega must equal the analytic overlap S(i,j) = Sum_R <chi_i|chi_j(R)>
// (both are integral chi_i^k* chi_j^k).  Tested at Rcut=0 by GPW.CollocationOverlapMatchesAnalytic but NOT with
// images -- if it diverges at Rcut>0 the density is inconsistent with the orbitals => every energy term balloons.
TEST(GPW_SCF, DISABLED_CollocationVsAnalyticOverlapWithImages)
{
    using BasisSet::Lattice_3D::GPW_IBS;
    using BasisSet::Lattice_3D::GPW_Evaluator;
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    for (double dE : {8.0, 12.0, 20.0, 30.0, 50.0})
    for (double rc : {0.0, 1.5*a})
    {
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), MakeBasisSR(cell), dE, rc, /*collRcut*/0.0);
        const GPW_Evaluator& ev = gpw;
        const BasisSet::Complex_OIBS& g = gpw;
        chmat_t S = g.Overlap();                 // analytic Bloch overlap
        G_ERI3 ov = ev.Overlap3CTensor();        // collocation weights
        int c0=-1; for (size_t c=0;c<ov.columns.size();c++)
            if (ov.columns[c].dm.x==0&&ov.columns[c].dm.y==0&&ov.columns[c].dm.z==0) c0=int(c);
        size_t n=S.rows(); double num=0,den=0, trW=0, trS=0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
        {
            double w=std::real(dcmplx(ov.weights[c0](i,j)))*ov.volume, s=std::real(dcmplx(S(i,j)));
            double d=w-s; num+=d*d; den+=s*s;
        }
        for (size_t i=0;i<n;i++){ trW+=std::real(dcmplx(ov.weights[c0](i,i)))*ov.volume; trS+=std::real(dcmplx(S(i,i))); }
        std::printf("dEcut=%4.0f Rcut=%.2f  ||W0*Om - S||/||S|| = %.4f   trace(W0*Om)=%.4f trace(S)=%.4f\n",
                    dE, rc/a, std::sqrt(num/den), trW, trS);
    }
}

// DIAGNOSTIC (disabled): which 1E/PP MATRIX balloons at Rcut>0?  Print trace + Frobenius norm of the kinetic,
// local-PP and separable-PP matrices at Rcut=0 vs 1.5a (no SCF).  A term whose matrix grows unphysically with
// images is the double-count.
TEST(GPW_SCF, DISABLED_PPMatrixTraceProbe)
{
    using BasisSet::Lattice_3D::GPW_IBS;
    auto frob=[](const chmat_t& M){ double s=0; for (size_t i=0;i<M.rows();i++) for (size_t j=0;j<M.rows();j++) s+=std::norm(dcmplx(M(i,j))); return std::sqrt(s); };
    auto tr  =[](const chmat_t& M){ double s=0; for (size_t i=0;i<M.rows();i++) s+=std::real(dcmplx(M(i,i))); return s; };
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto st=lat.GetStructure();
    Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH("Si","LDA",4);
    for (double rc : {0.0, 1.5*a})
    {
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), MakeBasisSR(cell), /*densityEcut*/12.0, rc, 0.0);
        const BasisSet::Complex_OIBS& g = gpw;
        chmat_t Kin =g.Kinetic();
        chmat_t Vloc=gpw.MakeLocalPotential   (st.get(), pp.local);
        chmat_t Vnl =gpw.MakeSeparablePotential(st.get(), pp.nonlocal);
        std::printf("Rcut=%.2f  Kin[tr=%.3f fro=%.3f]  Vloc[tr=%.3f fro=%.3f]  Vnl[tr=%.3f fro=%.3f]\n",
                    rc/a, tr(Kin),frob(Kin), tr(Vloc),frob(Vloc), tr(Vnl),frob(Vnl));
    }
}

// DIAGNOSTIC (disabled): TRANSLATION INVARIANCE.  Shifting ALL atoms by a constant must not change the total
// energy.  If it does at Rcut>0, the collocation/grid handling of atom position (esp. the corner atom at the
// origin, whose orbital straddles the cell boundary) is the bug.
TEST(GPW_SCF, DISABLED_SR_TranslationInvariance)
{
    const double a=10.26;
    auto run=[&](double shift, const char* lbl)->double
    {
        FCCUnitCell cell(a);
        cell.AddAtom(14, {0.0+shift, 0.0+shift, 0.0+shift});
        cell.AddAtom(14, {0.25+shift,0.25+shift,0.25+shift});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        GpwResult R=RunGPW(lat, MakeBasisSR(cell), 12.0, /*Rcut*/1.5*a, 8, "Si", lbl, false, 20);
        return R.E.GetTotalEnergy();
    };
    double e0=run(0.0,  "shift=0.00");
    double e1=run(0.125,"shift=0.125");   // move the origin atom off the cell corner
    std::printf("SR Gamma Rcut=1.5a: Etot(shift=0)=%.5f  Etot(shift=0.125)=%.5f  |diff|=%.5f\n",
                e0, e1, std::fabs(e0-e1));
}

// DIAGNOSTIC (disabled): is the corner-atom translation-invariance violation a grid-RESOLUTION effect?  Run the
// corner atom (shift=0) at rising densityEcut -- if Etot -> the off-corner value (~-8.4) the sharp Gaussian
// density was under-resolved (fix = higher densityEcut / smoother basis); if it stays ~-15 it is a real bug.
TEST(GPW_SCF, DISABLED_SR_CornerAtomVsDensityEcut)
{
    const double a=10.26;
    for (double dE : {12.0, 20.0, 30.0, 45.0})
    {
        FCCUnitCell cell(a);
        cell.AddAtom(14, {0,0,0});
        cell.AddAtom(14, {0.25,0.25,0.25});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        char lbl[48]; std::snprintf(lbl,sizeof lbl,"corner dE=%.0f",dE);
        RunGPW(lat, MakeBasisSR(cell), dE, /*Rcut*/1.5*a, 8, "Si", lbl, false, 20);
    }
}

// FAST DIAGNOSTIC (disabled): single-k (Gamma) minimal repro of the Rcut>0 over-binding -- compare charge +
// Etot for SR at Rcut=0 (home cell) vs Rcut=1.5a (images on).  If images DOUBLE the charge -> an occupation/
// normalisation double-count; if charge stays 8 but Etot ~doubles -> an energy-term (external-PP image) bug.
TEST(GPW_SCF, DISABLED_SR_GammaRcutChargeProbe)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));                 // Gamma only -- single k, fast
    GpwResult R0=RunGPW(lat, MakeBasisSR(cell), 8.0, /*Rcut*/0.0,   8, "Si", "SR Gamma Rcut=0",  false, 30);
    GpwResult R1=RunGPW(lat, MakeBasisSR(cell), 8.0, /*Rcut*/1.5*a, 8, "Si", "SR Gamma Rcut=1.5a",false, 30);
    std::cout << "SR Gamma: charge Rcut0="<<R0.charge<<" Rcut1.5a="<<R1.charge
              << " | Etot Rcut0="<<R0.E.GetTotalEnergy()<<" Rcut1.5a="<<R1.E.GetTotalEnergy() << std::endl;
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

// (1b) BULK Si: fold in the lattice images (Rcut>0) so the Gaussians hop between cells (true crystal, not
// Si2-in-a-box).  The diffuse SIPP Gaussians then go linearly dependent -> the Bloch overlap S is singular
// and the default Cholesky fails; canonical orthogonalisation (Eigen with a near-null eigenvalue cutoff)
// drops the redundant directions and the SCF runs.  The inter-cell screening relaxes the Rcut=0 over-binding
// toward the plane-wave bulk -7.2273.  (DISABLED while the tolerance/energy are tuned; see the sweep below.)
TEST(GPW_SCF, DISABLED_SiliconBulkOrtho)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The analytic single-sum Bloch overlap is INDEFINITE at small Rcut (a truncated single sum is not the
    // Gram matrix of the truncated Bloch functions) and only turns positive-definite once the image sphere is
    // large enough to converge it: min eig(S) = -0.12(1.5a) -> -0.0016(2a) -> +5.9e-5(3a) -> +5.9e-5(>=3a,
    // converged).  Overlap integrals are cheap so a large Rcut is affordable for the 1E matrices -- BUT the
    // density COLLOCATION re-sums all images (~450 cells at Rcut=3a) at every grid point, which is the cost.
    //
    // CONCLUSION (why this is not the bulk path): a single Gamma point -- even with a huge Rcut -- is one point
    // of the Brillouin zone, not the bulk (the bulk total energy is a BZ integral).  Real bulk needs MULTI-K
    // sampling (general-k Bloch phases e^{ik.R} in the GPW lattice sums + collocation), reduced to the
    // IRREDUCIBLE BZ wedge for efficiency.  That -- not Gamma + large Rcut -- is the next increment for solids.
    namespace L3=BasisSet::Lattice_3D;
    for (double rc : {1.5, 2.0, 3.0, 4.0})
    {
        std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasis(cell), /*densityEcut*/0.0, rc*a));
        chmat_t S;
        for (auto b : bs->Iterate<qchem::BasisSet::Complex_OIBS>()) S = b->Overlap();
        rvec_t w; mat_t<dcmplx> U; blazem::eigen(S, w, U);
        double lo=w[0]; for (size_t i=0;i<w.size();i++) lo=std::min(lo,w[i]);
        std::cout << "Rcut=" << rc << "a  min(eig S_analytic)=" << lo << (lo<0?"  (indefinite)":"  (PSD)") << std::endl;
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

// DIAGNOSTIC (disabled): does the atom-in-box TOTAL converge to the finite pseudo-atom energy as the box
// grows?  If yes -> the energy expression (G=0/alignment) is right and the crystal gap is purely Rcut=0-not-
// bulk.  If it retains a ~1/L tail / plateaus away from Efin -> the long-range G=0 bookkeeping is broken (the
// erf/erfc split is needed).  Run with --gtest_also_run_disabled_tests --gtest_filter=*BoxSizeSweep*.
TEST(GPW_SCF, DISABLED_SiAtomBoxSizeSweep)
{
    AtomCalculation cFin(14, 14-4, {.type=AtomType::Slater, .accuracy=BasisSetAccuracy::High, .pseudopotential=true});
    std::cout << "[finite target] Efin="<<cFin.Energy()<<std::endl;
    for (double a : {11.0, 15.0, 20.0})
    {
        UnitCell cell(a);
        cell.AddAtom(14, {0.5,0.5,0.5});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        char lbl[64]; std::snprintf(lbl, sizeof lbl, "a=%.0f", a);
        RunGPW(lat, MakeBasis(cell), /*densityEcut*/10.0, /*Rcut*/0.0, /*Nelec*/4, "Si", lbl, /*verbose*/false, /*nmax*/40);
    }
}
