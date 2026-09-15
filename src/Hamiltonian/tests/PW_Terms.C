// File: src/Hamiltonian/tests/PW_Terms.C  The plane-wave HAMILTONIAN TERM tests -- the PW pseudo/Hartree/XC
// terms against the basis's direct-grid oracles, the term cache, the Ewald ion-ion term (moved out of
// IntegrationTests/PlaneWaveDFTUT.C on 2026-09-15, doc/TestSuitePlan.md phase 5).  No SCF.
#include <map>
#include <set>
#include <memory>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include "gtest/gtest.h"

import qchem.Tests.PW_Fields;   // the direct-grid oracles + PWFixture
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
import qchem.BasisSet.PlaneWave.Evaluators;   // PW_Grid_Evaluator -- a UNIT TEST may reach the internal
import qchem.BasisSet.Lattice.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
import qchem.ScalarFunction;                 // ScalarFunction<double> -- arg of the moved-here field oracles
import qchem.Mesh;                           // qcMesh::MeshParams (Vee_Hartree's fit-basis factory arg; ignored)
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
import qchem.Blaze;          // mat_t<dcmplx>
import qchem.Math;           // Pi
import qchem.Hamiltonian.Internal.ExFunctional;     // the validated LDA functional interface
import qchem.Hamiltonian.Internal.SlaterExchange;   // Dirac exchange (alpha=2/3), eps_x = 3/4 v_x
import qchem.Hamiltonian.Internal.VWN_Correlation;  // VWN5 correlation (validated vs libxc)
import qchem.Hamiltonian.Internal.PWTerms;          // Ven_PP_Short/Long, Vee_Hartree (dcmplx Hamiltonian terms)
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
import qchem.Hamiltonian;                           // cStatic_HT / cDynamic_HT aliases (public term interfaces)
import qchem.Hamiltonian.Internal.Hamiltonian;      // cHamiltonianImp (the dcmplx Hamiltonian = sum of terms)
import qchem.Hamiltonian.Internal.Hamiltonians;     // Ham_PW_DFT (the assembled plane-wave LDA KS Hamiltonian)
import qchem.Hamiltonian.Internal.IonIon;           // IonIon<double> (ion-ion term: pair sum / Ewald via isFinite)
import qchem.Hamiltonian.Internal.Kinetic;          // Kinetic<dcmplx> (the shared kinetic term)
import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
import qchem.Energy;                                // EnergyBreakdown
import qchem.ChargeDensity.Imp.IrrepCD;             // PeriodicIrrepCD<dcmplx> (the periodic leaf -- 3c-2b)
import qchem.ChargeDensity.SeedCD;                  // SeedCD + PolarizedSeedCD (the spin-SAD staggering gate)
import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
import qchem.BasisSet.GMap;                             // ΔG_Map (the G-space coefficient map)
import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
import qchem.Symmetry.Irrep;                        // Irrep
import qchem.LASolver;                              // complex Hermitian eigensolver
import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

using namespace qchem;
using namespace qchem::tests::pw;

using BasisSet::PlaneWave::PlaneWave_IBS;
using BasisSet::PlaneWave::PW_Grid_Evaluator;   // internal grid evaluator (unit test may reach it directly)
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

namespace
{
// Vee_Hartree/Vxc_Quadrature now take their fit basis (from the basis's own factory) at construction, like
// FittedVee/FittedVxc.  These low-level term tests build one straight from the plane-wave basis at hand.
qchem::Hamiltonian::Vee_Hartree* NewPWHartree(const PlaneWave_IBS& pw)
{
    // Pure V_H[rho]: the Hartree term takes ONLY its fit basis now -- the long-range core-charge fold that
    // used to need a structure + local model here is its own term, Ven_PP_Long (doc/GPWPlan.md 0e-PP).
    return new qchem::Hamiltonian::Vee_Hartree(
        qchem::Hamiltonian::Vee_Hartree::fbs_t(pw.CreateCDFitBasisSet(nullptr, qcMesh::MeshParams{})));
}
// The XC term takes a QUADRATURE, and MakeDensitySampler reads off the strategy the fit basis supports --
// here a plane-wave (raster) basis, so the pair/collocation one.
qchem::Hamiltonian::Vxc_Quadrature* NewPWXC(const PlaneWave_IBS& pw, const qchem::Hamiltonian::Vxc_Quadrature::xc_t& xc)
{
    return new qchem::Hamiltonian::Vxc_Quadrature(xc, qchem::ChargeDensity::MakeDensitySampler(
        qchem::ChargeDensity::fitbasis_t(pw.CreateVxcFitBasisSet(nullptr, qcMesh::MeshParams{}))), SpinGroup::UnPolarized);
}


// The generalized IonIon Hamiltonian term: for a periodic Structure (isFinite()==false) it routes the
// ion-ion energy through the Ewald lattice sum (charges = the cell atoms' itsZ = the ion/valence charge),
// rather than the conditionally-convergent direct pair sum used for finite molecules.
TEST(PW_Terms, VnnPeriodicUsesEwald)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);   // FCC primitive
    auto cell=std::make_shared<UnitCell>(A);
    cell->AddAtom(4,rvec3_t(0.0,0.0,0.0));             // itsZ = Zion = 4 (Si valence)
    cell->AddAtom(4,rvec3_t(0.25,0.25,0.25));
    std::shared_ptr<const Structure> st=cell;
    EXPECT_FALSE(st->isFinite());

    qchem::Hamiltonian::IonIon<double> vnn(st);
    EnergyBreakdown eb;
    vnn.GetEnergy(eb, nullptr);                        // periodic branch ignores the density
    double ref=EwaldEnergy(*cell, rvec_t{4.0,4.0});
    EXPECT_NEAR(eb["Enn"], ref, 1e-9);                    // routes through Ewald
    EXPECT_NEAR(eb["Enn"], -8.40046, 1e-4);              // == the Si ion-ion Madelung energy
}


// The Ven_PP_Short Hamiltonian term (cStatic_HT<dcmplx>) routes through the framework's GetMatrix path
// and the abstract Integrals_Pseudo<dcmplx> dynamic_cast -- its matrix must equal the basis's assembly.
// This is the dependency inversion working end-to-end through a real Hamiltonian term.
TEST(PW_Terms, PWPseudoTermMatchesBasis)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> Amat(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(Amat);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    auto si=std::make_shared<Molecule>();
    si->Insert(new Atom(14, rvec3_t(0,0,0)));
    si->Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;

    // The TERMS own the models now and ask the basis to assemble them (the pseudo-wall).  The external
    // potential is THREE terms post-split (doc/GPWPlan.md 0e-PP): SHORT-range local + KB nonlocal here,
    // plus the LONG-range local in Ven_PP_Long.  Their matrices must sum to what the basis assembles --
    // the Hamiltonian adds them term by term, so check the same sum.
    qchem::Hamiltonian::Ven_PP_Short    extS(si, &loc);
    qchem::Hamiltonian::Ven_PP_NonLocal extN(si, &nl);
    qchem::Hamiltonian::cStatic_HT*     termS=&extS;    // the public term interface (as the Hamiltonian holds it)
    qchem::Hamiltonian::cStatic_HT*     termN=&extN;
    chmat_t M  = termS->GetMatrix(&pw, Spin::None) + termN->GetMatrix(&pw, Spin::None);
    chmat_t ref= pw.MakeSpeciesFieldMatrix(si.get(), loc, qchem::BasisSet::FieldRange::Short) + pw.MakeProjectorMatrix(si.get(),nl);

    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(M(i,j))-dcmplx(ref(i,j))), 0.0, 1e-12);
}


// The CP2K local-PP split invariant (doc/GPWPlan.md 0e-PP): the LONG (softened-Coulomb) + SHORT (poly-Gaussian)
// local-PP matrices sum EXACTLY to the full V_loc matrix -- the split only RELOCATES energy (into the Hartree
// term), it never changes the assembled electron-ion potential.  Pinned on the analytic plane-wave path (a
// pure form-factor callback swap), where the identity is machine-exact.
TEST(PW_Terms, LocalPPLongPlusShortEqualsFull)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> Amat(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell      cell(Amat);
    Lattice_3D    lat(cell, ivec3_t(1,1,1));
    PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    auto si=std::make_shared<Molecule>();
    si->Insert(new Atom(14, rvec3_t(0,0,0)));
    si->Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    const HGH_LocalPotential& loc=GetGTH("Si","LDA",4).local;

    chmat_t full =pw.MakeSpeciesFieldMatrix(si.get(), loc, qchem::BasisSet::FieldRange::Full);
    chmat_t split=pw.MakeSpeciesFieldMatrix(si.get(), loc, qchem::BasisSet::FieldRange::Long) + pw.MakeSpeciesFieldMatrix(si.get(), loc, qchem::BasisSet::FieldRange::Short);
    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(full(i,j))-dcmplx(split(i,j))), 0.0, 1e-12);
}


// The Vee_Hartree and Vxc_Quadrature dynamic terms route a complex density (PeriodicIrrepCD<dcmplx>) through the framework
// and must reproduce the basis's Repulsion / Overlap -- the inversion working for the
// density-dependent terms.  We feed a hand-built Hermitian density matrix (2 electrons in (e0+e1)/sqrt2).
TEST(PW_Terms, PWDynamicTermsMatchBasis)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();

    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);     // Hermitian density matrix
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;            // sets D(1,0)=conj(D(0,1)) -> a non-uniform density
    Irrep irr=F.pw.GetIrrep(Spin::None);
    qchem::ChargeDensity::PeriodicIrrepCD<dcmplx> cd(D, &F.pw, irr);

    // Hartree term matrix == basis Repulsion of the same density.
    std::unique_ptr<qchem::Hamiltonian::Vee_Hartree> h(NewPWHartree(F.pw));
    qchem::Hamiltonian::cDynamic_HT* ht=h.get();
    const chmat_t& Mh = ht->GetMatrix(&F.pw, Spin::None, &cd);
    chmat_t refh = RepulsionField(F.pw, cd);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(Mh(i,j))-dcmplx(refh(i,j))), 0.0, 1e-10);

    // XC term (Dirac exchange) matrix == the basis's FFT route: rho(r) via inverse FFT of rho-tilde,
    // v_xc applied pointwise on the grid, forward FFT to the matrix (what the term itself does).
    auto dirac=std::make_shared<qchem::Hamiltonian::SlaterExchange>(2.0/3.0);
    std::unique_ptr<qchem::Hamiltonian::Vxc_Quadrature> xc(NewPWXC(F.pw, dirac));
    qchem::Hamiltonian::cDynamic_HT* xt=xc.get();
    const chmat_t& Mx = xt->GetMatrix(&F.pw, Spin::None, &cd);
    rvec_t rho=GridOf(F.pw).RhoOnGrid(RhoTilde(F.pw, D));
    rvec_t vxc(rho.size());
    for (size_t q=0;q<rho.size();q++) vxc[q]=dirac->GetVxc(rho[q]);
    chmat_t refx = OverlapOnGrid(F.pw, vxc);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(Mx(i,j))-dcmplx(refx(i,j))), 0.0, 1e-10);
}


// Item K: the XC quadrature grid now comes from the FIT basis, so relCutoff (the CP2K REL_CUTOFF density-
// cutoff idea) is the genuine fit-accuracy lever.  (1) a larger relCutoff yields a DENSER fit {G}; (2) it
// feeds through to the Vxc matrix (relCutoff=1 -> the orbital grid, bit-identical to the FFT route; >1 -> a
// finer, less-aliased quadrature); (3) the matrix GRID-CONVERGES (successive refinements shrink).  The fit is
// non-variational, so acceptance is convergence of the matrix, NEVER an energy anchor.
TEST(PW_Terms, ItemK_RelCutoffDensifiesAndConvergesVxc)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;                       // non-uniform => nonlinear v_xc aliases on a coarse grid
    Irrep irr=F.pw.GetIrrep(Spin::None);
    qchem::ChargeDensity::PeriodicIrrepCD<dcmplx> cd(D, &F.pw, irr);
    auto dirac=std::make_shared<qchem::Hamiltonian::SlaterExchange>(2.0/3.0);

    auto vxcAt=[&](double relCutoff, size_t& nGfit)
    {
        qcMesh::MeshParams mp; mp.relCutoff=relCutoff;
        auto fb=qchem::ChargeDensity::fitbasis_t(F.pw.CreateVxcFitBasisSet(nullptr, mp));
        nGfit=fb->GetNumFunctions();
        std::unique_ptr<qchem::Hamiltonian::Vxc_Quadrature> xc(
            new qchem::Hamiltonian::Vxc_Quadrature(dirac, qchem::ChargeDensity::MakeDensitySampler(fb), SpinGroup::UnPolarized));
        return chmat_t(static_cast<qchem::Hamiltonian::cDynamic_HT*>(xc.get())->GetMatrix(&F.pw, Spin::None, &cd));
    };
    auto froDiff=[&](const chmat_t& A, const chmat_t& B)
    {
        double s=0.0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) s+=std::norm(dcmplx(A(i,j))-dcmplx(B(i,j)));
        return std::sqrt(s);
    };

    size_t n1,n4,n16;
    chmat_t M1=vxcAt(1.0,n1), M4=vxcAt(4.0,n4), M16=vxcAt(16.0,n16);

    // (1) relCutoff densifies the fit {G}.
    EXPECT_GT(n4, n1);
    EXPECT_GT(n16, n4);

    // (2) relCutoff=1 reproduces the orbital-grid FFT route exactly (the fix is inert at Gamma).
    rvec_t rho=GridOf(F.pw).RhoOnGrid(RhoTilde(F.pw, D));
    rvec_t vxc(rho.size());
    for (size_t q=0;q<rho.size();q++) vxc[q]=dirac->GetVxc(rho[q]);
    chmat_t refx=OverlapOnGrid(F.pw, vxc);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
        EXPECT_NEAR(std::abs(dcmplx(M1(i,j))-dcmplx(refx(i,j))), 0.0, 1e-10);

    // (3) the denser grid CHANGES the matrix (relCutoff is live) and it GRID-CONVERGES (Cauchy: the later,
    //     finer refinement moves the matrix less than the earlier one).
    double d1=froDiff(M4,M1), d2=froDiff(M16,M4);
    EXPECT_GT(d1, 1e-9);      // the fit grid actually feeds through to the quadrature
    EXPECT_LT(d2, d1);        // convergence
}


// Regression guard for the Hamiltonian-framework cache bug: rDynamic_HT_Imp::GetMatrix must invalidate
// its Irrep-keyed cache when the density changes (else it returns a STALE matrix to any caller that
// didn't happen to call GetTotalEnergy(new cd) in between).  Two different densities (same Irrep, no
// GetEnergy between): the second GetMatrix must reflect the second density, not return the first's
// (uniform => zero) Hartree matrix.  Was a DISABLED known-failure; fixed by the cd-change cache clear.
TEST(PW_Terms, DynamicTermCacheFreshAcrossDensity)
{
    PWFixture F;
    size_t  n=F.pw.GetNumFunctions();
    Irrep   irr=F.pw.GetIrrep(Spin::None);

    hmat_t<dcmplx> D1=blazem::zeroH<dcmplx>(n);  D1(0,0)=2.0;                            // uniform density
    hmat_t<dcmplx> D2=blazem::zeroH<dcmplx>(n);  D2(0,0)=1.0; D2(0,1)=1.0; D2(1,1)=1.0;  // modulated density
    qchem::ChargeDensity::PeriodicIrrepCD<dcmplx> cd1(D1,&F.pw,irr), cd2(D2,&F.pw,irr);

    std::unique_ptr<qchem::Hamiltonian::Vee_Hartree> hart(NewPWHartree(F.pw));
    qchem::Hamiltonian::cDynamic_HT* ht=hart.get();
    ht->GetMatrix(&F.pw, Spin::None, &cd1);                       // populates the cache for cd1
    const chmat_t& M2 = ht->GetMatrix(&F.pw, Spin::None, &cd2);   // BUG: returns cd1's stale matrix

    chmat_t ref2=RepulsionField(F.pw, cd2);                            // the correct V_H for cd2
    double  diff=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            diff=std::max(diff, std::abs(dcmplx(M2(i,j))-dcmplx(ref2(i,j))));
    EXPECT_NEAR(diff, 0.0, 1e-10);   // fails today (M2 is cd1's matrix); passes after the cache-invalidation fix
}


// G-space Hartree path: V_H assembled directly from the density's Fourier coefficients rho-tilde(dm)
// (PlaneWave_IBS::MakeFourierDensity -> Repulsion(ΔG_Map)) must equal the real-space route
// (sample rho(r) on the grid -> ForwardDFT -> V_H).  This is the O(n^2) G-space replacement for the
// O(Npts*n^2) pointwise sampling -- the foundation of the FFT speed-up.  Fast: no SCF, one block.
TEST(PW_Terms, HartreeFromFourierMatchesPointwise)
{
    const double a=10.26;
    FCCUnitCell    cell(a);
    Lattice_3D     lat(cell, ivec3_t(1,1,1));
    PlaneWave_IBS  pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    size_t n=pw.GetNumFunctions();
    ASSERT_GT(n, 3u);
    Irrep irr=pw.GetIrrep(Spin::None);

    // A non-trivial Hermitian density matrix: a uniform base plus some off-diagonal (cross-G) structure.
    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++) D(i,i)=8.0/double(n);
    D(0,1)=dcmplx(0.10,0.05);
    D(0,2)=dcmplx(-0.07,0.02);
    D(1,2)=dcmplx(0.04,-0.06);
    qchem::ChargeDensity::PeriodicIrrepCD<dcmplx> cd(D, &pw, irr);   // IS-A ScalarFunction rho(r)=phi^H D phi

    chmat_t VA=RepulsionField(pw, cd);                             // real-space: sample + ForwardDFT
    chmat_t VB=HartreeFromRhoTilde(pw, RhoTilde(pw, D));     // G-space: direct from D

    double maxd=0;                                           // the matrices agree elementwise => so does any derived E_H
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) maxd=std::max(maxd, std::abs(VA(i,j)-VB(i,j)));
    EXPECT_LT(maxd, 1e-9);
}

} //namespace
