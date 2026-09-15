// File: src/Hamiltonian/tests/GPW_XC_FD.C  The E_xc / H_xc FINITE-DIFFERENCE CONSISTENCY probes on the GPW
// evaluator's own XC route (moved out of the GPW basis unit tests on 2026-09-15, doc/TestSuitePlan.md phase 5:
// they are the only GPW tests that need the XC FUNCTIONAL OBJECTS, and the Lattice basis test exe is
// deliberately forbidden from linking the Hamiltonian).

// File GPW_UT.C  GPW at Gamma: periodic Gaussian 1-electron integrals + the DFT-tier collocation primitive.
//
// GPW puts GAUSSIAN orbitals on a lattice.  Its one-electron matrices are lattice sums of the ordinary
// (finite) two-centre integrals,  M_ij = Sum_R <chi_i | O | chi_j(.-R)>, computed by the molecular Gaussian
// basis (Gaussian::LatticeSum1E) and delegated to by GPW_Evaluator; GPW_IBS is the thin Orbital_1E_IBS on top.
//
// 1E validation mirrors L_PP: the SAME Si valence Gaussian basis gives the SAME overlap / kinetic (<p^2>) /
// nuclear matrices whether the atom is a finite Molecule or centred in a large periodic UnitCell.  Two teeth:
//   (1) home cell only (R={0}): GPW reproduces the finite matrices EXACTLY (same analytic M&D kernels);
//   (2) with the periodic images summed in: GPW matches to the (tiny) image tail, shrinking as the cell grows.
//
// DFT tier: GPW's genuinely-new primitive is COLLOCATION (rho=sum D chi chi on a grid -> FFT -> Poisson;
// integrate a grid potential back against the Gaussians).  Two grid-convergent checks isolate collocate +
// integrate-back against the analytic overlap, without a full SCF: the G=0 collocation weight is the grid-
// quadrature overlap, and a constant potential integrates back to V0*<i|j>.  (A physically rigorous periodic
// nuclear attraction (Ewald), general-k Bloch phases, and the full periodic SCF energy are later increments.)
#include "gtest/gtest.h"
#include <memory>
#include <cmath>
#include <complex>
#include <chrono>
#include <cstdlib>   // getenv/atof (the ill-conditioned charge probe's Ecut knob)
#include <iostream>

import qchem.Structure;                         // Molecule, Atom
import qchem.UnitCell;                          // UnitCell
import qchem.BasisSet;                          // Real_BS
import qchem.BasisSet.Orbital_1E_IBS;           // Real_OIBS / Complex_OIBS + cached Overlap()/Kinetic()/Nuclear()
import qchem.BasisSet.Gaussian.Point.Factory;         // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Gaussian.PG_Cart;         // direct PG_Cart construction (the diffuse-d V_long oracle gate)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;       // GPW_IBS (the basis under test)
import qchem.Pseudopotential.SeparablePotential; // HGH_SeparablePotential + the _R / _Gaussian faces (KB gate)
import qchem.Pseudopotential.GTH_Potentials;     // GetGTH (the Si GTH-LDA-q4 projector data)
import qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator; // GPW_Evaluator (tests may cheat-import internals) -- DFT tier
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;     // Gaussian::LatticeSum1E::CollocateDensity (analytic collocation)
import qchem.LASolver;                       // LASolver<dcmplx> (the k=1/4 spectrum gate)
import qchem.Symmetry.Factory;               // BlochFactory (arbitrary-shift k for the k=1/4 continuity gate)
import qchem.BasisSet.DeltaFit_IBS;          // DeltaFit_IBS (the delta representation's OverlapDiagonal gate)
import qchem.Fitting.FitOperations;          // OrthogonalFit -- the projection/metric invariant this gate pins
import qchem.Mesh.Quadrature;                // qcMesh::Mesh (the delta basis's quadrature)
import qchem.Symmetry.Lattice_3D.SpaceGroup;     // SpaceGroup::Detect + DirectOp (the T3 stream-fold unit gates)
import qchem.BasisSet.GMap;            // Projector3<dcmplx> / ΔG_Map (the collocation tensor + rho-tilde)
import qchem.Hamiltonian.Internal.ExFunctional;   // ExFunctional (the v_xc/eps_xc face; XC-consistency probe)
import qchem.Hamiltonian.Internal.SlaterExchange; // SlaterExchange (Dirac exchange -- the SCF's own X term)
import qchem.Hamiltonian.Internal.VWN_Correlation;// VWN_Correlation (VWN5 -- the SCF's own C term)
import qchem.Blaze;                             // hmat_t element access / rows()
import qchem.Math;                              // Pi (the Bloch-phase e^{2 pi i k.n})
import qchem.Vector3D;                          // rvec3_t arithmetic (r + R0)
import qchem.Types;

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;
using BasisSet::Complex_OIBS;
using BasisSet::Gaussian::GPW_IBS;
using BasisSet::Gaussian::GPW_Evaluator;
using qchem::BasisSet::Gaussian::BasisSetData;

namespace
{
// The valence Si Gaussian basis (SIPP, Cartesian) on ANY structure -- the L_PP builder.  The Engine argument
// is the integral-engine switch point (see Gaussian::LatticeSum1E): Engine::MnD here (its AtCenter + analytic
// 2C kernels make the periodic sum exact/trivial); Engine::LibCint would be the faster path once PG_LibCint
// realises LatticeSum1E -- GPW itself is unchanged either way.
std::unique_ptr<Real_BS> MakeBasis(const Structure& st)
{
    return std::unique_ptr<Real_BS>(
        BasisSet::Gaussian::Factory(BasisSetData::SIPP, &st,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
}
} //anon

// XC POTENTIAL-CONSISTENCY PROBE (doc/GPWPlan.md 0b instrument).  Question under test: is the assembled
// H_xc the EXACT D-derivative of the DISCRETE energy E_xc(D) = Sum_q w [eps_x+eps_c](rho_q) rho_q, where
// rho_q is the ball-limited grid density of the SCF's own chain?  The probe replicates the pair quadrature's
// route verbatim at the evaluator level (collocate -> nested {G_L} combine -> RhoOnGrid; v_xc pointwise ->
// raster ForwardFFT -> per-level restriction -> analytic IntegratePotential) and compares the central
// finite difference  [E_xc(D+h dD) - E_xc(D-h dD)]/2h  against  Re Tr(H_xc(D) dD).
//   - The HARTREE control (bilinear, kernel baked) isolates harness error: it must agree to FD accuracy.
//   - Probe 1: a PSD-like density (rho_q > 0 everywhere) -- the smooth-functional regime.
//   - Probe 2: an INDEFINITE D (rho_q < 0 over part of the grid) -- the Kerker-mixed-field regime that
//     exercises the functionals' rho<=0 guards (E integrand and v_xc must be a consistent kink).
// If both probes agree to FD accuracy, the E_xc/H_xc representation fork hypothesized for the NaF fine-grid
// attractor is FALSIFIED at this seam and the inconsistency must live elsewhere; if not, the disagreement
// magnitude localizes it.  Two step sizes separate FD truncation error from a genuine fork.
TEST(GPW_XC, PotentialConsistencyFD)
{
    const double a=10.26;
    FCCUnitCell cell(a);                                 // the production shape: cross-cell pairs + ladder
    cell.AddAtom(14,{0,0,0});
    cell.AddAtom(14,{0.25,0.25,0.25});
    std::shared_ptr<const Real_BS> mol=MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/6.0, BasisSet::Gaussian::CellImages::Periodic, 2.0,
                BasisSet::PlaneWave::RasterPolicy::AliasFree);   // EXACT-QUADRATURE gate (production default = BallOnly)
    const GPW_Evaluator& ev=gpw;
    const auto& grid=ev.DensityGrid();
    const size_t n=static_cast<const Complex_OIBS&>(gpw).GetNumFunctions();

    // The SCF's own functionals (Ham_PW_DFT::BuildTerms builds exactly these two PWFittedVxc terms).
    qchem::Hamiltonian::SlaterExchange  exch(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation corr;
    const qchem::Hamiltonian::ExFunctional* xcs[2]={&exch,&corr};

    Projector3<dcmplx> ov =ev.Overlap3CTensor();                     // rho-tilde (no kernel) -- the pair-quadrature route
    Projector3<dcmplx> cou=ev.Repulsion3CTensor();                   // V_H (Coulomb kernel baked) -- the control

    auto rhoOf=[&](const chmat_t& D)->rvec_t { return grid.RhoOnGrid(Contract(ov,D)); };
    auto Exc=[&](const rvec_t& rho)->double              // == Vxc_Quadrature::GetEnergy (both terms)
    {
        rvec_t e(rho.size());
        for (size_t q=0;q<rho.size();q++)
        {
            double s=0.0; for (auto xc : xcs) s+=xc->GetEpsXc(rho[q]);
            e[q]=s*rho[q];
        }
        return grid.Integral(e);
    };
    auto Hxc=[&](const rvec_t& rho)->chmat_t             // == Vxc_Quadrature::MakeMatrix (both terms summed)
    {
        rvec_t v(rho.size());
        for (size_t q=0;q<rho.size();q++)
        {
            double s=0.0; for (auto xc : xcs) s+=xc->GetVxc(rho[q]);
            v[q]=s;
        }
        cvec_t vt=grid.ForwardFFT(v);                    // full raster (the OrthoScalarFitter route)
        return ev.OverlapMatrix([&](const ivec3_t& dm)->dcmplx { return grid.GridCoeff(vt,dm); });
    };
    auto EH=[&](const chmat_t& D)->double                // == PW_Hartree::GetEnergy: 1/2 Tr(D H_H(D))
    {
        ΔG_Map VH=Contract(cou,D);
        chmat_t HH=ev.OverlapMatrix([&](const ivec3_t& dm)->dcmplx
            { auto it=VH.find(dm); return it==VH.end()?dcmplx(0.0):it->second; });
        dcmplx tr(0.0);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) tr+=dcmplx(D(i,j))*dcmplx(HH(j,i));
        return 0.5*std::real(tr);
    };
    auto trace=[&](const chmat_t& H, const chmat_t& dD)->double   // Re Tr(H dD), both Hermitian
    {
        dcmplx tr(0.0);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) tr+=dcmplx(H(i,j))*dcmplx(dD(j,i));
        return std::real(tr);
    };
    auto shifted=[&](const chmat_t& D, const chmat_t& dD, double h)->chmat_t
    {
        chmat_t Dh(D);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) Dh(i,j)=dcmplx(D(i,j))+h*dcmplx(dD(i,j));
        return Dh;
    };

    // Deterministic Hermitian (real-symmetric at Gamma) perturbation direction, O(0.1) entries.
    chmat_t dD(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++) dD(i,j)=dcmplx(0.1*std::sin(1.0+double(i)+2.0*double(j)),0.0);

    auto probe=[&](const chmat_t& D, const char* label)->double   // returns the XC rel error at h=1e-3
    {
        // Hartree control at h=1e-3 (bilinear: FD error only).
        const double hc=1e-3;
        double fdH=(EH(shifted(D,dD,+hc))-EH(shifted(D,dD,-hc)))/(2.0*hc);
        ΔG_Map VH=Contract(cou,D);
        rvec_t rho0=rhoOf(D);                            // ALSO leaves the colloc memo (screenD) at D
        chmat_t HH=ev.OverlapMatrix([&](const ivec3_t& dm)->dcmplx
            { auto it=VH.find(dm); return it==VH.end()?dcmplx(0.0):it->second; });
        double anH=trace(HH,dD);
        double relH=std::fabs(fdH-anH)/std::max(std::fabs(anH),1e-30);
        // XC probe at two step sizes (separates FD truncation from a genuine E/H fork).
        double relXc=0.0;
        for (double h : {1e-3, 1e-4})
        {
            double fd=(Exc(rhoOf(shifted(D,dD,+h)))-Exc(rhoOf(shifted(D,dD,-h))))/(2.0*h);
            rvec_t rho=rhoOf(D);                         // reset the memo/screenD to D before the H build
            double an=trace(Hxc(rho),dD);
            double rel=std::fabs(fd-an)/std::max(std::fabs(an),1e-30);
            if (h==1e-3) relXc=rel;
            std::cout << "[xc-consistency " << label << "] h=" << h << "  dE_fd=" << fd
                      << "  Tr(Hxc dD)=" << an << "  rel=" << rel
                      << "   (Hartree control rel=" << relH << ")" << std::endl;
        }
        double rmin=1e300, rmax=-1e300;
        for (size_t q=0;q<rho0.size();q++) { rmin=std::min(rmin,rho0[q]); rmax=std::max(rmax,rho0[q]); }
        std::cout << "[xc-consistency " << label << "] rho range on grid: [" << rmin << ", " << rmax << "]" << std::endl;
        EXPECT_LT(relH, 1e-6) << "Hartree control (bilinear) must agree to FD accuracy (" << label << ")";
        return relXc;
    };

    // Probe 1: PSD-like density -- rho_q > 0 (the converged-SCF regime).
    chmat_t D1(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D1(i,j)=(i==j)?dcmplx(1.0):dcmplx(0.3);
    double rel1=probe(D1,"positive");
    EXPECT_LT(rel1, 1e-5) << "H_xc must be the exact derivative of the discrete E_xc (smooth regime)";

    // Probe 2: indefinite D -- rho_q < 0 over part of the grid (the mixed-density regime; guard consistency).
    chmat_t D2(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++) D2(i,j)=(i==j)?dcmplx((i%2)?-0.6:0.8):dcmplx(0.3);
    double rel2=probe(D2,"indefinite");
    std::cout << "[xc-consistency] indefinite-D rel error = " << rel2
              << " (informational: FD across the rho=0 guard kink is not smooth)" << std::endl;
}

// 0.5(f2) RAW-COLLOCATION XC FEED (doc/GPWPlan).  The raw pair on Overlap3CTensor: applyRaw = rho_DM(r) on
// the integration raster (finest level RAW, others transferred spectrally -- NO ball restriction anywhere),
// applyRawAdjoint = its exact transpose (band-truncate per level -> analytic gather).  Three teeth:
//   (1) CHARGE: Integral(rho_raw) == Tr(D S^G) -- the spectral combine preserves the G=0 content exactly;
//   (2) POSITIVITY: for PSD D, rho_DM = phi^T D phi >= 0 pointwise to screening precision -- THE property
//       the ball-projected rho lacks (Gibbs lobes; the C=8 calibration driver) and the reason f2 exists;
//   (3) ADJOINT: H = applyRawAdjoint(v_xc(rho_raw)) is the FD-exact derivative of the raw discrete E_xc.
TEST(GPW_XC, RawConsistencyFD)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14,{0,0,0});
    cell.AddAtom(14,{0.25,0.25,0.25});
    std::shared_ptr<const Real_BS> mol=MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/6.0, BasisSet::Gaussian::CellImages::Periodic, 2.0,
                BasisSet::PlaneWave::RasterPolicy::AliasFree);   // EXACT-QUADRATURE gate (production default = BallOnly)
    const GPW_Evaluator& ev=gpw;
    const auto& grid=ev.DensityGrid();
    const size_t n=static_cast<const Complex_OIBS&>(gpw).GetNumFunctions();

    qchem::Hamiltonian::SlaterExchange  exch(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation corr;
    const qchem::Hamiltonian::ExFunctional* xcs[2]={&exch,&corr};

    Projector3<dcmplx> ov=ev.Overlap3CTensor();
    ASSERT_TRUE(bool(ov.applyRaw))        << "GPW must realise the raw-raster forward";
    ASSERT_TRUE(bool(ov.applyRawAdjoint)) << "GPW must realise the raw-raster adjoint";

    auto Exc=[&](const rvec_t& rho)->double
    {
        rvec_t e(rho.size());
        for (size_t q=0;q<rho.size();q++)
        {
            double s=0.0; for (auto xc : xcs) s+=xc->GetEpsXc(rho[q]);
            e[q]=s*rho[q];
        }
        return grid.Integral(e);
    };
    auto Hxc=[&](const rvec_t& rho)->chmat_t
    {
        rvec_t v(rho.size());
        for (size_t q=0;q<rho.size();q++)
        {
            double s=0.0; for (auto xc : xcs) s+=xc->GetVxc(rho[q]);
            v[q]=s;
        }
        return ov.applyRawAdjoint(v);
    };
    auto trace=[&](const chmat_t& H, const chmat_t& dD)->double
    {
        dcmplx tr(0.0);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) tr+=dcmplx(H(i,j))*dcmplx(dD(j,i));
        return std::real(tr);
    };
    auto shifted=[&](const chmat_t& D, const chmat_t& dD, double h)->chmat_t
    {
        chmat_t Dh(D);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) Dh(i,j)=dcmplx(D(i,j))+h*dcmplx(dD(i,j));
        return Dh;
    };

    // (1)+(2) on a PSD density (0.7 I + 0.3 J): charge and pointwise non-negativity of the raw feed.
    chmat_t D1(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D1(i,j)=(i==j)?dcmplx(1.0):dcmplx(0.3);
    rvec_t rhoRaw=ov.applyRaw(D1);
    const auto& S=static_cast<const Complex_OIBS&>(gpw).Overlap();
    dcmplx trDSc(0.0);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) trDSc+=dcmplx(D1(i,j))*dcmplx(S(j,i));
    const double trDS=std::real(trDSc);
    EXPECT_NEAR(grid.Integral(rhoRaw), trDS, 5e-3*std::fabs(trDS)) << "raw-feed charge == Tr(D S^G)";
    double rmin=1e300, rmax=-1e300;
    for (size_t q=0;q<rhoRaw.size();q++) { rmin=std::min(rmin,rhoRaw[q]); rmax=std::max(rmax,rhoRaw[q]); }
    rvec_t rhoBall=grid.RhoOnGrid(Contract(ov,D1));
    double bmin=1e300;
    for (size_t q=0;q<rhoBall.size();q++) bmin=std::min(bmin,rhoBall[q]);
    std::cout << "[raw-xc] rho_raw range [" << rmin << ", " << rmax << "]   (ball-path min " << bmin
              << " -- the Gibbs lobes the raw feed removes)" << std::endl;
    EXPECT_GT(rmin, -1e-6*rmax) << "PSD D must give pointwise-non-negative rho_DM (screening-eps only)";

    // (3) FD consistency of the raw pair (mirrors XCPotentialConsistencyFD probe 1).
    chmat_t dD(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++) dD(i,j)=dcmplx(0.1*std::sin(1.0+double(i)+2.0*double(j)),0.0);
    double rel1=0.0;
    for (double h : {1e-3, 1e-4})
    {
        double fd=(Exc(ov.applyRaw(shifted(D1,dD,+h)))-Exc(ov.applyRaw(shifted(D1,dD,-h))))/(2.0*h);
        rvec_t rho=ov.applyRaw(D1);                      // reset the colloc memo/screenD to D1 before H
        double an=trace(Hxc(rho),dD);
        double rel=std::fabs(fd-an)/std::max(std::fabs(an),1e-30);
        if (h==1e-3) rel1=rel;
        std::cout << "[raw-xc] h=" << h << "  dE_fd=" << fd << "  Tr(Hxc dD)=" << an
                  << "  rel=" << rel << std::endl;
    }
    EXPECT_LT(rel1, 1e-5) << "raw H_xc must be the exact derivative of the raw discrete E_xc";
}
