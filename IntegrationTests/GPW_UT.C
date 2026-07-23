// File GPW_UT.C  GPW at Gamma: periodic Gaussian 1-electron integrals + the DFT-tier collocation primitive.
//
// GPW puts GAUSSIAN orbitals on a lattice.  Its one-electron matrices are lattice sums of the ordinary
// (finite) two-centre integrals,  M_ij = Sum_R <chi_i | O | chi_j(.-R)>, computed by the molecular Gaussian
// basis (Molecule::LatticeSum1E) and delegated to by GPW_Evaluator; GPW_IBS is the thin Orbital_1E_IBS on top.
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
import qchem.BasisSet.Molecule.Factory;         // Molecule::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Lattice_3D.GPW_IBS;       // GPW_IBS (the basis under test)
import qchem.Pseudopotential.SeparablePotential; // HGH_SeparablePotential + the _R / _Gaussian faces (KB gate)
import qchem.Pseudopotential.GTH_Potentials;     // GetGTH (the Si GTH-LDA-q4 projector data)
import qchem.BasisSet.Lattice_3D.Evaluators.GPW; // GPW_Evaluator (tests may cheat-import internals) -- DFT tier
import qchem.BasisSet.Molecule.LatticeSum1E;     // Molecule::LatticeSum1E::CollocateDensity (analytic collocation)
import qchem.BasisSet.Internal.GMap;            // G_ERI3 / ΔG_Map (the collocation tensor + rho-tilde)
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
using BasisSet::Lattice_3D::GPW_IBS;
using BasisSet::Lattice_3D::GPW_Evaluator;
using qchem::BasisSet::Molecule::BasisSetData;

namespace
{
// The valence Si Gaussian basis (SIPP, Cartesian) on ANY structure -- the L_PP builder.  The Engine argument
// is the integral-engine switch point (see Molecule::LatticeSum1E): Engine::MnD here (its AtCenter + analytic
// 2C kernels make the periodic sum exact/trivial); Engine::LibCint would be the faster path once PG_LibCint
// realises LatticeSum1E -- GPW itself is unchanged either way.
std::unique_ptr<Real_BS> MakeBasis(const Structure& st)
{
    return std::unique_ptr<Real_BS>(
        BasisSet::Molecule::Factory(BasisSetData::SIPP, &st,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
}

// The single orbital block of a raw (no-SALC) single-atom basis (real or complex flavour).
template <class OIBS> const OIBS* OrbitalBlock(const Real_BS& bs)
{
    const OIBS* only=nullptr; int count=0;
    for (auto ibs : bs.Iterate<OIBS>()) { only=ibs; count++; }
    EXPECT_EQ(count,1) << "expected one orbital block for the raw single-atom basis";
    return only;
}

// Relative Frobenius distance  ||Re(A) - B||_F / ||B||_F  (A complex GPW matrix, B the real finite matrix).
template <class CM, class RM> double RelDiff(const CM& A, const RM& B)
{
    const size_t n=B.rows();
    double num=0.0, den=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            double d=std::real(A(i,j))-B(i,j);
            num+=d*d; den+=B(i,j)*B(i,j);
        }
    return std::sqrt(num/den);
}

// Largest |Im| over a complex matrix (must be ~0 at Gamma).
template <class CM> double MaxImag(const CM& A)
{
    const size_t n=A.rows();
    double m=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++) m=std::max(m, std::fabs(std::imag(A(i,j))));
    return m;
}

// The finite-Si reference matrices (atom at the origin): overlap, <p^2>, nuclear.  Translation-invariant, so
// they equal the GPW home-cell (R=0) matrices for the same atom centred anywhere in the cell.
struct FiniteRef
{
    std::shared_ptr<Molecule>   mol;
    std::unique_ptr<Real_BS>    bs;
    const Real_OIBS*            orb;
    FiniteRef()
        : mol(std::make_shared<Molecule>())
    {
        mol->Insert(new Atom(14,{0,0,0}));
        bs = MakeBasis(*mol);
        orb= OrbitalBlock<Real_OIBS>(*bs);
    }
};
} //anon

// Tooth (1): home cell only -- GPW reproduces the finite matrices EXACTLY (R=0 IS the finite integral).
TEST(GPW, HomeCellMatchesFiniteExactly)
{
    FiniteRef fin;

    const double a=20.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);

    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/0.0,
                BasisSet::Lattice_3D::CellImages::HomeCellOnly); // the finite-molecule mode
    const Complex_OIBS& g = gpw;

    ASSERT_EQ(g.GetNumFunctions(), fin.orb->GetNumFunctions());
    EXPECT_LT(MaxImag(g.Overlap()),        1e-14);
    EXPECT_LT(MaxImag(g.Kinetic()),        1e-14);
    EXPECT_LT(MaxImag(g.Nuclear(&cell)),   1e-13);
    EXPECT_LT(RelDiff(g.Overlap(),      fin.orb->Overlap()),        1e-12);
    EXPECT_LT(RelDiff(g.Kinetic(),      fin.orb->Kinetic()),        1e-12);
    EXPECT_LT(RelDiff(g.Nuclear(&cell), fin.orb->Nuclear(fin.mol.get())), 1e-12);
}

// Tooth (2): sum the periodic images in.  GPW matches the finite matrices to the image tail, and the tail
// shrinks as the cell grows (the large-cell limit).  Overlap / kinetic are compact, so the tail is tiny.
TEST(GPW, LatticeSumConvergesToFiniteAsCellGrows)
{
    FiniteRef fin;

    auto overlapResidual = [&](double a)->double
    {
        UnitCell cell(a);
        cell.AddAtom(14,{0.5,0.5,0.5});
        std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/0.0); // periodic (images internal)
        const Complex_OIBS& g = gpw;
        EXPECT_LT(MaxImag(g.Overlap()), 1e-13);          // still real at Gamma
        return RelDiff(g.Overlap(), fin.orb->Overlap());
    };

    // Clean exponential large-cell convergence (measured: ~two orders of magnitude per +4 a.u.):
    //   a=14 -> 1.5e-2,  a=22 -> 7.9e-6,  a=30 -> 5.9e-11.
    const double r14 = overlapResidual(14.0);
    const double r22 = overlapResidual(22.0);
    const double r30 = overlapResidual(30.0);

    EXPECT_GT(r14, r22);          // the image tail shrinks monotonically as the cell grows
    EXPECT_GT(r22, r30);
    EXPECT_LT(r30, 1e-9);         // and vanishes in the large-cell limit (GPW -> the finite molecule)
}

// === DFT tier: the collocation primitive =============================================================
// GPW's genuinely-new machinery is COLLOCATION: build rho(r)=sum D chi chi on a real grid, FFT to rho-tilde,
// Poisson for Hartree, integrate a grid potential back against the Gaussians for the KS matrix.  Two clean
// grid-convergent checks isolate the collocate + integrate-back primitives against the analytic overlap,
// without needing a full SCF / pseudopotential / G=0 background.
namespace
{
// The G=0 column of a collocation tensor (the fit function G_c = 0).
int G0Column(const G_ERI3& t)
{
    for (size_t c=0;c<t.columns.size();c++)
    {
        const ivec3_t& dm=t.columns[c].dm;
        if (dm.x==0 && dm.y==0 && dm.z==0) return int(c);
    }
    return -1;
}
} //anon

// Collocate (through the SCF seam): the G=0 component of the tensor's matrix-free `apply` map is the collocated
// charge, apply(D)[0]*Omega = Integral rho = Tr(D S), which must match the analytic 1E overlap trace as the
// density grid resolves the Gaussian products.  This exercises exactly what ContractG_ERI3 runs in the SCF.
TEST(GPW, CollocationOverlapMatchesAnalytic)
{
    const double a=10.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/30.0);
    // REFERENCE: the analytic collocation always sums the SCREENED cross-cell pair offsets, so the collocated
    // charge is Tr(D S^Bloch) -- the screened-complete Bloch overlap (generous Rcut enumeration; SIPP's diffuse
    // alpha=0.06 reaches several cells even in this box), NOT the home-only overlap.
    GPW_IBS gpwRef(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/0.0);

    const GPW_Evaluator& ev = gpw;
    G_ERI3 ov = ev.Overlap3CTensor();
    ASSERT_GE(G0Column(ov), 0);
    ASSERT_TRUE(bool(ov.apply)) << "the GPW tensor must carry the matrix-free analytic-collocation map";

    const auto& S = static_cast<const Complex_OIBS&>(gpwRef).Overlap();   // Bloch overlap S^Bloch (real at Gamma)
    const size_t n = S.rows();
    chmat_t D(n);                    // D = identity -> Integral rho = Tr(S^Bloch)
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=(i==j)?dcmplx(1.0):dcmplx(0.0);
    ΔG_Map rho = ContractG_ERI3(ov, D);                      // the SCF's own contraction (dispatches to apply)
    const double integral = std::real(rho[ivec3_t(0,0,0)]) * ov.volume;
    double trS=0.0; for (size_t i=0;i<n;i++) trS += std::real(dcmplx(S(i,i)));
    EXPECT_NEAR(integral, trS, 6e-2*std::fabs(trS));   // collocated charge == Bloch overlap trace to grid accuracy
}

// SHARPEST-PAIR charge conservation (doc/GPWPlan.md 0b').  MEASURED FINDING this test PINS: the collocated
// CHARGE of even the sharpest pair (the alpha_max+alpha_max product, whose pair->level requirement exceeds
// the reference grid) is exact to ~1e-9 WITH OR WITHOUT the top completion rung -- the G=0 coefficient
// survives the per-level BALL truncation by construction (only G>ball content is discarded), and the
// pow2-padded rasters keep the box SAMPLING error at ~e^{-50}.  So ball truncation is an ENERGY-tail
// effect (e^{-ecut/2p} coupling to V_H/v_xc), NEVER a charge leak -- and any e-scale grid-charge loss in
// an SCF must come from elsewhere (see DISABLED_IllConditionedChargeProbe: enumeration-scheme mismatch).
// The gate loads exactly that pair (unit D on the sharpest function, found via the kinetic diagonal
// <p^2> ~ alpha -- no exponent crosses the interface) and pins its collocated charge tight.
TEST(GPW, SharpestPairChargeConservation)
{
    const double a=10.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut AUTO*/-1.0);
    const GPW_Evaluator& ev = gpw;
    const Complex_OIBS&  g  = gpw;
    const size_t n=g.GetNumFunctions();

    const auto& T=g.Kinetic();                 // <p^2>_ii ~ alpha_i: the max diagonal marks the sharpest function
    size_t is=0;
    for (size_t i=1;i<n;i++) if (std::real(dcmplx(T(i,i)))>std::real(dcmplx(T(is,is)))) is=i;

    chmat_t D(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=dcmplx(0.0);
    D(is,is)=dcmplx(1.0);                      // unit load on the sharpest product -- the top rung's customer

    G_ERI3 ov=ev.Overlap3CTensor();
    ASSERT_TRUE(bool(ov.apply));
    ΔG_Map rho=ContractG_ERI3(ov, D);
    const double integral=std::real(rho[ivec3_t(0,0,0)])*ov.volume;
    const auto& S=g.Overlap();
    const double ref=std::real(dcmplx(S(is,is)));   // = 1 (normalized; the tight function's images are ~e^{-alpha a^2})
    std::cout << "[sharpest pair] i*=" << is << "  <p^2>_ii=" << std::real(dcmplx(T(is,is)))
              << "  Integral rho=" << integral << "  S_ii=" << ref
              << "  rel=" << std::fabs(integral-ref)/ref << std::endl;
    EXPECT_NEAR(integral, ref, 1e-4*ref) << "sharpest-pair collocated charge (the 0b' top-rung gate)";
}

// ILL-CONDITIONED-LOAD charge probe (doc/GPWPlan.md 0b' investigation) -- THE instrument that root-caused
// the NaF iteration-1 grid-charge loss.  D = S^-1 is the clean probe: PSD, entries ~1/lambda_min (the
// loading a mid-slosh SCF on a near-singular basis produces), and Tr(D S) = n EXACTLY, so
// [Integral rho - n] measures the |D|-amplified collocation error directly.  MEASURED (2026-07-16):
//   - Rcut=2a (the NaF SCF's setting): err = -2.247 e at max|D|=450 -- GRID-INDEPENDENT (identical at
//     Ecut=40 and auto=160, across different fp32 tiering) => an ANALYTIC mismatch, not noise: the
//     collocation enumerates its cross-cell offsets INTERNALLY to the complete magnitude screen (pair
//     reach ~33 au for VALENCE_LOWQ_SR), while S is built over the caller's Rcut=2a=17.5 au images --
//     the "two self-consistent schemes" pin violated by the config.  Mid-slosh D loads exactly the
//     near-null (diffuse) directions where the truncated S is most wrong -> the e-scale charge swings.
//   - Rcut=AUTO (complete): err/|D| drops ~15000x; the complete-enumeration S is however genuinely
//     near-singular (lambda_min ~ 1e-6, |S^-1| ~ 1e6 -- the 2a truncation was doubling as a conditioning
//     crutch).  Residual at |D|~1e6: -0.361 pure-fp64 vs -0.355 with fp32 streams => the fp32 tier
//     contributes ~7e-3 and the kScreenEps screening tails ~0.36 EVEN at million-scale loading -- the
//     precision machinery is vindicated; per-term floors hold.
// Knobs: GPW_ILLCOND_ECUT (-1 = production auto), GPW_ILLCOND_RCUT (-1 = AUTO; default 2a = the SCF's).
// DISABLED: minutes-long stream builds; run explicitly when investigating conditioning/enumeration.
TEST(GPW, DISABLED_IllConditionedChargeProbe)
{
    const double a=8.73;                                     // the NaF rocksalt cell (GPW_SCF_UT's)
    FCCUnitCell cell(a);
    cell.AddAtom(11,{0,0,0});
    cell.AddAtom(9, {0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> mol(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ_SR, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    const char* e=std::getenv("GPW_ILLCOND_ECUT");
    const double ecut=e?std::atof(e):40.0;
    // NOTE (banish-Rcut): the historical Rcut=2a leg that MEASURED the -2.247 e scheme mismatch is now
    // UNREPRESENTABLE -- enumeration lives inside the seam and no truncated configuration can be built.
    // The probe now reports the complete-enumeration error only (precision floors under extreme loading).
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, ecut);
    const GPW_Evaluator& ev=gpw;
    const Complex_OIBS&  g =gpw;
    const size_t n=g.GetNumFunctions();

    // D = S^-1 (real at Gamma): PSD with ~1/lambda_min entries; Tr(D S) = n exactly.
    const auto& S=g.Overlap();
    rmat_t Sr(n,n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) Sr(i,j)=std::real(dcmplx(S(i,j)));
    rmat_t Si=blazem::inv(Sr);
    double dmax=0.0;
    chmat_t D(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++) { D(i,j)=dcmplx(0.5*(Si(i,j)+Si(j,i)),0.0); dmax=std::max(dmax,std::fabs(Si(i,j))); }

    G_ERI3 ov=ev.Overlap3CTensor();
    ASSERT_TRUE(bool(ov.apply));
    ΔG_Map rho=ContractG_ERI3(ov, D);
    const double integral=std::real(rho[ivec3_t(0,0,0)])*ov.volume;
    std::cout << "[illcond probe] Ecut=" << ecut << "  n=" << n << "  max|D|=" << dmax
              << "  Integral rho=" << integral << "  Tr(D S)=" << double(n)
              << "  err=" << integral-double(n) << std::endl;
}

// Integrate-back: a CONSTANT potential V(r)=V0 (Vtilde nonzero only at G=0) must give <chi_i|V0|chi_j> = V0 S_ij
// (grid-quadrature).  This validates OverlapMatrix (the RhoOnGrid inverse-FFT + the grid quadrature adjoint).
TEST(GPW, OverlapWithConstantFieldEqualsV0Overlap)
{
    const double a=10.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/30.0);
    // REFERENCE: the analytic integrate-back sums the screened cross-cell offsets, so a constant field gives
    // V0 * S^Bloch (screened-complete Bloch overlap), not V0 * S_home -- see CollocationOverlapMatchesAnalytic.
    GPW_IBS gpwRef(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/0.0);

    const GPW_Evaluator& ev = gpw;
    const double V0 = 0.7;
    const chmat_t M = ev.OverlapMatrix([V0](const ivec3_t& dm)->dcmplx
        { return (dm.x==0 && dm.y==0 && dm.z==0) ? dcmplx(V0) : dcmplx(0.0); });

    const auto& S = static_cast<const Complex_OIBS&>(gpwRef).Overlap();
    const size_t n = S.rows();
    double num=0.0, den=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            double m = std::real(M(i,j));
            double s = V0*std::real(S(i,j));
            double d=m-s; num+=d*d; den+=s*s;
        }
    EXPECT_LT(std::sqrt(num/den), 6e-2);   // <i|V0|j> == V0<i|j> to grid accuracy (== the collocation residual)
}

// ANALYTIC COLLOCATION charge conservation (GPWPlan.md S0 Increment A -- the CP2K rewrite).
// LatticeSum1E::CollocateDensity collocates rho = Sum_ij D_ij chi_i chi_j analytically per pair on compact
// exp-tail boxes, modulo-wrapped onto the grid (NO image sum, NO Rcut).  The defining invariant: Integral of
// the collocated rho over the cell equals Tr(D S) (the density integrates to the electron count) -- to grid
// tolerance.  Tested with the atom in the INTERIOR (no wrap) AND at the CORNER {0,0,0} (the box wraps around
// every face -> exercises the modulo-wrap-IS-the-image-sum mechanism; charge must be identical, translation-
// invariant, with no ringing).  Kernel-level, K=1 (single grid); the multi-grid ladder is exercised by the
// seam-level gates below.
TEST(GPW, AnalyticCollocationConservesCharge)
{
    auto probe=[](const rvec3_t& frac, const char* where)->double
    {
        const double a=12.0;
        UnitCell cell(a);
        cell.AddAtom(14, frac);                                     // Si
        std::shared_ptr<const Real_BS> mol = MakeBasis(cell);       // SIPP Si
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/12.0);
        // REFERENCE: Tr(D S^Bloch) -- the collocation always includes the screened cross-cell pair offsets
        // (SIPP's diffuse alpha=0.06 reaches neighbour cells even at a=12), so the home-only Tr(D S) is ~3% off.
        GPW_IBS gpwRef(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/0.0);
        const GPW_Evaluator& ev=gpw;
        const auto* lat=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(OrbitalBlock<Real_OIBS>(*mol));
        EXPECT_TRUE(lat) << "orbital block must realise LatticeSum1E";
        const ivec3_t N=ev.DensityGrid().FFTGrid();
        const size_t  n=ev.size();
        chmat_t D(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=(i==j)?dcmplx(1.0):dcmplx(0.0);
        auto gamma=[](const ivec3_t&)->dcmplx { return dcmplx(1.0); };           // Gamma: every offset phase 1
        rvec_t rho=lat->CollocateDensity(D, gamma, cell, {N}, {ev.DensityGrid().Ecut()})[0];   // K=1
        double integral=blazem::sum(rho)*cell.GetCellVolume()/double(rho.size());
        const auto& S=static_cast<const Complex_OIBS&>(gpwRef).Overlap();        // Bloch overlap S^Bloch
        double trDS=0.0; for (size_t i=0;i<n;i++) trDS+=std::real(dcmplx(S(i,i)));
        std::cout << "[collocate " << where << "] Integral rho=" << integral << "  Tr(D S^Bloch)=" << trDS
                  << "  rel=" << std::fabs(integral-trDS)/std::fabs(trDS) << std::endl;
        EXPECT_NEAR(integral, trDS, 5e-3*std::fabs(trDS)) << "collocated charge vs Tr(D S^Bloch) (" << where << ")";
        return integral;
    };
    double cInterior=probe({0.5,0.5,0.5}, "interior");
    double cCorner  =probe({0.0,0.0,0.0}, "corner-wrapped");
    EXPECT_NEAR(cInterior, cCorner, 1e-6) << "collocated charge must be translation-invariant (wrap == interior)";
}

namespace
{
// Scoped env override with restore (the stream-budget knobs are read per EnsureStreams call).
struct EnvGuard
{
    std::string name, old; bool had;
    EnvGuard(const char* n, const std::string& v) : name(n)
    {   const char* o=std::getenv(n); had=(o!=nullptr); if (o) old=o; setenv(n, v.c_str(), 1); }
    void Set(const std::string& v) { setenv(name.c_str(), v.c_str(), 1); }
    ~EnvGuard() { if (had) setenv(name.c_str(), old.c_str(), 1); else unsetenv(name.c_str()); }
};
// The last "[stream cache]" build readout in a captured stderr blob.  built=false means EnsureStreams HIT an
// existing cache (no build ran) -- itself an assertable outcome (rebuild stability).
struct StreamReadout { bool built=false; size_t pts64=0, pts32=0, dropped=0; };
StreamReadout ParseStreamReadout(const std::string& err)
{
    StreamReadout r;
    size_t p=err.rfind("[stream cache]");
    if (p==std::string::npos) return r;
    r.built=true;
    auto ptsInParens=[&](const std::string& tag)->size_t     // "fp64 N (P pts)" -> P
    {
        size_t b=err.find(tag, p);
        if (b==std::string::npos) { ADD_FAILURE() << "readout tag missing: " << tag; return 0; }
        return std::stoull(err.substr(err.find('(', b)+1));
    };
    r.pts64=ptsInParens("fp64 "); r.pts32=ptsInParens("fp32 ");
    size_t d=err.find("dropped ", p);
    if (d==std::string::npos) { ADD_FAILURE() << "readout tag missing: dropped"; return r; }
    r.dropped=std::stoull(err.substr(d+8));
    return r;
}
} //anon

// STREAM-CACHE RESIDENCY (doc/GPWPlan.md 0.5(b)).  The pair-box streams live on the SHARED molecular
// evaluator keyed by ladder shape, with a GLOBAL point budget -- so in a grid-continuation run the RESIDENT
// coarse-stage caches starve the fine stage to ~0% coverage (measured: the 8.45-h NaF full-SR diagnostic,
// billions of points re-evaluated per iteration).  The fix under test, both halves:
//   (1) RELEASE: destroying a GPW block (bsC.reset() in the SCF test) hands its ladder's streams back to the
//       budget (GPW_Evaluator dtor -> LatticeSum1E::ReleaseStreams);
//   (2) SELF-HEAL: a cache built STARVED (the fine shape is built during the seed handoff, while the coarse
//       stage still squats) rebuilds when the headroom grows -- and ONLY then (a complete cache never
//       rebuilds: bit-stable replay; an unchanged starved cache never churns).
TEST(GPW, StreamCacheReleaseUnstarvesLaterGrid)
{
    const double a=12.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});                            // Si, interior
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);      // ONE molecular basis shared by both grids
    const double ecutC=16.0, ecutF=32.0;                       // coarse/fine density grids (>= the 8*alpha_max floor)

    // One OverlapMatrix(V) call drives IntegratePotential -> EnsureStreams for the block's ladder shape.
    // A FRESH constant per call defeats the static-field IntegrateMemo (same shape + same field replays the
    // memoised reductions and never consults the stream cache), so every trigger reaches EnsureStreams.
    auto trigger=[call=0](const GPW_Evaluator& ev) mutable -> StreamReadout
    {
        const dcmplx v(double(++call));
        testing::internal::CaptureStderr();
        ev.OverlapMatrix([v](const ivec3_t&)->dcmplx { return v; });
        return ParseStreamReadout(testing::internal::GetCapturedStderr());
    };

    // MEASURE the two shapes' demand under an ample single-tier budget (fp32 tier off: one number per shape).
    EnvGuard b64("GPW_STREAM_BUDGET_PTS",     "2000000000");
    EnvGuard b32("GPW_STREAM_BUDGET_PTS_F32", "0");
    size_t ptsC=0, ptsF=0;
    {
        GPW_IBS c(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, ecutC);
        auto r=trigger(c);
        ASSERT_TRUE(r.built); ASSERT_EQ(r.dropped,0u); ptsC=r.pts64+r.pts32;
    }   // dtor releases the coarse shape
    {
        GPW_IBS f(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, ecutF);
        auto r=trigger(f);
        ASSERT_TRUE(r.built) << "coarse dtor must have released its shape (fresh build expected)";
        ASSERT_EQ(r.dropped,0u); ptsF=r.pts64+r.pts32;
    }
    ASSERT_GE(ptsC,2u);
    ASSERT_GT(ptsF,ptsC) << "the finer grid must demand more stream points";

    // THE GRID-CONTINUATION CONFIGURATION under a budget that fits EITHER shape but not BOTH:
    b64.Set(std::to_string(ptsF + ptsC/2));
    auto coarse=std::make_unique<GPW_IBS>(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, ecutC);
    auto r=trigger(*coarse);
    ASSERT_TRUE(r.built); EXPECT_EQ(r.dropped,0u);             // coarse fits alone
    GPW_IBS fine(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, ecutF);
    r=trigger(fine);
    ASSERT_TRUE(r.built); EXPECT_GT(r.dropped,0u)              // STARVED by the resident coarse caches
        << "expected the fine shape to be starved while the coarse caches squat on the budget";
    r=trigger(fine);
    EXPECT_FALSE(r.built) << "starved cache with UNCHANGED headroom must not churn/rebuild";
    coarse.reset();                                            // the bsC.reset() fix: refund the coarse points
    r=trigger(fine);
    ASSERT_TRUE(r.built); EXPECT_EQ(r.dropped,0u)              // SELF-HEAL: rebuilt into the refunded budget
        << "released coarse budget must reach the starved fine shape";
    r=trigger(fine);
    EXPECT_FALSE(r.built) << "complete cache must never rebuild (bit-stable replay)";
}

// ANALYTIC COLLOCATION on a CRYSTAL (cross-cell pairs), through the SCF SEAM (Overlap3CTensor's matrix-free
// `apply` -> the REL_CUTOFF multi-grid ladder).  The periodic Gamma density is a product of BLOCH orbitals,
// chi_i^G chi_j^G = Sum_R'' chi_i^0 chi_j^R'' -- so collocation must sum the screened CROSS-CELL offsets, not
// just the home pair (R''=0).  Invariant: Omega x the G=0 component of apply(D) == Integral rho == Tr(D S^G)
// with S^G the Bloch overlap (its own screened image sum, via a generous Rcut).  A 2-atom Si crystal with real
// inter-cell overlap (which the single-atom gate above cannot exercise).
TEST(GPW, AnalyticCollocationCrystalChargeConservation)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);              // SIPP Si (one orbital block over both atoms)
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/6.0);   // periodic: S^G by construction
    const GPW_Evaluator& ev=gpw;
    const Complex_OIBS& g=gpw;
    const size_t n=g.GetNumFunctions();
    chmat_t D(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=(i==j)?dcmplx(1.0):dcmplx(0.0);
    G_ERI3 ov=ev.Overlap3CTensor();
    ASSERT_TRUE(bool(ov.apply));
    ΔG_Map rho=ContractG_ERI3(ov, D);                                  // the SCF's own contraction (multi-grid)
    const double integral=std::real(rho[ivec3_t(0,0,0)])*ov.volume;
    const auto& S=g.Overlap();                                         // Bloch overlap (screened images) = S^G
    double trDS=0.0; for (size_t i=0;i<n;i++) trDS+=std::real(dcmplx(S(i,i)));
    std::cout << "[collocate crystal] Integral rho=" << integral << "  Tr(D S^G)=" << trDS
              << "  rel=" << std::fabs(integral-trDS)/std::fabs(trDS) << std::endl;
    EXPECT_NEAR(integral, trDS, 5e-2*std::fabs(trDS)) << "crystal collocated charge vs Tr(D S^G)";
}

// ANALYTIC INTEGRATE-BACK (GPWPlan.md S0 Increment B): IntegratePotential is the exact adjoint of
// CollocateDensity (same screened offsets, boxes, wrap, level assignment).  Two gates:
//   (1) KERNEL-level (K=1): <collocate(D),V> == Tr(D h) to machine precision (variational -- the KS matrix is
//       the exact gradient of the grid energy);
//   (2) SEAM-level, MULTI-GRID: Tr(D OverlapMatrix(Vtilde)) == Omega Sum_G apply(D)(G) Vtilde(G) -- the same
//       identity through the level ladder: per level the grid quadrature is Parseval-exact against the
//       spectrally-restricted V_l, and the nested G-space combine matches the per-pair level assignment, so
//       the multigrid density and multigrid KS matrix are exact adjoints too.
TEST(GPW, AnalyticIntegrateBackAdjoint)
{
    const double a=12.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/12.0);
    const GPW_Evaluator& ev=gpw;
    const auto* lat=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(OrbitalBlock<Real_OIBS>(*mol));
    EXPECT_TRUE(lat);
    const ivec3_t N=ev.DensityGrid().FFTGrid();
    const size_t  n=ev.size();
    auto gamma=[](const ivec3_t&)->dcmplx { return dcmplx(1.0); };
    auto field=[](const ivec3_t& dm)->dcmplx        // smooth all-G field (Hartree/XC-like), symmetric in G->-G
        { double g2=double(dm.x*dm.x+dm.y*dm.y+dm.z*dm.z); return dcmplx(1.0/(1.0+g2),0.0); };
    chmat_t D(n);                                   // Hermitian test density
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=(i==j)?dcmplx(1.0):dcmplx(0.3);

    // (1) kernel-level adjoint at K=1: <collocate(D),V> == Tr(D h) to machine precision
    ΔG_Map vmap; for (const ivec3_t& dm : ev.DensityGrid().Gs()) vmap[dm]=field(dm);
    rvec_t V=ev.DensityGrid().RhoOnGrid(vmap);                 // V(r) on the fine grid
    const double ecut=ev.DensityGrid().Ecut();
    chmat_t h=lat->IntegratePotential({V}, gamma, cell, {N}, {ecut});
    rvec_t rho=lat->CollocateDensity(D, gamma, cell, {N}, {ecut})[0];
    const double w=cell.GetCellVolume()/double(V.size());
    double lhs=0.0; for (size_t p=0;p<V.size();p++) lhs+=rho[p]*V[p]*w;
    dcmplx rhs(0.0); for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) rhs+=dcmplx(D(i,j))*dcmplx(h(j,i));
    EXPECT_NEAR(lhs, std::real(rhs), 1e-8*std::fabs(lhs)) << "adjoint <collocate(D),V> == Tr(D h)";
    EXPECT_LT(std::fabs(std::imag(rhs)), 1e-10) << "Tr(D h) must be real (both Hermitian)";

    // (2) seam-level, MULTI-GRID: Tr(D OverlapMatrix(field)) == Omega Sum_G apply(D)(G) field(G)
    chmat_t H=ev.OverlapMatrix(field);                         // the SCF KS bridge (multi-grid analytic)
    dcmplx trDH(0.0); for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) trDH+=dcmplx(D(i,j))*dcmplx(H(j,i));
    G_ERI3 ovt=ev.Overlap3CTensor();
    ΔG_Map rhoT=ContractG_ERI3(ovt, D);                        // the SCF density map (multi-grid, nested)
    dcmplx eG(0.0); for (const auto& [dm,c] : rhoT) eG+=c*field(dm);
    eG*=ovt.volume;
    std::cout << "[integrate-back] kernel adjoint lhs=" << lhs << " rhs=" << std::real(rhs)
              << "   seam adjoint Tr(D H)=" << std::real(trDH) << " Omega Sum rho.V=" << std::real(eG) << std::endl;
    EXPECT_NEAR(std::real(trDH), std::real(eG), 1e-8*std::fabs(std::real(trDH)))
        << "multi-grid seam adjoint: Tr(D MakeOverlap(V)) == <apply(D),V>";
}

// XC POTENTIAL-CONSISTENCY PROBE (doc/GPWPlan.md 0b instrument).  Question under test: is the assembled
// H_xc the EXACT D-derivative of the DISCRETE energy E_xc(D) = Sum_q w [eps_x+eps_c](rho_q) rho_q, where
// rho_q is the ball-limited grid density of the SCF's own chain?  The probe replicates the PW_XC term's
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
TEST(GPW, XCPotentialConsistencyFD)
{
    const double a=10.26;
    FCCUnitCell cell(a);                                 // the production shape: cross-cell pairs + ladder
    cell.AddAtom(14,{0,0,0});
    cell.AddAtom(14,{0.25,0.25,0.25});
    std::shared_ptr<const Real_BS> mol=MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/6.0);
    const GPW_Evaluator& ev=gpw;
    const auto& grid=ev.DensityGrid();
    const size_t n=static_cast<const Complex_OIBS&>(gpw).GetNumFunctions();

    // The SCF's own functionals (Ham_PW_DFT::BuildTerms builds exactly these two PW_XC terms).
    qchem::Hamiltonian::SlaterExchange  exch(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation corr;
    const qchem::Hamiltonian::ExFunctional* xcs[2]={&exch,&corr};

    G_ERI3 ov =ev.Overlap3CTensor();                     // rho-tilde (no kernel) -- the PW_XC route
    G_ERI3 cou=ev.Repulsion3CTensor();                   // V_H (Coulomb kernel baked) -- the control

    auto rhoOf=[&](const chmat_t& D)->rvec_t { return grid.RhoOnGrid(ContractG_ERI3(ov,D)); };
    auto Exc=[&](const rvec_t& rho)->double              // == PW_XC::GetEnergy (both terms)
    {
        rvec_t e(rho.size());
        for (size_t q=0;q<rho.size();q++)
        {
            double s=0.0; for (auto xc : xcs) s+=xc->GetEpsXc(rho[q]);
            e[q]=s*rho[q];
        }
        return grid.Integral(e);
    };
    auto Hxc=[&](const rvec_t& rho)->chmat_t             // == PW_XC::CalcMatrix (both terms summed)
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
        ΔG_Map VH=ContractG_ERI3(cou,D);
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
        ΔG_Map VH=ContractG_ERI3(cou,D);
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

// === General-k (Step 1): the Bloch phase e^{ik.R} enters the lattice sums ============================
// These isolate the general-k machinery at the matrix level (no SCF): the phase is inert at the home cell,
// LIVE once images are summed, obeys the Bloch translation law, and conjugates under k -> -k.  (The full
// general-k DFT/PP path is exercised by the multi-k bulk SCF in GPW_SCF_UT.)

// (a) Home cell only (Rcut=0): the phase multiplies only the origin (=1), so the matrices are k-INDEPENDENT
// and real -- identical to the finite molecule at ANY k (the Gamma-reduction invariant, lifted to k!=0).
TEST(GPW, GeneralK_HomeCellIsKInvariantAndReal)
{
    FiniteRef fin;
    const double a=20.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);

    GPW_IBS gpw(cell, ivec3_t(4,4,4), ivec3_t(1,0,0), molCell, /*densityEcut=*/0.0,
                BasisSet::Lattice_3D::CellImages::HomeCellOnly); // k=(1/4,0,0), finite mode
    const Complex_OIBS& g = gpw;
    EXPECT_LT(MaxImag(g.Overlap()), 1e-14);                        // phase inert at R=0 -> real
    EXPECT_LT(MaxImag(g.Kinetic()), 1e-14);
    EXPECT_LT(RelDiff(g.Overlap(), fin.orb->Overlap()), 1e-12);    // == finite (k-independent at Rcut=0)
    EXPECT_LT(RelDiff(g.Kinetic(), fin.orb->Kinetic()), 1e-12);
}

// (b) With images summed in (Rcut>0), the phase is LIVE at k!=0: the Bloch matrices acquire a genuine
// imaginary part (they are real only at Gamma or with no images).  Diagonals stay real & positive.
TEST(GPW, GeneralK_PhaseIsLiveWithImages)
{
    const double a=8.0;                          // small cell so the image overlaps are non-negligible
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);

    GPW_IBS gpw(cell, ivec3_t(4,4,4), ivec3_t(1,0,0), molCell, /*densityEcut=*/0.0);
    const Complex_OIBS& g = gpw;
    const auto& S = g.Overlap();
    EXPECT_GT(MaxImag(S),           1e-4) << "k!=0 with images must give a genuinely complex overlap";
    EXPECT_GT(MaxImag(g.Kinetic()), 1e-4);
    for (size_t i=0;i<S.rows();i++)
    {
        EXPECT_LT(std::fabs(std::imag(S(i,i))), 1e-13);   // Hermitian diagonal real
        EXPECT_GT(std::real(S(i,i)),            0.0);     // overlap diagonal positive
    }
}

// (c) Bloch translation law: chi^k(r + R0) = e^{ik.R0} chi^k(r) for a lattice vector R0 (exact in the
// infinite sum; the truncation error decays with Rcut, so a generous sphere makes it tight).
TEST(GPW, GeneralK_BlochTranslationCondition)
{
    const double a=8.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);

    const ivec3_t N(4,4,4), ik(1,0,0);
    GPW_IBS gpw(cell, N, ik, molCell, /*densityEcut=*/0.0);
    const GPW_Evaluator& ev = gpw;

    const rvec3_t r  = cell.ToCartesian(rvec3_t(0.5,0.5,0.5));               // cell centre
    const ivec3_t n0(1,0,0);
    const rvec3_t R0 = cell.ToCartesian(rvec3_t(double(n0.x),double(n0.y),double(n0.z)));
    const double  kR0 = 2.0*Pi*(double(ik.x)/N.x*n0.x + double(ik.y)/N.y*n0.y + double(ik.z)/N.z*n0.z);
    const dcmplx  phase = std::exp(dcmplx(0.0,kR0));

    cvec_t f0 = ev.Eval(r);
    cvec_t f1 = ev.Eval(r+R0);
    double num=0.0, den=0.0;
    for (size_t i=0;i<f0.size();i++)
    {
        dcmplx d = f1[i] - phase*f0[i];
        num += std::norm(d); den += std::norm(f0[i]);
    }
    EXPECT_LT(std::sqrt(num/den), 1e-6) << "Bloch translation law chi^k(r+R0)=e^{ik.R0}chi^k(r)";
}

// (d) Time reversal at the matrix level: S(-k) = conj(S(k)) elementwise -- the phases conjugate under
// k -> -k while the underlying real 2-centre integrals are k-independent.  Exact (same image set, no tol).
TEST(GPW, GeneralK_ConjugateUnderKtoMinusK)
{
    const double a=8.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);

    const ivec3_t N(4,4,4);
    GPW_IBS gk (cell, N, ivec3_t(1,0,0), molCell, 0.0);   //  k = (1/4,0,0)
    GPW_IBS gmk(cell, N, ivec3_t(3,0,0), molCell, 0.0);   // -k = (3/4,0,0) == -(1/4,0,0) mod 1
    const Complex_OIBS& Sk  = gk;
    const Complex_OIBS& Smk = gmk;
    const auto& A = Sk.Overlap();
    const auto& B = Smk.Overlap();
    const size_t n=A.rows();
    double m=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
            m=std::max(m, std::abs(B(i,j)-std::conj(A(i,j))));
    EXPECT_LT(m, 1e-12) << "S(-k) must equal conj(S(k))";
}

// ANALYTIC KB == MESH KB.  The GTH projector is polynomial x Gaussian, so when the model exposes its closed
// Gaussian form (SeparablePotential_Gaussian) MakeSeparablePP assembles <chi_i^k|beta Y_lm> ANALYTICALLY via
// the molecular <chi|g> lattice-sum seam -- no mesh, exact.  This gate pins the analytic path against the
// legacy mesh quadrature it replaced: a wrapper hiding the Gaussian face forces the mesh path on the SAME
// model, and the two matrices must agree to the MESH's own quadrature error (the analytic one is exact).
namespace
{
class MeshOnlyKB : public Pseudopotential::SeparablePotential, public virtual Pseudopotential::SeparablePotential_R
{
    const Pseudopotential::HGH_SeparablePotential& h;
public:
    explicit MeshOnlyKB(const Pseudopotential::HGH_SeparablePotential& h_) : h(h_) {}
    virtual size_t NumProjectors  (int Z)           const override {return h.NumProjectors(Z);}
    virtual double Coefficient    (int Z, size_t p) const override {return h.Coefficient(Z,p);}
    virtual int    AngularMomentum(int Z, size_t p) const override {return h.AngularMomentum(Z,p);}
    virtual double Projector      (int Z, size_t p, double q) const override {return h.Projector(Z,p,q);}
    virtual double BetaR          (int Z, size_t p, double r) const override {return h.BetaR(Z,p,r);}
};
} //anon
TEST(GPW, AnalyticSeparablePPMatchesMesh)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});                    // the corner atom = the wrap-sensitive case
    cell.AddAtom(14, {0.25,0.25,0.25});
    // SR basis + AUTO Rcut (the production SCF configuration): the comparison needs a SCREENED-COMPLETE image
    // set.  At an under-enumerated Rcut the two paths truncate DIFFERENTLY -- the mesh's Bloch orbital catches
    // chi-image x beta-image separations up to 2 Rcut while the analytic single sum stops at Rcut (the "two
    // self-consistent schemes" pin, doc/GPWPlan.md) -- measured 9.3e-2 for diffuse SIPP at Rcut=1.5a vs
    // agreement at AUTO.  Complete enumeration is the production setting; the gate pins THAT.
    std::shared_ptr<const Real_BS> mol(
        BasisSet::Molecule::Factory(BasisSetData::SIPP_SR, &cell,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);

    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);
    MeshOnlyKB meshOnly(gth.nonlocal);
    // Call the EVALUATOR directly: the IBS's MakeSeparablePotential is framework-cached by BasisSetID (which
    // does not key the potential model), so a second call through the IBS would replay the first matrix.
    const GPW_Evaluator& ev = gpw;
    auto Va = ev.MakeSeparablePP(&cell, gth.nonlocal);   // closed-Gaussian face -> ANALYTIC
    auto Vm = ev.MakeSeparablePP(&cell, meshOnly);       // face hidden -> legacy mesh quadrature

    ASSERT_EQ(Va.rows(), Vm.rows());
    double num=0.0, den=0.0, imax=0.0;
    for (size_t i=0;i<Va.rows();i++)
        for (size_t j=0;j<Va.columns();j++)
        {
            const dcmplx va=Va(i,j), vm=Vm(i,j);
            num += std::norm(va-vm);
            den += std::norm(vm);
            imax = std::max(imax, std::fabs(va.imag()));
        }
    double rel=std::sqrt(num/den);
    std::cout << "[analytic KB] ||Va-Vm||_F/||Vm||_F = " << rel << "  max|Im(Va)| = " << imax << std::endl;
    EXPECT_LT(imax, 1e-12);                       // Gamma: analytic KB matrix is real
    EXPECT_LT(rel,  1e-8) << "analytic KB must match the mesh quadrature to the mesh's own error (pinned 4.6e-11)";
}

// The local-PP sweep's ABSOLUTE pair->level rule is STANDALONE-exact (doc/GPWPlan.md 0e-PP step (a)):
// req = kappa*(alpha_i+alpha_j) bounds every pair's spectral tail by e^{-kappa/2} independent of the
// field's sharpness, so doubling kappa (e^{-15} -> e^{-30}) must leave <i|V|j> unchanged to tolerance --
// for the FULL V_loc AND for the split long / short pieces SEPARATELY (each piece must stand alone with
// no cancellation partner, the property step (b)'s analytic short requires).  Also pins the split
// linearity Long + Short == Full (the sweep is linear in the form factor).
TEST(GPW, LocalPPKappaSelfConverged)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    std::shared_ptr<const Real_BS> mol(
        BasisSet::Molecule::Factory(BasisSetData::SIPP_SR, &cell,
                                    BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    // densityEcut=10 keeps the sweeps cheap: the absolute rule's tail bound e^{-kappa/2} references NO grid
    // (that is its point), so kappa-independence is testable at any ladder.  (First measured at Ecut=20:
    // Full 7.6e-9 / Long 1.4e-9 / Short 1.6e-8 between kappa=30 and 60 -- the e^{-15} class on the nose.)
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/10.0);
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);
    const GPW_Evaluator& ev = gpw;   // evaluator directly: bypass the framework cache (not model-keyed)

    auto relDiff=[](const chmat_t& A, const chmat_t& B)
    {
        double num=0.0, den=0.0;
        for (size_t i=0;i<A.rows();i++)
            for (size_t j=0;j<A.columns();j++) { num+=std::norm(A(i,j)-B(i,j)); den+=std::norm(B(i,j)); }
        return std::sqrt(num/den);
    };
    using LP=GPW_Evaluator::LocalPart;
    const chmat_t Vf=ev.MakeLocalPP(&cell, gth.local, LP::Full);           // kappa=30 (the default)
    const chmat_t Vl=ev.MakeLocalPP(&cell, gth.local, LP::Long);
    const chmat_t Vs=ev.MakeLocalPP(&cell, gth.local, LP::Short);
    setenv("GPW_LOCALPP_RELCUTOFF","60",1);                                // e^{-30}: the converged reference
    const chmat_t Vl60=ev.MakeLocalPP(&cell, gth.local, LP::Long);
    const chmat_t Vs60=ev.MakeLocalPP(&cell, gth.local, LP::Short);
    unsetenv("GPW_LOCALPP_RELCUTOFF");
    const double relL=relDiff(Vl,Vl60), relS=relDiff(Vs,Vs60);
    std::cout << "[localPP kappa] long ||V30-V60||/||V60|| = "<<relL<<"  short = "<<relS<<std::endl;
    EXPECT_LT(relL, 1e-5) << "LONG: kappa=30 not converged (tail bound violated)";
    EXPECT_LT(relS, 1e-5) << "SHORT: kappa=30 not converged (tail bound violated)";
    chmat_t Vls=Vl; Vls+=Vs;                                               // Full's convergence follows by linearity
    EXPECT_LT(relDiff(Vls,Vf), 1e-12) << "split linearity: Long + Short == Full (same sweep, same levels)";

    // Step (b) cross-validation: the ANALYTIC short (exact 3-centre Gaussian lattice sums, the production
    // path since 2026-07-22; periodic G=0 mean subtracted -- the 5.7% double-count bug this gate caught) vs
    // the kappa-ruled GRID short on the PERIODIC crystal.  NOTE the tolerance: the GRID side is the
    // approximate one here -- its SHARPEST pairs saturate at this cheap ladder's top (Ecut=10 -> rung 20;
    // kappa*p unreachable), an error the kappa SELF-convergence above cannot see (both kappas hit the same
    // top).  Measured 3.6e-3 at Ecut=10 (the saturation tail e^{-ecut_top/2p} * V-tilde(ball edge) class);
    // the tolerance 1e-2 pins CONVENTION-class bugs (G=0, phases, normalisation), and the Si/NaF SCF anchors
    // carry the mHa-precision verification of the analytic path.
    const chmat_t Va=ev.MakeLocalPPShort(&cell, gth.local);                // closed-Gaussian face -> ANALYTIC
    const double relA=relDiff(Va,Vs);
    std::cout << "[localPP short] ||analytic-grid||/||grid|| = "<<relA<<std::endl;
    EXPECT_LT(relA, 1e-2) << "analytic short vs kappa-ruled grid short: convention-class disagreement";
}
