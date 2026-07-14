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
#include <iostream>

import qchem.Structure;                         // Molecule, Atom
import qchem.UnitCell;                          // UnitCell
import qchem.BasisSet;                          // Real_BS
import qchem.BasisSet.Orbital_1E_IBS;           // Real_OIBS / Complex_OIBS + cached Overlap()/Kinetic()/Nuclear()
import qchem.BasisSet.Molecule.Factory;         // Molecule::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Lattice_3D.GPW_IBS;       // GPW_IBS (the basis under test)
import qchem.BasisSet.Lattice_3D.Evaluators.GPW; // GPW_Evaluator (tests may cheat-import internals) -- DFT tier
import qchem.BasisSet.Molecule.LatticeSum1E;     // Molecule::LatticeSum1E::CollocateDensity (analytic collocation)
import qchem.BasisSet.Internal.GMap;            // G_ERI3 / ΔG_Map (the collocation tensor + rho-tilde)
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

    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/0.0, /*Rcut=*/0.0); // home cell only
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
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/0.0, /*Rcut=*/1.5*a); // images
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

// Collocate: the G=0 collocation weight is the grid-quadrature overlap, W_0(i,j)*Omega = integral chi_i chi_j,
// which must converge to the analytic 1E overlap as the density grid resolves the Gaussian products.
TEST(GPW, CollocationOverlapMatchesAnalytic)
{
    const double a=10.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> molCell = MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), molCell, /*densityEcut=*/30.0);

    const GPW_Evaluator& ev = gpw;
    const Complex_OIBS&   g  = gpw;
    G_ERI3 ov = ev.Overlap3CTensor();
    int c0 = G0Column(ov);
    ASSERT_GE(c0, 0);

    const auto& S = g.Overlap();     // analytic 1E overlap (chmat, real at Gamma)
    const size_t n = S.rows();
    double num=0.0, den=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            double w = std::real(ov.weights[c0](i,j)) * ov.volume;   // grid-quadrature overlap
            double s = std::real(S(i,j));
            double d=w-s; num+=d*d; den+=s*s;
        }
    EXPECT_LT(std::sqrt(num/den), 6e-2);   // collocated == analytic to grid-quadrature accuracy (~4% at Ecut=30)
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

    const GPW_Evaluator& ev = gpw;
    const Complex_OIBS&   g  = gpw;
    const double V0 = 0.7;
    const chmat_t M = ev.OverlapMatrix([V0](const ivec3_t& dm)->dcmplx
        { return (dm.x==0 && dm.y==0 && dm.z==0) ? dcmplx(V0) : dcmplx(0.0); });

    const auto& S = g.Overlap();
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

// PERF micro-benchmark (DISABLED): isolate OverlapMatrix -- the SCF integrate-back that is ~44% of the NaF
// profile -- at TRUE NaF scale (FCC Na+F, VALENCE_LOWQ_SR n=32, densityEcut AUTO->160) but with a FAST setup:
// Rcut=0 makes PhiOnGrid a single-image build (seconds), and OverlapMatrix's cost depends only on (n, Npts),
// which Rcut does not change.  Times K back-to-back calls on a representative all-G field (Hartree/XC-like),
// so a fixed setup is amortised -- an A/B lever for the OverlapMatrix contraction (scalar loop vs GEMM etc.).
TEST(GPW, DISABLED_BenchOverlapMatrix)
{
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11,{0,0,0});          // Na
    cell.AddAtom(9,{0.5,0.5,0.5});     // F  (tight 40-a.u. exponent -> densityEcut floor 160)
    auto mol=std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ_SR, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut AUTO*/-1.0, /*Rcut*/0.0);
    const GPW_Evaluator& ev=gpw;
    auto field=[](const ivec3_t& dm)->dcmplx                     // smooth, all-G (realistic Vtilde)
        { double g2=double(dm.x*dm.x+dm.y*dm.y+dm.z*dm.z); return dcmplx(1.0/(1.0+g2),0.0); };
    chmat_t warm=ev.DenseOverlapMatrix(field);                    // warm the PhiOnGrid cache (excluded from timing)
    chmat_t warmMG=ev.MultiGridOverlapMatrix(field);              // warm the per-level column caches too
    const int K=20;
    auto t0=std::chrono::steady_clock::now();
    chmat_t M;
    for (int q=0;q<K;q++) M=ev.DenseOverlapMatrix(field);
    auto t1=std::chrono::steady_clock::now();
    for (int q=0;q<K;q++) M=ev.MultiGridOverlapMatrix(field);
    auto t2=std::chrono::steady_clock::now();
    double msDense=std::chrono::duration<double,std::milli>(t1-t0).count();
    double msMG   =std::chrono::duration<double,std::milli>(t2-t1).count();
    std::cout << "[bench integrate-back] n=" << warm.rows() << " calls=" << K
              << "  dense=" << msDense/K << " ms/call   multi-grid=" << msMG/K << " ms/call   speedup="
              << msDense/msMG << "x" << std::endl;
    SUCCEED();
}

// PATCH-SPARSITY PROBE (DISABLED): measure-before-optimize for the GPW_Plan.md S0 patch-collocation rewrite.
// The per-iteration hotspot is the integrate-back  M_ij = w Sum_p conj(Phi_pi) V_p Phi_pj  over the FULL grid
// (n^2 x Npts, ~44% of the NaF profile).  The plan proposes replacing the dense contraction with LOCAL patches
// (points where the orbital / orbital-product is non-negligible), turning O(n^2 Npts) into O(n_pairs x patch).
// The achievable speedup depends entirely on how sparse Phi actually is -- diffuse Gaussians reach far, so a
// per-ORBITAL patch may be large while the per-PAIR PRODUCT patch (bounded by the tighter product Gaussian,
// exponent alpha_i+alpha_j) is compact.  This probe reconstructs Phi on the density grid (via the public Eval /
// DensityGrid surface, no internals) at TRUE NaF scale and prints, for a few screening thresholds, the exact
// contraction cost of each candidate strategy vs the dense baseline, so the granularity choice is data-driven:
//   dense           : n^2 * Npts                              (what the zgemm does today)
//   per-point scatter: Sum_p k_p^2, k_p=#{i:|Phi_pi|>eps}     (skip ~0 orbitals per point; loses vectorization)
//   per-pair product : Sum_{i,j} #{p:|Phi_pi Phi_pj|>eps}     (CP2K-style: contract each pair on its own patch)
TEST(GPW, DISABLED_PatchSparsityProbe)
{
    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11,{0,0,0});          // Na
    cell.AddAtom(9,{0.5,0.5,0.5});     // F  (tight 40-a.u. exponent -> densityEcut floor 160)
    auto mol=std::shared_ptr<const Real_BS>(BasisSet::Molecule::Factory(
        BasisSetData::VALENCE_LOWQ_SR, &cell, BasisSet::Molecule::Engine::MnD, BasisSet::Molecule::Angular::Cartesian));
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut AUTO*/-1.0, /*Rcut*/0.0);
    const GPW_Evaluator& ev=gpw;

    const rvec3vec_t& pts=ev.DensityGrid().GridPoints();
    const size_t Npts=pts.size(), n=ev.size();
    std::vector<double> Phi(Npts*n);                 // |Phi_pi| (magnitudes; real at Gamma) row-major [p][i]
    double phiMax=0.0;
    for (size_t p=0;p<Npts;p++)
    {
        cvec_t v=ev.Eval(pts[p]);
        for (size_t i=0;i<n;i++) { double m=std::abs(v[i]); Phi[p*n+i]=m; phiMax=std::max(phiMax,m); }
    }
    const double dense=double(n)*double(n)*double(Npts);
    std::cout << "[patch probe] n=" << n << " Npts=" << Npts << " max|Phi|=" << phiMax
              << "  dense n^2*Npts=" << dense << " complex MACs/call\n";

    for (double ratio : {1e-6, 1e-8, 1e-10})
    {
        const double epsOrb =ratio*phiMax;               // per-orbital screen
        const double epsProd=ratio*phiMax*phiMax;        // per-product screen (max product = phiMax^2)
        double scatter=0.0; size_t nz=0, unionPts=0;
        double pairCost=0.0; size_t pairsNZ=0;
        for (size_t p=0;p<Npts;p++)
        {
            const double* row=&Phi[p*n];
            size_t kp=0; for (size_t i=0;i<n;i++) if (row[i]>epsOrb) kp++;
            scatter += double(kp)*double(kp); nz += kp; if (kp) unionPts++;
        }
        // per-pair product patch sizes (O(n^2 Npts) one-off; fine for a disabled probe)
        std::vector<size_t> pairPts(n*n,0);
        for (size_t p=0;p<Npts;p++)
        {
            const double* row=&Phi[p*n];
            for (size_t i=0;i<n;i++) { double vi=row[i]; if (vi<=0) continue;
                for (size_t j=0;j<n;j++) if (vi*row[j]>epsProd) pairPts[i*n+j]++; }
        }
        for (size_t ij=0;ij<n*n;ij++) { pairCost+=double(pairPts[ij]); if (pairPts[ij]) pairsNZ++; }
        std::cout << "  ratio=" << ratio
                  << " | scatter Sum k_p^2=" << scatter << " (" << scatter/dense << "x dense),"
                  << " nz=" << nz << " unionPts=" << unionPts << "/" << Npts
                  << " | pair-product cost=" << pairCost << " (" << pairCost/dense << "x dense),"
                  << " nz-pairs=" << pairsNZ << "/" << n*n << "\n";
    }
    SUCCEED();
}

// PATCHED INTEGRATE-BACK bit-consistency (GPW_Plan.md S0, Increment 1).  The molecular-side patched
// collocation adjoint (LatticeSum1E::MakePotentialMatrix -- per-orbital Gaussian-support patches, contract
// each pair on its support overlap) must reproduce the dense zgemm OverlapMatrix(Vtilde) to the screening
// tolerance.  This is the SEAM gate for the multi-grid rewrite; the dense path stays the default (byte-
// identical anchors) until the grid levels make the patched path win.  Exercised at Gamma (real orbitals)
// AND at a genuinely-complex k with images (the conj(bra) slot + the Bloch image sum in EnsureSupports).
TEST(GPW, PatchedIntegrateBackMatchesDense)
{
    auto relDiff=[](const chmat_t& A, const chmat_t& B)
    {
        double num=0.0, den=0.0;
        for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.rows();j++)
        { dcmplx d=A(i,j)-B(i,j); num+=std::norm(d); den+=std::norm(A(i,j)); }
        return std::sqrt(num/den);
    };
    auto field=[](const ivec3_t& dm)->dcmplx                     // smooth all-G field (Hartree/XC-like)
        { double g2=double(dm.x*dm.x+dm.y*dm.y+dm.z*dm.z); return dcmplx(1.0/(1.0+g2),0.0); };

    // (a) Gamma, home cell: real orbitals, single image.
    {
        const double a=12.0;
        UnitCell cell(a);
        cell.AddAtom(14,{0.5,0.5,0.5});
        std::shared_ptr<const Real_BS> mol = MakeBasis(cell);           // SIPP Si
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/12.0, /*Rcut*/0.0);
        const GPW_Evaluator& ev=gpw;
        EXPECT_LT(relDiff(ev.OverlapMatrix(field), ev.PatchedOverlapMatrix(field)), 1e-8) << "Gamma";
    }
    // (b) complex k=(1/4,0,0) with images: the Bloch sum is live (conj slot + image support build).
    {
        const double a=10.0;                                            // small cell -> non-negligible images
        UnitCell cell(a);
        cell.AddAtom(14,{0.5,0.5,0.5});
        std::shared_ptr<const Real_BS> mol = MakeBasis(cell);
        GPW_IBS gpw(cell, ivec3_t(4,4,4), ivec3_t(1,0,0), mol, /*densityEcut*/12.0, /*Rcut*/1.5*a);
        const GPW_Evaluator& ev=gpw;
        chmat_t Md=ev.OverlapMatrix(field);
        EXPECT_GT(MaxImag(Md), 1e-6) << "k!=0 with images must give a complex integrate-back";
        EXPECT_LT(relDiff(Md, ev.PatchedOverlapMatrix(field)), 1e-8) << "complex k";
    }
}

// MULTI-GRID integrate-back MECHANICS (GPW_Plan.md S0 Increment 2).  MultiGridOverlapMatrix maps each orbital
// pair to the coarsest grid level resolving its product (from the internal exponents), restricts V per level,
// and contracts each pair on its assigned level.  STRUCTURAL invariant (grid points / weights / columns /
// level assignment / V restriction all sound): for a CONSTANT potential the integrate-back is V0*<i|j> on ANY
// grid, so multi-grid must equal the fine grid to the coarse-grid quadrature error of the orbital norms.
// The FIELD-dependent coarsening error (dropping V's high-G coupling to diffuse pairs) is MEASURED but NOT
// gated: it is large for a PEAKED V (e.g. 1/|G|^2, the local PP), which is why diffuse pairs that overlap a
// sharp external potential cannot be coarsened independently of it -- doing Increment 2 right needs a
// CONSISTENT whole-multigrid (density collocation AND integrate-back per level) + PP-aware assignment
// (an integrate-back-only multigrid gave Si Gamma Etot -21.4 vs -8.25; see doc/GPWPlan.md).
TEST(GPW, MultiGridIntegrateBackMechanics)
{
    auto relDiff=[](const chmat_t& A, const chmat_t& B)
    {
        double num=0.0, den=0.0;
        for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.rows();j++)
        { dcmplx d=A(i,j)-B(i,j); num+=std::norm(d); den+=std::norm(A(i,j)); }
        return std::sqrt(num/den);
    };
    const double a=12.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);           // SIPP Si
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/12.0, /*Rcut*/0.0);
    const GPW_Evaluator& ev=gpw;

    auto cfield=[](const ivec3_t& dm)->dcmplx { return (dm.x==0&&dm.y==0&&dm.z==0)?dcmplx(1.0):dcmplx(0.0); };
    double rdConst=relDiff(ev.PatchedOverlapMatrix(cfield), ev.MultiGridOverlapMatrix(cfield));
    EXPECT_LT(rdConst, 5e-3) << "multi-grid must match the fine grid for a constant potential (structural)";

    auto pfield=[](const ivec3_t& dm)->dcmplx        // peaked field -> the coarsening error (informational)
        { double g2=double(dm.x*dm.x+dm.y*dm.y+dm.z*dm.z); return dcmplx(1.0/(1.0+g2),0.0); };
    std::cout << "[multi-grid] const-field relDiff=" << rdConst << "  peaked-field relDiff="
              << relDiff(ev.PatchedOverlapMatrix(pfield), ev.MultiGridOverlapMatrix(pfield)) << std::endl;
}

// ANALYTIC COLLOCATION charge conservation (GPWPlan.md S0 Increment A -- the CP2K rewrite).  The new
// LatticeSum1E::CollocateDensity collocates rho = Sum_ij D_ij chi_i chi_j analytically per pair on compact
// exp-tail boxes, modulo-wrapped onto the grid (NO image sum, NO Rcut).  The defining invariant: Integral of
// the collocated rho over the cell equals Tr(D S) (the density integrates to the electron count) -- to grid
// tolerance.  Tested with the atom in the INTERIOR (no wrap) AND at the CORNER {0,0,0} (the box wraps around
// every face -> exercises the modulo-wrap-IS-the-image-sum mechanism; charge must be identical, translation-
// invariant, with no ringing).
TEST(GPW, AnalyticCollocationConservesCharge)
{
    auto probe=[](const rvec3_t& frac, const char* where)->double
    {
        const double a=12.0;
        UnitCell cell(a);
        cell.AddAtom(14, frac);                                     // Si
        std::shared_ptr<const Real_BS> mol = MakeBasis(cell);       // SIPP Si
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/12.0, /*Rcut*/0.0);
        const GPW_Evaluator& ev=gpw;
        const Complex_OIBS& g=gpw;
        const auto* lat=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(OrbitalBlock<Real_OIBS>(*mol));
        EXPECT_TRUE(lat) << "orbital block must realise LatticeSum1E";
        const ivec3_t N=ev.DensityGrid().FFTGrid();
        const size_t  n=g.GetNumFunctions();
        rmat_t D(n,n,0.0); for (size_t i=0;i<n;i++) D(i,i)=1.0;      // D = identity -> Integral rho = Tr(S)
        rvec_t rho=lat->CollocateDensity(D, cell, N);
        double integral=blazem::sum(rho)*cell.GetCellVolume()/double(rho.size());
        const auto& S=g.Overlap();
        double trDS=0.0; for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) trDS+=D(i,j)*std::real(dcmplx(S(i,j)));
        std::cout << "[collocate " << where << "] Integral rho=" << integral << "  Tr(D S)=" << trDS
                  << "  rel=" << std::fabs(integral-trDS)/std::fabs(trDS) << std::endl;
        EXPECT_NEAR(integral, trDS, 5e-2*std::fabs(trDS)) << "collocated charge vs Tr(D S) (" << where << ")";
        return integral;
    };
    double cInterior=probe({0.5,0.5,0.5}, "interior");
    double cCorner  =probe({0.0,0.0,0.0}, "corner-wrapped");
    EXPECT_NEAR(cInterior, cCorner, 1e-6) << "collocated charge must be translation-invariant (wrap == interior)";
}

// ANALYTIC COLLOCATION on a CRYSTAL (cross-cell pairs).  The periodic Gamma density is a product of BLOCH
// orbitals, chi_i^G chi_j^G = Sum_R'' chi_i^0 chi_j^R'' -- so collocation must sum the screened CROSS-CELL
// offsets, not just the home pair (R''=0).  Invariant: Integral of the collocated rho == Tr(D S^G) with S^G
// the Bloch overlap (its own screened image sum, via a generous Rcut).  A 2-atom Si crystal with real
// inter-cell overlap (which the single-atom AnalyticCollocationConservesCharge test could not exercise).
TEST(GPW, DISABLED_AnalyticCollocationCrystalChargeConservation)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);              // SIPP Si (one orbital block over both atoms)
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/6.0, /*Rcut*/4.0*a);  // large Rcut -> S^G
    const GPW_Evaluator& ev=gpw;
    const Complex_OIBS& g=gpw;
    const auto* lat=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(OrbitalBlock<Real_OIBS>(*mol));
    EXPECT_TRUE(lat);
    const ivec3_t N=ev.DensityGrid().FFTGrid();
    const size_t  n=g.GetNumFunctions();
    rmat_t D(n,n,0.0); for (size_t i=0;i<n;i++) D(i,i)=1.0;            // D = identity -> Integral rho = Tr(S^G)
    rvec_t rho=lat->CollocateDensity(D, cell, N);
    double integral=blazem::sum(rho)*cell.GetCellVolume()/double(rho.size());
    const auto& S=g.Overlap();                                        // Bloch overlap (screened images) = S^G
    double trDS=0.0; for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) trDS+=D(i,j)*std::real(dcmplx(S(i,j)));
    std::cout << "[collocate crystal] Integral rho=" << integral << "  Tr(D S^G)=" << trDS
              << "  rel=" << std::fabs(integral-trDS)/std::fabs(trDS) << std::endl;
    EXPECT_NEAR(integral, trDS, 5e-2*std::fabs(trDS)) << "crystal collocated charge vs Tr(D S^G)";
}

// ANALYTIC INTEGRATE-BACK (GPWPlan.md S0 Increment B): LatticeSum1E::IntegratePotential is the exact adjoint
// of CollocateDensity (same box + wrap).  Two gates: (1) ADJOINT consistency <collocate(D),V> == <D,integrate(V)>
// to machine precision (variational -- the KS matrix is the exact gradient of the grid energy); (2) it
// reproduces the dense sampled integrate-back (DenseOverlapMatrix) to grid tolerance.
TEST(GPW, DISABLED_AnalyticIntegrateBackAdjointAndDense)
{
    const double a=12.0;
    UnitCell cell(a);
    cell.AddAtom(14,{0.5,0.5,0.5});
    std::shared_ptr<const Real_BS> mol = MakeBasis(cell);
    GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/12.0, /*Rcut*/0.0);
    const GPW_Evaluator& ev=gpw;
    const auto* lat=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(OrbitalBlock<Real_OIBS>(*mol));
    EXPECT_TRUE(lat);
    const ivec3_t N=ev.DensityGrid().FFTGrid();
    const size_t  n=ev.size();

    auto field=[](const ivec3_t& dm)->dcmplx
        { double g2=double(dm.x*dm.x+dm.y*dm.y+dm.z*dm.z); return dcmplx(1.0/(1.0+g2),0.0); };
    ΔG_Map vmap; for (const ivec3_t& dm : ev.DensityGrid().Gs()) vmap[dm]=field(dm);
    rvec_t V=ev.DensityGrid().RhoOnGrid(vmap);                 // V(r) on the grid
    rmat_t h=lat->IntegratePotential(V, cell, N);

    // (1) adjoint: <collocate(D),V> == Tr(D integrate(V)) to machine precision (same box + wrap on both sides)
    rmat_t D(n,n,0.0); for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) D(i,j)=(i==j)?1.0:0.3;  // symmetric
    rvec_t rho=lat->CollocateDensity(D, cell, N);
    const double w=cell.GetCellVolume()/double(V.size());
    double lhs=0.0; for (size_t g=0;g<V.size();g++) lhs+=rho[g]*V[g]*w;
    double rhs=0.0; for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) rhs+=D(i,j)*h(i,j);
    EXPECT_NEAR(lhs, rhs, 1e-8*std::fabs(rhs)) << "adjoint <collocate(D),V> == <D,integrate(V)>";

    // (2) vs the dense sampled integrate-back, to grid tolerance
    chmat_t hd=ev.DenseOverlapMatrix(field);
    double num=0.0, den=0.0;
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
    { double d=h(i,j)-std::real(dcmplx(hd(i,j))); num+=d*d; den+=std::norm(dcmplx(hd(i,j))); }
    std::cout << "[integrate-back] adjoint lhs=" << lhs << " rhs=" << rhs
              << "  vs-dense relDiff=" << std::sqrt(num/den) << std::endl;
    EXPECT_LT(std::sqrt(num/den), 5e-2) << "analytic integrate-back vs dense to grid tolerance";
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

    GPW_IBS gpw(cell, ivec3_t(4,4,4), ivec3_t(1,0,0), molCell, /*densityEcut=*/0.0, /*Rcut=*/0.0); // k=(1/4,0,0)
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

    GPW_IBS gpw(cell, ivec3_t(4,4,4), ivec3_t(1,0,0), molCell, /*densityEcut=*/0.0, /*Rcut=*/1.5*a);
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
    GPW_IBS gpw(cell, N, ik, molCell, /*densityEcut=*/0.0, /*Rcut=*/2.5*a);
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
    GPW_IBS gk (cell, N, ivec3_t(1,0,0), molCell, 0.0, 1.5*a);   //  k = (1/4,0,0)
    GPW_IBS gmk(cell, N, ivec3_t(3,0,0), molCell, 0.0, 1.5*a);   // -k = (3/4,0,0) == -(1/4,0,0) mod 1
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
