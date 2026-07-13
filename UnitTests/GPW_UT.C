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
    chmat_t warm=ev.OverlapMatrix(field);                         // warm the PhiOnGrid cache (excluded from timing)
    const int K=20;
    auto t0=std::chrono::steady_clock::now();
    chmat_t M;
    for (int q=0;q<K;q++) M=ev.OverlapMatrix(field);
    auto t1=std::chrono::steady_clock::now();
    double ms=std::chrono::duration<double,std::milli>(t1-t0).count();
    std::cout << "[bench OverlapMatrix] n=" << warm.rows() << " calls=" << K
              << " total=" << ms << " ms  per-call=" << ms/double(K) << " ms" << std::endl;
    SUCCEED();
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
