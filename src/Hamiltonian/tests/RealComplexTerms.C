// File: RealComplexTermsUT.C  Step-3c-1 gates: the periodic terms serve a REAL TRIM block (doc/RealComplexPlan.md).
//
// Step 3a made the real TRIM block (tGPW_IBS<double>) exist; this increment makes the TERM STACK able to
// answer it.  Each periodic term now carries the real block faces tStatic_CC<double> / tDynamic_CC<double> (the
// V1.6/V1.7 cross-cast idiom), implemented by ONE scalar-generic assembly body per term -- so a real
// block's matrix is the SAME arithmetic narrowed at the end, and must equal the complex block's real
// part BITWISE (EXPECT_EQ on doubles; Step 0's exact +/-1 phases are what make zero tolerance possible).
//
// Covered here: the static set (Kinetic, Ven_PP_Short/Long/NonLocal) and the density-dependent
// Vee_Hartree + the Becke-route cVxc_DFT.  Each scalar arm drives its OWN production-shaped
// density (a real child on the real block, a complex child on the complex block -- the 3c-3 mixed
// shape): the GPW integrate-back's D-aware screen is INSTANCE-PAIRED state (the collocation and its
// adjoint must go through the same block instance -- the 2026-08-18 cross-run-pollution fix scoped the
// 3C tensors per instance, see tGPW_IBS::Repulsion3C), so a cross-asked arrangement (density collocated
// on cx, matrix asked of re) is no longer screen-consistent and was only ever bitwise by the two twins
// accidentally sharing one DBCache'd closure.  The pair-quadrature RAW-route real path needs the full
// collocation feed and is gated with the 3c-2 SCF wiring.
#include "gtest/gtest.h"
#include <algorithm>
#include <cstdio>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>

import qchem.Structure;
import qchem.Tests.PW_Fields;   // OneBlockBS (a one-block whole basis for the refresh)
import qchem.UnitCell;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;          // tGPW_IBS<T> (both block alternatives)
import qchem.BasisSet.Gaussian.Lattice.RasterLadder;     // MakeRasterLadder (one ladder shared by the k-fold gate's blocks)
import qchem.Hamiltonian;                          // the terms and their real block faces
import qchem.Hamiltonian.Internal.Kinetic;         // Kinetic<T> (tests may cheat-import internals)
import qchem.Hamiltonian.Internal.Hamiltonian;     // tHamiltonianImp<dcmplx> (the 3c-2 assembly gate)
import qchem.Hamiltonian.Internal.Terms;            // Vee_Fit<dcmplx> (the periodic Hartree term)
import qchem.Fitting.CoulombOperator;               // MakeGSpaceDensityFitter
import qchem.Hamiltonian.Internal.PWTerms;         // the periodic term set + the two DensitySampler strategies
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
                                                  // since 2026-09-08; .Internal. modules are
                                                  // never re-exported, so name it directly)
import qchem.BasisSet.DeltaFit_IBS;                // DeltaFit_IBS -- the delta basis the singles strategy runs on
import qchem.Symmetry.Factory;                     // BlochFactory (its Gamma irrep)
import qchem.Hamiltonian.Internal.SlaterExchange;  // SlaterExchange (the Dirac-exchange functional)
import qchem.BasisSet.Orbital_PP_IBS;              // the real child's own PP pullbacks (the force-walk gate)
import qchem.Hamiltonian.Internal.GeometryGradient;   // tStatic_PC (the force-walk gate)
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH (Si GTH-LDA-q4)
import qchem.Pseudopotential.LocalPotential;       // GaussianSmearedNucleus (a local field with no closed-Gaussian form)
import qchem.ChargeDensity.Imp.IrrepCD;            // PeriodicIrrepCD (the periodic leaf -- 3c-2b)
import qchem.CompositeCD;                          // tComposite_CD<dcmplx> (the run-shaped density)
import qchem.Mesh;                                 // qcMesh::MeshParams
import qchem.Mesh.XCPolicy;                        // BeckeXCParams (the Becke fold gate)
import qchem.Fitting.FunctionFitter;               // DeltaScalarFit
import qchem.BasisSet.GSpaceDensityProjector;      // ρ̃ through the one projection face (V-ONE-OF-N 3c)
import qchem.BasisSet.RasterDensityProjector;      // raw ρ_DM through the same face (R1.0n/o)
import qchem.Types;
import qchem.Reporting;                            // the factored leaf's health record (E-c (b))
import qchem.Structure.MeshFactory;   // MakeMeshMaker (the GPW block demands its cell's mesh maker, V-MESH-FACTORY §4.4)

using namespace qchem;
using namespace qchem::Hamiltonian;
using BasisSet::Real_BS;
using BasisSet::Gaussian::tGPW_IBS;

namespace
{
struct Rig
{
    UnitCell cell;
    std::shared_ptr<const Structure> st;
    std::shared_ptr<const Real_BS>   mol;
    std::unique_ptr<tGPW_IBS<double>> re;   // the real Gamma TRIM block (Step 3a)
    std::unique_ptr<tGPW_IBS<dcmplx>> cx;   // the same physical block, complex
    Rig() : cell(8.0)
    {
        cell.AddAtom(14,{0.5,0.5,0.5});
        st=cell.Clone();
        namespace M = qchem::BasisSet::Gaussian;
        mol.reset(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
        re=std::make_unique<tGPW_IBS<double>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
        cx=std::make_unique<tGPW_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    }
    //! A one-block complex density in the run's shape (tComposite_CD<dcmplx>): D(0,0)=2 -- two electrons
    //! in the first orbital.  Enough to make V_H / rho genuinely nonzero; the terms never see past the
    //! composite face.
    std::unique_ptr<ChargeDensity::cDM_CD> MakeDensity() const
    {
        using namespace qchem::ChargeDensity;
        chmat_t D(cx->GetNumFunctions());
        D(0,0)=2.0;
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(
        new PeriodicIrrepCD<dcmplx>(D, cx.get(), cx->GetIrrep(Spin::None))), cx->GetIrrep(Spin::None));
        return cd;
    }
    //! The REAL twin of \c MakeDensity: the same D(0,0)=2 as a REAL child on the REAL block inside the
    //! complex-faced composite -- the 3c-3 mixed shape.  The real arm of every density-dependent gate
    //! drives THIS density, so its collocation and integrate-back pair on the real instance exactly as
    //! a production mixed run pairs them (the screen-consistency note in the file header).
    std::unique_ptr<ChargeDensity::cDM_CD> MakeRealDensity() const
    {
        using namespace qchem::ChargeDensity;
        rsmat_t D(re->GetNumFunctions());
        D(0,0)=2.0;
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<tDM_CD<double>>(
        new PeriodicIrrepCD<double>(D, re.get(), re->GetIrrep(Spin::None))), re->GetIrrep(Spin::None));
        return cd;
    }
};

void ExpectBitwiseEqual(const rsmat_t& R, const chmat_t& C, const char* what)
{
    ASSERT_EQ(R.rows(), C.rows()) << what;
    for (size_t i=0;i<R.rows();i++)
        for (size_t j=0;j<R.columns();j++)
            EXPECT_EQ(R(i,j), std::real(C(i,j))) << what << " (" << i << "," << j << ")";
}
// The Becke quadrature GEMM comparison CANNOT be bitwise, and that is the Step-3 win showing itself:
// the real block's Phi table drives blaze's REAL dgemm-shaped kernel, whose accumulation order differs
// from the complex kernel's -- same exact integrand, last-ulp summation differences (~3e-15 measured).
// Everything ELEMENTWISE or same-order (S/T/V/KB, the Hartree ContractAdjoint) stays EXPECT_EQ above.
void ExpectMachineEqual(const rsmat_t& R, const chmat_t& C, const char* what)
{
    ASSERT_EQ(R.rows(), C.rows()) << what;
    for (size_t i=0;i<R.rows();i++)
        for (size_t j=0;j<R.columns();j++)
            EXPECT_NEAR(R(i,j), std::real(C(i,j)), 1e-13) << what << " (" << i << "," << j << ")";
}
} //anon

// ---- R1.0c: the FACTORED forward contraction, on BOTH scalars ------------------------------------------
//
// One integral, two representations of D: applyRaw(D) is the full quadratic form [Phi D Phi^dag]_gg and
// applyRawFactored(L) is Sum_m |[Phi L]_gm|^2 for a caller holding D = L L^dag.  They must agree, and
// NOTHING WAS CHECKING THAT.  Measured 2026-08-23 by instrumenting DeltaFit_IBS::ForwardFactoredT over the
// whole suite: the factored branch fires exactly ONCE, on the `double` instantiation, so
// ForwardFactoredT<dcmplx> -- live code on the default route -- was exercised by no enabled test at all.
//
// It matters more than a normal coverage gap because a wrong factored rho is NOT LOUD: it is a slightly
// wrong density feeding a pointwise-nonlinear functional, which reads as a degraded SCF rather than a bug.
// A UNIT gate, deliberately, not an SCF one: build D from L so the two inputs are the same object by
// construction, and compare the two contractions of the ONE tensor.
namespace
{
inline rvec3_t rig_pt(double a,double b,double c) {return rvec3_t(a,b,c);}

//! A delta fit basis PLUS the fitter over it, held as its projection face -- what a density projects
//! itself onto since 2026-08-24 (a density no longer reaches the fit basis directly).  Built through the
//! same factory the sampler uses, so a probe cannot get a fitter production would not.
struct DeltaProjector
{
    std::shared_ptr<const BasisSet::cFIT_SF_ABS>    basis;
    std::unique_ptr<Fitting::DeltaScalarFit<dcmplx>> fitter;
    explicit DeltaProjector(qchem::BasisSet::FitQuadrature q)
        : basis(std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), 
                    std::move(q), Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0))))
        , fitter(Fitting::MakeDeltaScalarFitter(basis)) {}
    const BasisSet::PairProjector<dcmplx>& operator()() const
        {return *fitter;}
    const BasisSet::DeltaFit_IBS<dcmplx>&   Basis() const
        {return dynamic_cast<const BasisSet::DeltaFit_IBS<dcmplx>&>(*basis);}
};

//! A deterministic thin factor: no RNG, so a failure is reproducible from the test name alone.
template <class U> U FactorEntry(size_t i, size_t m)
{
    const double a=std::sin(1.0+0.7*double(i)+1.3*double(m))/(1.0+0.1*double(i));
    if constexpr (std::is_same_v<U,dcmplx>) return dcmplx(a, std::cos(0.4+0.9*double(i)-0.6*double(m))*0.5);
    else                                    return a;
}

//! \f$\rho\f$ from \f$L\f$ two ways over ONE \c PairProjection: the factored contraction against the full
//! quadratic form of \f$D=LL^\dagger\f$.  PSD by construction, so the factored route is always applicable.
template <class U>
void ExpectFactoredMatchesFull(const BasisSet::DeltaFit_IBS<dcmplx>& dfit,
                               const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb,
                               const char* what)
{
    const size_t n=orb.GetNumFunctions(), r=3;
    ASSERT_GT(n, r) << what << ": the factor must be thin for this to be the interesting case";
    mat_t<U> L(n,r);
    for (size_t i=0;i<n;i++) for (size_t m=0;m<r;m++) L(i,m)=FactorEntry<U>(i,m);
    // D = L L^dag, built ELEMENTWISE so hmat_t's Hermitian invariant holds exactly rather than to roundoff.
    hmat_t<U> D(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
    {
        U s=U(0.0);
        // std::conj, NOT blazem::conj: blaze's conj() conjugates a MATRIX or VECTOR but is the IDENTITY
        // on a std::complex SCALAR (measured 2026-08-24).  That trap is what this gate first tripped over,
        // and then found in production -- see R2.21 in doc/CleanupCandidates.md.
        for (size_t m=0;m<r;m++)
        {
            if constexpr (std::is_same_v<U,dcmplx>) s += L(i,m)*std::conj(L(j,m));
            else                                    s += L(i,m)*L(j,m);
        }
        if (i==j) s=U(std::real(std::complex<double>(s)));   // a Hermitian diagonal is real by construction
        D(i,j)=s;
    }

    const PairProjection<U>& o3=static_cast<const BasisSet::Integrals_Overlap3C<U,dcmplx>&>(dfit).Overlap3C(orb);   // through the face: the real-block overload lives on the periodic base
    ASSERT_TRUE(bool(o3.applyRaw))         << what << ": the delta basis must realise the raw forward";
    ASSERT_TRUE(bool(o3.applyRawFactored)) << what << ": ...and its factored sibling";
    const rvec_t full=o3.applyRaw(D), fact=o3.applyRawFactored(L);
    ASSERT_EQ(full.size(), fact.size()) << what;

    // Not bitwise, and the reason is the point of having two routes: Sum_m |[Phi L]_gm|^2 accumulates over
    // the RANK where [Phi D Phi^dag]_gg accumulates over the BASIS, so the summation orders differ.
    double scale=0.0;
    for (size_t g=0; g<full.size(); g++) scale=std::max(scale, std::fabs(full[g]));
    EXPECT_GT(scale, 1e-6) << what << ": a vacuous (all-zero) rho would pass anything";
    for (size_t g=0; g<full.size(); g++)
        EXPECT_NEAR(fact[g], full[g], 1e-12*std::max(scale,1.0)) << what << " at point " << g;
}
} //anon

namespace {
//! THE CONTRACT (doc/KAxisPlan.md P2; V1.41 stage 3): a term's state is refreshed for a density, then its blocks are built from it,
//! before GetMatrix -- a pure read -- is asked; GetMatrix THROWS otherwise.  A test that drives a term directly does what the eager
//! phase does, on the one block \a b it will read.
template <class Term, class CD, class B> void RefreshTerm(const Term& t, const CD* cd, const B* b)
{
    t.RefreshForDensity(cd);
    using TRun=std::remove_cvref_t<decltype(*cd)>;
    using S=std::conditional_t<std::is_base_of_v<qchem::ChargeDensity::tChargeDensity<dcmplx>,TRun>, dcmplx, double>;
    qchem::tests::pw::OneBlockBS<S> bs(b);
    t.BuildBlocks(&bs, cd);
}
}

TEST(RealComplexTerms, FactoredRhoMatchesFullQuadraticFormBothScalars)
{
    Rig rig;
    // Three arbitrary interior points -- this is about the CONTRACTION algebra, not the quadrature, so a
    // tiny mesh keeps it a unit gate.  Weight 1 each: applyRaw/applyRawFactored return rho itself.
    rvec3vec_t pts(3);
    pts[0]=rig.cell.ToCartesian(rvec3_t(0.3,0.4,0.7));
    pts[1]=rig.cell.ToCartesian(rvec3_t(0.5,0.5,0.5));
    pts[2]=rig.cell.ToCartesian(rvec3_t(0.1,0.9,0.2));
    DeltaProjector dp(BasisSet::FitQuadrature{std::make_shared<const qcMesh::Mesh>(pts, rvec_t(pts.size(),1.0)), {}});
    const BasisSet::DeltaFit_IBS<dcmplx>& dfit=dp.Basis();

    ExpectFactoredMatchesFull<double>(dfit, *rig.re, "ForwardFactoredT<double> (the real TRIM block)");
    ExpectFactoredMatchesFull<dcmplx>(dfit, *rig.cx, "ForwardFactoredT<dcmplx> (the Bloch block)");
}

// ...and the SECOND half of R1.0c, which is the half that bites: the gate above supplies its OWN factor,
// so it pins the CONTRACTION pair but never reaches LowRankFactor -- the code that BUILDS L from D on the
// production default route (QCHEM_DM_LOWRANK is on unless explicitly disabled).  So drive the two leaves
// against the same D and demand the same rho: FactoredRho (pivoted Cholesky, then the thin contraction)
// against the plain periodic leaf (the full quadratic form).
//
// It is deliberately END-TO-END rather than a direct L L^dag == D check, because LowRankFactor lives in an
// implementation unit and, more to the point, rho is the quantity anything downstream actually consumes.
// NB the existing PSD guard inside LowRankFactor cannot stand in for this: it compares Tr(L L^H) to Tr(D)
// through |L|^2, and a conjugation error does not move a modulus.
TEST(RealComplexTerms, FactoredLeafMatchesDirectLeafBothScalars)
{
    using namespace qchem::ChargeDensity;
    rvec3vec_t pts(4);
    pts[0]=rig_pt(0.3,0.4,0.7); pts[1]=rig_pt(0.5,0.5,0.5);
    pts[2]=rig_pt(0.1,0.9,0.2); pts[3]=rig_pt(0.62,0.15,0.44);
    Rig rig;
    for (size_t g=0; g<pts.size(); g++) pts[g]=rig.cell.ToCartesian(pts[g]);
    DeltaProjector dp(BasisSet::FitQuadrature{std::make_shared<const qcMesh::Mesh>(pts, rvec_t(pts.size(),1.0)), {}});

    auto arm=[&](auto tag, const auto& orb, const char* what)
    {
        using U = decltype(tag);
        const size_t n=orb.GetNumFunctions(), r=3;
        ASSERT_GT(n, 2*r) << what << ": D must be low-rank enough that FactoredRho does not fall back";
        mat_t<U> L(n,r);
        for (size_t i=0;i<n;i++) for (size_t m=0;m<r;m++) L(i,m)=FactorEntry<U>(i,m);
        hmat_t<U> D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
        {
            U s=U(0.0);
            for (size_t m=0;m<r;m++)
            {
                if constexpr (std::is_same_v<U,dcmplx>) s += L(i,m)*std::conj(L(j,m));
                else                                    s += L(i,m)*L(j,m);
            }
            if (i==j) s=U(std::real(std::complex<double>(s)));
            D(i,j)=s;
        }
        FactoredRho<PeriodicIrrepCD<U>> fac(D, &orb, orb.GetIrrep(Spin::None));
        PeriodicIrrepCD<U>              dir(D, &orb, orb.GetIrrep(Spin::None));
        // The factored leaf reports its HEALTH where the caller's context is (E-c (b), 2026-11-24): one density, factored, rank r of n.
        const std::string run=std::string("factored health ")+what;
        rvec_t rf;
        { qchem::report::Run filed(run); rf=fac.Project(dp()).Vector(); }
        const rvec_t rd=dir.Project(dp()).Vector();
        ASSERT_EQ(rf.size(), rd.size()) << what;
        {
            const auto& G=qchem::report::GlobalReport();
            const qchem::report::json* last=nullptr;
            for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(run+"@",0)==0) last=&it.value();
            ASSERT_TRUE(last && last->contains("densityFactor")) << what << ": the factored leaf files a densityFactor record";
            const auto& h=(*last)["densityFactor"];
            EXPECT_EQ(h.value("densities", 0.0), 1.0)              << what;
            EXPECT_EQ(h.value("factored", 0.0), 1.0)               << what << ": rank 3 of n pays";
            EXPECT_EQ(h.value("fallbackRank", 0.0), 0.0)           << what;
            // The ACHIEVED rank is D's true rank: the real arm's columns sin(x_i + 1.3m)/(1+0.1i) span only {sin x_i, cos x_i}/(1+0.1i)
            // (sin(x+c) = sin x cos c + cos x sin c), so its D has rank 2; the complex arm's cos imaginary part lifts it to r = 3.
            const double trueRank = std::is_same_v<U,double> ? 2.0 : double(r);
            EXPECT_EQ(h.value("rankSum", 0.0), trueRank)           << what << ": the achieved pivoted-Cholesky rank";
            EXPECT_EQ(h.value("functionSum", 0.0), double(n))      << what;
        }
        double scale=0.0;
        for (size_t g=0; g<rd.size(); g++) scale=std::max(scale, std::fabs(rd[g]));
        EXPECT_GT(scale, 1e-6) << what << ": a vacuous rho would pass anything";
        for (size_t g=0; g<rd.size(); g++)
            EXPECT_NEAR(rf[g], rd[g], 1e-11*std::max(scale,1.0)) << what << " at point " << g;
    };
    arm(double(0.0), *rig.re, "FactoredRho<PeriodicIrrepCD<double>> (the real TRIM block)");
    arm(dcmplx(0.0), *rig.cx, "FactoredRho<PeriodicIrrepCD<dcmplx>> (Bloch at GAMMA -- Phi is real)");

    // ...AND A BLOCK WHERE Phi IS GENUINELY COMPLEX, which is the only arm that discriminates.  Every k on
    // a 2x2x2 mesh is TRIM (2k is a reciprocal lattice vector), so rig.cx carries REAL Phi stored in a
    // complex matrix -- and with Im(Phi)=0 a transposed D gives the SAME rho, because
    //     rho' - rho = sum_{i<j} 4 Im(D_ij) Im(Phi_gi conj(Phi_gj)).
    // k=(1/4,0,0) on a 4x4x4 mesh is not TRIM, so the second factor is nonzero and a conjugation error in
    // the factor becomes visible.  This is the arm that caught R2.21.
    tGPW_IBS<dcmplx> kb(std::make_shared<::qchem::BasisSet::BasisCache>(), rig.cell, MakeMeshMaker(rig.cell), ivec3_t(4,4,4), ivec3_t(1,0,0), rig.mol, /*densityEcut*/20.0);
    arm(dcmplx(0.0), kb, "FactoredRho<PeriodicIrrepCD<dcmplx>> (Bloch at k=(1/4,0,0) -- Phi complex)");
}

// THE NOT-PSD CANARY IS PER BLOCK (D-MUTABLE-STATICS, 2026-11-23).  An indefinite D makes FactoredRho fall back to the full-rank
// route and print one warning.  It was a process-wide cap of 5 warnings, so a test binary or a geometry series went SILENT after its
// fifth non-PSD block.  Now a block warns once when it leaves the cone (not again for a later serial that is still outside), and
// every other block keeps its own canary.
TEST(RealComplexTerms, TheNotPSDCanaryWarnsOncePerBlockNotOncePerProcess)
{
    using namespace qchem::ChargeDensity;
    rvec3vec_t pts(2);
    pts[0]=rig_pt(0.3,0.4,0.7); pts[1]=rig_pt(0.5,0.5,0.5);
    Rig rig;
    for (size_t g=0; g<pts.size(); g++) pts[g]=rig.cell.ToCartesian(pts[g]);
    DeltaProjector dp(BasisSet::FitQuadrature{std::make_shared<const qcMesh::Mesh>(pts, rvec_t(pts.size(),1.0)), {}});
    const auto& orb=*rig.re;
    const size_t n=orb.GetNumFunctions(), r=3;
    ASSERT_GT(n, r+1);
    hmat_t<double> D(n);   // rank 3 PSD, minus a small negative weight on the last function: outside the cone
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
    { double s=0.0; for (size_t m=0;m<r;m++) s+=FactorEntry<double>(i,m)*FactorEntry<double>(j,m); D(i,j)=s; }
    D(n-1,n-1)-=1e-3;
    const std::string what="D is not positive semi-definite";
    auto count=[&](const std::string& s){ size_t c=0; for (size_t p=s.find(what); p!=std::string::npos; p=s.find(what,p+1)) c++; return c; };

    FactoredRho<PeriodicIrrepCD<double>> one(D, &orb, orb.GetIrrep(Spin::None));
    testing::internal::CaptureStdout();
    (void)one.Project(dp());
    one.ReScale(2.0);              // a new serial of the same block (1.0 is an exact no-op), still outside the cone
    (void)one.Project(dp());
    for (int b=0; b<6; b++)        // six more blocks: the old process-wide cap (5) would have silenced the last ones
    {
        FactoredRho<PeriodicIrrepCD<double>> other(D, &orb, orb.GetIrrep(Spin::None));
        (void)other.Project(dp());
    }
    const std::string out=testing::internal::GetCapturedStdout();
    EXPECT_EQ(count(out), 7u) << "one warning for the first block over its two serials, one for each of the six others:\n" << out;
}

// The STATIC set: kinetic + the three PP terms.  One capability cross-cast each, then bitwise equality.
TEST(RealComplexTerms, StaticTermsServeTheRealBlockBitwise)
{
    Rig rig;
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    Kinetic<dcmplx> kin;
    Ven_PP_Short    shrt(rig.st, &gth.local, 4.0*rig.st->GetNumAtoms());
    Ven_PP_Long     lng (rig.st, &gth.local, 4.0*rig.st->GetNumAtoms());
    Ven_PP_NonLocal nl  (rig.st, &gth.nonlocal);

    const std::vector<std::pair<const cStatic_HT*,const char*>> terms =
        {{&kin,"Kinetic"},{&shrt,"Ven_PP_Short"},{&lng,"Ven_PP_Long"},{&nl,"Ven_PP_NonLocal"}};
    for (auto [t,name] : terms)
    {
        auto* rb=dynamic_cast<const tStatic_CC<double>*>(t);
        ASSERT_NE(rb,nullptr) << name << " must carry the real-block capability (Step 3c)";
        ExpectBitwiseEqual(rb->GetMatrix(rig.re.get(), Spin::None),
                           t ->GetMatrix(rig.cx.get(), Spin::None), name);
    }
}

// The DYNAMIC pair reachable without the full SCF: Hartree (G-space Poisson off the shared V_H(G)) and
// the Becke-route cVxc_DFT (the engine's rho raster is block-independent; the real block gets its
// own REAL Phi table, so its quadrature GEMM runs in real arithmetic -- the first realized Step-3 win).
TEST(RealComplexTerms, HartreeAndBeckeXcServeTheRealBlockBitwise)
{
    Rig rig;
    auto cd  = rig.MakeDensity();       // the complex arm's density (complex child on cx)
    auto cdr = rig.MakeRealDensity();   // the real arm's TWIN (real child on re -- the 3c-3 mixed shape)

    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(rig.cx->CreateCDFitBasisSet(rig.st.get()));
    Vee_Fit<dcmplx> vh(Fitting::MakeGSpaceDensityFitter(fb), 0.0, fb);
    {
        auto* rb=dynamic_cast<const tDynamic_CC<double>*>(&vh);
        ASSERT_NE(rb,nullptr) << "the periodic Hartree term must carry the real-block capability (Step 3c)";
        // Concrete-typed calls are ambiguous across the two GetMatrix base classes; go through the faces.
        RefreshTerm(static_cast<const cDynamic_HT&>(vh), cd.get(), rig.cx.get());    // each arm refreshed for ITS density (KAxisPlan P2)
        const chmat_t Vc=static_cast<const cDynamic_HT&>(vh).GetMatrix(rig.cx.get(), Spin::None, DensityID(cd.get()));
        RefreshTerm(static_cast<const cDynamic_HT&>(vh), cdr.get(), rig.re.get());
        const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, DensityID(cdr.get()));
        ExpectBitwiseEqual(Vr, Vc, "Vee_Fit<dcmplx>");
    }

    auto q = MakeMeshMaker(rig.cell)->Make(qcMesh::MeshParams{}, qcMesh::NoSymmetry());   // a bare block: no imposition, the plain mesh
    auto dfb=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), q, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
    // Through the FACTORY: a delta fit basis can only resolve to SINGLES, so naming the strategy bought
    // nothing but an Internal import (user, 2026-09-10).
    auto engine=ChargeDensity::MakeDensitySampler(dfb, q);   // ONE bundle to both collaborators
    cVxc_DFT vxc(std::make_shared<SlaterExchange>(2.0/3.0), engine, SpinGroup::UnPolarized);
    {
        auto* rb=dynamic_cast<const tDynamic_CC<double>*>(&vxc);
        ASSERT_NE(rb,nullptr) << "cVxc_DFT must carry the real-block capability (Step 3c)";
        RefreshTerm(static_cast<const cDynamic_HT&>(vxc), cd.get(), rig.cx.get());    // each arm refreshed for ITS density (KAxisPlan P2)
        const chmat_t Vc=static_cast<const cDynamic_HT&>(vxc).GetMatrix(rig.cx.get(), Spin::None, DensityID(cd.get()));
        RefreshTerm(static_cast<const cDynamic_HT&>(vxc), cdr.get(), rig.re.get());
        const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, DensityID(cdr.get()));
        ExpectMachineEqual(Vr, Vc, "cVxc_DFT");
    }
}


// PRODUCTION-SHAPED rig (Step 3c-3 divergence hunt): the FCC-diamond 2-atom cell the acceptance runs,
// where the single-atom cubic Rig above cannot reach multi-centre / oblique-cell defects.
struct FccRig
{
    FCCUnitCell cell;
    std::shared_ptr<const Structure> st;
    std::shared_ptr<const Real_BS>   mol;
    std::unique_ptr<tGPW_IBS<double>> re;
    std::unique_ptr<tGPW_IBS<dcmplx>> cx;
    FccRig() : cell(10.26)
    {
        cell.AddAtom(14,{0,0,0});
        cell.AddAtom(14,{0.25,0.25,0.25});
        st=cell.Clone();
        namespace M = qchem::BasisSet::Gaussian;
        mol.reset(M::Factory(M::BasisSetData::SIPP_SR, &cell, M::Engine::MnD, M::Angular::Cartesian));
        re=std::make_unique<tGPW_IBS<double>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
        cx=std::make_unique<tGPW_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    }
};

TEST(RealComplexTerms, StaticTermsServeTheRealBlockBitwise_FccDiamond)
{
    FccRig rig;
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    Kinetic<dcmplx> kin;
    Ven_PP_Short    shrt(rig.st, &gth.local, 4.0*rig.st->GetNumAtoms());
    Ven_PP_Long     lng (rig.st, &gth.local, 4.0*rig.st->GetNumAtoms());
    Ven_PP_NonLocal nl  (rig.st, &gth.nonlocal);

    const std::vector<std::pair<const cStatic_HT*,const char*>> terms =
        {{&kin,"Kinetic"},{&shrt,"Ven_PP_Short"},{&lng,"Ven_PP_Long"},{&nl,"Ven_PP_NonLocal"}};
    for (auto [t,name] : terms)
    {
        auto* rb=dynamic_cast<const tStatic_CC<double>*>(t);
        ASSERT_NE(rb,nullptr) << name;
        ExpectBitwiseEqual(rb->GetMatrix(rig.re.get(), Spin::None),
                           t ->GetMatrix(rig.cx.get(), Spin::None), name);
    }
}

// The SEED-SHAPED density on the production cell (3c-3 hunt): the run's two arms feed the SAME uniform
// D = (N/n)·I through DIFFERENT leaves (real vs complex).  Their collocated rho and V_H(G) must agree
// bitwise, or the very first Fock differs.
TEST(RealComplexTerms, SeedDensityTrioMatches_FccDiamond)
{
    using namespace qchem::ChargeDensity;
    FccRig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); chmat_t Dc(n);
    for (size_t i=0;i<n;i++) { Dr(i,i)=8.0/double(n); Dc(i,i)=8.0/double(n); }
    PeriodicIrrepCD<double> cdr(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None));
    PeriodicIrrepCD<dcmplx> cdc(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None));
    EXPECT_NEAR(cdr.GetTotalCharge(), cdc.GetTotalCharge(), 1e-13);

    std::shared_ptr<const BasisSet::cFIT_SF_ABS> sf(rig.cx->CreateVxcFitBasisSet(rig.st.get(), nullptr, 1.0));
    // ρ̃(dm) through the one projection face (V-ONE-OF-N 3c; V_H = k·ρ̃ is the Hartree fitter's, so this was the V_H claim too).
    const auto proj=BasisSet::MakeGSpaceDensityProjector(sf);
    const ΔG_Map rr=cdr.Project(*proj).GMap(), rc=cdc.Project(*proj).GMap();
    ASSERT_EQ(rr.size(), rc.size());
    for (const auto& [dm,v] : rr)
    {
        auto it=rc.find(dm);
        ASSERT_NE(it, rc.end());
        EXPECT_EQ(v, it->second) << "rho-tilde(dm) mismatch";
    }
    const auto raster=BasisSet::MakeRasterDensityProjector(sf);   // the raw ρ_DM through the same face (R1.0n/o)
    const rvec_t gr=cdr.Project(*raster).Vector(), gc=cdc.Project(*raster).Vector();
    ASSERT_EQ(gr.size(), gc.size());
    for (size_t g=0; g<gr.size(); g++) EXPECT_EQ(gr[g], gc[g]) << "raw rho_DM mismatch at g="<<g;

    // ...and the SAME through the run-shaped COMPOSITES (the exact seed objects the two arms feed the
    // raw XC route): a mixed composite's rho visit must reach the real child identically.
    tComposite_CD<dcmplx> mixed;  mixed.Insert(std::unique_ptr<tDM_CD<double>>(
        new PeriodicIrrepCD<double>(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None))), rig.re->GetIrrep(Spin::None));
    tComposite_CD<dcmplx> cplx;   cplx.Insert(std::unique_ptr<tDM_CD<dcmplx>>(
        new PeriodicIrrepCD<dcmplx>(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None))), rig.cx->GetIrrep(Spin::None));
    const rvec_t cgr=mixed.Project(*raster).Vector(), cgc=cplx.Project(*raster).Vector();
    ASSERT_EQ(cgr.size(), cgc.size());
    for (size_t g=0; g<cgr.size(); g++) EXPECT_EQ(cgr[g], cgc[g]) << "composite raw rho mismatch at g="<<g;
}

// R1.0n/o (2026-11-19): the raw ρ_DM on the fit's raster, asked through the one projection face (MakeRasterDensityProjector), is the
// block's own raw collocation (applyRaw on the block's D) BIT FOR BIT -- for the real TRIM leaf, the Bloch leaf, and the run-shaped
// composites over each (one block, no group: the sum is that block's).  It replaced the interim RasterDensity::GetRhoOnGrid, which
// was exactly this call; a seed answers none (empty).
TEST(RealComplexTerms, TheRasterProjectorIsTheBlocksRawCollocationBitForBit)
{
    using namespace qchem::ChargeDensity;
    FccRig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); chmat_t Dc(n);
    for (size_t i=0;i<n;i++) { Dr(i,i)=8.0/double(n); Dc(i,i)=8.0/double(n); }
    Dr(0,1)=Dr(1,0)=0.25; Dc(0,1)=Dc(1,0)=0.25;   // an off-diagonal pair, so the collocation sums pairs
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> sf(rig.cx->CreateVxcFitBasisSet(rig.st.get(), nullptr, 1.0));
    const auto raster=BasisSet::MakeRasterDensityProjector(sf);
    EXPECT_NE(dynamic_cast<const BasisSet::RasterPairCollocation*>(rig.cx.get()), nullptr) << "a GPW block has the raw route";

    auto same=[](const rvec_t& a, const rvec_t& b, const char* what)
    {
        ASSERT_EQ(a.size(), b.size()) << what;
        ASSERT_GT(a.size(), 0u) << what;
        for (size_t g=0; g<a.size(); g++) ASSERT_EQ(a[g], b[g]) << what << " at g=" << g;
    };
    PeriodicIrrepCD<double> cdr(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None));
    PeriodicIrrepCD<dcmplx> cdc(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None));
    const rvec_t rawR=static_cast<const BasisSet::Orbital_DFT_IBS<double,dcmplx>&>(*rig.re).Overlap3C(*sf).applyRaw(Dr);
    const rvec_t rawC=static_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(*rig.cx).Overlap3C(*sf).applyRaw(Dc);
    same(cdr.Project(*raster).Vector(), rawR, "real TRIM leaf");
    same(cdc.Project(*raster).Vector(), rawC, "Bloch leaf");

    tComposite_CD<dcmplx> mixed;  mixed.Insert(std::unique_ptr<tDM_CD<double>>(
        new PeriodicIrrepCD<double>(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None))), rig.re->GetIrrep(Spin::None));
    tComposite_CD<dcmplx> cplx;   cplx.Insert(std::unique_ptr<tDM_CD<dcmplx>>(
        new PeriodicIrrepCD<dcmplx>(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None))), rig.cx->GetIrrep(Spin::None));
    same(mixed.Project(*raster).Vector(), rawR, "mixed composite");
    same(cplx .Project(*raster).Vector(), rawC, "complex composite");
}

// The RAW-route pair-quadrature gate (Step 3c-3): the production Uniform-mesh XC.  The collocated rho_DM
// feed is block-independent; the real block's matrix comes back through the SAME raw adjoint tensor
// (TFit==dcmplx either way) narrowed at the end, so it must equal the complex block's real part
// BITWISE.  (This was the one term route 3c-1 could not gate without the full collocation feed.)
TEST(RealComplexTerms, RawRouteXcServesTheRealBlockBitwise)
{
    Rig rig;
    auto cd  = rig.MakeDensity();       // the complex arm's density
    auto cdr = rig.MakeRealDensity();   // the real arm's twin (screen-consistent pairing; file header)
    ChargeDensity::fitbasis_t fb(rig.cx->CreateVxcFitBasisSet(rig.st.get(), nullptr, 1.0));
    cVxc_DFT vxc(std::make_shared<SlaterExchange>(2.0/3.0), ChargeDensity::MakeDensitySampler(fb), SpinGroup::UnPolarized);
    auto* rb=dynamic_cast<const tDynamic_CC<double>*>(&vxc);
    ASSERT_NE(rb,nullptr) << "cVxc_DFT must carry the real-block capability (Step 3c)";
    RefreshTerm(static_cast<const cDynamic_HT&>(vxc), cd.get(), rig.cx.get());    // each arm refreshed for ITS density (KAxisPlan P2)
    const chmat_t Vc=static_cast<const cDynamic_HT&>(vxc).GetMatrix(rig.cx.get(), Spin::None, DensityID(cd.get()));
    RefreshTerm(static_cast<const cDynamic_HT&>(vxc), cdr.get(), rig.re.get());
    const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, DensityID(cdr.get()));
    ExpectBitwiseEqual(Vr, Vc, "cVxc_DFT (pair) raw");
}

// The ASSEMBLY gate (Step 3c-2): the Hamiltonian's Ham_RealBlock fold over the term faces must equal the
// native complex assembly's real part BITWISE -- each term's block already does (the 3c-1 gates), and the
// fold accumulates them elementwise in the same order.  This is the exact matrix a tIrrepWF<double> child
// receives from its CalculateH inside a complex run.
namespace { struct MixedHam : Hamiltonian::tHamiltonianImp<dcmplx> { MixedHam() : tHamiltonianImp<dcmplx>(SpinGroup::UnPolarized) {} }; }

TEST(RealComplexTerms, HamiltonianAssemblyServesTheRealBlockBitwise)
{
    Rig rig;
    auto cd  = rig.MakeDensity();       // the complex arm's density
    auto cdr = rig.MakeRealDensity();   // the real arm's twin (screen-consistent pairing; file header)
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    MixedHam ham;
    ham.Add(new Kinetic<dcmplx>);
    ham.Add(new Ven_PP_Short   (rig.st, &gth.local, 4.0*rig.st->GetNumAtoms()));
    ham.Add(new Ven_PP_Long    (rig.st, &gth.local, 4.0*rig.st->GetNumAtoms()));
    ham.Add(new Ven_PP_NonLocal(rig.st, &gth.nonlocal));
    {
        std::shared_ptr<const BasisSet::cFIT_CD_ABS> hfb(rig.cx->CreateCDFitBasisSet(rig.st.get()));
        ham.Add(new Vee_Fit<dcmplx>(Fitting::MakeGSpaceDensityFitter(hfb), 0.0, hfb));
    }

    const qchem::tests::pw::OneBlockBS<dcmplx> cxBS(rig.cx.get()), reBS(rig.re.get());   // the block each arm reads: BuildBlocks builds it (V1.41 stage 3)
    ham.RefreshForDensity(&cxBS, cd.get());    // the contract (KAxisPlan P2): each arm refreshed for ITS density
    const chmat_t Hc = ham.GetMatrix(rig.cx.get(), Spin::None, cd.get(), nullptr);
    auto& rb = dynamic_cast<Hamiltonian::Ham_RealBlock&>(ham);   // the face a real WF child drives
    ham.RefreshForDensity(&reBS, cdr.get());
    const rsmat_t Hr = rb.GetMatrix(rig.re.get(), Spin::None, cdr.get(), nullptr);
    ExpectBitwiseEqual(Hr, Hc, "Ham_RealBlock fold");
}


// ===== Step 3c-2b: the PERIODIC LEAF SPLIT + the mixed-composite energy/density arms ==================
// PeriodicIrrepCD<double> is the real TRIM block's density: real D, real DM-side arithmetic, and the
// SAME reciprocal trio as its complex twin -- the lineage-as-class split (V1.7's scalar⇔lineage
// conflation retired).  The trio must match the complex twin BITWISE: the contraction runs the same
// tensor in the same order, only D's (exactly zero) imaginary parts differ.
TEST(RealComplexTerms, RealPeriodicLeafFourierTrioMatchesComplexBitwise)
{
    using namespace qchem::ChargeDensity;
    Rig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); Dr(0,0)=2.0;
    chmat_t Dc(n); Dc(0,0)=2.0;
    PeriodicIrrepCD<double> cdr(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None));
    PeriodicIrrepCD<dcmplx> cdc(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None));

    // NEITHER leaf carries the finite (HF) faces -- the honest probe.
    EXPECT_EQ(dynamic_cast<const tHF_Pair_CD<double>*>(&cdr), nullptr);

    std::shared_ptr<const BasisSet::cFIT_SF_ABS> sf(rig.cx->CreateVxcFitBasisSet(rig.st.get(), nullptr, 1.0));
    // ρ̃(dm) through the one projection face (V-ONE-OF-N 3c; V_H = k·ρ̃ is the Hartree fitter's, so this was the V_H claim too).
    const auto proj=BasisSet::MakeGSpaceDensityProjector(sf);
    const ΔG_Map rr=cdr.Project(*proj).GMap(), rc=cdc.Project(*proj).GMap();
    ASSERT_EQ(rr.size(), rc.size());
    for (const auto& [dm,v] : rr)
    {
        auto it=rc.find(dm);
        ASSERT_NE(it, rc.end());
        EXPECT_EQ(v, it->second) << "rho-tilde(dm) mismatch";
    }
    const auto raster=BasisSet::MakeRasterDensityProjector(sf);   // the raw ρ_DM (R1.0n/o)
    const rvec_t gr=cdr.Project(*raster).Vector(), gc=cdc.Project(*raster).Vector();
    ASSERT_EQ(gr.size(), gc.size());
    ASSERT_GT(gr.size(), 0u);
    for (size_t g=0; g<gr.size(); g++) EXPECT_EQ(gr[g], gc[g]) << "raw rho_DM mismatch at g="<<g;
}

// The MIXED composite's energy/density arms: ONE real child inside a complex-faced composite must give
// the same static/dynamic energy contractions and the same rho(points) as the all-complex twin.
TEST(RealComplexTerms, MixedCompositeEnergyAndRhoMatchComplex)
{
    using namespace qchem::ChargeDensity;
    Rig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); Dr(0,0)=2.0;
    chmat_t Dc(n); Dc(0,0)=2.0;

    tComposite_CD<dcmplx> mixed;   // complex face, REAL child (the 3c-3 shape)
    mixed.Insert(std::unique_ptr<tDM_CD<double>>(
        new PeriodicIrrepCD<double>(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None))), rig.re->GetIrrep(Spin::None));
    tComposite_CD<dcmplx> cplx;    // the all-complex twin
    cplx.Insert(std::unique_ptr<tDM_CD<dcmplx>>(
        new PeriodicIrrepCD<dcmplx>(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None))), rig.cx->GetIrrep(Spin::None));

    // STATIC arm: the term's real-block face IS a tStatic_CC<double>, so the real child contracts natively.
    Kinetic<dcmplx> kin;
    EXPECT_EQ(mixed.DM_Contract(static_cast<const cStatic_CC*>(&kin)),
              cplx .DM_Contract(static_cast<const cStatic_CC*>(&kin)));

    // DYNAMIC arm: the real block face tDynamic_CC<double>::GetEMatrix, which the periodic Hartree term BUILDS from its refreshed fit (the
    // energy pass builds no blocks, V1.41 stage 3), against the complex one.
    auto cdrun = rig.MakeDensity();   // the run's density (defines V_H)
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> vfb(rig.cx->CreateCDFitBasisSet(rig.st.get()));
    Vee_Fit<dcmplx> vh(Fitting::MakeGSpaceDensityFitter(vfb), 0.0, vfb);
    static_cast<const cDynamic_HT&>(vh).RefreshForDensity(cdrun.get());   // the contract (KAxisPlan P2): the energy pass's refresh
    EXPECT_EQ(mixed.DM_Contract(static_cast<const cDynamic_CC*>(&vh), cdrun.get()),
              cplx .DM_Contract(static_cast<const cDynamic_CC*>(&vh), cdrun.get()));

    // RHO arm: same rho(r) whether the run is mixed-scalar or all-complex.  Since 2026-08-22 each child
    // asks the QUADRATURE for its own typed table (the real child a real one, its complex siblings complex
    // ones) instead of the caller threading two table maps -- so this arm now exercises the GEMM path on
    // both sides rather than the retired pointwise fallback, which is the stronger test of the same claim.
    rvec3vec_t pts(3);
    pts[0]=rig.cell.ToCartesian(rvec3_t(0.3,0.4,0.7));
    pts[1]=rig.cell.ToCartesian(rvec3_t(0.5,0.5,0.5));
    pts[2]=rig.cell.ToCartesian(rvec3_t(0.1,0.9,0.2));
    DeltaProjector qpts(BasisSet::FitQuadrature{std::make_shared<const qcMesh::Mesh>(pts, rvec_t(pts.size(),1.0)), {}});
    const rvec_t rm=mixed.Project(qpts()).Vector();
    const rvec_t rc=cplx .Project(qpts()).Vector();
    for (size_t g=0; g<pts.size(); g++) EXPECT_EQ(rm[g], rc[g]) << "rho mismatch at point "<<g;

    // And the composite's ρ̃ (the generic visit, through the one projection face) reaches the REAL child too.
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> sf2(rig.cx->CreateVxcFitBasisSet(rig.st.get(), nullptr, 1.0));
    const auto proj2=BasisSet::MakeGSpaceDensityProjector(sf2);
    const ΔG_Map vm=mixed.Project(*proj2).GMap(), vc=cplx.Project(*proj2).GMap();
    ASSERT_EQ(vm.size(), vc.size());
    for (const auto& [dm,v] : vm) { auto it=vc.find(dm); ASSERT_NE(it,vc.end()); EXPECT_EQ(v,it->second); }
}


// ===== THE FORCE WALK ON A MIXED RUN (GeometryOptimizationPlan §13.2; lane B 2026-11-09) =====================================================
// A real TRIM block inside a complex run is walked through the term's REAL-BLOCK pullback client, tStatic_PC<double>, as DM_Contract reaches
// tStatic_CC<double>.  The walk's arithmetic is gated on a real block whose pullbacks the molecular gates hold (the GPW real TRIM
// block has its own since 2026-11-10, ARealTRIMGPWBlockGivesKineticAndPulayForceAndVirialAndNamesTheRest): a molecular Gaussian block (MnD, whose pullbacks ride libcint) posing as the real child of a complex-faced basis.  The claim is the WIRING: the mixed walk
// must equal the same block's walk in an all-real run (whose Kinetic and Pulay pullbacks the molecular gates hold to the finite difference),
// bit for bit, since it is the same body on the same scalar.
namespace
{
struct RealChildRig
{
    UnitCell cell{8.0};
    std::shared_ptr<const Real_BS> mol;
    const BasisSet::Orbital_1E_IBS<double>* blk=nullptr;
    rsmat_t D;
    RealChildRig()
    {
        cell.AddAtom(14,{0.40,0.50,0.50});
        cell.AddAtom(14,{0.62,0.55,0.47});   // two atoms, off any symmetric arrangement: every per-atom pullback is nonzero
        namespace M = qchem::BasisSet::Gaussian;
        mol.reset(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));   // MnD: its pullbacks ride libcint
        blk=(*mol)[0];
        const size_t n=blk->GetNumFunctions();
        D=rsmat_t(n);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) D(i,j)=D(j,i)=std::sin(1.0+i+2.0*j)/(1.0+i+j);   // symmetric, every entry nonzero
    }
    Irrep Q() const {return blk->GetIrrep(Spin::None);}
    std::unique_ptr<ChargeDensity::tComposite_CD<dcmplx>> Mixed() const
    {
        auto cd=std::make_unique<ChargeDensity::tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<ChargeDensity::tDM_CD<double>>(new ChargeDensity::PeriodicIrrepCD<double>(D, blk, Q())), Q());
        return cd;
    }
    std::unique_ptr<ChargeDensity::tComposite_CD<double>> Real() const
    {
        auto cd=std::make_unique<ChargeDensity::tComposite_CD<double>>();
        cd->Insert(std::unique_ptr<ChargeDensity::tDM_CD<double>>(new ChargeDensity::PeriodicIrrepCD<double>(D, blk, Q())), Q());
        return cd;
    }
};
struct RealHam : Hamiltonian::tHamiltonianImp<double> { RealHam() : tHamiltonianImp<double>(SpinGroup::UnPolarized) {} };
void ExpectSameForces(const rvec3vec_t& a, const rvec3vec_t& b, const char* what)
{
    ASSERT_EQ(a.size(), b.size()) << what;
    double big=0.0;
    for (size_t A=0;A<a.size();A++)
    {
        EXPECT_EQ(a[A].x, b[A].x) << what << " atom " << A;
        EXPECT_EQ(a[A].y, b[A].y) << what << " atom " << A;
        EXPECT_EQ(a[A].z, b[A].z) << what << " atom " << A;
        big=std::max({big, std::fabs(b[A].x), std::fabs(b[A].y), std::fabs(b[A].z)});
    }
    EXPECT_GT(big, 1e-3) << what << ": a vanishing force cannot show the walk is wired";
}
} // namespace

TEST(RealComplexTerms, TheForceWalksARealTRIMBlockThroughItsRealPullbackFace)
{
    const RealChildRig rig;
    const auto mixed=rig.Mixed();
    const auto real =rig.Real();
    const size_t nA=rig.cell.GetNumAtoms();

    // The term: Kinetic<dcmplx> on the mixed composite against Kinetic<double> on the all-real one.
    Kinetic<dcmplx> kc;
    Kinetic<double> kr;
    rvec3vec_t Fc(nA, rvec3_t(0,0,0)), Fr(nA, rvec3_t(0,0,0));
    kc.AddForces(Fc, mixed.get(), &rig.cell);
    kr.AddForces(Fr, real.get(),  &rig.cell);
    ExpectSameForces(Fc, Fr, "Kinetic");

    // The Hamiltonian: MissingGradients walks the real child (no throw, nothing missing), and MakeForces' Pulay sum reaches it (W = D).
    const qchem::tests::pw::OneBlockBS<dcmplx> mixedBS(rig.blk);
    const qchem::tests::pw::OneBlockBS<double> realBS(rig.blk);
    MixedHam hc; hc.Add(new Kinetic<dcmplx>);
    RealHam  hr; hr.Add(new Kinetic<double>);
    EXPECT_TRUE(hc.MissingGradients(&mixedBS, &rig.cell).empty());
    ExpectSameForces(hc.MakeForces(&mixedBS, mixed.get(), mixed.get(), &rig.cell),
                     hr.MakeForces(&realBS,  real.get(),  real.get(),  &rig.cell), "Kinetic + Pulay");

    // The PW PP terms' real face: the walk hands the real child to the block's own species-field / projector pullback.
    const auto gth=Pseudopotential::GetGTH("Si","LDA",4);
    const auto st=rig.cell.Clone();
    const auto* pp=dynamic_cast<const BasisSet::Orbital_PP_IBS<double>*>(rig.blk);
    ASSERT_NE(pp, nullptr);
    ASSERT_TRUE(pp->CanPullbackSpeciesField(gth.local, BasisSet::FieldRange::Short));
    ASSERT_TRUE(pp->CanPullbackProjectors(gth.nonlocal));
    Ven_PP_Short    vs(st, &gth.local, 4.0*st->GetNumAtoms());
    Ven_PP_NonLocal vn(st, &gth.nonlocal);
    EXPECT_TRUE(vs.MissingGradient(&mixedBS, st.get()).empty());
    EXPECT_TRUE(vn.MissingGradient(&mixedBS, st.get()).empty());
    ExpectSameForces(mixed->DM_Pullback(static_cast<const tStatic_PC<dcmplx>*>(&vs), st.get()),
                     *pp->MakeSpeciesFieldPullback(st.get(), gth.local, BasisSet::FieldRange::Short, rig.D), "Ven_PP_Short");
    ExpectSameForces(mixed->DM_Pullback(static_cast<const tStatic_PC<dcmplx>*>(&vn), st.get()),
                     *pp->MakeProjectorPullback(st.get(), gth.nonlocal, rig.D), "Ven_PP_NonLocal");
}

// THE GPW ONE-ELECTRON FORCE AND VIRIAL THROUGH BOTH BLOCK SCALARS (GeometryOptimizationPlan §13.2, doc/Records/CellHalf_Design.md §5;
// lane B 2026-11-10).  The real TRIM block a GPW run builds (tGPW_IBS<double>) now gives the kinetic and overlap pullbacks and virials
// (Enzyme over the MnD kernel per image, gated against the energy's finite difference in UTGaussian_BS M_LatticePullback), so Kinetic and the Pulay sum have
// their force and virial on it, and on the Bloch block, IDENTICALLY (the real block widens its D onto the same walk at a real phase).  Since
// 2026-11-19 the three GPW PP terms do too (lattice Enzyme walks and the long range's gathers, gated against their energies' finite
// differences in UTLattice_BS GPW_Force), so they join the identity; a local field WITHOUT the closed-Gaussian form still refuses by name.
TEST(RealComplexTerms, ARealTRIMGPWBlockGivesTheOneElectronAndPPForcesAndVirialsAndNamesTheRest)
{
    UnitCell cell(8.0);
    cell.AddAtom(14,{0.40,0.50,0.50});
    cell.AddAtom(14,{0.62,0.55,0.47});
    const auto st=cell.Clone();
    namespace M = qchem::BasisSet::Gaussian;
    std::shared_ptr<const Real_BS> mol(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
    tGPW_IBS<double> re(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    const size_t n=cx.GetNumFunctions();
    rsmat_t Dr(n); chmat_t Dc(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) { const double d=std::sin(1.0+i+2.0*j)/(1.0+i+j); Dr(i,j)=Dr(j,i)=d; Dc(j,i)=dcmplx(d,0.0); }
    using namespace qchem::ChargeDensity;
    auto cdc=std::make_unique<tComposite_CD<dcmplx>>();
    cdc->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(Dc, &cx, cx.GetIrrep(Spin::None))), cx.GetIrrep(Spin::None));
    auto cdr=std::make_unique<tComposite_CD<dcmplx>>();
    cdr->Insert(std::unique_ptr<tDM_CD<double>>(new PeriodicIrrepCD<double>(Dr, &re, re.GetIrrep(Spin::None))), re.GetIrrep(Spin::None));
    const qchem::tests::pw::OneBlockBS<dcmplx> cxBS(&cx), reBS(&re);

    const auto gth=Pseudopotential::GetGTH("Si","LDA",4);
    const double N=4.0*st->GetNumAtoms();
    Kinetic<dcmplx> kin;
    Ven_PP_Short    vs(st, &gth.local, N);
    Ven_PP_Long     vl(st, &gth.local, N);
    Ven_PP_NonLocal vn(st, &gth.nonlocal);
    for (const auto* bs : {&cxBS, &reBS})
    {
        EXPECT_EQ(kin.MissingGradient(bs, st.get()), "") << "a GPW block has the kinetic pullback";
        EXPECT_EQ(kin.MissingVirial  (bs, st.get()), "") << "a GPW block has the kinetic virial";
        for (const auto* t : {static_cast<const tGeometryGradient_HT<dcmplx>*>(&vs), static_cast<const tGeometryGradient_HT<dcmplx>*>(&vl),
                              static_cast<const tGeometryGradient_HT<dcmplx>*>(&vn)})
            EXPECT_EQ(t->MissingGradient(bs, st.get()), "") << "a GPW block has the PP pullbacks";
        for (const auto* t : {static_cast<const tStrainGradient_HT<dcmplx>*>(&vs), static_cast<const tStrainGradient_HT<dcmplx>*>(&vl),
                              static_cast<const tStrainGradient_HT<dcmplx>*>(&vn)})
            EXPECT_EQ(t->MissingVirial(bs, st.get()), "") << "a GPW block has the PP virials";
    }
    MixedHam ham; ham.Add(new Kinetic<dcmplx>); ham.Add(new Ven_PP_Short(st, &gth.local, N)); ham.Add(new Ven_PP_Long(st, &gth.local, N));
    ham.Add(new Ven_PP_NonLocal(st, &gth.nonlocal));
    EXPECT_TRUE(ham.MissingGradients(&cxBS, st.get()).empty()) << "Kinetic + PP + Pulay: nothing missing on a Bloch GPW block";
    EXPECT_TRUE(ham.MissingGradients(&reBS, st.get()).empty()) << "Kinetic + PP + Pulay: nothing missing on a real TRIM GPW block";
    EXPECT_TRUE(ham.MissingVirials  (&cxBS, st.get()).empty()) << "Kinetic + PP + Pulay virials: nothing missing on a Bloch GPW block";
    EXPECT_TRUE(ham.MissingVirials  (&reBS, st.get()).empty()) << "Kinetic + PP + Pulay virials: nothing missing on a real TRIM GPW block";

    const rvec3vec_t Fc=ham.MakeForces(&cxBS, cdc.get(), cdc.get(), st.get()), Fr=ham.MakeForces(&reBS, cdr.get(), cdr.get(), st.get());
    const rmat3d_t   Vc=ham.MakeVirial(&cxBS, cdc.get(), cdc.get(), st.get()), Vr=ham.MakeVirial(&reBS, cdr.get(), cdr.get(), st.get());
    double bigF=0.0, dF=0.0, bigV=0.0, dV=0.0;
    for (size_t A=0;A<2;A++)
    {
        bigF=std::max({bigF, std::fabs(Fc[A].x), std::fabs(Fc[A].y), std::fabs(Fc[A].z)});
        dF  =std::max({dF, std::fabs(Fc[A].x-Fr[A].x), std::fabs(Fc[A].y-Fr[A].y), std::fabs(Fc[A].z-Fr[A].z)});
    }
    for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) { bigV=std::max(bigV, std::fabs(Vc(a,b))); dV=std::max(dV, std::fabs(Vc(a,b)-Vr(a,b))); }
    std::printf("  [GPW kinetic + PP + Pulay] max|F| %.3e |F_c - F_r| %.2e   max|V| %.3e |V_c - V_r| %.2e\n", bigF, dF, bigV, dV);
    ASSERT_GT(bigF, 1e-3) << "a vanishing force cannot show the walk is wired";
    ASSERT_GT(bigV, 1e-3) << "a vanishing virial cannot show the walk is wired";
    EXPECT_LT(dF, 1e-12*bigF) << "the real TRIM block and the Bloch block give different kinetic + PP + Pulay forces";
    EXPECT_LT(dV, 1e-12*bigV) << "the real TRIM block and the Bloch block give different kinetic + PP + Pulay virials";

    // The STRESS is the virial over the cell volume, unsymmetrised (CellHalf_Design §5); with a term that has no virial it REFUSES, naming it.
    const rmat3d_t σ=ham.At(&cxBS, cdc.get()).Stress(cdc.get(), cell);
    const double Ω=cell.GetCellVolume();
    for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) EXPECT_DOUBLE_EQ(σ(a,b), Vc(a,b)/Ω) << "Stress = MakeVirial / Ω, component " << a << b;
    // A local field WITHOUT the closed-Gaussian form (a Gaussian-smeared nucleus) takes GPW's grid sweep, which has no derivative: the term
    // HAS the face, and its GPW block refuses ("...gives no species-field virial").
    const Pseudopotential::GaussianSmearedNucleus smeared(0.5);
    MixedHam withPP; withPP.Add(new Kinetic<dcmplx>); withPP.Add(new Ven_PP_Short(st, &smeared, N));
    const auto mv=withPP.MissingVirials(&cxBS, st.get());
    ASSERT_EQ(mv.size(), 1u) << "exactly the grid-swept short PP has no virial";
    EXPECT_NE(mv[0].find("Ven_PP_Short"), std::string::npos) << mv[0];
    EXPECT_NE(mv[0].find("virial"), std::string::npos) << mv[0];
    try { (void)withPP.At(&cxBS, cdc.get()).Stress(cdc.get(), cell); FAIL() << "a Stress with a missing virial must throw"; }
    catch (const std::logic_error& e) { EXPECT_NE(std::string(e.what()).find("MissingVirials"), std::string::npos) << e.what(); }
}

TEST(RealComplexTerms, AClientWithoutARealPullbackFaceThrowsByNameOnARealChild)
{
    struct ComplexOnly : tStatic_PC<dcmplx>
    {
        rvec3vec_t Pullback(const cobs_t*, const Spin&, const chmat_t&, const Structure*) const override {return {};}
    };
    const RealChildRig rig;
    const ComplexOnly client;
    try { rig.Mixed()->DM_Pullback(&client, &rig.cell); FAIL() << "a real child walked by a complex-only client must throw"; }
    catch (const std::logic_error& e) { EXPECT_NE(std::string(e.what()).find("tStatic_PC<double>"), std::string::npos) << e.what(); }
}

// THE GPW HARTREE FORCE AT THE TERM (GeometryOptimizationPlan §13.2; lane B 2026-11-09): Vee_Fit<dcmplx> on the plane-wave fit has its
// force face, MissingGradient is empty for a Bloch AND a real TRIM GPW block, and the two arms give the same force (the real child walked
// through the term's tStatic_PC<double> face, its D promoted onto the same complex Coulomb tensor).  The force is the tensors' pullback,
// gated against the energy's finite difference in UTLattice_BS (GPW_Force); here the PHYSICS check is translation: E_H at fixed D does
// not move when every atom moves together, so the per-atom forces sum to the raster's egg-box while each one does not.
TEST(RealComplexTerms, TheGPWHartreeForceWalksBothBlockScalarsAndSumsToTheEggBox)
{
    UnitCell cell(8.0);
    cell.AddAtom(14,{0.40,0.50,0.50});
    cell.AddAtom(14,{0.62,0.55,0.47});
    namespace M = qchem::BasisSet::Gaussian;
    std::shared_ptr<const Real_BS> mol(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
    tGPW_IBS<double> re(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    tGPW_IBS<dcmplx> cx(std::make_shared<::qchem::BasisSet::BasisCache>(), cell, MakeMeshMaker(cell), ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    const size_t n=cx.GetNumFunctions();
    rsmat_t Dr(n); chmat_t Dc(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) { const double d=std::sin(1.0+i+2.0*j)/(1.0+i+j); Dr(i,j)=Dr(j,i)=d; Dc(j,i)=dcmplx(d,0.0); }
    using namespace qchem::ChargeDensity;
    auto cdc=std::make_unique<tComposite_CD<dcmplx>>();
    cdc->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(Dc, &cx, cx.GetIrrep(Spin::None))), cx.GetIrrep(Spin::None));
    auto cdr=std::make_unique<tComposite_CD<dcmplx>>();
    cdr->Insert(std::unique_ptr<tDM_CD<double>>(new PeriodicIrrepCD<double>(Dr, &re, re.GetIrrep(Spin::None))), re.GetIrrep(Spin::None));

    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(cx.CreateCDFitBasisSet(&cell));
    Vee_Fit<dcmplx> vh(Fitting::MakeGSpaceDensityFitter(fb), 0.0, fb);
    const qchem::tests::pw::OneBlockBS<dcmplx> cxBS(&cx), reBS(&re);
    EXPECT_EQ(vh.MissingGradient(&cxBS, &cell), "") << "a Bloch GPW block has the Hartree force";
    EXPECT_EQ(vh.MissingGradient(&reBS, &cell), "") << "a real TRIM GPW block has the Hartree force";

    rvec3vec_t Fc(2, rvec3_t(0,0,0)), Fr(2, rvec3_t(0,0,0));
    RefreshTerm(static_cast<const cDynamic_HT&>(vh), cdc.get(), &cx);
    vh.AddForces(Fc, cdc.get(), &cell);
    RefreshTerm(static_cast<const cDynamic_HT&>(vh), cdr.get(), &re);
    vh.AddForces(Fr, cdr.get(), &cell);
    double big=0.0, d=0.0;
    rvec3_t sum(0,0,0);
    for (size_t A=0;A<2;A++)
    {
        big=std::max({big, std::fabs(Fc[A].x), std::fabs(Fc[A].y), std::fabs(Fc[A].z)});
        d  =std::max({d, std::fabs(Fc[A].x-Fr[A].x), std::fabs(Fc[A].y-Fr[A].y), std::fabs(Fc[A].z-Fr[A].z)});
        sum=sum+Fc[A];
        std::printf("  [GPW Hartree force] atom %zu  (% .10e, % .10e, % .10e)\n", A, Fc[A].x, Fc[A].y, Fc[A].z);
    }
    const double s=std::max({std::fabs(sum.x), std::fabs(sum.y), std::fabs(sum.z)});
    std::printf("  [GPW Hartree force] max|F| %.3e  |F_complex - F_real| %.2e  |Σ F| %.3e\n", big, d, s);
    ASSERT_GT(big, 1e-3) << "a vanishing force cannot show the term is wired";
    EXPECT_LT(d, 1e-12*big) << "the real TRIM block and the Bloch block give different Hartree forces";
    // MEASURED 2026-11-09 (|Σ F| / max|F|): Ecut 10 -> 1.0e-1, 20 -> 1.3e-4, 40 -> 5.3e-5.  It falls with the grid, so it is the DISCRETE
    // energy's egg-box (the force is that energy's exact derivative: GPW_Force, M_CollocationPullback), floored by the ladder's relative
    // level rule, not by this code.  The bound catches a missing or doubled side (which would leave Σ F of order max|F|).
    EXPECT_LT(s, 1e-3*big) << "a rigid translation at fixed D must cost the Hartree energy only the raster's egg-box";
}

// THE k-FOLD OF THE BECKE XC PULLBACK AND VIRIAL (user ruling 2026-11-22 "yes"; the controller's gate: "each overriding term's folded
// pullback equals its per-block sum").  Three GPW blocks of ONE molecular basis on ONE raster ladder -- a real TRIM block at Γ and two Bloch
// blocks at general k -- share a fold domain, so the singles sampler takes them through the first one's folded point-set pullback (the
// functions evaluated once per point for all three).  Against the per-block sum (the base class's default, called non-virtually) the
// forces and the virial agree to roundoff: the same terms in a different order.
TEST(RealComplexTerms, TheBeckeXcPullbackFoldsTheKBlocksOfOneLatticeBasisToThePerBlockSum)
{
    UnitCell cell(8.0);
    cell.AddAtom(14,{0.40,0.50,0.50});
    cell.AddAtom(14,{0.62,0.55,0.47});
    const auto st=cell.Clone();
    namespace M = qchem::BasisSet::Gaussian;
    std::shared_ptr<const Real_BS> mol(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
    const auto ladder=M::MakeRasterLadder(cell, *mol, {.densityEcut=20.0});
    auto cache=std::make_shared<::qchem::BasisSet::BasisCache>();
    tGPW_IBS<double> re (cache, cell, MakeMeshMaker(cell), Symmetry::BlochFactory(ivec3_t(3,3,3), ivec3_t(0,0,0)), mol, ladder);
    tGPW_IBS<dcmplx> cx1(cache, cell, MakeMeshMaker(cell), Symmetry::BlochFactory(ivec3_t(3,3,3), ivec3_t(1,0,0)), mol, ladder);
    tGPW_IBS<dcmplx> cx2(cache, cell, MakeMeshMaker(cell), Symmetry::BlochFactory(ivec3_t(3,3,3), ivec3_t(0,1,1)), mol, ladder);
    ASSERT_NE(re.FoldDomain()[0], 0u);
    ASSERT_TRUE(re.FoldDomain()==cx1.FoldDomain() && cx1.FoldDomain()==cx2.FoldDomain()) << "one basis and one ladder: one fold domain";
    const size_t n=re.GetNumFunctions();
    rsmat_t Dr(n); chmat_t D1(n), D2(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++)
    {
        const double d=std::sin(1.0+i+2.0*j)/(1.0+i+j);
        Dr(i,j)=Dr(j,i)=d;
        D1(j,i)=dcmplx(0.8*d, i==j ? 0.0 : 0.3*std::cos(0.5+i+j)/(1.0+i+j));
        D2(j,i)=dcmplx(1.1*d, i==j ? 0.0 : -0.2*std::sin(2.0+i-j)/(1.0+i+j));
    }
    auto q = MakeMeshMaker(cell)->Make(qcMesh::BeckeXCParams(30, -1.0, 11), qcMesh::NoSymmetry());   // a small atom-centred mesh: owners and a weight response
    auto dfb=std::make_shared<const BasisSet::DeltaFit_IBS<dcmplx>>(std::make_shared<::qchem::BasisSet::BasisCache>(), q, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
    auto engine=ChargeDensity::MakeDensitySampler(dfb, q);
    {
        const auto gap=engine->MakeQuadraturePullback(rvec_t(dfb->GetNumFunctions(), 0.0), st.get());
        ASSERT_TRUE(gap) << "the Becke mesh must give the derivatives: " << gap.Error().why;
    }
    ASSERT_TRUE(engine->CanPullbackMatrix(&cx1)) << "the GPW block must give the point-set pullback";
    rvec_t v(dfb->GetNumFunctions());
    for (size_t g=0; g<v.size(); g++) v[g]=0.5+std::sin(0.37*double(g));   // any field: the identity is linear in it
    ChargeDensity::SampledBlocks<dcmplx> b;
    b.blocks.push_back({&cx1, &v, &D1});
    b.blocks.push_back({&cx2, &v, &D2});
    b.real  .push_back({&re,  &v, &Dr});
    const ChargeDensity::tDensitySampler<dcmplx>& face=*engine;
    const auto Ff=face.MakeMatrixPullbackBlocks(b, st.get());
    const auto Fp=face.ChargeDensity::tDensitySampler<dcmplx>::MakeMatrixPullbackBlocks(b, st.get());   // the per-block default
    ASSERT_TRUE(Ff) << Ff.Error().why;
    ASSERT_TRUE(Fp) << Fp.Error().why;
    double big=0.0, d=0.0;
    for (size_t A=0;A<Ff->size();A++)
    {
        const rvec3_t a=(*Ff)[A], c=(*Fp)[A];
        big=std::max({big, std::fabs(c.x), std::fabs(c.y), std::fabs(c.z)});
        d  =std::max({d, std::fabs(a.x-c.x), std::fabs(a.y-c.y), std::fabs(a.z-c.z)});
    }
    std::printf("  [Becke XC fold] max|F| %.3e  |F_folded - F_per-block| %.2e\n", big, d);
    ASSERT_GT(big, 1e-4) << "a vanishing pullback cannot show the fold is wired";
    EXPECT_LT(d, 1e-12*big) << "the folded Becke XC pullback differs from the per-block sum";
    const auto Vf=face.MakeMatrixVirialBlocks(b, st.get());
    const auto Vp=face.ChargeDensity::tDensitySampler<dcmplx>::MakeMatrixVirialBlocks(b, st.get());
    ASSERT_TRUE(Vf) << Vf.Error().why;
    ASSERT_TRUE(Vp) << Vp.Error().why;
    double vbig=0.0, vd=0.0;
    for (size_t a=1;a<=3;a++) for (size_t c=1;c<=3;c++)
    {
        vbig=std::max(vbig, std::fabs((*Vp)(a,c)));
        vd  =std::max(vd,   std::fabs((*Vf)(a,c)-(*Vp)(a,c)));
    }
    std::printf("  [Becke XC fold] max|V| %.3e  |V_folded - V_per-block| %.2e\n", vbig, vd);
    ASSERT_GT(vbig, 1e-4);
    EXPECT_LT(vd, 1e-12*vbig) << "the folded Becke XC virial differs from the per-block sum";
}

// THE k-FOLD OF THE GPW HARTREE PULLBACK AND VIRIAL (user ruling 2026-11-22 "yes"; the controller's gate: the folded pullback equals the
// per-block sum).  The same three blocks of one basis on one ladder as the Becke gate, in ONE composite density; Vee_Fit<dcmplx> refreshed
// for it, then its PullbackBlocks / VirialBlocks (one jet walk through the first block's Coulomb tensor) against the base class's per-block
// default, called non-virtually, on the same list.
TEST(RealComplexTerms, TheGPWHartreePullbackFoldsTheKBlocksOfOneLatticeBasisToThePerBlockSum)
{
    UnitCell cell(8.0);
    cell.AddAtom(14,{0.40,0.50,0.50});
    cell.AddAtom(14,{0.62,0.55,0.47});
    const auto st=cell.Clone();
    namespace M = qchem::BasisSet::Gaussian;
    std::shared_ptr<const Real_BS> mol(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
    const auto ladder=M::MakeRasterLadder(cell, *mol, {.densityEcut=20.0});
    auto cache=std::make_shared<::qchem::BasisSet::BasisCache>();
    tGPW_IBS<double> re (cache, cell, MakeMeshMaker(cell), Symmetry::BlochFactory(ivec3_t(3,3,3), ivec3_t(0,0,0)), mol, ladder);
    tGPW_IBS<dcmplx> cx1(cache, cell, MakeMeshMaker(cell), Symmetry::BlochFactory(ivec3_t(3,3,3), ivec3_t(1,0,0)), mol, ladder);
    tGPW_IBS<dcmplx> cx2(cache, cell, MakeMeshMaker(cell), Symmetry::BlochFactory(ivec3_t(3,3,3), ivec3_t(0,1,1)), mol, ladder);
    const size_t n=re.GetNumFunctions();
    rsmat_t Dr(n); chmat_t D1(n), D2(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++)
    {
        const double d=std::sin(1.0+i+2.0*j)/(1.0+i+j);
        Dr(i,j)=Dr(j,i)=d;
        D1(j,i)=dcmplx(0.8*d, i==j ? 0.0 : 0.3*std::cos(0.5+i+j)/(1.0+i+j));
        D2(j,i)=dcmplx(1.1*d, i==j ? 0.0 : -0.2*std::sin(2.0+i-j)/(1.0+i+j));
    }
    using namespace qchem::ChargeDensity;
    auto cd=std::make_unique<tComposite_CD<dcmplx>>();
    cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D1, &cx1, cx1.GetIrrep(Spin::None))), cx1.GetIrrep(Spin::None));
    cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(new PeriodicIrrepCD<dcmplx>(D2, &cx2, cx2.GetIrrep(Spin::None))), cx2.GetIrrep(Spin::None));
    cd->Insert(std::unique_ptr<tDM_CD<double>>(new PeriodicIrrepCD<double>(Dr, &re,  re .GetIrrep(Spin::None))), re .GetIrrep(Spin::None));
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(cx1.CreateCDFitBasisSet(&cell));
    Vee_Fit<dcmplx> vh(Fitting::MakeGSpaceDensityFitter(fb), 0.0, fb);
    RefreshTerm(static_cast<const cDynamic_HT&>(vh), cd.get(), &cx1);
    DensityBlocks<dcmplx> b;
    b.blocks.push_back({&cx1, Spin::None, &D1});
    b.blocks.push_back({&cx2, Spin::None, &D2});
    b.real  .push_back({&re,  Spin::None, &Dr});
    const tStatic_PC<dcmplx>& pc=vh;
    const rvec3vec_t Ff=pc.PullbackBlocks(b, st.get());
    const rvec3vec_t Fp=pc.tStatic_PC<dcmplx>::PullbackBlocks(b, st.get());   // the per-block default
    double big=0.0, d=0.0;
    for (size_t A=0;A<Ff.size();A++)
    {
        big=std::max({big, std::fabs(Fp[A].x), std::fabs(Fp[A].y), std::fabs(Fp[A].z)});
        d  =std::max({d, std::fabs(Ff[A].x-Fp[A].x), std::fabs(Ff[A].y-Fp[A].y), std::fabs(Ff[A].z-Fp[A].z)});
    }
    std::printf("  [GPW Hartree fold] max|F| %.3e  |F_folded - F_per-block| %.2e\n", big, d);
    ASSERT_GT(big, 1e-4) << "a vanishing pullback cannot show the fold is wired";
    EXPECT_LT(d, 1e-12*big) << "the folded Hartree pullback differs from the per-block sum";
    const tStatic_VC<dcmplx>& vc=vh;
    const rmat3d_t Vf=vc.VirialBlocks(b, st.get());
    const rmat3d_t Vp=vc.tStatic_VC<dcmplx>::VirialBlocks(b, st.get());
    double vbig=0.0, vd=0.0;
    for (size_t a=1;a<=3;a++) for (size_t c=1;c<=3;c++) { vbig=std::max(vbig, std::fabs(Vp(a,c))); vd=std::max(vd, std::fabs(Vf(a,c)-Vp(a,c))); }
    std::printf("  [GPW Hartree fold] max|V| %.3e  |V_folded - V_per-block| %.2e\n", vbig, vd);
    ASSERT_GT(vbig, 1e-4);
    EXPECT_LT(vd, 1e-12*vbig) << "the folded Hartree virial differs from the per-block sum";
}
