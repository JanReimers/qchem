// File: BasisSet/Gaussian/tests/M_GeometryDerivatives.C  The derivative methods on the integral faces (GeometryOptimizationPlan §13.2),
// one unit gate per derivative method per family that implements it (the user, 2026-11-08: "we will need to make sure these are all
// correct").  FD tolerances are stated against the central difference's own truncation error, h²/6 |f‴|, at each test.
// - PointJets, filled with analytic s Gaussians through its test friend (no family's JetsOn exists yet): each function's jet, ρ's jet and
//   the density pullback against central finite differences, each point riding with its owner and each function with its centre.
// - PG_Cart's overlap / kinetic / nuclear pullbacks, asked of the faces the block already is: libcint's gradient matrices (gated in
//   M_LibCint.C) contracted with the adjoint; translation invariance; a mismatched adjoint throws.
// - PG_Cart's Coulomb-fit (Dunlap) pullbacks, MakeRepulsion3CPullback and its fit's MakeRepulsionPullback, against the central difference
//   of the same faces' contractions at displaced atoms; translation invariance; a mismatched adjoint throws.
// - The SALC decorator's pullbacks (1E, and the ERI pair and whole-system ones) against the central difference of its OWN blocks' matrices,
//   the raw basis rebuilt at displaced atoms under the SAME SALC columns.
// - A family with no derivative yet answers Can...() false and refuses with NoGradient.
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>
#include "forward.H"
import qchem.StructureData;
import qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_HF_IBS;                     // the ERI4 contractions + their pullbacks (the HF pieces)
import qchem.BasisSet.Orbital_DFT_IBS;                    // Repulsion3C + the Coulomb-fit faces (the Dunlap pieces)
import qchem.Mesh;                                        // qcMesh::MeshParams (the fit basis's own mesh, unused by the Coulomb metric)
import qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Point.PG_Spherical;
import qchem.BasisSet.Gaussian.Point.PG_LibCint;               // its own hand Gradient(r) (the hand-gradient audit)
import qchem.BasisSet.Gaussian.Point.BasisFiles;          // BasisFile (the contracted dzvp)
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;
import qchem.BasisSet.Gaussian.PG_Cart.Symmetry;               // ExtractAoShells, StructureToSymPoints (the SALC gate)
import qchem.BasisSet.Internal.SymmetryAdapted_IBS;
import qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet;
import qchem.Symmetry.Molecule.SALC;                            // BuildAbelianGroup, BuildSALCs
import qchem.Structure;
import qchem.Types;

using namespace qchem;
using BasisSet::Gaussian::PG_Cart::Orbital_IBS;

namespace
{
// THREE s Gaussians, one per atom of a bent triatomic: φ_i(r) = exp(-α_i |r - R_i|²), ∇_r φ_i = -2 α_i (r - R_i) φ_i.
const std::vector<rvec3_t> kR={rvec3_t(0,0,0.2), rvec3_t(0,1.4,-0.9), rvec3_t(0.3,-1.3,-0.8)};
const double kAlpha[3]={0.9, 0.5, 0.7};
double Phi(const std::vector<rvec3_t>& R, size_t i, const rvec3_t& r) {const rvec3_t d=r-R[i]; return std::exp(-kAlpha[i]*(d*d));}
const std::vector<rvec3_t> kOffsets={rvec3_t(0.3,-0.2,0.1), rvec3_t(-0.4,0.5,0.2), rvec3_t(0.1,0.1,-0.6)};
rmat_t SomeD() { rmat_t D(3,3); const double v[9]={0.9,0.2,-0.1, 0.2,0.6,0.05, -0.1,0.05,0.4}; for (size_t k=0;k<9;k++) D(k/3,k%3)=v[k]; return D; }
double Rho(const std::vector<rvec3_t>& R, const rvec3_t& r, const rmat_t& D)
{ double s=0.0; for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) s+=Phi(R,i,r)*D(i,j)*Phi(R,j,r); return s; }
double& C(rvec3_t& v, int k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
double  C(const rvec3_t& v, int k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
}

//! The friend (forward.H): fills a PointJets with the analytic functions above, at the points of one owner's slab.
class PointJetsTests
{
public:
    static BasisSet::PointJets<double> Analytic(const std::vector<rvec3_t>& R, size_t owner, const rvec3vec_t& pts)
    {
        BasisSet::PointJets<double> J;
        auto t=J.Fill(owner, pts.size(), {0,1,2});
        for (size_t q=0;q<pts.size();q++) for (size_t i=0;i<3;i++)
        {
            const double f=Phi(R,i,pts[q]); const rvec3_t g=(-2.0*kAlpha[i]*f)*(pts[q]-R[i]);
            t.phi(q,i)=f; t.grad[0](q,i)=g.x; t.grad[1](q,i)=g.y; t.grad[2](q,i)=g.z;
        }
        return J;
    }
};

TEST(PointJets, AFunctionsAndADensitysJetsAreTheFiniteDifferenceWithThePointRidingWithItsOwner)
{
    const rmat_t D=SomeD();
    const double h=1e-6;
    for (size_t owner=0;owner<3;owner++)
    {
        rvec3vec_t pts(kOffsets.size()); for (size_t q=0;q<pts.size();q++) pts[q]=kR[owner]+kOffsets[q];
        const auto J=PointJetsTests::Analytic(kR, owner, pts);
        ASSERT_EQ(J.NumPoints(), pts.size());
        for (size_t q=0;q<pts.size();q++)
        {
            const auto rho=J.Density(q, D);
            EXPECT_NEAR(rho.Value(), Rho(kR, pts[q], D), 1e-15);
            for (size_t A=0;A<3;A++) for (int k=0;k<3;k++)
            {
                auto p=kR, m=kR; C(p[A],k)+=h; C(m[A],k)-=h;
                const rvec3_t rp=p[owner]+kOffsets[q], rm=m[owner]+kOffsets[q];   // the point rides with its owner
                EXPECT_NEAR(C(rho.D(A),k), (Rho(p,rp,D)-Rho(m,rm,D))/(2*h), 1e-8) << "dρ/dR owner " << owner << " point " << q << " atom " << A << " axis " << k;
                for (size_t i=0;i<3;i++)
                    EXPECT_NEAR(C(J.Jet(q,i).D(A),k), (Phi(p,i,rp)-Phi(m,i,rm))/(2*h), 1e-8) << "dφ_" << i << "/dR atom " << A << " axis " << k;
            }
            EXPECT_LT(norm(J.Jet(q,owner).D(owner)), 1e-15) << "a function riding with its own point does not change";
        }
    }
}

TEST(PointJets, TheDensityPullbackIsTheWeightedSumOfTheDensityJets)
{
    const rmat_t D=SomeD();
    rvec3vec_t pts(kOffsets.size()); for (size_t q=0;q<pts.size();q++) pts[q]=kR[1]+kOffsets[q];
    const auto J=PointJetsTests::Analytic(kR, 1, pts);
    rvec_t a(pts.size()); for (size_t q=0;q<a.size();q++) a[q]=0.3+0.2*q;
    rvec3vec_t F(3, rvec3_t(0,0,0)), G(3, rvec3_t(0,0,0));
    J.AddDensityPullback(D, a, F);
    for (size_t q=0;q<pts.size();q++) { const auto rho=J.Density(q,D); for (size_t A=0;A<3;A++) G[A]+=a[q]*rho.D(A); }
    for (size_t A=0;A<3;A++) EXPECT_LT(norm(F[A]-G[A]), 1e-14) << "atom " << A;
    rvec3_t sum(0,0,0); for (size_t A=0;A<F.size();A++) sum+=F[A];   // an index loop: Blaze iterators lack op!= here (CLAUDE.md)
    EXPECT_LT(norm(sum), 1e-14) << "translation invariance";
    rvec3vec_t small(2);
    EXPECT_THROW(J.AddDensityPullback(D, a, small), std::out_of_range);
    EXPECT_THROW(J.AddDensityPullback(D, rvec_t(2), F), std::invalid_argument);
}

TEST(GeometryDerivatives, PG_CartsPullbacksAreLibcintsGradientsContractedAndTranslationInvariant)
{
    std::unique_ptr<Molecule> h2o(new Molecule(StructureData::GetMolecule("H2O")));
    const Orbital_IBS ibs(rvec_t{1.0, 0.25}, 2, h2o.get());
    const BasisSet::Orbital_1E_IBS<double>& block=ibs;   // a term holds this face; it asks it, no cast
    ASSERT_TRUE(block.CanPullbackOverlap() && block.CanPullbackKinetic() && block.CanPullbackNuclear());
    const size_t n=block.GetNumFunctions();
    std::mt19937 rng(11); std::uniform_real_distribution<double> U(-1,1);
    hmat_t<double> A(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) A(i,j)=U(rng);
    const BasisSet::Gaussian::Evaluators::PG_LibCint::NR_Evaluator lc(dynamic_cast<const BasisSet::Gaussian::Evaluators::PG_Cart_MnD::PGData&>(ibs), h2o.get());
    const auto dS=lc.OverlapGradient(), dK=lc.KineticGradient(), dV=lc.NuclearGradient();
    const auto FS=block.MakeOverlapPullback(A,h2o.get()), FK=block.MakeKineticPullback(A,h2o.get()), FV=block.MakeNuclearPullback(A,h2o.get());
    ASSERT_TRUE(FS && FK && FV);
    ASSERT_EQ(FS->size(), h2o->GetNumAtoms());
    rvec3_t sS(0,0,0), sK(0,0,0), sV(0,0,0);
    for (size_t at=0;at<FS->size();at++)
    {
        for (int k=0;k<3;k++)
        {
            double s=0, kk=0, v=0;
            for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) { s+=A(i,j)*dS[at][k](i,j); kk+=A(i,j)*dK[at][k](i,j); v+=A(i,j)*dV[at][k](i,j); }
            EXPECT_DOUBLE_EQ(C((*FS)[at],k), s); EXPECT_DOUBLE_EQ(C((*FK)[at],k), kk); EXPECT_DOUBLE_EQ(C((*FV)[at],k), v);
        }
        sS+=(*FS)[at]; sK+=(*FK)[at]; sV+=(*FV)[at];
    }
    EXPECT_LT(norm(sS), 1e-11); EXPECT_LT(norm(sK), 1e-10); EXPECT_LT(norm(sV), 1e-10);
    EXPECT_THROW((void)block.MakeOverlapPullback(hmat_t<double>(n+1), h2o.get()), std::invalid_argument) << "a mismatched adjoint throws";
}

namespace
{
std::unique_ptr<Molecule> Water() {return std::make_unique<Molecule>(StructureData::GetMolecule("H2O"));}
std::vector<rvec3_t> Sites(const Structure& st) {std::vector<rvec3_t> R; st.ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); }); return R;}
//! The contracted DZVP block of \a st (s, p, d shells, contracted Cartesian Gaussians).
std::unique_ptr<Orbital_IBS> Dzvp(const Structure* st)
{
    BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile("dzvp.bsd"));
    return std::make_unique<Orbital_IBS>(&reader, st);
}
double TrAM(const hmat_t<double>& A, const hmat_t<double>& M)
{ double t=0.0; for (size_t i=0;i<A.rows();i++) for (size_t j=0;j<A.rows();j++) t+=A(i,j)*M(i,j); return t; }
}

// GATE 1, PG_Cart's pullbacks: Σ A ∂M/∂R_A against the central difference of Tr(A M) of the SAME matrices (our M&D engine) at displaced
// atoms, on the contracted DZVP block.  Step h = 1e-4: the truncation error is h²/6 |∂³Tr(AM)/∂R³| ~ 2e-9 |∂³|, and with a random ±1 A
// over ~25² entries |∂³| stays below ~1e2, so 2e-6 (relative above 1) is ~10× the estimate.  Translation invariance: S and K (only the
// functions move), and V with the nuclei moving too.
TEST(GeometryDerivatives, PG_CartsPullbacksAreTheFiniteDifferenceOfItsOwnMatricesAndTranslationInvariant)
{
    const auto h2o=Water();
    const auto ibs=Dzvp(h2o.get());
    const BasisSet::Orbital_1E_IBS<double>& block=*ibs;
    const size_t n=block.GetNumFunctions();
    std::mt19937 rng(7); std::uniform_real_distribution<double> U(-1,1);
    hmat_t<double> A(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) A(i,j)=U(rng);
    const auto FS=block.MakeOverlapPullback(A,h2o.get()), FK=block.MakeKineticPullback(A,h2o.get()), FV=block.MakeNuclearPullback(A,h2o.get());
    ASSERT_TRUE(FS && FK && FV);
    const double h=1e-4;
    const auto R=Sites(*h2o);
    rvec3_t sS(0,0,0), sK(0,0,0), sV(0,0,0);
    for (size_t at=0;at<R.size();at++)
    {
        sS+=(*FS)[at]; sK+=(*FK)[at]; sV+=(*FV)[at];
        for (int k=0;k<3;k++)
        {
            auto at_=[&](double s){ auto Rs=R; C(Rs[at],k)+=s; std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(Rs);
                                     auto b=Dzvp(st.get()); const BasisSet::Orbital_1E_IBS<double>& b1=*b;
                                     return std::array<double,3>{TrAM(A,b1.Overlap()), TrAM(A,b1.Kinetic()), TrAM(A,b1.Nuclear(st.get()))}; };
            const auto p=at_(h), m=at_(-h);
            const double fS=(p[0]-m[0])/(2*h), fK=(p[1]-m[1])/(2*h), fV=(p[2]-m[2])/(2*h);
            EXPECT_NEAR(C((*FS)[at],k), fS, 2e-6*std::max(1.0,std::fabs(fS))) << "overlap atom " << at << " axis " << k;
            EXPECT_NEAR(C((*FK)[at],k), fK, 2e-6*std::max(1.0,std::fabs(fK))) << "kinetic atom " << at << " axis " << k;
            EXPECT_NEAR(C((*FV)[at],k), fV, 2e-6*std::max(1.0,std::fabs(fV))) << "nuclear atom " << at << " axis " << k;
        }
    }
    EXPECT_LT(norm(sS), 1e-10); EXPECT_LT(norm(sK), 1e-9); EXPECT_LT(norm(sV), 1e-9);
}

// GATE 2, PG_Cart's hand Gradient(r), load-bearing for the jets: against the central difference of operator()(r), every component of
// every shell, for s through g (uncontracted, two exponents) and the contracted DZVP, at points near (0.3 bohr) and far (3 bohr) from the
// centres.  Step h = 1e-5: the error h²/6 |φ‴| is about 2e-11 |φ‴|, and |φ‴| ≲ (2α)^{3/2} for the exponents probed (≲ 1e3 for the DZVP
// cores at 0.3 bohr), so 1e-7 absolute (relative above 1) is ≥ 5× the estimate.  The magnitude screen (|φ| < 1e-10 → 0) never fires here.
TEST(GeometryDerivatives, PG_CartsHandGradientIsTheFiniteDifferenceOfItsValuesForEveryShellType)
{
    const auto h2o=Water();
    const auto R=Sites(*h2o);
    std::vector<std::unique_ptr<Orbital_IBS>> blocks;
    blocks.push_back(std::make_unique<Orbital_IBS>(rvec_t{1.3, 0.35}, 4, h2o.get()));   // s, p, d, f, g
    blocks.push_back(Dzvp(h2o.get()));                                                 // contracted s, p, d
    const double h=1e-5;
    for (const auto& b : blocks)
        for (const rvec3_t& base : R)
            for (const rvec3_t& off : {rvec3_t(0.17,-0.21,0.12), rvec3_t(1.9,1.2,-1.8)})
            {
                const rvec3_t r=base+off;
                const rvec3vec_t g=b->Gradient(r);
                for (int k=0;k<3;k++)
                {
                    rvec3_t p=r, m=r; C(p,k)+=h; C(m,k)-=h;
                    const rvec_t vp=(*b)(p), vm=(*b)(m);
                    for (size_t i=0;i<b->GetNumFunctions();i++)
                    {
                        const double fd=(vp[i]-vm[i])/(2*h);
                        EXPECT_NEAR(C(g[i],k), fd, 1e-7*std::max(1.0,std::fabs(fd))) << "function " << i << " axis " << k << " at " << r;
                    }
                }
            }
}

// GATE 3, PG_Cart's jets: JetsOn's φ jets, ρ's jet and the density pullback against central differences in the ATOM positions, the block
// rebuilt at each displaced geometry and each point riding with its owner.  Step h = 1e-5, the error as GATE 2: 1e-7 (relative above 1).
TEST(GeometryDerivatives, PG_CartsJetsAreTheFiniteDifferenceInTheAtomPositionsWithThePointRidingItsOwner)
{
    const auto h2o=Water();
    const auto R=Sites(*h2o);
    const auto ibs=Dzvp(h2o.get());
    const BasisSet::Evaluatable_IBS<double>& eval=*ibs;
    ASSERT_TRUE(eval.CanGiveJets());
    const size_t n=eval.GetNumFunctions();
    rmat_t D(n,n); std::mt19937 rng(3); std::uniform_real_distribution<double> U(-0.5,0.5);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=D(j,i)=U(rng);
    const std::vector<rvec3_t> offsets={rvec3_t(0.2,-0.3,0.1), rvec3_t(-0.6,0.5,0.4), rvec3_t(1.1,0.2,-0.7)};
    const double h=1e-5;
    for (size_t owner=0;owner<R.size();owner++)
    {
        rvec3vec_t pts(offsets.size()); for (size_t q=0;q<pts.size();q++) pts[q]=R[owner]+offsets[q];
        BasisSet::PointJets<double> J;
        const auto filled=eval.JetsOn(pts, owner, h2o.get(), J);
        ASSERT_TRUE(filled); ASSERT_EQ(*filled, pts.size());
        rvec_t a(pts.size()); for (size_t q=0;q<a.size();q++) a[q]=0.4+0.3*q;
        rvec3vec_t F(R.size(), rvec3_t(0,0,0));
        J.AddDensityPullback(D, a, F);
        for (size_t A=0;A<R.size();A++) for (int k=0;k<3;k++)
        {
            auto rebuilt=[&](double s){ auto Rs=R; C(Rs[A],k)+=s; std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(Rs);
                                         auto b=Dzvp(st.get()); std::vector<rvec_t> v;
                                         for (const auto& off : offsets) v.push_back((*b)(Rs[owner]+off));   // the point rides with its owner
                                         return v; };
            const auto vp=rebuilt(h), vm=rebuilt(-h);
            double sumFD=0.0;
            for (size_t q=0;q<pts.size();q++)
            {
                double rp=0.0, rm=0.0;
                for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) { rp+=vp[q][i]*D(i,j)*vp[q][j]; rm+=vm[q][i]*D(i,j)*vm[q][j]; }
                const double fdRho=(rp-rm)/(2*h);
                EXPECT_NEAR(C(J.Density(q,D).D(A),k), fdRho, 1e-7*std::max(1.0,std::fabs(fdRho))) << "dρ/dR owner " << owner << " point " << q << " atom " << A << " axis " << k;
                sumFD+=a[q]*fdRho;
                for (size_t i=0;i<n;i++)
                {
                    const double fd=(vp[q][i]-vm[q][i])/(2*h);
                    EXPECT_NEAR(C(J.Jet(q,i).D(A),k), fd, 1e-7*std::max(1.0,std::fabs(fd))) << "dφ_" << i << "/dR owner " << owner << " atom " << A << " axis " << k;
                }
            }
            EXPECT_NEAR(C(F[A],k), sumFD, 1e-7*std::max(1.0,std::fabs(sumFD))) << "pullback owner " << owner << " atom " << A << " axis " << k;
        }
    }
}

// GATE 4, the refusals: a family with no derivative yet answers Can...() false AND its methods refuse with a reason (the Outcome); a family
// that gives them answers Can...() true AND its methods answer.  PG_Spherical refuses all four in increment 1; PG_Cart gives all four.
TEST(GeometryDerivatives, EachCanQueryAgreesWithWhetherItsMethodAnswers)
{
    const auto h2o=Water();
    BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile("dzvp.bsd"));
    const BasisSet::Gaussian::PG_Spherical::Orbital_IBS sph(&reader, h2o.get());
    const auto cart=Dzvp(h2o.get());
    for (const BasisSet::Orbital_1E_IBS<double>* b : {static_cast<const BasisSet::Orbital_1E_IBS<double>*>(&sph), static_cast<const BasisSet::Orbital_1E_IBS<double>*>(cart.get())})
    {
        const hmat_t<double> A(b->GetNumFunctions());
        BasisSet::PointJets<double> J;
        const auto s=b->MakeOverlapPullback(A,h2o.get()), k=b->MakeKineticPullback(A,h2o.get()), v=b->MakeNuclearPullback(A,h2o.get());
        const auto j=b->JetsOn(rvec3vec_t(1, rvec3_t(0.1,0.2,0.3)), 0, h2o.get(), J);
        EXPECT_EQ(b->CanPullbackOverlap(), bool(s)); EXPECT_EQ(b->CanPullbackKinetic(), bool(k));
        EXPECT_EQ(b->CanPullbackNuclear(), bool(v)); EXPECT_EQ(b->CanGiveJets(), bool(j));
        if (!s) EXPECT_FALSE(s.Error().why.empty());
        if (!j) EXPECT_FALSE(j.Error().why.empty());
    }
    EXPECT_FALSE(sph.CanPullbackOverlap()) << "PG_Spherical has no geometry derivative in increment 1";
    EXPECT_TRUE(cart->CanPullbackOverlap() && cart->CanGiveJets()) << "PG_Cart gives them";
}

// GATE 5, the Dunlap pieces (§13.2): the DZVP block's MakeRepulsion3CPullback(fit, D, a) = Σ D_ij a_k ∂(ij|k)/∂R and its Coulomb fit's
// MakeRepulsionPullback(A) = Σ A_kl ∂(k|l)/∂R, against the central difference of the SAME faces' contractions (our M&D engine:
// Repulsion3C and the fit's Repulsion) at displaced atoms, the block, the fit and their centres all rebuilt.  Step h = 1e-4: the
// truncation error h²/6 |f‴| is ~2e-9 |f‴|, and with random ±1 adjoints (D over ~25², a over the fit, A over its square) |f‴| stays below
// ~1e2 relative to |f′|, so 2e-6 (relative above 1) is ~10× the estimate.  Translation invariance: every centre is an atom's.
TEST(GeometryDerivatives, PG_CartsCoulombFitPullbacksAreTheFiniteDifferenceOfItsOwnIntegralsAndTranslationInvariant)
{
    const auto h2o=Water();
    const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3};
    auto build=[&](const Structure* st)
    {
        auto b=Dzvp(st);
        std::unique_ptr<BasisSet::rFIT_CD_ABS> f(b->CreateCDFitBasisSet(st, mp));
        return std::pair{std::move(b), std::move(f)};
    };
    auto [ibs, fit]=build(h2o.get());
    const BasisSet::Orbital_DFT_IBS<double>& dft=*ibs;
    const auto& metric=dynamic_cast<const BasisSet::FIT_CD_NonOrtho&>(*fit);
    const size_t n=dft.GetNumFunctions(), nf=fit->GetNumFunctions();
    std::mt19937 rng(11); std::uniform_real_distribution<double> U(-1,1);
    hmat_t<double> D(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=U(rng);
    rvec_t a(nf); for (size_t k=0;k<nf;k++) a[k]=U(rng);
    rsmat_t A(nf); for (size_t k=0;k<nf;k++) for (size_t l=k;l<nf;l++) A(k,l)=U(rng);
    ASSERT_TRUE(dft.CanPullbackRepulsion3C(*fit)); ASSERT_TRUE(metric.CanPullbackRepulsion());
    const auto F3=dft.MakeRepulsion3CPullback(*fit, D, a, h2o.get());
    const auto F2=metric.MakeRepulsionPullback(A, h2o.get());
    ASSERT_TRUE(F3 && F2);
    auto contractions=[&](const Structure* st)
    {
        auto [b, f]=build(st);
        const BasisSet::Orbital_DFT_IBS<double>& d=*b;
        const auto& P=d.Repulsion3C(*f);
        double s3=0.0; for (size_t k=0;k<nf;k++) s3+=a[k]*TrAM(D, P.dense[k]);
        const auto& G=dynamic_cast<const BasisSet::FIT_CD_NonOrtho&>(*f).Repulsion();
        double s2=0.0; for (size_t k=0;k<nf;k++) for (size_t l=0;l<nf;l++) s2+=A(k,l)*G(k,l);
        return std::array<double,2>{s3, s2};
    };
    const double h=1e-4;
    const auto R=Sites(*h2o);
    rvec3_t s3(0,0,0), s2(0,0,0);
    double largest=0.0;
    for (size_t at=0;at<R.size();at++)
    {
        s3+=(*F3)[at]; s2+=(*F2)[at];
        for (int k=0;k<3;k++)
        {
            auto at_=[&](double s){ auto Rs=R; C(Rs[at],k)+=s; std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(Rs); return contractions(st.get()); };
            const auto p=at_(h), m=at_(-h);
            const double f3=(p[0]-m[0])/(2*h), f2=(p[1]-m[1])/(2*h);
            EXPECT_NEAR(C((*F3)[at],k), f3, 2e-6*std::max(1.0,std::fabs(f3))) << "Σ D a ∂(ij|k)/∂R, atom " << at << " axis " << k;
            EXPECT_NEAR(C((*F2)[at],k), f2, 2e-6*std::max(1.0,std::fabs(f2))) << "Σ A ∂(k|l)/∂R, atom " << at << " axis " << k;
            largest=std::max({largest, std::fabs(f3), std::fabs(f2)});
        }
    }
    EXPECT_GT(largest, 0.1) << "the contractions move with the atoms";
    EXPECT_LT(norm(s3), 1e-9); EXPECT_LT(norm(s2), 1e-9);
    EXPECT_THROW((void)dft.MakeRepulsion3CPullback(*fit, hmat_t<double>(n+1), a, h2o.get()), std::invalid_argument) << "a mismatched D";
    EXPECT_THROW((void)metric.MakeRepulsionPullback(rsmat_t(nf+1), h2o.get()), std::invalid_argument) << "a mismatched adjoint";
}

// The Coulomb-fit Can queries agree with their methods: PG_Cart with its own fit answers both; PG_Spherical, and its own fit, refuse both
// with a reason.
TEST(GeometryDerivatives, EachCoulombFitCanQueryAgreesWithWhetherItsMethodAnswers)
{
    const auto h2o=Water();
    const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3};
    BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile("dzvp.bsd"));
    const BasisSet::Gaussian::PG_Spherical::Orbital_IBS sph(&reader, h2o.get());
    const auto cart=Dzvp(h2o.get());
    for (const BasisSet::Orbital_DFT_IBS<double>* b : {static_cast<const BasisSet::Orbital_DFT_IBS<double>*>(&sph), static_cast<const BasisSet::Orbital_DFT_IBS<double>*>(cart.get())})
    {
        std::unique_ptr<BasisSet::rFIT_CD_ABS> f(b->CreateCDFitBasisSet(h2o.get(), mp));
        const auto& metric=dynamic_cast<const BasisSet::FIT_CD_NonOrtho&>(*f);
        const auto r3=b->MakeRepulsion3CPullback(*f, hmat_t<double>(b->GetNumFunctions()), rvec_t(f->GetNumFunctions(), 0.0), h2o.get());
        const auto r2=metric.MakeRepulsionPullback(rsmat_t(f->GetNumFunctions()), h2o.get());
        EXPECT_EQ(b->CanPullbackRepulsion3C(*f), bool(r3)); EXPECT_EQ(metric.CanPullbackRepulsion(), bool(r2));
        if (!r3) EXPECT_FALSE(r3.Error().why.empty());
        if (!r2) EXPECT_FALSE(r2.Error().why.empty());
    }
    std::unique_ptr<BasisSet::rFIT_CD_ABS> fc(cart->CreateCDFitBasisSet(h2o.get(), mp)), fs(sph.CreateCDFitBasisSet(h2o.get(), mp));
    EXPECT_TRUE(cart->CanPullbackRepulsion3C(*fc)) << "PG_Cart with its own fit gives it";
    EXPECT_FALSE(sph.CanPullbackRepulsion3C(*fs)) << "PG_Spherical has no geometry derivative yet";
    EXPECT_FALSE(cart->CanPullbackRepulsion3C(*fs)) << "a fit of another family: refused, not cast";
}

// GATE 6, the HF pieces (§13.2): PG_Cart's MakeDirectPullback / MakeExchangePullback(Da, Db, other) against the central difference of the
// SAME block's own contractions, Σ Da (J or K)[Db] through AccumulateDirect / AccumulateExchange (our M&D ERI4), at displaced atoms, all
// four centres rebuilt.  Two pairings: the DZVP block with itself, and the DZVP block against an even-tempered s+p block (the cross
// pair: two different shell sets).  Step h = 1e-4: the truncation error h²/6 |f‴| is ~2e-9 |f‴|, and with random ±1 densities |f‴| stays
// below ~1e2 relative to |f′|, so 2e-6 (relative above 1) is ~10× the estimate.  Translation invariance; the Can query agrees.
TEST(GeometryDerivatives, PG_CartsERI4PullbacksAreTheFiniteDifferenceOfItsOwnContractionsAndTranslationInvariant)
{
    const auto h2o=Water();
    const auto R=Sites(*h2o);
    auto even=[](const Structure* st){ return std::make_unique<Orbital_IBS>(rvec_t{1.0,0.25}, 1, st); };
    for (bool cross : {false, true})
    {
        const auto a=Dzvp(h2o.get());
        const auto bOwn=even(h2o.get());
        const BasisSet::Orbital_HF_IBS<double>& A=*a;
        const BasisSet::Orbital_HF_IBS<double>& B=cross ? static_cast<const BasisSet::Orbital_HF_IBS<double>&>(*bOwn) : A;
        const size_t na=a->GetNumFunctions(), nb=cross ? bOwn->GetNumFunctions() : na;
        std::mt19937 rng(cross ? 17 : 13); std::uniform_real_distribution<double> U(-1,1);
        rsmat_t Da(na), Db(nb);
        for (size_t i=0;i<na;i++) for (size_t j=i;j<na;j++) Da(i,j)=U(rng);
        for (size_t i=0;i<nb;i++) for (size_t j=i;j<nb;j++) Db(i,j)=U(rng);
        ASSERT_TRUE(A.CanPullbackERI4(&B));
        const auto FJ=A.MakeDirectPullback(Da, Db, &B, h2o.get()), FK=A.MakeExchangePullback(Da, Db, &B, h2o.get());
        ASSERT_TRUE(FJ && FK);
        auto contractions=[&](const Structure* st)
        {
            const auto a1=Dzvp(st); const auto b1=even(st);
            const BasisSet::Orbital_HF_IBS<double>& A1=*a1;
            const BasisSet::Orbital_HF_IBS<double>& B1=cross ? static_cast<const BasisSet::Orbital_HF_IBS<double>&>(*b1) : A1;
            rsmat_t J(na), K(na);
            A1.AccumulateDirect(J, Db, &B1); A1.AccumulateExchange(K, Db, &B1);
            return std::array<double,2>{TrAM(Da, J), TrAM(Da, K)};
        };
        const double h=1e-4;
        rvec3_t sJ(0,0,0), sK(0,0,0); double largest=0.0;
        for (size_t at=0;at<R.size();at++)
        {
            sJ+=(*FJ)[at]; sK+=(*FK)[at];
            for (int k=0;k<3;k++)
            {
                auto at_=[&](double s){ auto Rs=R; C(Rs[at],k)+=s; std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(Rs); return contractions(st.get()); };
                const auto p=at_(h), m=at_(-h);
                const double fJ=(p[0]-m[0])/(2*h), fK=(p[1]-m[1])/(2*h);
                EXPECT_NEAR(C((*FJ)[at],k), fJ, 2e-6*std::max(1.0,std::fabs(fJ))) << (cross ? "cross " : "self ") << "Coulomb, atom " << at << " axis " << k;
                EXPECT_NEAR(C((*FK)[at],k), fK, 2e-6*std::max(1.0,std::fabs(fK))) << (cross ? "cross " : "self ") << "exchange, atom " << at << " axis " << k;
                largest=std::max({largest, std::fabs(fJ), std::fabs(fK)});
            }
        }
        EXPECT_GT(largest, 0.1) << "the contractions move with the atoms";
        EXPECT_LT(norm(sJ), 1e-8); EXPECT_LT(norm(sK), 1e-8);
    }
    BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile("dzvp.bsd"));
    const BasisSet::Gaussian::PG_Spherical::Orbital_IBS sph(&reader, h2o.get());
    const BasisSet::Orbital_HF_IBS<double>& S=sph;
    EXPECT_FALSE(S.CanPullbackERI4(&S)) << "PG_Spherical has no geometry derivative yet";
    const auto r=S.MakeDirectPullback(rsmat_t(sph.GetNumFunctions()), rsmat_t(sph.GetNumFunctions()), &S, h2o.get());
    EXPECT_FALSE(bool(r)); if (!r) EXPECT_FALSE(r.Error().why.empty());
}

// GATE 6, the SALC decorator (§13.2; lane B 2026-11-09): every C2v block of water's DZVP gets a random symmetric density D_Γ, and
//   Σ_Γ MakeXPullback(D_Γ)                         (X = overlap, kinetic, nuclear),
//   Σ_ΓΔ MakeDirect/ExchangePullback(D_Γ, D_Δ, Δ)  (every ordered pair), and the whole-system MakeAOPullback(Σ_Γ O_Γ D_Γ O_Γᵀ),
// against the central difference of Σ_Γ Tr(D_Γ M_Γ) and Σ_ΓΔ Tr(D_Γ X_Γ[D_Δ]) of the SAME blocks' own matrices and contractions, the raw
// basis rebuilt at displaced atoms under the SAME SALC columns (O is held: it is made of symmetry-operation images, not of R).  Step and
// tolerance as gates 1 and 5; the whole-system pullback is the pair sum to roundoff (the pair pullback is bilinear).
TEST(GeometryDerivatives, SALCPullbacksAreTheFiniteDifferenceOfTheirOwnBlocksWithTheColumnsHeld)
{
    using SALCBlock=BasisSet::SymmetryAdapted_IBS;
    const auto h2o=Water();
    const auto R=Sites(*h2o);
    const auto raw=Dzvp(h2o.get());
    const auto pts=BasisSet::Gaussian::PG_Cart::StructureToSymPoints(*h2o);
    const auto grp=Symmetry::Molecule::BuildAbelianGroup(pts, 1e-4);
    const auto salc=Symmetry::Molecule::BuildSALCs(BasisSet::Gaussian::PG_Cart::ExtractAoShells(*raw), grp, Symmetry::Molecule::Centroid(pts), 1e-4);
    auto blocks=[&](const BasisSet::Gaussian::SymmetryAdaptedBasisSet& sab)
    {
        std::vector<const SALCBlock*> b;
        for (auto oi : sab.Iterate<BasisSet::Real_OIBS>()) b.push_back(dynamic_cast<const SALCBlock*>(oi));   // a test may cheat
        return b;
    };
    const BasisSet::Gaussian::SymmetryAdaptedBasisSet sab(raw.get(), salc);
    const auto B=blocks(sab);
    ASSERT_GE(B.size(), 2u);
    std::mt19937 rng(23); std::uniform_real_distribution<double> U(-1,1);
    std::vector<rsmat_t> D;
    rsmat_t Dao(raw->GetNumFunctions());
    for (const SALCBlock* b : B)
    {
        ASSERT_NE(b, nullptr);
        ASSERT_TRUE(b->CanPullbackOverlap() && b->CanPullbackKinetic() && b->CanPullbackNuclear());
        const size_t n=b->GetNumFunctions();
        rsmat_t X(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) X(i,j)=U(rng);
        b->AddAODensity(Dao, X);
        D.push_back(X);
    }
    const size_t nat=R.size();
    rvec3vec_t FS(nat, rvec3_t(0,0,0)), FK=FS, FV=FS, FJ=FS, FX=FS;
    auto add=[](rvec3vec_t& F, const Outcome<rvec3vec_t,AD::NoGradient>& g){ ASSERT_TRUE(bool(g)) << (g ? "" : g.Error().why); for (size_t A=0;A<F.size();A++) F[A]+=(*g)[A]; };
    for (size_t k=0;k<B.size();k++)
    {
        add(FS, B[k]->MakeOverlapPullback(D[k], h2o.get())); add(FK, B[k]->MakeKineticPullback(D[k], h2o.get()));
        add(FV, B[k]->MakeNuclearPullback(D[k], h2o.get()));
        for (size_t l=0;l<B.size();l++)
        {
            ASSERT_TRUE(B[k]->CanPullbackERI4(B[l]));
            add(FJ, B[k]->MakeDirectPullback(D[k], D[l], B[l], h2o.get())); add(FX, B[k]->MakeExchangePullback(D[k], D[l], B[l], h2o.get()));
        }
    }
    const auto WJ=B[0]->MakeAOPullback(Dao, false, h2o.get()), WK=B[0]->MakeAOPullback(Dao, true, h2o.get());
    ASSERT_TRUE(WJ && WK);
    for (size_t A=0;A<nat;A++) { EXPECT_LT(norm((*WJ)[A]-FJ[A]), 1e-10) << "whole-system Coulomb, atom " << A;
                                 EXPECT_LT(norm((*WK)[A]-FX[A]), 1e-10) << "whole-system exchange, atom " << A; }
    auto contractions=[&](const Structure* st)
    {
        const auto r1=Dzvp(st);
        const BasisSet::Gaussian::SymmetryAdaptedBasisSet s1(r1.get(), salc);   // the SAME columns
        const auto B1=blocks(s1);
        std::array<double,5> t{0,0,0,0,0};
        for (size_t k=0;k<B1.size();k++)
        {
            t[0]+=TrAM(D[k], B1[k]->Overlap()); t[1]+=TrAM(D[k], B1[k]->Kinetic()); t[2]+=TrAM(D[k], B1[k]->Nuclear(st));
            for (size_t l=0;l<B1.size();l++)
            {
                rsmat_t J(D[k].rows()), K(D[k].rows());
                B1[k]->AccumulateDirect(J, D[l], B1[l]); B1[k]->AccumulateExchange(K, D[l], B1[l]);
                t[3]+=TrAM(D[k], J); t[4]+=TrAM(D[k], K);
            }
        }
        return t;
    };
    const std::array<const rvec3vec_t*,5> F{&FS, &FK, &FV, &FJ, &FX};
    const char* name[5]={"overlap", "kinetic", "nuclear", "Coulomb", "exchange"};
    const double h=1e-4;
    std::array<rvec3_t,5> sum; sum.fill(rvec3_t(0,0,0));
    double largest=0.0;
    for (size_t at=0;at<nat;at++)
    {
        for (size_t q=0;q<5;q++) sum[q]+=(*F[q])[at];
        for (int k=0;k<3;k++)
        {
            auto at_=[&](double s){ auto Rs=R; C(Rs[at],k)+=s; std::shared_ptr<Structure> st=h2o->Clone(); st->MoveSites(Rs); return contractions(st.get()); };
            const auto p=at_(h), m=at_(-h);
            for (size_t q=0;q<5;q++)
            {
                const double fd=(p[q]-m[q])/(2*h);
                EXPECT_NEAR(C((*F[q])[at],k), fd, 2e-6*std::max(1.0,std::fabs(fd))) << name[q] << ", atom " << at << " axis " << k;
                largest=std::max(largest, std::fabs(fd));
            }
        }
    }
    EXPECT_GT(largest, 0.1) << "the contractions move with the atoms";
    for (size_t q=0;q<5;q++) EXPECT_LT(norm(sum[q]), 1e-8) << name[q] << ": translation invariance";
    // A partner with no whole-system face (the raw AO block itself) cannot fold into this block's AO space: refused by name, Can agreeing.
    const BasisSet::Orbital_HF_IBS<double>& rawHF=*raw;
    EXPECT_FALSE(B[0]->CanPullbackERI4(&rawHF));
    const auto r=B[0]->MakeDirectPullback(D[0], rsmat_t(raw->GetNumFunctions()), &rawHF, h2o.get());
    EXPECT_FALSE(bool(r)); if (!r) EXPECT_FALSE(r.Error().why.empty());
}

// GATE 7, THE HAND-GRADIENT AUDIT (lane B 2026-11-09, ahead of the user's Enzyme decision): every other point-Gaussian family's hand
// Gradient(r) against the central difference of its own operator()(r), every function, at points near (≈0.3 bohr) and far (≈3 bohr) from
// each atom.  PG_Spherical (s through g uncontracted, and the contracted DZVP), PG_LibCint Cartesian and spherical (the same two), and a
// SALC block of water's DZVP (Oᵀ applied to the raw gradients).  Step and tolerance as GATE 2; the largest derivative probed must be
// 0.1 or more (each family's is O(1) except the SALC A2 block, d functions only, ~0.36), so a sign or a missing factor would show.
TEST(GeometryDerivatives, EveryPointGaussianFamilysHandGradientIsTheFiniteDifferenceOfItsValues)
{
    const auto h2o=Water();
    const auto R=Sites(*h2o);
    BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile("dzvp.bsd"));
    std::vector<std::pair<std::string, std::unique_ptr<BasisSet::Orbital_1E_IBS<double>>>> blocks;
    blocks.emplace_back("PG_Spherical spdfg", std::make_unique<BasisSet::Gaussian::PG_Spherical::Orbital_IBS>(rvec_t{1.3, 0.35}, 4, h2o.get()));
    blocks.emplace_back("PG_Spherical DZVP",  std::make_unique<BasisSet::Gaussian::PG_Spherical::Orbital_IBS>(&reader, h2o.get()));
    for (bool sph : {false, true})
    {
        const std::string tag=sph ? "PG_LibCint spherical " : "PG_LibCint Cartesian ";
        blocks.emplace_back(tag+"spdfg", std::make_unique<BasisSet::Gaussian::PG_LibCint::Orbital_IBS>(rvec_t{1.3, 0.35}, 4, h2o.get(), sph));
        blocks.emplace_back(tag+"DZVP",  std::make_unique<BasisSet::Gaussian::PG_LibCint::Orbital_IBS>(&reader, h2o.get(), sph));
    }
    const auto raw=Dzvp(h2o.get());
    const auto pts=BasisSet::Gaussian::PG_Cart::StructureToSymPoints(*h2o);
    const auto grp=Symmetry::Molecule::BuildAbelianGroup(pts, 1e-4);
    const auto salc=Symmetry::Molecule::BuildSALCs(BasisSet::Gaussian::PG_Cart::ExtractAoShells(*raw), grp, Symmetry::Molecule::Centroid(pts), 1e-4);
    const BasisSet::Gaussian::SymmetryAdaptedBasisSet sab(raw.get(), salc);
    std::vector<std::pair<std::string, const BasisSet::Orbital_1E_IBS<double>*>> all;
    for (const auto& [n,b] : blocks) all.emplace_back(n, b.get());
    for (auto oi : sab.Iterate<BasisSet::Real_OIBS>()) all.emplace_back("SALC "+oi->Name(), dynamic_cast<const BasisSet::Orbital_1E_IBS<double>*>(oi));
    const double h=1e-5;
    for (const auto& [name, b] : all)
    {
        ASSERT_NE(b, nullptr) << name;
        double largest=0.0;
        for (const rvec3_t& base : R)
            for (const rvec3_t& off : {rvec3_t(0.17,-0.21,0.12), rvec3_t(1.9,1.2,-1.8)})
            {
                const rvec3_t r=base+off;
                const rvec3vec_t g=b->Gradient(r);
                ASSERT_EQ(g.size(), b->GetNumFunctions()) << name;
                for (int k=0;k<3;k++)
                {
                    rvec3_t p=r, m=r; C(p,k)+=h; C(m,k)-=h;
                    const rvec_t vp=(*b)(p), vm=(*b)(m);
                    for (size_t i=0;i<b->GetNumFunctions();i++)
                    {
                        const double fd=(vp[i]-vm[i])/(2*h);
                        EXPECT_NEAR(C(g[i],k), fd, 1e-7*std::max(1.0,std::fabs(fd))) << name << ": function " << i << " axis " << k << " at " << r;
                        largest=std::max(largest, std::fabs(fd));
                    }
                }
            }
        EXPECT_GT(largest, 0.1) << name << ": the derivatives probed are not trivially small";
    }
}
