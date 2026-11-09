// File: M_PP_ClosedForm.C  The closed-form molecular PP integrals (Gaussian::EOrbital_PP_IBS) against an
// independent route: the same field quadratured on a fine Becke mesh.  The quadrature converges to the
// closed form, and a wrong sign, norm or Y_lm convention misses by far more than the tolerance.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <memory>
#include <vector>
import qchem.BasisSet;
import qchem.BasisSet.Orbital_PP_IBS;
import qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Point.PG_Spherical;
import qchem.Math;
import qchem.Math.Gaussian;
import qchem.Math.SolidHarmonics;
import qchem.Mesh;
import qchem.Mesh.Quadrature;
import qchem.Mesh.XCPolicy;
import qchem.ScalarFunction;
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using BasisSet::FieldRange;

namespace
{
//! An HGH-shaped field, with every short-range term present so the r^6 monomials are reached.
struct TestField : public virtual BasisSet::SpeciesRadialField_Gaussian, public virtual BasisSet::SpeciesRadialField_GaussianCharge
{
    static double Rloc(int Z) {return Z==9 ? 0.40 : 0.55;}
    static double Zion(int Z) {return Z==9 ? 7.0  : 4.0;}
    double ValueR(int Z, double r, FieldRange range) const override
    {
        double v=0;
        if (range!=FieldRange::Long)
            for (const auto& t : AsGaussians(Z,FieldRange::Short)) v+=t.c*std::pow(r,2*t.n)*std::exp(-t.alpha*r*r);
        if (range!=FieldRange::Short)
        {
            const auto q=AsGaussianCharge(Z);
            v-= r<1e-10 ? q.q*2*std::sqrt(q.alpha/Pi) : q.q*std::erf(std::sqrt(q.alpha)*r)/r;
        }
        return v;
    }
    double ValueQ(int, double, FieldRange) const override {return 0.0;}
    double CellMeanQ(int, FieldRange) const override {return 0.0;}
    std::vector<Math::Gaussian> AsGaussians(int Z, FieldRange range) const override
    {
        if (range!=FieldRange::Short) return {};
        const double a=0.5/(Rloc(Z)*Rloc(Z));
        return {{-3.1,0,a},{0.9,1,a},{-0.35,2,a},{0.12,3,a}};
    }
    BasisSet::GaussianCharge AsGaussianCharge(int Z) const override {return {Zion(Z), 0.5/(Rloc(Z)*Rloc(Z))};}
};

//! One projector per l of 0..2, the l=1 one with a two-term radial.
struct TestProjectors : public virtual BasisSet::SpeciesProjectorSet_R, public virtual BasisSet::SpeciesProjectorSet_Gaussian
{
    size_t Count(int) const override {return 3;}
    int    L(int, size_t p) const override {return int(p);}
    double Weight(int, size_t p) const override {return 1.3-0.8*p;}
    double RadialQ(int, size_t, double) const override {return 0.0;}
    double SharpnessR(int, size_t) const override {return 1.2;}
    std::vector<Math::Gaussian> AsGaussians(int, size_t p) const override
    {
        if (p==1) return {{0.7,0,1.2},{-0.4,1,1.2}};
        return {{1.1,0,1.2}};
    }
    double RadialR(int Z, size_t p, double r) const override
    {
        double s=0;
        for (const auto& t : AsGaussians(Z,p)) s+=t.c*std::pow(r,2*t.n)*std::exp(-t.alpha*r*r);
        return s*std::pow(r,int(p));
    }
};

struct FieldOnMesh : public ScalarFunction<double>
{
    FieldOnMesh(const Structure& s, const BasisSet::SpeciesRadialField& f, FieldRange r) : st(s), v(f), range(r) {}
    double operator()(const rvec3_t& r) const override
    {
        double s=0;
        for (const Site& a : st.Sites()) s+=v.ValueR(a.Z, norm(r-a.R), range);
        return s;
    }
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}
    const Structure& st; const BasisSet::SpeciesRadialField& v; FieldRange range;
};

struct BetaOnMesh : public ScalarFunction<double>
{
    BetaOnMesh(const rvec3_t& c, const TestProjectors& ps, int z, size_t pp, int mm) : R(c), v(ps), Z(z), p(pp), m(mm) {}
    double operator()(const rvec3_t& r) const override
    {
        const rvec3_t d=r-R; const double rr=norm(d); const int l=v.L(Z,p);
        if (rr<1e-12) return l==0 ? v.RadialR(Z,p,0.0)*Math::RealYlm(0,0,0,0,0) : 0.0;
        return v.RadialR(Z,p,rr)*Math::RealYlm(l,m,d.x/rr,d.y/rr,d.z/rr);
    }
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}
    rvec3_t R; const TestProjectors& v; int Z; size_t p; int m;
};

std::shared_ptr<Molecule> TwoSpecies()
{
    auto m=std::make_shared<Molecule>();
    m->Insert(Site{9, rvec3_t(0,0,0)});
    m->Insert(Site{14,rvec3_t(0.6,1.1,2.0)});   // no symmetry axis: every (l,m) channel is reached
    return m;
}

double MaxDiff(const rsmat_t& a, const rsmat_t& b)
{
    double d=0;
    for (size_t i=0;i<a.rows();i++) for (size_t j=0;j<a.columns();j++) d=std::max(d,std::fabs(a(i,j)-b(i,j)));
    return d;
}

template <class IBS> void CheckAgainstMesh()
{
    const auto st=TwoSpecies();
    const IBS bs(rvec_t{0.45,1.6}, 2, st.get());
    const qcMesh::Mesh mesh=st->CreateIntegrationMesh(qcMesh::BeckeXCParams(200,-1.0,41));
    const TestField field;
    ASSERT_TRUE(bs.CanAssembleField(field));
    for (FieldRange range : {FieldRange::Short, FieldRange::Long, FieldRange::Full})
        EXPECT_LT(MaxDiff(bs.MakeSpeciesFieldMatrix(st.get(),field,range),
                          qcMesh::MatrixOverlap(mesh,bs,FieldOnMesh(*st,field,range))), 2e-6);

    const TestProjectors ps;
    ASSERT_TRUE(bs.CanAssembleProjectors(ps));
    const size_t n=bs.GetVectorSize();
    rmat_t V(n,n,0.0);
    for (const Site& a : st->Sites())
        for (size_t p=0;p<ps.Count(a.Z);p++)
            for (int m=-ps.L(a.Z,p); m<=ps.L(a.Z,p); m++)
            {
                const rvec_t b=qcMesh::Overlap(mesh,bs,BetaOnMesh(a.R,ps,a.Z,p,m));
                for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) V(i,j)+=ps.Weight(a.Z,p)*b[i]*b[j];
            }
    EXPECT_LT(MaxDiff(bs.MakeProjectorMatrix(st.get(),ps), rsmat_t(V)), 2e-6);
}
} //anon

TEST(PP_ClosedForm, CartesianEqualsFineMesh) {CheckAgainstMesh<BasisSet::Gaussian::PG_Cart::Orbital_IBS>();}
TEST(PP_ClosedForm, SphericalEqualsFineMesh) {CheckAgainstMesh<BasisSet::Gaussian::PG_Spherical::Orbital_IBS>();}

//! A field without the closed-form faces stays on the term's quadrature route.
TEST(PP_ClosedForm, RefusesAFieldWithNoClosedForm)
{
    struct Plain : public virtual BasisSet::SpeciesRadialField
    {
        double ValueR(int, double, FieldRange) const override {return 0.0;}
        double ValueQ(int, double, FieldRange) const override {return 0.0;}
        double CellMeanQ(int, FieldRange) const override {return 0.0;}
    };
    const auto st=TwoSpecies();
    const BasisSet::Gaussian::PG_Cart::Orbital_IBS bs(rvec_t{0.45}, 0, st.get());
    EXPECT_FALSE(bs.CanAssembleField(Plain()));
}

// THE PULLBACKS (GeometryOptimizationPlan §13.2; lane B 2026-11-08): PG_Cart's MakeSpeciesFieldPullback (each range) and
// MakeProjectorPullback, libcint's derivative integrals with each field Gaussian an auxiliary shell, against the central difference of
// Tr(D M) of the SAME closed-form matrices (our M&D engine) at displaced atoms, the block rebuilt on the moved sites.  Every term of the
// test field is reached (r^0..r^6 short-range monomials, the erf long part, l = 0..2 projectors, a two-term radial), on a two-species
// pair with no symmetry axis.  Step h = 1e-4: the truncation error h²/6 |f‴| is ~2e-9 |f‴|, and with a random ±1 D over the block
// |f‴| stays below ~1e2 relative to |f′|, so 2e-6 (relative above 1) is ~10× the estimate.  Translation invariance: the basis and the
// species' sites all ride on the two atoms.
TEST(PP_ClosedForm, CartesianPullbacksAreTheFiniteDifferenceOfItsOwnMatricesAndTranslationInvariant)
{
    const auto st=TwoSpecies();
    using IBS=BasisSet::Gaussian::PG_Cart::Orbital_IBS;
    const IBS bs(rvec_t{0.45,1.6}, 2, st.get());
    const size_t n=bs.GetVectorSize();
    std::mt19937 rng(5); std::uniform_real_distribution<double> U(-1,1);
    hmat_t<double> D(n); for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=U(rng);
    const TestField field; const TestProjectors ps;
    ASSERT_TRUE(bs.CanPullbackSpeciesField(field)); ASSERT_TRUE(bs.CanPullbackProjectors(ps));
    auto TrDM=[&](const rsmat_t& M){ double t=0.0; for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) t+=D(i,j)*M(i,j); return t; };
    std::vector<rvec3_t> R; st->ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
    const double h=1e-4;
    auto check=[&](const rvec3vec_t& F, auto value, const char* what)
    {
        rvec3_t sum(0,0,0); double largest=0.0;
        for (size_t A=0;A<R.size();A++)
        {
            sum+=F[A];
            for (int k=0;k<3;k++)
            {
                auto at=[&](double s){ auto Rs=R; (k==0 ? Rs[A].x : k==1 ? Rs[A].y : Rs[A].z)+=s; std::shared_ptr<Structure> m=st->Clone(); m->MoveSites(Rs);
                                       const IBS b(rvec_t{0.45,1.6}, 2, m.get()); return value(b, m.get()); };
                const double fd=(at(h)-at(-h))/(2*h), an = k==0 ? F[A].x : k==1 ? F[A].y : F[A].z;
                EXPECT_NEAR(an, fd, 2e-6*std::max(1.0,std::fabs(fd))) << what << ", atom " << A << " axis " << k;
                largest=std::max(largest, std::fabs(fd));
            }
        }
        EXPECT_GT(largest, 1e-2) << what << " moves with the atoms";
        EXPECT_LT(norm(sum), 1e-9) << what << ": translation invariance";
    };
    for (FieldRange range : {FieldRange::Short, FieldRange::Long, FieldRange::Full})
    {
        const auto F=bs.MakeSpeciesFieldPullback(st.get(), field, range, D);
        ASSERT_TRUE(F);
        check(*F, [&](const IBS& b, const Structure* m){ return TrDM(b.MakeSpeciesFieldMatrix(m, field, range)); },
              range==FieldRange::Short ? "local short" : range==FieldRange::Long ? "local long" : "local full");
    }
    const auto P=bs.MakeProjectorPullback(st.get(), ps, D);
    ASSERT_TRUE(P);
    check(*P, [&](const IBS& b, const Structure* m){ return TrDM(b.MakeProjectorMatrix(m, ps)); }, "projectors");
    EXPECT_THROW((void)bs.MakeSpeciesFieldPullback(st.get(), field, FieldRange::Full, hmat_t<double>(n+1)), std::invalid_argument) << "a mismatched D";
    EXPECT_THROW((void)bs.MakeProjectorPullback(st.get(), ps, hmat_t<double>(n+1)), std::invalid_argument) << "a mismatched D";
}

//! The pullbacks' Can queries agree with their methods: PG_Cart answers both; PG_Spherical refuses both with a reason.
TEST(PP_ClosedForm, EachPullbackCanQueryAgreesWithWhetherItsMethodAnswers)
{
    const auto st=TwoSpecies();
    const BasisSet::Gaussian::PG_Cart::Orbital_IBS cart(rvec_t{0.45,1.6}, 2, st.get());
    const BasisSet::Gaussian::PG_Spherical::Orbital_IBS sph(rvec_t{0.45,1.6}, 2, st.get());
    const TestField field; const TestProjectors ps;
    for (const BasisSet::Orbital_PP_IBS<double>* b : {static_cast<const BasisSet::Orbital_PP_IBS<double>*>(&cart), static_cast<const BasisSet::Orbital_PP_IBS<double>*>(&sph)})
    {
        const size_t n=dynamic_cast<const BasisSet::Orbital_1E_IBS<double>&>(*b).GetNumFunctions();
        const auto f=b->MakeSpeciesFieldPullback(st.get(), field, FieldRange::Full, hmat_t<double>(n));
        const auto p=b->MakeProjectorPullback(st.get(), ps, hmat_t<double>(n));
        EXPECT_EQ(b->CanPullbackSpeciesField(field), bool(f)); EXPECT_EQ(b->CanPullbackProjectors(ps), bool(p));
        if (!f) EXPECT_FALSE(f.Error().why.empty());
        if (!p) EXPECT_FALSE(p.Error().why.empty());
    }
    EXPECT_TRUE(cart.CanPullbackSpeciesField(field) && cart.CanPullbackProjectors(ps)) << "PG_Cart gives them";
    EXPECT_FALSE(sph.CanPullbackSpeciesField(field) || sph.CanPullbackProjectors(ps)) << "PG_Spherical has no geometry derivative yet";
}
