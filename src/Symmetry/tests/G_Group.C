// File: Symmetry/tests/G_Group.C  The Group face (Group-face note, migration step 1): Contains, the continuous-family table,
// T(L1) ∩ T(L2), crystal ∩ O(3)_c = the site stabiliser, crystal ∩ crystal, centring and magnetic cosets.
//
// A group has no public order or element list (user 2026-11-19), so equality is asserted the way a client would: G ⊆ H and H ⊆ G.
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
import qchem.Symmetry.Group;
import qchem.Symmetry.Lattice_3D.SpaceGroup;
import qchem.Symmetry.Molecule.PointGroup;
using namespace qchem;
using namespace qchem::Symmetry;
using Mat=Matrix3D<double>;

namespace
{
const GroupTolerance tol{1e-8, 1e-10};

Mat Diag(double a, double b, double c) {return Mat(a,0,0, 0,b,0, 0,0,c);}
//! Columns a, b, c.
Mat Cols(const rvec3_t& a, const rvec3_t& b, const rvec3_t& c) {return Mat(a.x,b.x,c.x, a.y,b.y,c.y, a.z,b.z,c.z);}
Mat I3() {return Diag(1,1,1);}
//! An op about the point o: r -> R (r - o) + o.
SymOp About(const Mat& R, const rvec3_t& o, SpinAction s=SpinAction::None) {return {R, o-R*o, s};}
SymOp Translation(const rvec3_t& t, SpinAction s=SpinAction::None) {return {I3(), t, s};}
Mat C4(const rvec3_t& axis) {return Molecule::SymOp::Cn(axis, 4).Matrix();}
Mat C2(const rvec3_t& axis) {return Molecule::SymOp::Cn(axis, 2).Matrix();}
Mat MirrorNormal(const rvec3_t& n) {return Molecule::SymOp::Sigma(n).Matrix();}
Mat Inversion() {return Diag(-1,-1,-1);}

bool Has(const GroupPtr& g, const SymOp& op) {return g->Contains(Generator::Op(op));}
void ExpectEqual(const GroupPtr& a, const GroupPtr& b)
{
    EXPECT_TRUE(a->Contains(*b)) << a->Describe() << " should contain " << b->Describe();
    EXPECT_TRUE(b->Contains(*a)) << b->Describe() << " should contain " << a->Describe();
}
GroupPtr Meet(const GroupPtr& a, const GroupPtr& b)
{
    auto r=a->Intersect(*b);
    if (!r) ADD_FAILURE() << r.Error();
    auto s=b->Intersect(*a);                    // the table is symmetric: both orders give the same group
    if (!s) ADD_FAILURE() << s.Error();
    if (r && s) ExpectEqual(*r, *s);
    return r ? *r : GroupPtr();
}
//! Simple cubic, one atom at the origin: Pm-3m, order 48.
Lattice_3D::SpaceGroup SimpleCubic(double a) {return Lattice_3D::SpaceGroup::Detect(Diag(a,a,a), {{1, rvec3_t(0,0,0)}});}
} //anon

// ---------------------------------------------------------------------------------------------------------------------------
TEST(Group, AGeneratorReadsOnlyItsOwnRealization)
{
    const Generator g=Generator::FullRotations(rvec3_t(0,0,0));
    EXPECT_EQ(g.GetKind(), Generator::Kind::FullRotations);
    EXPECT_NO_THROW(g.GetCentre());
    EXPECT_THROW(g.GetOp(), std::logic_error);
    EXPECT_THROW(g.GetLattice(), std::logic_error);
    EXPECT_THROW(Generator::Translations(Diag(1,1,0)), std::invalid_argument);
}

TEST(Group, ContainsAGroupIsContainsEveryGenerator)
{
    const rvec3_t c(0.3,-0.2,0.7), n(1,2,2);
    const GroupPtr O3=PointRotationGroup(c, tol), line=LineRotationGroup(c, n, tol), lineElsewhere=LineRotationGroup(c+rvec3_t(1,0,0), n, tol);
    EXPECT_TRUE (O3->Contains(*line));
    EXPECT_FALSE(line->Contains(*O3));
    EXPECT_FALSE(O3->Contains(*lineElsewhere));
    EXPECT_TRUE (O3->Contains(Generator::AxialRotations(c-2.0*n, n)))    << "an axis through c, wherever it is stated";
    EXPECT_FALSE(O3->Contains(Generator::Translations(I3())))            << "a finite or rotation group holds no lattice";
    EXPECT_TRUE (O3->Contains(Generator::Op(About(Inversion(), c))))     << "O(3) holds the improper ops";
    EXPECT_FALSE(O3->Contains(Generator::Op(About(Inversion(), c, SpinAction::Flip)))) << "a spatial group holds no spin flip";
    EXPECT_EQ(O3->ActsOn(), Space::R);
}

// ---------------------------------------------------------------------------------------------------------------------------
//  The continuous-family table (note §2.3): two pointwise stabilisers meet in the stabiliser of the affine hull.
TEST(Group, TwoFullRotationGroupsMeetInTheAxialGroupOfTheirLine)
{
    const rvec3_t c(0,0,0), cp(0,0,2.5), z(0,0,1);
    ExpectEqual(Meet(PointRotationGroup(c, tol), PointRotationGroup(c, tol)), PointRotationGroup(c, tol));
    const GroupPtr m=Meet(PointRotationGroup(c, tol), PointRotationGroup(cp, tol));
    ASSERT_TRUE(m);
    ExpectEqual(m, LineRotationGroup(c, z, tol));
    EXPECT_TRUE (m->Contains(Generator::AxialRotations(cp, z)));
    EXPECT_TRUE (Has(m, About(MirrorNormal(rvec3_t(1,1,0)), c)))  << "a mirror containing the axis";
    EXPECT_FALSE(Has(m, About(C2(rvec3_t(1,0,0)), c)))             << "a C2 perpendicular to the axis moves c'";
    EXPECT_FALSE(m->Contains(Generator::FullRotations(c)));
}

TEST(Group, ALineAndAPointOffItMeetInTheirPlanesMirror)
{
    const GroupPtr m=Meet(LineRotationGroup(rvec3_t(0,0,0), rvec3_t(0,0,1), tol), PointRotationGroup(rvec3_t(1.5,0,0), tol));
    ASSERT_TRUE(m);
    EXPECT_TRUE (Has(m, About(MirrorNormal(rvec3_t(0,1,0)), rvec3_t(0,0,0)))) << "the mirror in the xz plane fixes both";
    EXPECT_FALSE(Has(m, About(MirrorNormal(rvec3_t(1,0,0)), rvec3_t(0,0,0))));
    EXPECT_FALSE(m->Contains(Generator::AxialRotations(rvec3_t(0,0,0), rvec3_t(0,0,1))));
    ExpectEqual(m, FiniteGroup(std::vector<SymOp>{About(MirrorNormal(rvec3_t(0,1,0)), rvec3_t(0,0,0))}, tol));
}

TEST(Group, TwoSkewLinesShareOnlyTheIdentity)
{
    const GroupPtr m=Meet(LineRotationGroup(rvec3_t(0,0,0), rvec3_t(0,0,1), tol), LineRotationGroup(rvec3_t(1,0,0), rvec3_t(0,1,0), tol));
    ASSERT_TRUE(m);
    ExpectEqual(m, FiniteGroup(std::vector<SymOp>{}, tol));
    // parallel distinct lines span a plane: its mirror survives
    const GroupPtr p=Meet(LineRotationGroup(rvec3_t(0,0,0), rvec3_t(0,0,1), tol), LineRotationGroup(rvec3_t(1,0,0), rvec3_t(0,0,-1), tol));
    ASSERT_TRUE(p);
    ExpectEqual(p, FiniteGroup(std::vector<SymOp>{About(MirrorNormal(rvec3_t(0,1,0)), rvec3_t(0,0,0))}, tol));
}

TEST(Group, ALatticeMeetsARotationGroupOnlyAtTheIdentity)
{
    const GroupPtr m=Meet(TranslationGroup(Diag(2,3,4), tol), PointRotationGroup(rvec3_t(0.5,0.5,0.5), tol));
    ASSERT_TRUE(m);
    ExpectEqual(m, FiniteGroup(std::vector<SymOp>{}, tol));
}

// ---------------------------------------------------------------------------------------------------------------------------
//  T(L1) ∩ T(L2) = T(L1 ∩ L2): exact, by integer arithmetic.
TEST(Group, TwoLatticesMeetInTheirCommonSublattice)
{
    struct Case {Mat L1, L2, common; const char* what;};
    const rvec3_t a(1,1,0), b(1,-1,0), z(0,0,1);
    const std::vector<Case> cases={
        {Diag(2,1,1),  Diag(1,3,1),     Diag(2,3,1),         "orthogonal: the lcm per axis"},
        {I3(),         Diag(1.5,1.5,1.5), Diag(3,3,3),       "a rational ratio 3/2"},
        {I3(),         Cols(a, b, z),   Cols(a, b, z),       "a rotated sublattice is its own intersection"},
        {Cols(a, b, z), Diag(1,1,2),    Cols(a, b, 2.0*z),   "neither nested: a+b even and c even"},
    };
    for (const Case& c:cases)
    {
        const GroupPtr m=Meet(TranslationGroup(c.L1, tol), TranslationGroup(c.L2, tol));
        ASSERT_TRUE(m) << c.what;
        ExpectEqual(m, TranslationGroup(c.common, tol));
        EXPECT_TRUE(m->Contains(Generator::Translations(c.common))) << c.what;
    }
}

TEST(Group, IncommensurateLatticesAreARefusalNotAnApproximation)
{
    const auto r=TranslationGroup(I3(), tol)->Intersect(*TranslationGroup(Diag(std::sqrt(2.0), 1, 1), tol));
    ASSERT_FALSE(r);
    EXPECT_NE(r.Error().find("incommensurate"), std::string::npos) << r.Error();
}

// ---------------------------------------------------------------------------------------------------------------------------
//  Crystals.
TEST(Group, ACrystalMeetsO3AtAPointInItsSiteStabiliser)
{
    const double a=2.0;
    const auto sg=SimpleCubic(a);
    ASSERT_EQ(sg.Order(), 48u);
    const GroupPtr G=CrystalGroup(sg, tol);
    for (const rvec3_t f:{rvec3_t(0,0,0), rvec3_t(0.5,0.5,0.5), rvec3_t(0.5,0,0), rvec3_t(0.1,0.2,0.3)})
    {
        const rvec3_t c=a*f;
        const GroupPtr m=Meet(G, PointRotationGroup(c, tol));
        ASSERT_TRUE(m);
        std::vector<SymOp> expected;
        for (const auto& op:sg.SiteStabilizer(f)) expected.push_back({op.W, a*op.tau, SpinAction::None});   // cubic: A W A^-1 = W
        // SiteStabilizer's ops map f to f mod 1; about c they are {W | c - W c}
        std::vector<SymOp> aboutC;
        for (const SymOp& g:expected) aboutC.push_back(About(g.W, c));
        ExpectEqual(m, FiniteGroup(aboutC, tol));
    }
    const GroupPtr edge=Meet(G, PointRotationGroup(rvec3_t(a/2,0,0), tol));
    EXPECT_TRUE (Has(edge, About(C4(rvec3_t(1,0,0)), rvec3_t(a/2,0,0)))) << "C4 about x fixes the edge centre";
    EXPECT_FALSE(Has(edge, About(C4(rvec3_t(0,0,1)), rvec3_t(a/2,0,0))));
}

TEST(Group, ACrystalMeetsASupercellsTranslationsInThoseTranslations)
{
    const double a=1.5;
    const GroupPtr G=CrystalGroup(SimpleCubic(a), tol);
    ExpectEqual(Meet(G, TranslationGroup(Diag(2*a,2*a,2*a), tol)), TranslationGroup(Diag(2*a,2*a,2*a), tol));
    ExpectEqual(Meet(G, TranslationGroup(Diag(a/2,a/2,a/2), tol)), TranslationGroup(Diag(a,a,a), tol));
    EXPECT_TRUE (G->Contains(Generator::Translations(Diag(2*a,a,3*a))));
    EXPECT_FALSE(G->Contains(Generator::Translations(Diag(a/2,a,a))));
}

TEST(Group, ACrystalMeetsItsSupercellsGroupInTheSupercellsGroup)
{
    const double a=1.0;
    const GroupPtr G=CrystalGroup(SimpleCubic(a), tol), G2=CrystalGroup(SimpleCubic(2*a), tol);
    ExpectEqual(Meet(G, G2), G2);
    EXPECT_TRUE (Has(G2, About(C4(rvec3_t(0,0,1)), rvec3_t(0,0,0))));
    EXPECT_FALSE(G2->Contains(Generator::Translations(Diag(a,a,a))));
    EXPECT_TRUE (G->Contains(*G2));
    EXPECT_FALSE(G2->Contains(*G));
    // shifted by half a cell: the supercell group about (a/2, 0, 0) keeps only the ops that fix the x = a/2 plane family
    const GroupPtr shifted=CrystalGroup(Lattice_3D::SpaceGroup::Detect(Diag(2*a,2*a,2*a), {{1, rvec3_t(0.25,0,0)}}), tol);
    const GroupPtr m=Meet(G, shifted);
    ASSERT_TRUE(m);
    EXPECT_TRUE (m->Contains(Generator::Translations(Diag(2*a,2*a,2*a))));
    EXPECT_TRUE (Has(m, About(C4(rvec3_t(1,0,0)), rvec3_t(a/2,0,0))));
    EXPECT_FALSE(Has(m, About(C4(rvec3_t(0,0,1)), rvec3_t(a/2,0,0)))) << "C4 about z through (a/2,0,0) moves the lattice points of G";
}

TEST(Group, ACentringTranslationJoinsTheLattice)
{
    // bcc in its conventional cubic cell: two identical atoms, every tau enumerated (ShubnikovOps on an undecorated basis)
    const Mat A=Diag(1,1,1);
    const std::vector<Lattice_3D::AtomSite> basis={{26, rvec3_t(0,0,0)}, {26, rvec3_t(0.5,0.5,0.5)}};
    const auto sg=Lattice_3D::SpaceGroup::Detect(A, basis);
    const GroupPtr G=CrystalGroup(sg.ShubnikovOps(basis), A, tol);
    const Mat bcc=Cols(rvec3_t(-0.5,0.5,0.5), rvec3_t(0.5,-0.5,0.5), rvec3_t(0.5,0.5,-0.5));
    EXPECT_TRUE(G->Contains(Generator::Translations(bcc)));
    EXPECT_TRUE(Has(G, Translation(rvec3_t(0.5,0.5,0.5))));
    EXPECT_EQ(G->ActsOn(), Space::R);
    ExpectEqual(Meet(G, TranslationGroup(I3(), tol)), TranslationGroup(I3(), tol));
}

TEST(Group, AMagneticGroupPairsTheSublatticeTranslationWithAFlip)
{
    // CsCl-type AFM: the two sublattices carry opposite moments, so {E|½½½} survives only with the flip (note §2.4)
    const Mat A=Diag(1,1,1);
    const std::vector<Lattice_3D::AtomSite> grey={{26, rvec3_t(0,0,0)}, {26, rvec3_t(0.5,0.5,0.5)}};
    const std::vector<Lattice_3D::AtomSite> afm ={{26, rvec3_t(0,0,0), +1}, {26, rvec3_t(0.5,0.5,0.5), -1}};
    const auto sg=Lattice_3D::SpaceGroup::Detect(A, grey);
    const GroupPtr G=CrystalGroup(sg.ShubnikovOps(grey), A, tol), M=CrystalGroup(sg.ShubnikovOps(afm), A, tol);
    EXPECT_EQ(M->ActsOn(), Space::RSpin);
    EXPECT_TRUE (Has(M, Translation(rvec3_t(0.5,0.5,0.5), SpinAction::Flip)));
    EXPECT_FALSE(Has(M, Translation(rvec3_t(0.5,0.5,0.5))));
    EXPECT_FALSE(M->Contains(Generator::Translations(Cols(rvec3_t(-0.5,0.5,0.5), rvec3_t(0.5,-0.5,0.5), rvec3_t(0.5,0.5,-0.5)))));
    EXPECT_TRUE (M->Contains(Generator::Translations(I3())));
    EXPECT_FALSE(G->Contains(*M)) << "the grey group's elements carry no flip, so it is not a supergroup of the magnetic one";
    // the magnetic group meets O(3) at a moment-carrying site in the site's spatial ops, all without a flip
    const GroupPtr site=Meet(M, PointRotationGroup(rvec3_t(0,0,0), tol));
    ASSERT_TRUE(site);
    EXPECT_TRUE(Has(site, About(Inversion(), rvec3_t(0,0,0))));
    EXPECT_TRUE(Has(site, About(C4(rvec3_t(0,0,1)), rvec3_t(0,0,0))));
    EXPECT_EQ(site->ActsOn(), Space::R);
    ExpectEqual(Meet(M, M), M);
}

TEST(Group, AMolecularPointGroupIsClosedAndMeetsTheRest)
{
    const rvec3_t o(0.2,0.1,-0.3);
    const GroupPtr C4v=FiniteGroup({Molecule::SymOp::Cn(rvec3_t(0,0,1), 4), Molecule::SymOp::Sigma(rvec3_t(1,0,0))}, o, tol);
    EXPECT_TRUE (Has(C4v, About(C2(rvec3_t(0,0,1)), o)))              << "closure: C4^2";
    EXPECT_TRUE (Has(C4v, About(MirrorNormal(rvec3_t(1,1,0)), o)))    << "closure: a diagonal mirror";
    EXPECT_FALSE(Has(C4v, About(Inversion(), o)));
    ExpectEqual(Meet(C4v, PointRotationGroup(o, tol)), C4v);
    const GroupPtr offAxis=Meet(C4v, PointRotationGroup(o+rvec3_t(1,0,0), tol));
    ASSERT_TRUE(offAxis);
    ExpectEqual(offAxis, FiniteGroup(std::vector<SymOp>{About(MirrorNormal(rvec3_t(0,1,0)), o)}, tol));
    ExpectEqual(Meet(C4v, TranslationGroup(I3(), tol)), FiniteGroup(std::vector<SymOp>{}, tol));
    ExpectEqual(Meet(C4v, LineRotationGroup(o, rvec3_t(0,0,1), tol)), C4v);
    EXPECT_THROW(FiniteGroup(std::vector<SymOp>{Translation(rvec3_t(1,0,0))}, tol), std::invalid_argument) << "not a finite group";
}

//! {e} (Group face step 3a, what a fit that forces nothing answers): it holds the identity and nothing else, every group contains it, and it
//! contains no other group but itself.
TEST(Group, TheTrivialGroupIsTheIdentityAlone)
{
    const auto e=TrivialGroup();
    const rvec3_t c(0.2,0.1,-0.4);
    EXPECT_TRUE (Has(e, Translation(rvec3_t(0,0,0))));
    EXPECT_FALSE(Has(e, About(Inversion(), c)));
    EXPECT_FALSE(e->Contains(Generator::Translations(I3())));
    EXPECT_FALSE(e->Contains(Generator::FullRotations(c)));
    for (const GroupPtr& g : {PointRotationGroup(c, tol), TranslationGroup(Diag(3,4,5), tol), CrystalGroup(SimpleCubic(5.0), tol), e})
        EXPECT_TRUE(g->Contains(*e)) << g->Describe();
    EXPECT_FALSE(e->Contains(*PointRotationGroup(c, tol)));
    EXPECT_TRUE (e->Contains(*e));
}
