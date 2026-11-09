// File: Common/tests/PositionJet.C  PositionJet<T> (qchem.PositionJet, qcMath): the first-order jet in the atom positions.  Its arithmetic is
// the chain rule, checked against analytic derivatives and central finite differences; it stays sparse in atoms; a translation-invariant
// quantity's differential sums to zero; the pullback accumulates a weighted differential into a force array.
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
import qchem.PositionJet;
import qchem.Vector3D;

using namespace qchem;
using Jet=PositionJet<double>;
using V3=Vector3D<double>;

namespace
{
const std::vector<V3> kR={V3(0.0,0.0,0.1), V3(1.4,0.2,-0.9), V3(-1.3,0.5,-0.8), V3(0.3,2.2,0.6)};
//! |R_a - R_b| as a jet: built from coordinate jets, so its differential is the chain rule's.
Jet Distance(const std::vector<V3>& R, size_t a, size_t b)
{
    Jet s(0.0);
    for (int k=0;k<3;k++)
    {
        const double ra = k==0 ? R[a].x : k==1 ? R[a].y : R[a].z, rb = k==0 ? R[b].x : k==1 ? R[b].y : R[b].z;
        const Jet d=Jet::Coordinate(a,k,ra)-Jet::Coordinate(b,k,rb);
        s+=d*d;
    }
    return sqrt(s);
}
//! A composite with every operation in it: a Becke-like ratio of smooth functions of three distances.
Jet Composite(const std::vector<V3>& R)
{
    const Jet r01=Distance(R,0,1), r02=Distance(R,0,2), r13=Distance(R,1,3);
    const Jet mu=(r01-r02)/r13;
    const Jet s=0.5*(1.0-mu*(3.0-mu*mu)*0.5);
    return s*exp(-0.3*r01)/(s+r02*0.1)+mu*mu;
}
double& C(V3& v, int k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
double  C(const V3& v, int k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
}

TEST(PositionJet, ADistancesDifferentialIsTheUnitVectorAndTouchesOnlyItsTwoAtoms)
{
    const Jet r=Distance(kR,1,3);
    EXPECT_NEAR(r.Value(), norm(kR[1]-kR[3]), 1e-15);
    ASSERT_EQ(r.Derivatives().size(), 2u) << "sparse: atoms 1 and 3 only";
    const V3 e=(1.0/norm(kR[1]-kR[3]))*(kR[1]-kR[3]);
    EXPECT_LT(norm(r.D(1)-e), 1e-15);
    EXPECT_LT(norm(r.D(3)+e), 1e-15);
    EXPECT_EQ(norm(r.D(0)), 0.0) << "an absent atom has an exactly zero derivative";
}

TEST(PositionJet, ACompositesDifferentialIsTheFiniteDifference)
{
    const Jet f=Composite(kR);
    const double h=1e-6;
    for (size_t A=0;A<kR.size();A++)
        for (int k=0;k<3;k++)
        {
            auto p=kR, m=kR; C(p[A],k)+=h; C(m[A],k)-=h;
            EXPECT_NEAR(C(f.D(A),k), (Composite(p).Value()-Composite(m).Value())/(2*h), 1e-9) << "atom " << A << " axis " << k;
        }
}

TEST(PositionJet, ATranslationInvariantQuantitysDifferentialSumsToZero)
{
    const Jet f=Composite(kR);
    V3 sum(0,0,0);
    for (const auto& e : f.Derivatives()) sum+=e.d;
    EXPECT_LT(norm(sum), 1e-14);
}

TEST(PositionJet, ThePullbackAccumulatesAWeightedDifferential)
{
    const Jet f=Distance(kR,0,2), g=Distance(kR,1,3);
    std::vector<V3> F(kR.size(), V3(0,0,0));
    f.AddWeightedTo(F, 2.0);
    g.AddWeightedTo(F, -0.5);
    const Jet sum=2.0*f-0.5*g;   // the same thing, by arithmetic
    for (size_t A=0;A<kR.size();A++) EXPECT_LT(norm(F[A]-sum.D(A)), 1e-15) << "atom " << A;
    std::vector<V3> small(2);
    EXPECT_THROW(g.AddWeightedTo(small, 1.0), std::out_of_range) << "atom 3 is past a two-atom force array";
}

TEST(PositionJet, MergingKeepsTheAtomListSortedWithNoDuplicates)
{
    const Jet a=Jet::Coordinate(3,0,1.0)+Jet::Coordinate(1,1,2.0), b=Jet::Coordinate(2,2,0.5)*Jet::Coordinate(1,0,3.0);
    const Jet c=a*b+a;
    size_t last=0; bool first=true;
    for (const auto& e : c.Derivatives()) { if (!first) EXPECT_LT(last, e.atom); last=e.atom; first=false; }
    EXPECT_EQ(c.Derivatives().size(), 3u) << "atoms 1, 2, 3";
    EXPECT_THROW(Jet::Coordinate(0,3,1.0), std::invalid_argument);
}
