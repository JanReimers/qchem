// File: UnitTests/BasisSet_Atom1.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <iostream>
#include <nlohmann/json.hpp>
#include <blaze/math/expressions/DMatDMatEqualExpr.h> //op== inside gtest header.

using std::cout;
using std::endl;

import qchem.BasisSet.Radial.Factory;
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;  // tests may cheat: the ERI4 substrate face
import qchem.BasisSet.Orbital_DFT_IBS;
using namespace qchem;

using namespace qchem::BasisSet::Radial;
using enum Type;
using enum BasisSetAccuracy;

class DBCach1Tests : public ::testing::Test
{
public:
    DBCach1Tests() 
        : N(3),Z(1)
        , cl_hydrogen    (new Atom(1,0.0,Vector3D(0,0,0)))
        , cl_hydrogen_100(new Atom(1,0.0,Vector3D(1,0,0)))
        , cl_helium      (new Atom(2,0.0,Vector3D(0,0,0)))
    {
    }
    ~DBCach1Tests()
    {
        delete cl_hydrogen;
        delete cl_hydrogen_100;
        delete cl_helium;
        delete bs1;
        delete bs2;
    }
    void Init(BasisSet::Radial::Type type)
    {
        bs1=Factory(Low,type,Z);
        bs2=Factory(Low,type,Z);
    }
  
    // THE BASIS'S OWN TABLES ARE THE BASIS'S (doc/IntegralCacheOwnershipPlan.md S3a): a block asked twice returns the SAME stored object
    // (cached in its parent's cache), and a second, independently built basis computes the EQUAL table into its OWN cache -- no sharing
    // across owners (the user: sharing between two equal, independently built bases "should be very low priority").
    template <class M> static void SameOwnerSharesOtherOwnerRecomputes(const M& a1, const M& a1again, const M& a2)
    {
        EXPECT_EQ(&a1,&a1again);   // one owner: one stored object
        EXPECT_EQ(a1,a2);          // another owner: the same numbers...
        EXPECT_NE(&a1,&a2);        // ...in its own cache
    }
    void TestOverlap() const
    {
        using OIBS=BasisSet::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            SameOwnerSharesOtherOwnerRecomputes(ibs1->Overlap(), ibs1->Overlap(), (*ibs2)->Overlap());
            ++ibs2;
        }
    }
    void TestKinetic() const
    {
        using OIBS=BasisSet::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            // Kinetic() here is just a representative cached 1e matrix (the <p^2> block); its physics is incidental.
            SameOwnerSharesOtherOwnerRecomputes(ibs1->Kinetic(), ibs1->Kinetic(), (*ibs2)->Kinetic());
            ++ibs2;
        }
    }
    // Nuclear is still in the process-wide cache (until S4 moves it to the Ven term): equal keys share one object.
    void TestNuclear() const
    {
        using OIBS=BasisSet::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Nuclear(cl_hydrogen);
            auto& S2=(*ibs2)->Nuclear(cl_hydrogen);
            auto& S3=(*ibs2)->Nuclear(cl_hydrogen_100);
            auto& S4=(*ibs2)->Nuclear(cl_helium);
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(S1,S3);
            EXPECT_NE(S1,S4);
            EXPECT_EQ(&S1,&S2);
            EXPECT_NE(&S1,&S3);
            EXPECT_NE(&S1,&S4);
            ++ibs2;
        }
    }
    // ⟨ab|c⟩ against ONE fit basis: keyed by (the block, the fit's instance serial) in the block's parent's cache.
    void TestOverlap3C(double eps) const
    {
        using BasisSet::Real_DFT_OIBS;
        auto ibs2=bs2->Iterate<Real_DFT_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<Real_DFT_OIBS>())
        {
            auto ff=ibs1->CreateVxcFitBasisSet(cl_hydrogen,qcMesh::MeshParams{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=2, .mhl_alpha=2.0, .angular=qcMesh::AngularKind::Lebedev, .angularDegree=3, .beckeOrder=3});
            const PairProjection<double>& E1=ibs1->Overlap3C(*ff);
            const PairProjection<double>& E2=(*ibs2)->Overlap3C(*ff);
            EXPECT_EQ(&E1,&ibs1->Overlap3C(*ff));
            EXPECT_EQ(E1.dense,E2.dense);
            EXPECT_NE(&E1,&E2);
            ++ibs2;
            delete ff;
        }
    }
    void TestRepulsion3C(double eps) const
    {
        using BasisSet::Real_DFT_OIBS;
        auto ibs2=bs2->Iterate<Real_DFT_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<Real_DFT_OIBS>())
        {
            auto ff=ibs1->CreateCDFitBasisSet(cl_hydrogen);
            const PairProjection<double>& E1=ibs1->Repulsion3C(*ff);
            const PairProjection<double>& E2=(*ibs2)->Repulsion3C(*ff);
            EXPECT_EQ(&E1,&ibs1->Repulsion3C(*ff));
            EXPECT_EQ(E1.dense,E2.dense);
            EXPECT_NE(&E1,&E2);
            ++ibs2;
            delete ff;
        }
    }
    // The ERI4 store is canonical-only (doc/ERI4Rework.md §5.2): only the a<=b block of an unordered pair is built and stored, a by the
    // block's index in its parent's cache, i.e. the parent's insertion (= iteration) order.  So (1) fetch only the canonical block and check
    // it is stored once per owner, (2) verify the bra-ket symmetry J(a,b)=J(b,a)^T against the UNCACHED MakeDirect of the partner, (3) assert
    // the non-canonical request is forbidden, and (4) a pair across two parents is refused (it has no single owner).
    void TestDirect(double eps) const
    {
        cout << *bs1 << endl;
        using BasisSet::Real_ERI4_OIBS;
        auto ibs21=bs2->Iterate<Real_ERI4_OIBS>().begin();
        size_t i=0;
        for (auto ibs11:bs1->Iterate<Real_ERI4_OIBS>())
        {
            auto ibs22=bs2->Iterate<Real_ERI4_OIBS>().begin();
            size_t j=0;
            for (auto ibs12:bs1->Iterate<Real_ERI4_OIBS>())
            {
                if (i<=j)
                {
                    const ERI4& J1=ibs11->Direct(*ibs12);
                    SameOwnerSharesOtherOwnerRecomputes(J1, ibs11->Direct(*ibs12), (*ibs21)->Direct(**ibs22));
                    ERI4 J1ba=ibs12->MakeDirect(*ibs11);          // partner built UNCACHED (the store forbids it)
                    EXPECT_NEAR(fnorm(J1,J1ba.Transpose()),0.0,eps);
                }
                else
                    EXPECT_ANY_THROW(ibs11->Direct(*ibs12));      // non-canonical request is forbidden
                EXPECT_ANY_THROW(ibs11->Direct(**ibs22));         // two parents: refused
                ++ibs22; ++j;
            }
            ++ibs21; ++i;
        }
    }
    void TestExchange(double eps) const
    {
        using BasisSet::Real_ERI4_OIBS;
        auto ibs21=bs2->Iterate<Real_ERI4_OIBS>().begin();
        size_t i=0;
        for (auto ibs11:bs1->Iterate<Real_ERI4_OIBS>())
        {
            auto ibs22=bs2->Iterate<Real_ERI4_OIBS>().begin();
            size_t j=0;
            for (auto ibs12:bs1->Iterate<Real_ERI4_OIBS>())
            {
                if (i<=j)
                {
                    const ERI4& K1=ibs11->Exchange(*ibs12);
                    SameOwnerSharesOtherOwnerRecomputes(K1, ibs11->Exchange(*ibs12), (*ibs21)->Exchange(**ibs22));
                    ERI4 K1ba=ibs12->MakeExchange(*ibs11);            // partner built UNCACHED
                    EXPECT_NEAR(fnorm(K1,K1ba.Transpose()),0.0,eps);
                }
                else
                    EXPECT_ANY_THROW(ibs11->Exchange(*ibs12));        // non-canonical request is forbidden
                ++ibs22; ++j;
            }
            ++ibs21; ++i;
        }
    }
    size_t N,Z;
    Structure *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    BasisSet::Real_BS *bs1,*bs2;
};




TEST_F(DBCach1Tests,BSplineOverlap)
{
    Init(BSpline6);
    TestOverlap();
}
TEST_F(DBCach1Tests,GaussianOverlap)
{
    Init(Gaussian);
    TestOverlap();
}
TEST_F(DBCach1Tests,SlaterOverlap)
{
    Init(Slater);
    TestOverlap();
}


TEST_F(DBCach1Tests,BSplineKinetic)
{
    Init(BSpline6);
    TestKinetic();
}
TEST_F(DBCach1Tests,GaussianKinetic)
{
    Init(Gaussian);
    TestKinetic();
}
TEST_F(DBCach1Tests,SlaterKinetic)
{
    Init(Slater);
    TestKinetic();
}
TEST_F(DBCach1Tests,BSplineNuclear)
{
    Init(BSpline6);
    TestNuclear();
}
TEST_F(DBCach1Tests,GaussianNuclear)
{
    Init(Gaussian);
    TestNuclear();
}
TEST_F(DBCach1Tests,SlaterNuclear)
{
    Init(Slater);
    TestNuclear();
}
TEST_F(DBCach1Tests,GaussianOverlap3C)
{
    Init(Gaussian);
    TestOverlap3C(1e-14);
}
TEST_F(DBCach1Tests,SlaterOverlap3C)
{
    Init(Slater);
    TestOverlap3C(1e-14);
}
TEST_F(DBCach1Tests,GaussianRepulsion3C)
{
    Init(Gaussian);
    TestRepulsion3C(1e-14);
}
TEST_F(DBCach1Tests,SlaterRepulsion3C)
{
    Init(Slater);
    TestRepulsion3C(1e-14);
}


//
//  For efficiency legacy calculations do Jcdab=Jabcd.Transpose().   
//  The BasisSet version recalculates Jcdab directly and does not yet do
//  the Transpose call.  For Guassian and Slater this results in a minor error <1e-15.  For BSpline near the origin
//  we get much bigger errors.  Perhaps this related to problems encountered when trying to reproduce Saitos results with
//  rmin=0.0001.
//
TEST_F(DBCach1Tests,BSplineDirect)
{
    Init(BSpline6);
    // ⚠ 6e-13 -> 8e-13 (2026-09-04): the residual is ||J - J^T||, mathematically EXACTLY zero, so the
    // number here is a roundoff floor and not a physics bar.  It was fitted just above what one -O3 build
    // produced (1.05x headroom), and `-march=native`'s FMA contraction moved it to 6.30e-13.  Widened to
    // keep ~1.3x headroom over the contracted value; a real bra-ket asymmetry is O(1e-3), not O(1e-13).
    TestDirect(8e-13);
}
TEST_F(DBCach1Tests,BSplinerDirect)
{
    Init(BSpliner6);
    TestDirect(1.2e-13);
}
TEST_F(DBCach1Tests,GaussianDirect)
{
    Init(Gaussian);
    TestDirect(8e-14);
}
TEST_F(DBCach1Tests,SlaterDirect)
{
    Init(Slater);
    TestDirect(4e-15);
}
TEST_F(DBCach1Tests,BSplineExchange)
{
    Init(BSpline6);
    TestExchange(2e-13);   // 1.4e-13 -> 2e-13, same roundoff-floor reason as TestDirect above (FMA: 1.63e-13)
}
TEST_F(DBCach1Tests,BSplinerExchange)
{
    Init(BSpliner6);
    TestExchange(6e-14);
}
TEST_F(DBCach1Tests,GaussianExchange)
{
    Init(Gaussian);
    TestExchange(6e-14);
}
TEST_F(DBCach1Tests,SlaterExchange)
{
    Init(Slater);
    TestExchange(4e-15);
}
