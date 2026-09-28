// File: SymQN.C  Unite tests for symmetrys and QN classes

#include "gtest/gtest.h"
#include <set>
#include <iostream>
import qchem.Symmetry.Orbital; 
import qchem.Streamable;
import qchem.Symmetry.Factory;
import qchem.Symmetry.Atom.Spherical;
import qchem.Blaze;
import qchem.Symmetry.Unit;                 // UnitQN (a C1 block, for the Invariant selection rule)
import qchem.Symmetry.Lattice_3D.BlochQN;   // MeshShift / IsShiftOf / CommensurateShifts (LinearResponsePlan S1)
using namespace qchem;

using std::cout; 
using std::endl;
using namespace qchem::Symmetry;
class SymQNTests : public ::testing::Test
{
public:
    SymQNTests() 
    : LMax(4) 
    , κ_max(4)
    , n_max(20)
    , quiet(true)
    {};

    sym_t Y(int l) const {return YFactory(l);}
    sym_t Y(int l, const ivec_t& mls) const {return YFactory(l,mls);}
    sym_t Ω(int   κ,rvec_t mjs={}) const {return ΩFactory(κ,mjs);}

    static Spin makespin(int ms)
    {
        if (ms==0)
            return Spin::Down;
        else if (ms==1)
            return Spin::None;
        return Spin::Up;
    }
    static ivec_t make_mls(int ml0, int ml1)
    {
        assert(ml1>=ml0);
        size_t N=ml1-ml0+1;
        return blazem::linspace(N,ml0,ml1);
    }
    static rvec_t make_mjs(double mj0, double mj1)
    {
        assert(mj1>=mj0);
        size_t N=mj1-mj0+1;
        return blazem::linspace(N,mj0,mj1);
    }
    size_t LMax;
    int κ_max;
    int n_max;
    bool quiet;
};

TEST_F(SymQNTests, Yl_SequenceIndex)
{
    for (size_t l1=0;l1<=LMax;l1++)
    {
        sym_t yl1=Y(l1);
        for (size_t l2=0;l2<l1;l2++)
        {
            // if (!quiet) cout << "{l1,l2}={" << l1 << "," << l2 << "}" << endl;
            sym_t yl2=Y(l2);
            EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
        }
    }
}
TEST_F(SymQNTests, Ylm_SequenceIndex)
{
    for (size_t l1=0;l1<=LMax;l1++)
    for (int m1=-(int)l1;m1<=(int)l1;m1++)
    {
        sym_t yl1=Y(l1,make_mls(-(int)l1,m1));
        if (!quiet) cout << "{l1,m1,sn}={" << l1 << "," << m1 << "," << yl1->SequenceIndex() << "}" << endl;
        for (size_t l2=0;l2<=l1;l2++)
        for (int m2=-(int)l2;m2<=(int)l2;m2++)
        {
            // if (!quiet) cout << "{l1,l2,m1,m2}={" << l1 << "," << l2 << "," << m1 << "," << m2 << "}" << endl;
            sym_t yl2=Y(l2,make_mls(-(int)l2,m2));
            if (l1!=l2 || m1!=m2)
            {
                EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
            }
        }
    }
}
TEST_F(SymQNTests, Yl_Ylm_CrossSequenceIndex)
{
    for (size_t l1=0;l1<=LMax;l1++)
    {
        sym_t yl1=Y(l1);
        // if (!quiet) cout << "{l1,m1,sn}={" << l1 << "," << m1 << "," << yl1->SequenceIndex() << "}" << endl;
        for (size_t l2=0;l2<=l1;l2++)
        for (int m2=-(int)l2;m2<=(int)l2;m2++)
        {
            // if (!quiet) cout << "{l1,l2,m1,m2}={" << l1 << "," << l2 << "," << m1 << "," << m2 << "}" << endl;
            sym_t yl2=Y(l2,make_mls(-(int)l2,m2));
            EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
        }
    }
}

TEST_F(SymQNTests, Ωκ_SequenceIndex)
{
    for (int κ1=-κ_max;κ1<κ_max;κ1++) //Leave out the uppermost κ, it corresponds to LMAX+1
    {
        sym_t Ol1=Ω(κ1);
        if (!quiet) cout << "{k1,sn}={" << κ1 << "," << Ol1->SequenceIndex() << "}" << endl;
        for (int κ2=-κ_max;κ2<κ1;κ2++)
        {
            sym_t Ol2=Ω(κ2);
            EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
        }
    }
}
TEST_F(SymQNTests, Ωκmj_SequenceIndex)
{
    for (int κ1=-κ_max;κ1<κ_max;κ1++)
    {
        double j1=::qchem::Symmetry::Atom::SphericalSpinor::j(κ1);
        for (double mj1=-j1;mj1<=j1;mj1++)
        {
            sym_t Ol1=Ω(κ1,make_mjs(-j1,mj1));
            if (!quiet) cout << "{k1,mj1,sn}={" << κ1 << "," << mj1 << "," << Ol1->SequenceIndex() << "}" << endl;
            for (int κ2=-κ_max;κ2<κ_max;κ2++)
            {
                double j2=::qchem::Symmetry::Atom::SphericalSpinor::j(κ2);
                for (double mj2=-j2;mj2<=j2;mj2++)   
                {
                    sym_t Ol2=Ω(κ2,make_mjs(-j2,mj2));
                    // if (!quiet) cout << "{k1,k2,mj1,mj2,sn}={" << κ1 << "," << κ2 << "," << mj1 << "," << mj2 << "," << Ol2->SequenceIndex() << "}" << endl;
                    if (κ1!=κ2 || mj1!=mj2)
                    {
                        EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
                    }
                }
            }
        }
    }
}
TEST_F(SymQNTests, Omega_k_kmj_CrossSequenceIndex)
{
    for (int κ1=-κ_max;κ1<κ_max;κ1++)
    {
        {
            sym_t Ol1=Ω(κ1);
            if (!quiet) cout << "{k1,sn}={" << κ1 << "," << Ol1->SequenceIndex() << "}" << endl;
            for (int κ2=-κ_max;κ2<κ_max;κ2++)
            {
                double j2=::qchem::Symmetry::Atom::SphericalSpinor::j(κ2);
                for (double mj2=-j2;mj2<=j2;mj2++)   
                {
                    sym_t Ol2=Ω(κ2,make_mjs(-j2,mj2));
                    // if (!quiet) cout << "{k1,k2,mj1,mj2,sn}={" << κ1 << "," << κ2 << "," << mj1 << "," << mj2 << "," << Ol2->SequenceIndex() << "}" << endl;
                    EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
                }
            }
        }
    }
}
TEST_F(SymQNTests, Orbital_QNs_Yl_SequenceIndex)
{
    for (int n1=1;n1<=n_max;n1++)
    for (int n2=1;n2<=n1;n2++)
    for (int ms1=0;ms1<=2;ms1++)
    for (int ms2=0;ms2<=2;ms2++)
    for (size_t l1=0;l1<=LMax;l1++)
    {
        Spin s1=makespin(ms1);
        Spin s2=makespin(ms2);
        auto yl1=Y(l1);
        Orbital_QNs oqn1(n1,s1,yl1);
        for (size_t l2=0;l2<l1;l2++)
        {
            // if (!quiet) cout << "{l1,l2}={" << l1 << "," << l2 << "}" << endl;
            auto yl2=Y(l2);
            Orbital_QNs oqn2(n2,s2,yl2);
            EXPECT_NE(oqn1.SequenceIndex(),oqn2.SequenceIndex());
        }
    }
}
TEST_F(SymQNTests, Orbital_QNs_set)
{
    std::set<Orbital_QNs> qns;

    for (int n1=1;n1<=n_max;n1++)
    for (int ms1=0;ms1<=2;ms1++)
    for (size_t l1=0;l1<=LMax;l1++)
    {
        Spin s1=makespin(ms1);
        auto yl1=Y(l1);
        qns.insert(Orbital_QNs(n1,s1,yl1));
        // delete yl1;
    }
    for (auto qn:qns) cout << qn << " ";
    cout << endl;
}


// doc/RealComplexPlan.md Step 1: the basis-type half of the realness rule.  Ordinary spatial symmetries
// default true; BlochQN answers the TRIM question by EXACT integer arithmetic (no float-k tolerance).
TEST_F(SymQNTests, IsReal)
{
    EXPECT_TRUE(Y(0)->IsReal());                     // atomic shells: real by default
    EXPECT_TRUE(Y(2,make_mls(-2,2))->IsReal());

    ivec3_t N(4,4,4);
    EXPECT_TRUE (BlochFactory(N,{0,0,0})->IsReal());     // Gamma
    EXPECT_TRUE (BlochFactory(N,{2,2,0})->IsReal());     // zone boundary k=(1/2,1/2,0)
    EXPECT_TRUE (BlochFactory(N,{-2,2,2})->IsReal());    // negative index, still TRIM
    EXPECT_FALSE(BlochFactory(N,{1,0,0})->IsReal());     // k=(1/4,0,0): not TRIM
    EXPECT_FALSE(BlochFactory(N,{2,2,3})->IsReal());     // one non-TRIM component poisons the block

    // Gamma-centred 2x2x2 is TRIM THROUGHOUT; the MP shift=1/2 mesh (k=+/-1/4) never is.
    ivec3_t N2(2,2,2);
    for (int x=0;x<2;x++) for (int y=0;y<2;y++) for (int z=0;z<2;z++)
    {
        EXPECT_TRUE (BlochFactory(N2,{x,y,z})->IsReal());
        EXPECT_FALSE(BlochFactory(N2,{x,y,z},1.0,{0.5,0.5,0.5})->IsReal());
    }
    EXPECT_TRUE(BlochFactory({1,1,1},{0,0,0},1.0,{0.5,0.5,0.5})->IsReal());   // N=1 shifted: k=(1/2,1/2,1/2) IS TRIM
}

TEST_F(SymQNTests, BlochQNs)
{
    ivec3_t N(5,6,7);
    ivec3_t k1;
    
    for (k1.x=-N.x;k1.x<=N.x;k1.x++)
    for (k1.y=-N.y;k1.y<=N.y;k1.y++)
    for (k1.z=-N.z;k1.z<=N.z;k1.z++)
    {
        auto bq1=BlochFactory(N,k1);
        // cout << k1 << " " << bq1 << " " << bq1.SequenceIndex() << endl;
        ivec3_t k2;
        for (k2.x=k1.x;k2.x<=N.x;k2.x++)
        for (k2.y=k1.y;k2.y<=N.y;k2.y++)
        for (k2.z=k1.z;k2.z<=N.z;k2.z++)
        {
            if (k1==k2) continue;
            auto bq2=BlochFactory(N,k2);
            EXPECT_NE(bq1->SequenceIndex(),bq2->SequenceIndex());
        }
    }
}
// ------------------------------------------------------------------ MeshShift (doc/LinearResponsePlan.md S1)
// The k -> k+q pairing a monochromatic linear response needs.  Pinned: every k on the mesh has EXACTLY ONE
// partner for every commensurate q (so the pairing is a permutation), q=0 pairs a point with itself, an
// incommensurate q-mesh is a FAILED Outcome (never a throw), and a shifted Monkhorst-Pack mesh pairs within
// itself -- q is a difference of mesh points, so k+q stays on the SHIFTED mesh.
namespace {
std::vector<sym_t> Mesh(ivec3_t N, rvec3_t shift={0,0,0})
{
    std::vector<sym_t> out;
    for (int x=0;x<N.x;x++) for (int y=0;y<N.y;y++) for (int z=0;z<N.z;z++)
        out.push_back(BlochFactory(N,ivec3_t(x,y,z),1.0/(N.x*N.y*N.z),shift));
    return out;
}
}

TEST_F(SymQNTests, MeshShift_PartnerIsAPermutation)
{
    using Lattice_3D::IsShiftOf;
    for (rvec3_t shift : {rvec3_t(0,0,0), rvec3_t(0.5,0.5,0.5)})
    {
        const ivec3_t N(4,2,6);
        auto mesh=Mesh(N,shift);
        auto qs=Lattice_3D::CommensurateShifts(*mesh[0], ivec3_t(2,2,3));
        ASSERT_TRUE(qs.IsOk());
        EXPECT_EQ(qs->size(), 12u);
        for (const auto& q : *qs)
        {
            std::set<size_t> hit;
            for (size_t b=0;b<mesh.size();b++)
            {
                size_t n=0, p=0;
                for (size_t c=0;c<mesh.size();c++) if (IsShiftOf(*mesh[c],*mesh[b],q)) {n++; p=c;}
                EXPECT_EQ(n,1u) << "q=" << q;
                hit.insert(p);
            }
            EXPECT_EQ(hit.size(), mesh.size()) << "q=" << q << ": the pairing is not a permutation";
        }
    }
}

TEST_F(SymQNTests, MeshShift_ZeroPairsWithItselfAndNegativeIndicesWork)
{
    const ivec3_t N(3,3,3);
    auto qs=Lattice_3D::CommensurateShifts(*BlochFactory(N,ivec3_t(0,0,0)), N);
    ASSERT_TRUE(qs.IsOk());
    EXPECT_TRUE((*qs)[0].IsZero());
    auto k=BlochFactory(N,ivec3_t(-1,2,0)), kk=BlochFactory(N,ivec3_t(2,-1,0));   // -1 == 2 (mod 3)
    EXPECT_TRUE (Lattice_3D::IsShiftOf(*k,*kk,(*qs)[0]));
    EXPECT_TRUE (Lattice_3D::IsShiftOf(*k,*k,(*qs)[0]));
    // q = (1,0,0) steps: (-1,2,0)+(1,0,0) = (0,2,0)
    auto q=(*qs)[9];                          // ix=1, iy=0, iz=0
    EXPECT_EQ(q.Steps(), ivec3_t(1,0,0));
    EXPECT_TRUE (Lattice_3D::IsShiftOf(*BlochFactory(N,ivec3_t(0,2,0)),*k,q));
    EXPECT_FALSE(Lattice_3D::IsShiftOf(*BlochFactory(N,ivec3_t(1,2,0)),*k,q));
}

TEST_F(SymQNTests, MeshShift_IncommensurateQMeshFails)
{
    auto qs=Lattice_3D::CommensurateShifts(*BlochFactory(ivec3_t(4,4,4),ivec3_t(0,0,0)), ivec3_t(3,2,2));
    EXPECT_FALSE(qs.IsOk());
    EXPECT_NE(qs.Error().find("not commensurate"), std::string::npos);
}

TEST_F(SymQNTests, MeshShift_DifferentMeshesNeverPair)
{
    auto a=BlochFactory(ivec3_t(4,4,4),ivec3_t(1,0,0));
    auto b=BlochFactory(ivec3_t(4,4,4),ivec3_t(1,0,0),1.0/64,rvec3_t(0.5,0.5,0.5));   // shifted MP twin
    auto qs=Lattice_3D::CommensurateShifts(*a, ivec3_t(1,1,1));
    ASSERT_TRUE(qs.IsOk());
    EXPECT_FALSE(Lattice_3D::IsShiftOf(*a,*b,(*qs)[0]));
}

//! S1 (doc/LinearResponsePlan.md §3c): a MeshShift seen through the SelectionRule face IS IsShiftOf, and the
//! Invariant rule (q = 0 / a totally symmetric perturbation) pairs every block with itself -- the same
//! pairing as the zero MeshShift, but for any symmetry (a molecule has no Bloch points).
TEST_F(SymQNTests, SelectionRule_MeshShiftAndInvariant)
{
    const ivec3_t N(2,3,2);
    auto mesh=Mesh(N,rvec3_t(0,0,0));
    auto qs=Lattice_3D::CommensurateShifts(*mesh[0], N);
    ASSERT_TRUE(qs.IsOk());
    const Invariant inv;
    for (const auto& q : *qs)
    {
        const SelectionRule& rule=q;
        for (size_t b=0;b<mesh.size();b++)
            for (size_t c=0;c<mesh.size();c++)
            {
                EXPECT_EQ(rule.Couples(*mesh[c],*mesh[b]), Lattice_3D::IsShiftOf(*mesh[c],*mesh[b],q));
                if (q.IsZero()) EXPECT_EQ(rule.Couples(*mesh[c],*mesh[b]), inv.Couples(*mesh[c],*mesh[b]));
            }
    }
    const UnitQN c1;                          // a C1 molecule's one block
    EXPECT_TRUE(inv.Couples(c1,c1));
}
