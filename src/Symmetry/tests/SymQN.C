// File: SymQN.C  Unite tests for symmetrys and QN classes

#include "gtest/gtest.h"
#include <set>
#include <iostream>
import qchem.Symmetry.Orbital; 
import qchem.Streamable;
import qchem.Symmetry.Factory;
import qchem.Symmetry.Spherical;
import qchem.Blaze;

using std::cout; 
using std::endl;
using namespace Symmetry;
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
        double j1=::Symmetry::SphericalSpinor::j(κ1);
        for (double mj1=-j1;mj1<=j1;mj1++)
        {
            sym_t Ol1=Ω(κ1,make_mjs(-j1,mj1));
            if (!quiet) cout << "{k1,mj1,sn}={" << κ1 << "," << mj1 << "," << Ol1->SequenceIndex() << "}" << endl;
            for (int κ2=-κ_max;κ2<κ_max;κ2++)
            {
                double j2=::Symmetry::SphericalSpinor::j(κ2);
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
                double j2=::Symmetry::SphericalSpinor::j(κ2);
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