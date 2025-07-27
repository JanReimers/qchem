// File: SymQN.C  Unite tests for symmetrys and QN classes

#include "gtest/gtest.h"
#include <set>
#include <iostream>
import qchem.Symmetry.Orbital;
import qchem.Symmetry.Yl;
import qchem.Symmetry.Ylm;
import qchem.Symmetry.Okmj;
import qchem.Symmetry.BlochQN;
import qchem.Streamable;

using std::cout;
using std::endl;

class SymQNTests : public ::testing::Test
{
public:
    typedef Irrep_QNs::sym_t sym_t;
    SymQNTests() 
    : LMax(4) 
    , kappa_max(4)
    , n_max(20)
    , quiet(true)
    {};

    static Spin makespin(int ms)
    {
        if (ms==0)
            return Spin::Down;
        else if (ms==1)
            return Spin::None;
        return Spin::Up;
    }
    static std::vector<int> make_mls(int m0, int m1)
    {
        std::vector<int> ms;
        for (int m=m0;m<=m1;m++) ms.push_back(m);
        return ms;
    }
    size_t LMax;
    int kappa_max;
    int n_max;
    bool quiet;
};

TEST_F(SymQNTests, Yl_Sym)
{
    for (size_t l1=0;l1<=LMax;l1++)
    {
        Symmetry* yl1=new Yl_Sym(l1);
        for (size_t l2=0;l2<l1;l2++)
        {
            if (!quiet) cout << "{l1,l2}={" << l1 << "," << l2 << "}" << endl;
            Symmetry* yl2=new Yl_Sym(l2);
            EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
            delete yl2;
        }
        delete yl1;
    }
}
TEST_F(SymQNTests, Omega_k_Sym)
{
    for (int kappa1=-kappa_max;kappa1<=kappa_max;kappa1++)
    {
        Symmetry* Ol1=new Omega_k_Sym(kappa1);
        for (int kappa2=-kappa_max;kappa2<kappa1;kappa2++)
        {
            Symmetry* Ol2=new Omega_k_Sym(kappa2);
            EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
            delete Ol2;
        }
        delete Ol1;
    }
}
TEST_F(SymQNTests, Ylm_Sym_multi)
{
    for (size_t l1=0;l1<=LMax;l1++)
    for (int m1=-(int)l1;m1<=(int)l1;m1++)
    {
        Symmetry* yl1=new Ylm_Sym(l1,make_mls(-(int)l1,m1));
        for (size_t l2=0;l2<=l1;l2++)
        for (int m2=-(int)l2;m2<=(int)l2;m2++)
        {
            if (!quiet) cout << "{l1,l2,m1,m2}={" << l1 << "," << l2 << "," << m1 << "," << m2 << "}" << endl;
            Symmetry* yl2=new Ylm_Sym(l2,make_mls(-(int)l2,m2));
            if (m1!=m2)
            {
                EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
            }
            delete yl2;
        }
        delete yl1;
    }
}
TEST_F(SymQNTests, Omega_kmj_Sym)
{
    for (int kappa1=-kappa_max;kappa1<=kappa_max;kappa1++)
    {
        double j1=Omega_k_Sym::j(kappa1);
        for (double mj1=-j1;mj1<=j1;mj1++)
        {
            Symmetry* Ol1=new Omega_kmj_Sym(kappa1,mj1);
            for (int kappa2=-kappa_max;kappa2<=kappa1;kappa2++)
            {
                double j2=Omega_k_Sym::j(kappa2);
                for (double mj2=-j2;mj2<=j2;mj2++)   
                {
                    Symmetry* Ol2=new Omega_kmj_Sym(kappa2,mj2);
                    if (!quiet) cout << "{k1,k2,mj1,mj2,sn}={" << kappa1 << "," << kappa2 << "," << mj1 << "," << mj2 << "," << Ol2->SequenceIndex() << "}" << endl;
                    if (mj1!=mj2)
                    {
                        EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
                    }
                    delete Ol2;
                }
            }
            delete Ol1;
        }
    }
}
TEST_F(SymQNTests, Orbital_QNs_YlQN)
{
    for (int n1=1;n1<=n_max;n1++)
    for (int n2=1;n2<=n1;n2++)
    for (int ms1=0;ms1<=2;ms1++)
    for (int ms2=0;ms2<=2;ms2++)
    for (size_t l1=0;l1<=LMax;l1++)
    {
        Spin s1=makespin(ms1);
        Spin s2=makespin(ms2);
        auto yl1=sym_t(new Yl_Sym(l1));
        Orbital_QNs oqn1(n1,s1,yl1);
        for (size_t l2=0;l2<l1;l2++)
        {
            if (!quiet) cout << "{l1,l2}={" << l1 << "," << l2 << "}" << endl;
            auto yl2=sym_t(new Yl_Sym(l2));
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
        auto yl1=sym_t(new Yl_Sym(l1));
        qns.insert(Orbital_QNs(n1,s1,yl1));
        // delete yl1;
    }
    for (auto qn:qns) cout << qn << " ";
    cout << endl;
}


TEST_F(SymQNTests, BlochQNs)
{
    IVec3 N(5,6,7);
    IVec3 k1;
    
    for (k1.x=-N.x;k1.x<=N.x;k1.x++)
    for (k1.y=-N.y;k1.y<=N.y;k1.y++)
    for (k1.z=-N.z;k1.z<=N.z;k1.z++)
    {
        BlochQN bq1(N,k1);
        // cout << k1 << " " << bq1 << " " << bq1.SequenceIndex() << endl;
        IVec3 k2;
        for (k2.x=k1.x;k2.x<=N.x;k2.x++)
        for (k2.y=k1.y;k2.y<=N.y;k2.y++)
        for (k2.z=k1.z;k2.z<=N.z;k2.z++)
        {
            if (k1==k2) continue;
            BlochQN bq2(N,k2);
            EXPECT_NE(bq1.SequenceIndex(),bq2.SequenceIndex());
        }
    }
}