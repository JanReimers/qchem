// File: SymQN.C  Unite tests for symmetrys and QN classes

#include "gtest/gtest.h"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Symmetry/YlmQN.H"
#include "Imp/Symmetry/OkmjQN.H"
#include <Orbital_QNs.H>
#include <set>
#include <iostream>

using std::cout;
using std::endl;

class SymQNTests : public ::testing::Test
{
public:
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
    int LMax;
    int kappa_max;
    int n_max;
    bool quiet;
};

TEST_F(SymQNTests, YlQN)
{
    for (size_t l1=0;l1<=LMax;l1++)
    {
        QNs* yl1=new YlQN(l1);
        for (size_t l2=0;l2<l1;l2++)
        {
            if (!quiet) cout << "{l1,l2}={" << l1 << "," << l2 << "}" << endl;
            QNs* yl2=new YlQN(l2);
            EXPECT_FALSE(*yl1==*yl2);
            EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
            delete yl2;
        }
        delete yl1;
    }
}
TEST_F(SymQNTests, Omega_kQN)
{
    for (int kappa1=-kappa_max;kappa1<=kappa_max;kappa1++)
    {
        QNs* Ol1=new Omega_kQN(kappa1);
        for (size_t kappa2=-kappa_max;kappa2<kappa1;kappa2++)
        {
            QNs* Ol2=new Omega_kQN(kappa2);
            EXPECT_FALSE(*Ol1==*Ol2);
            EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
            delete Ol2;
        }
        delete Ol1;
    }
}
TEST_F(SymQNTests, YlmQN)
{
    for (int l1=0;l1<=LMax;l1++)
    for (int m1=-l1;m1<=l1;m1++)
    {
        QNs* yl1=new YlmQN(l1,m1);
        for (int l2=0;l2<=l1;l2++)
        for (int m2=-l2;m2<=l2;m2++)
        {
            if (!quiet) cout << "{l1,l2,m1,m2}={" << l1 << "," << l2 << "," << m1 << "," << m2 << "}" << endl;
            QNs* yl2=new YlmQN(l2,m2);
            if (!(*yl1==*yl2))
                EXPECT_NE(yl1->SequenceIndex(),yl2->SequenceIndex());
            delete yl2;
        }
        delete yl1;
    }
}
TEST_F(SymQNTests, Omega_kmjQN)
{
    for (int kappa1=-kappa_max;kappa1<=kappa_max;kappa1++)
    {
        double j1=Omega_kQN::j(kappa1);
        for (double mj1=-j1;mj1<=j1;mj1++)
        {
            QNs* Ol1=new Omega_kmjQN(kappa1,mj1);
            for (int kappa2=-kappa_max;kappa2<=kappa1;kappa2++)
            {
                double j2=Omega_kQN::j(kappa2);
                for (double mj2=-j2;mj2<=j2;mj2++)   
                {
                    QNs* Ol2=new Omega_kmjQN(kappa2,mj2);
                    if (!quiet) cout << "{k1,k2,mj1,mj2,sn}={" << kappa1 << "," << kappa2 << "," << mj1 << "," << mj2 << "," << Ol2->SequenceIndex() << "}" << endl;
                    if (!(*Ol1==*Ol2))
                        EXPECT_NE(Ol1->SequenceIndex(),Ol2->SequenceIndex());
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
        QNs* yl1=new YlQN(l1);
        Orbital_QNs oqn1(n1,s1,yl1);
        for (size_t l2=0;l2<l1;l2++)
        {
            if (!quiet) cout << "{l1,l2}={" << l1 << "," << l2 << "}" << endl;
            QNs* yl2=new YlQN(l2);
            Orbital_QNs oqn2(n2,s2,yl2);
            EXPECT_FALSE(oqn1==oqn2);
            EXPECT_NE(oqn1.SequenceIndex(),oqn2.SequenceIndex());
            delete yl2;
        }
        delete yl1;
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
        QNs* yl1=new YlQN(l1);
        qns.insert(Orbital_QNs(n1,s1,yl1));
        // delete yl1;
    }

    for (auto qn:qns) cout << qn << endl;
}
