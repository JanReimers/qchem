// File: UnitTests/BasisSet_Atom1.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <iostream>
// #include <cmath>
#include <blaze/Math.h>
using std::cout;
using std::endl;

import qchem.BasisSet.DB_Cache1;
import qchem.BasisSet.Atom.BSpline.NR.BS1;
import BasisSet.Atom.BSpline.NR.BS_Evaluator;
import qchem.Symmetry.Yl;

bool operator==(const ERI4& a, const ERI4& b); //Defined in UnitTests/BasisSet_Atom.C


class DBCach1Tests : public ::testing::Test
{
public:
    DBCach1Tests() 
        : cl_hydrogen    (new Atom(1,0.0,Vector3D(0,0,0)))
        , cl_hydrogen_100(new Atom(1,0.0,Vector3D(1,0,0)))
        , cl_helium      (new Atom(2,0.0,Vector3D(0,0,0)))
        , yl(new Yl_Sym(0))
        ,  bs1(new AtomBS::BSpline1::BasisSet<6,BSpline_r_BS>(3,0.1,10.0,Atom_EC(86)))
        ,  bs2(new AtomBS::BSpline1::BasisSet<6,BSpline_r_BS>(3,0.1,10.0,Atom_EC(86)))
    {
        theGlobalCache=new IntegralsCache_RAM<double>();
    }
    ~DBCach1Tests()
    {
        delete cl_hydrogen;
        delete cl_hydrogen_100;
        delete cl_helium;
        delete bs1;
        delete bs2;
        delete theGlobalCache;
    }
    
    Cluster *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    Irrep_QNs::sym_t yl;
    ::BasisSet1 *bs1,*bs2;
};


TEST_F(DBCach1Tests,BSplineOverlap)
{
    auto ibs2=bs2->Iterate<Real_OIBS1>().begin();
    for (auto ibs1:bs1->Iterate<Real_OIBS1>())
    {
        auto& S1=ibs1->Overlap();
        auto& S2=(*ibs2)->Overlap();
        EXPECT_EQ(S1,S2);
        EXPECT_EQ(&S1,&S2);
        ++ibs2;
    }
}
TEST_F(DBCach1Tests,BSplineKinetic)
{
    auto ibs2=bs2->Iterate<Real_OIBS1>().begin();
    for (auto ibs1:bs1->Iterate<Real_OIBS1>())
    {
        auto& S1=ibs1->Kinetic();
        auto& S2=(*ibs2)->Kinetic();
        EXPECT_EQ(S1,S2);
        EXPECT_EQ(&S1,&S2);
        ++ibs2;
    }
}
TEST_F(DBCach1Tests,BSplineNuclear)
{
    auto ibs2=bs2->Iterate<Real_OIBS1>().begin();
    for (auto ibs1:bs1->Iterate<Real_OIBS1>())
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

TEST_F(DBCach1Tests,BSplineDirect)
{
    auto ibs21=bs2->Iterate<Orbital_HF_IBS1<double>>().begin();
    for (auto ibs11:bs1->Iterate<Orbital_HF_IBS1<double>>())
    {
        auto ibs22=bs2->Iterate<Orbital_HF_IBS1<double>>().begin();
        for (auto ibs12:bs1->Iterate<Orbital_HF_IBS1<double>>())
        {
            const ERI4& J1=ibs11->Direct(*ibs12);
            const ERI4& J2=(*ibs21)->Direct(**ibs22);
            EXPECT_EQ(J1,J2);
            EXPECT_EQ(&J1,&J2);
            ++ibs22;
        }
        ++ibs21;
    }
}
TEST_F(DBCach1Tests,BSplineExchange)
{
    auto ibs21=bs2->Iterate<Orbital_HF_IBS1<double>>().begin();
    for (auto ibs11:bs1->Iterate<Orbital_HF_IBS1<double>>())
    {
        auto ibs22=bs2->Iterate<Orbital_HF_IBS1<double>>().begin();
        for (auto ibs12:bs1->Iterate<Orbital_HF_IBS1<double>>())
        {
            const ERI4& K1=ibs11->Exchange(*ibs12);
            const ERI4& K2=(*ibs21)->Exchange(**ibs22);
            EXPECT_EQ(K1,K2);
            EXPECT_EQ(&K1,&K2);
            ++ibs22;
        }
        ++ibs21;
    }
}
