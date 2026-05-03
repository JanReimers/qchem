// File: UnitTests/BasisSet_Atom.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <iostream>
// #include <cmath>
#include <blaze/Math.h>
using std::cout;
using std::endl;

import qchem.BasisSet.DB_Cache1;
import BasisSet.Atom.BSpline.NR.BS_Evaluator;
import BasisSet.Atom.BSpline.NR.IBS_Evaluator;

import qchem.IrrepBasisSet1;
import qchem.BasisSet.Atom.IBS1;
import qchem.BasisSet.Atom.IE1;
import qchem.BasisSet.Internal.Common1;

import qchem.Orbital_1E_IBS1;
import qchem.BasisSet1;

import qchem.Symmetry.Yl;
import qchem.Symmetry.AtomEC;
import qchem.Cluster;
import qchem.Types;

bool operator==(const ERI4& a, const ERI4& b);
// {
//     static double eps=5e-16;
//     if (a.size()!=b.size()) return false;
//     for (size_t i=0;i<a.Nab();i++)
//         for (size_t j=0;j<a.Nab();j++)
//             if (norm(a(i,j)-b(i,j))>eps) 
//             {
//                 std::cout << "a(" << i << "," << j << ")=" << a(i,j);
//                 std::cout << "b(" << i << "," << j << ")=" << b(i,j);
//                 std::cout << "[a-b](" << i << "," << j << ")=" << a(i,j)-b(i,j);
//                 std::cout << "norm(a(i,j)-b(i,j))=" << norm(a(i,j)-b(i,j)) << std::endl;
//                 return false;
//             }
//     return true;
// }

namespace AtomBS
{

namespace BSpline1
{
template <size_t K,class Evaluator> class Orbital_IBS
    : public AtomBS::Orbital_HF_IBS1
    , private Evaluator
{
public:
    Orbital_IBS(BS_Evaluator* bse,size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : AtomBS::Orbital_HF_IBS1(bse,yl)
    , Evaluator(N,rmin,rmax,yl)
    {};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet1*,const Cluster*) const {return 0;}
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet1*,const Cluster*) const {return 0;}
    virtual size_t GetNumFunctions() const {return Evaluator::size();}
    virtual const IBS_Evaluator* GetEvaluator() const {return this;}
    virtual       IBS_Evaluator* GetEvaluator()       {return this;}
    // virtual size_t  size           () const {return BSpline_IBS<K>::size();}

};


template <size_t K,template<size_t> class Evaluator> class BasisSet
    : public virtual ::BasisSet1
    , public ::BS_Common1
    , public Evaluator<K> 
{
    using oibs_t=Orbital_IBS<K,typename Evaluator<K>::IBS_Evaluator_t>; //Corresponding Orbital IBS type
public:
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(this,N,rmin,rmax,ir));  
     
        Evaluator<K>::BuildCache(LMax);
    }
private:
    void Insert(oibs_t* oibs)
    {
        ::BS_Common1::Insert(oibs);
        Evaluator<K>::Register(oibs->GetEvaluator());
    }
};


}} //namespace 


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
// TESTF(DBCach1Tests,Overlap)
// {
//     Init(86);
//     std::vector<rsmat_t*> S1s;
//     std::vector<rsmat_t*> S2s;

//     for (auto ibs:*bs1)
//         S1s.push_back(&ibs->Overlap());
//     for (auto ibs:*bs2)
//         S2s.push_back(&ibs->Overlap());

//     for (auto i:iv_t(0,S1s.size()))
//     {
//         EXPECT_EQ(S1s[i],S2s[i]);
//     }
// }
