// File: UnitTests/BasisSet_Atom.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <iostream>
// #include <cmath>
#include <blaze/Math.h>
using std::cout;
using std::endl;

import qchem.Cluster;
import qchem.BasisSet.DB_Cache1;
import BasisSet.Atom.BSpline.NR.IBS_Evaluator;
import qchem.Types;
import qchem.Orbital_1E_IBS;
import qchem.IrrepBasisSet1;
import qchem.BasisSet.Atom.IBS1;
import qchem.BasisSet.Atom.IE1;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.Orbital_1E_IBS1;

namespace AtomBS1
{
namespace BSpline1
{
template <size_t K> class Orbital_IBS
    : public Orbital_IBS1<double>
    , private BSpline_IBS<K>
    , public AtomBS::IrrepBasisSet1
    , public AtomIE_Overlap1<double>
    , public AtomIE_Kinetic1<double>
    , public AtomIE_Nuclear1<double>
{
public:
    Orbital_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_IBS1<double>(yl)
    , BSpline_IBS<K>(N,rmin,rmax,yl)
    , AtomBS::IrrepBasisSet1(this)
    , AtomIE_Overlap1<double>(this)
    , AtomIE_Kinetic1<double>(this)
    , AtomIE_Nuclear1<double>(this)
    {};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const {return 0;}
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const {return 0;}
    virtual size_t GetNumFunctions() const {return BSpline_IBS<K>::size();}
    // virtual size_t  size           () const {return BSpline_IBS<K>::size();}

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
        , ibs1(new AtomBS1::BSpline1::Orbital_IBS<6>(3,.5,2.0,yl))
        , ibs2(new AtomBS1::BSpline1::Orbital_IBS<6>(3,.5,2.0,yl))
        // , bs1(0), bs2(0)
    {
        theGlobalCache=new IntegralsCache_RAM<double>();
    }
    ~DBCach1Tests()
    {
        delete cl_hydrogen;
        delete cl_hydrogen_100;
        delete cl_helium;
        delete ibs1;
        delete ibs2;
        // delete bs1;
        // delete bs2;
        delete theGlobalCache;
    }
    void Init(size_t Z)
    {
        // bs1=new AtomBS::BSpline::BasisSet1<6>(3,0.1,10.0,Atom_EC(Z));
        // bs2=new AtomBS::BSpline::BasisSet1<6>(3,0.1,10.0,Atom_EC(Z));
    }
    
    Cluster *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    Irrep_QNs::sym_t yl;
     Real_OIBS1 *ibs1,*ibs2;
    // BasisSet *bs1,*bs2;
};


TEST_F(DBCach1Tests,BSplineOverlap)
{
    auto& S1=ibs1->Overlap();
    auto& S2=ibs2->Overlap();
    EXPECT_EQ(S1,S2);
    EXPECT_EQ(&S1,&S2);
}
TEST_F(DBCach1Tests,BSplineKinetic)
{
    auto& S1=ibs1->Kinetic();
    auto& S2=ibs2->Kinetic();
    EXPECT_EQ(S1,S2);
    EXPECT_EQ(&S1,&S2);
}
TEST_F(DBCach1Tests,BSplineNuclear)
{
    auto& S1=ibs1->Nuclear(cl_hydrogen);
    auto& S2=ibs2->Nuclear(cl_hydrogen);
    auto& S3=ibs2->Nuclear(cl_hydrogen_100);
    auto& S4=ibs2->Nuclear(cl_helium);
    EXPECT_EQ(S1,S2);
    EXPECT_EQ(S1,S3);
    EXPECT_NE(S1,S4);
    EXPECT_EQ(&S1,&S2);
    EXPECT_NE(&S1,&S3);
    EXPECT_NE(&S1,&S4);

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
