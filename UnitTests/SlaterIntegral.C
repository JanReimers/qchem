// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater/BasisSet.H"
//#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
//#include "Imp/Integrals/SlaterRadialIntegrals.H"
//#include "Imp/Integrals/Wigner3j.H"
#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"
#include "Mesh/MeshIntegrator.H"
#include "Misc/DFTDefines.H"
#include "oml/imp/ran250.h"
#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;


//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class SlaterRadialIntegralTests : public ::testing::Test
{
public:
    SlaterRadialIntegralTests()
    : lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , ie(new Slater::IntegralEngine())
    , db(new HeapDB<double>())
    , ibs(new Slater::IrrepBasisSet(lap,db,5,.01,100.0,0))
    , bs(new Slater::BasisSet(lap,5,.01,100.0,4))
    , mesh(0)
    {
        StreamableObject::SetToPretty();
        RadialMesh*  rm=new MHLRadialMesh(200,3U,2.0); //mem leak
        AngularMesh* am=new GaussAngularMesh(1);      //mem leak
        mesh=new AtomMesh(*rm,*am); 
    }
    
    LAParams lap;
    AnalyticIE<double>* ie;
    IntegralDataBase<double>* db;
    Slater::IrrepBasisSet* ibs;
    Slater::BasisSet* bs;
    Mesh* mesh;
};

TEST_F(SlaterRadialIntegralTests, Overlap)
{
    MeshIntegrator<double> mintegrator(mesh);
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> S=ie->MakeOverlap(*i);
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        cout << S << endl;
        SMatrix<double> Snum = mintegrator.Overlap(**i);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);

    }
}

