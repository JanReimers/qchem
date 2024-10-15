// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
//#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
//#include "Imp/Integrals/SlaterRadialIntegrals.H"
//#include "Imp/Integrals/Wigner3j.H"
#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"
#include "Mesh/MeshIntegrator.H"
#include "Misc/DFTDefines.H"
//#include "Cluster/Atom.H"
#include "Cluster/Molecule.H"
#include "Cluster.H"
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
class GaussianRadialIntegralTests : public ::testing::Test
{
public:
    GaussianRadialIntegralTests()
    : Lmax(4    )
    , Z(1)
    , lap({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , ie(new SphericalGaussian::IntegralEngine())
    , db(new HeapDB<double>())
    , ibs(new SphericalGaussian::IrrepBasisSet(lap,db,5,.01,100.0,0))
    , bs(new SphericalGaussian::BasisSet(lap,5,.01,100.0,Lmax))
    , mesh(0)
    , cl(new Molecule())
    , mintegrator()
    {
        StreamableObject::SetToPretty();
        RadialMesh*  rm=new MHLRadialMesh(200,3U,2.0); //mem leak
        AngularMesh* am=new GaussAngularMesh(1);      //mem leak
        mesh=new AtomMesh(*rm,*am); 
        mintegrator=new MeshIntegrator<double>(mesh);
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
    }
    
    int Lmax, Z;
    LAParams lap;
    AnalyticIE<double>* ie;
    IntegralDataBase<double>* db;
    SphericalGaussian::IrrepBasisSet* ibs;
    SphericalGaussian::BasisSet* bs;
    Mesh* mesh;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};

TEST_F(GaussianRadialIntegralTests, Overlap)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> S=ie->MakeOverlap(*i);
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        SMatrix<double> Snum = mintegrator->Overlap(**i);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);

    }
}

TEST_F(GaussianRadialIntegralTests, Nuclear)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> Hn=ie->MakeNuclear(*i,*cl);
        //cout << S << endl;
        SMatrix<double> Hnnum = -1*mintegrator->Nuclear(**i);
        EXPECT_NEAR(Max(fabs(Hn-Hnnum)),0.0,1e-8);

    }
}

TEST_F(GaussianRadialIntegralTests, Kinetic)
{
    
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> K=ie->MakeKinetic(*i);
        //cout << S << endl;
        SMatrix<double> Knum = 0.5*mintegrator->Grad(**i); //This give the wrong answer for l>0

        // We need to add the l*(l+1) term that comes from the angular integrals.
        // Lost of dynamic cast just to get at L!
        const QuantumNumber& qn=i->GetQuantumNumber();
        const SphericalSymmetryQN& sqn=dynamic_cast<const SphericalSymmetryQN& >(qn);
        int l=sqn.GetL();
        const SphericalGaussian::IrrepBasisSet* sg=dynamic_cast<const SphericalGaussian::IrrepBasisSet*>(*i);
        assert(sg);
        for (auto i:Knum.rows())
            for (auto j:Knum.cols(i))
                Knum(i,j)+=0.5*((l)*(l+1))*GaussianIntegral(sg->es(i)+sg->es(j),2*l-2)*sg->ns(i)*sg->ns(j);
            
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-8);

    }
}

