// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"
#include "Imp/BasisSet/Atom/l/Gaussian_IE.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_IBS.H"
#include "Symmetry/Yl.H"
#include "Mesh/MeshIntegrator.H"
#include "Common/DFTDefines.H"
#include "Cluster/Molecule.H"
#include "Cluster/Atom.H"
#include <Mesh/MeshParams.H>
#include <Cluster.H>
#include "oml/smatrix.h"
#include "oml/matrix.h"
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
    , bs(new Atoml::Gaussian::BasisSet(5,.01,100.0,Lmax))
    , cl(new Molecule())
    , mintegrator()
    {
        bs->Set(lap);
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        
    }
    
    int Lmax, Z;
    LAParams lap;
    Atoml::Gaussian::BasisSet* bs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};




TEST_F(GaussianRadialIntegralTests, Overlap)
{
    for (auto oi:bs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> S=oi->Overlap();

        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        SMatrix<double> Snum = mintegrator->Overlap(*oi);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);
       
    }
}

TEST_F(GaussianRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> Hn=oi->Nuclear(cl);
        //cout << S << endl;
        SMatrix<double> Hnnum = -1*mintegrator->Inv_r1(*oi);
        EXPECT_NEAR(Max(fabs(Hn-Hnnum)),0.0,1e-8);

    }
}

TEST_F(GaussianRadialIntegralTests, Kinetic)
{
    
    for (auto oi:bs->Iterate<TOrbital_IBS<double> >())
    {
        SMatrix<double> K=oi->Kinetic();
        //cout << S << endl;
        int l=dynamic_cast<const Angular_Sym* >(oi->GetSymmetry().get())->GetL();
        SMatrix<double> Knum = mintegrator->Grad(*oi) + l*(l+1)*mintegrator->Inv_r2(*oi);
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-12);
        
    }
}

