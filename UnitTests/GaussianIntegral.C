// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
import qchem.LAParams;
#include <iostream>
#include <fstream>
#include <cmath>
#include <blaze/Math.h>

import qchem.Factory;
import Common.Constants;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Symmetry.Yl;

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
    , bs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    , mintegrator()
    {
        nlohmann::json js = {
        {"type",BasisSetAtomFactory::Type::Gaussian},
        {"N", 5}, {"emin", 0.01}, {"emax", 100.0},
        };
        bs=BasisSetAtomFactory::Factory(js,75);
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        
    }
    
    int Lmax, Z;
    Real_BS* bs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};




TEST_F(GaussianRadialIntegralTests, Overlap)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();

        for (auto d:blaze::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        rsmat_t Snum = mintegrator->Overlap(*oi);
        EXPECT_NEAR(max(abs(S-Snum)),0.0,1e-8);
       
    }
}

TEST_F(GaussianRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t Hn=oi->Nuclear(cl);
        //cout << S << endl;
        rsmat_t Hnnum = -1*mintegrator->Inv_r1(*oi);
        EXPECT_NEAR(max(abs(Hn-Hnnum)),0.0,1e-8);

    }
}

TEST_F(GaussianRadialIntegralTests, Kinetic)
{
    
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();
        //cout << S << endl;
        int l=oi->CastSymmetry<Angular_Sym>().GetL();
        rsmat_t Knum = mintegrator->Grad2(*oi) + l*(l+1)*mintegrator->Inv_r2(*oi);
        EXPECT_NEAR(max(abs(K-Knum)),0.0,1e-12);
        
    }
}

