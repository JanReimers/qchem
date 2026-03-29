// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
import qchem.LAParams;
#include <iostream>
#include <fstream>
#include <cmath>
#include <blaze/Math.h>

import qchem.Factory;
import qchem.BasisSet;
import qchem.IrrepBasisSet;
import Common.Constants;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Atom;
import qchem.Molecule;
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
    , lap({qchem::SVD,1e-6})
    , bs(0)
    , cl(new Molecule())
    , mintegrator()
    {
        nlohmann::json js = {
        {"type",BasisSetAtom::Type::Gaussian},
        {"N", 5}, {"emin", 0.01}, {"emax", 100.0},
        };
        bs=BasisSetAtom::Factory(js,75);
        bs->Set(lap);
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        
    }
    
    int Lmax, Z;
    LAParams lap;
    BasisSet* bs;
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
        int l=dynamic_cast<const Angular_Sym* >(oi->GetSymmetry().get())->GetL();
        rsmat_t Knum = mintegrator->Grad2(*oi) + l*(l+1)*mintegrator->Inv_r2(*oi);
        EXPECT_NEAR(max(abs(K-Knum)),0.0,1e-12);
        
    }
}

