// File: UnitTests/SlaterIntegral.C 


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet;
import qchem.Constants;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Symmetry.Spherical;
import qchem.Blaze;

using std::cout;
using std::endl;
using Real_BS  =BasisSet::Real_BS;
using Real_OIBS=BasisSet::Real_OIBS;
//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class SlaterRadialIntegralTests : public ::testing::Test
{
public:
    SlaterRadialIntegralTests()
    : Lmax(4    )
    , Z(1)
    , bs(0)
    , cl(new Atom(Z,0.0,Vector3D(0,0,0)))
    {
        nlohmann::json js = {
        {"type",BasisSet::Atom::Type::Slater},
        {"N", 6}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSet::Atom::Factory(js,75);
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        MeshParams rmp({qchem::MHL,200,3,2.0,qchem::Gauss,32,0,0,3});
        rmintegrator=new MeshIntegrator<double>(cl->CreateMesh(rmp));
        //cout << *bs << endl;
    }
    
    // bool   supported(const Evaluator&,const Evaluator&,int ia, int ib, int ic, int id) const;
    // double R0       (const Evaluator&,const Evaluator&,int ia, int ib, int ic, int id) const;
    
    
    int Lmax, Z;
    Real_BS* bs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
    MeshIntegrator<double>* rmintegrator;
};


TEST_F(SlaterRadialIntegralTests, Overlap)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        for (auto d:blazem::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        rsmat_t Snum = mintegrator->Overlap(*oi);
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,1e-8);
    }
}

TEST_F(SlaterRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t Hn=oi->Nuclear(cl);
        rsmat_t Hnnum = -1*mintegrator->Inv_r1(*oi);
        EXPECT_NEAR(blazem::max(blazem::abs(Hn-Hnnum)),0.0,1e-7);

    }
}

TEST_F(SlaterRadialIntegralTests, Kinetic)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();
        //cout << S << endl;
        int l=Getl(oi->GetSymmetry());;
        rsmat_t Knum = mintegrator->Grad2(*oi) + l*(l+1)*mintegrator->Inv_r2(*oi);
        EXPECT_NEAR(blazem::max(blazem::abs(K-Knum)),0.0,1e-10);
        
        // cout << "K=" << K << endl;
        // cout << "Knum=" << Knum << endl;
    }
}


