// File: ERIList.C  Test the DFT persistance classes
#include <blaze/math/expressions/DMatDMatEqualExpr.h>
#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>

import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS;
import qchem.BasisSet.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet;

import qchem.Constants;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Symmetry.Spherical;
import qchem.Symmetry.Factory;
import qchem.Blaze;

using std::cout;
using std::endl;
using Real_BS=BasisSet::Real_BS;
using Real_OIBS=BasisSet::Real_OIBS;

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
        {"type",BasisSet::Atom::Type::Gaussian},
        {"N", 5}, {"emin", 0.01}, {"emax", 100.0},
        };
        bs=BasisSet::Atom::Factory(js,75);
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

        for (auto d:blazem::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        rsmat_t Snum = mintegrator->Overlap(*oi);
        EXPECT_NEAR(blazem::max(blazem::abs(S-Snum)),0.0,1e-8);
       
    }
}

TEST_F(GaussianRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t Hn=oi->Nuclear(cl);
        //cout << S << endl;
        rsmat_t Hnnum = -1*mintegrator->Inv_r1(*oi);
        EXPECT_NEAR(blazem::max(blazem::abs(Hn-Hnnum)),0.0,1e-8);

    }
}

TEST_F(GaussianRadialIntegralTests, Kinetic)
{
    
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();
        //cout << S << endl;
        int l=Getl(oi->GetSymmetry());;
        rsmat_t Knum = mintegrator->Grad2(*oi) + l*(l+1)*mintegrator->Inv_r2(*oi);
        EXPECT_NEAR(blazem::max(blazem::abs(K-Knum)),0.0,1e-12);
        
    }
}

TEST_F(GaussianRadialIntegralTests,RkSymmetry_l0)
{
    using namespace BasisSet::Atom::Evaluators::Gaussian;
    typedef rvec11_t rvec11_t; 
    auto eval=new NR_Evaluator(15,.03,20.0,Symmetry::YFactory(0));
    auto cache4=eval->MakeCache4();
    cache4->Register(eval);
    auto ns=eval->Norm();

    rvec11_t Ak({1,1,1,1,1,1,1,1,1,1,1});

    {
        {
            auto rk0001=dynamic_cast<const Rk*>(cache4->Create(0,0,0,1));
            auto rk0100=dynamic_cast<const Rk*>(cache4->Create(0,1,0,0));
            auto rk0010=dynamic_cast<const Rk*>(cache4->Create(0,0,1,0));
            cout << " 0001=" << std::setprecision(12) << rk0001->DirectR0  (0,0   )*ns[0]*ns[0]*ns[0]*ns[1]*FourPi2 << endl;
            cout << " 0100=" << std::setprecision(12) << rk0100->DirectR0  (0,0   )*ns[0]*ns[0]*ns[0]*ns[1]*FourPi2 << endl;
            cout << " 0010=" << std::setprecision(12) << rk0010->DirectR0  (0,0   )*ns[0]*ns[0]*ns[1]*ns[0]*FourPi2 << endl;
            cout << "J0001=" << std::setprecision(12) << eval->direct(rk0001,0,0,Ak)*ns[0]*ns[0]*ns[0]*ns[1]*FourPi2 << endl;
            cout << "J0100=" << std::setprecision(12) << eval->direct(rk0100,0,0,Ak)*ns[0]*ns[1]*ns[0]*ns[0]*FourPi2 << endl;
            cout << "J0010=" << std::setprecision(12) << eval->direct(rk0010,0,0,Ak)*ns[0]*ns[0]*ns[1]*ns[0]*FourPi2 << endl;

            EXPECT_NEAR(rk0001->DirectR0  (0,0   ),rk0100->DirectR0  (0,0   ),1e-15);
            EXPECT_NEAR(rk0001->DirectR0  (0,0   ),rk0100->DirectR0  (0,0   ),1e-15);
            EXPECT_NEAR(rk0001->DirectRk  (0,0,Ak),rk0100->DirectRk  (0,0,Ak),1e-15);
            EXPECT_NEAR(rk0001->DirectRk  (0,0,Ak),rk0100->DirectRk  (0,0,Ak),1e-15);
            EXPECT_NEAR(rk0001->ExchangeRk(0,0,Ak),rk0100->ExchangeRk(0,0,Ak),1e-15);
            EXPECT_NEAR(rk0001->ExchangeRk(0,0,Ak),rk0100->ExchangeRk(0,0,Ak),1e-15);
        }
    }
}



