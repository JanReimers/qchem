// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <blaze/Math.h>

import qchem.LAParams;
import qchem.Factory;
import qchem.Orbital_HF_IBS;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.HeapDB;

import qchem.BasisSet;
import qchem.IrrepBasisSet;
import Common.Constants;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.Atom;
import qchem.Molecule;
import qchem.Symmetry.Angular;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;

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
    : Lmax(4    )
    , Z(1)
    , lap({qchem::SVD,1e-6})
    , bs(0)
    , cl(new Molecule())
    {
        nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", 6}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSetAtom::Factory(js,75);
        bs->Set(lap);
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        MeshParams rmp({qchem::MHL,200,3,2.0,qchem::Gauss,32,0,0,3});
        rmintegrator=new MeshIntegrator<double>(cl->CreateMesh(rmp));
        //cout << *bs << endl;
    }
    
    // bool   supported(const Slater_IBS&,const Slater_IBS&,int ia, int ib, int ic, int id) const;
    // double R0       (const Slater_IBS&,const Slater_IBS&,int ia, int ib, int ic, int id) const;
    
    
    int Lmax, Z;
    LAParams lap;
    BasisSet* bs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
    MeshIntegrator<double>* rmintegrator;
};


TEST_F(SlaterRadialIntegralTests, Overlap)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t S=oi->Overlap();
        for (auto d:blaze::diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        rsmat_t Snum = mintegrator->Overlap(*oi);
        EXPECT_NEAR(max(abs(S-Snum)),0.0,1e-8);
    }
}

TEST_F(SlaterRadialIntegralTests, Nuclear)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t Hn=oi->Nuclear(cl);
        rsmat_t Hnnum = -1*mintegrator->Inv_r1(*oi);
        EXPECT_NEAR(max(abs(Hn-Hnnum)),0.0,1e-7);

    }
}

TEST_F(SlaterRadialIntegralTests, Kinetic)
{
    for (auto oi:bs->Iterate<Real_OIBS >())
    {
        rsmat_t K=oi->Kinetic();
        //cout << S << endl;
        int l=dynamic_cast<const Angular_Sym* >(oi->GetSymmetry().get())->GetL();
        rsmat_t Knum = mintegrator->Grad2(*oi) + l*(l+1)*mintegrator->Inv_r2(*oi);
        EXPECT_NEAR(max(abs(K-Knum)),0.0,1e-10);
        
        // cout << "K=" << K << endl;
        // cout << "Knum=" << Knum << endl;
    }
}


