// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/Symmetry/YlQN.H"
#include "Imp/Integrals/MeshIntegrator.H"
#include "Imp/Misc/DFTDefines.H"
#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
#include <MeshParams.H>
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
    , ie(new SphericalGaussian::IntegralEngine())
    , bs(new SphericalGaussian::BasisSet(lap,5,.01,100.0,Lmax))
    , cl(new Molecule())
    , mintegrator()
    {
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(Z,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
        
    }
    
    int Lmax, Z;
    LAParams lap;
    AnalyticIE<double>* ie;
    SphericalGaussian::BasisSet* bs;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};

TEST_F(GaussianRadialIntegralTests, Overlap)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> S=ie->MakeOverlap(*i);
        auto oi=dynamic_cast<const Orbital_IBS<double>*>(*i);
        SMatrix<double> S1=oi->Integrals(qchem::Overlap1);

        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        //cout << S << endl;
        SMatrix<double> Snum = mintegrator->Overlap(**i);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-8);
        EXPECT_NEAR(Max(fabs(S-S1)),0.0,1e-14);

    }
}

TEST_F(GaussianRadialIntegralTests, Nuclear)
{
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> Hn=ie->MakeNuclear(*i,*cl);
        auto oi=dynamic_cast<const Orbital_IBS<double>*>(*i);
        SMatrix<double> Hn1=oi->Integrals(qchem::Nuclear1,cl);
        //cout << S << endl;
        SMatrix<double> Hnnum = -1*mintegrator->Nuclear(**i);
        EXPECT_NEAR(Max(fabs(Hn-Hnnum)),0.0,1e-8);
        EXPECT_NEAR(Max(fabs(Hn-Hn1)),0.0,1e-14);

    }
}

TEST_F(GaussianRadialIntegralTests, Kinetic)
{
    
    for (auto i=bs->beginT();i!=bs->end();i++)
    {
        SMatrix<double> K=ie->MakeKinetic(*i);
        auto oi=dynamic_cast<const Orbital_IBS<double>*>(*i);
        SMatrix<double> K1=oi->Integrals(qchem::Kinetic1);
        //cout << S << endl;
        SMatrix<double> Knum = 0.5*mintegrator->Grad(**i); //This give the wrong answer for l>0

        // We need to add the l*(l+1) term that comes from the angular integrals.
        // Lost of dynamic cast just to get at L!
        const QuantumNumber& qn=i->GetQuantumNumber();
        const YlQN& sqn=dynamic_cast<const YlQN& >(qn);
        int l=sqn.GetL();
        const SphericalGaussian::IrrepBasisSet* sg=dynamic_cast<const SphericalGaussian::IrrepBasisSet*>(*i);
        assert(sg);
        for (auto i:Knum.rows())
            for (auto j:Knum.cols(i))
                Knum(i,j)+=0.5*((l)*(l+1))*GaussianIntegral(sg->es(i)+sg->es(j),2*l-2)*sg->ns(i)*sg->ns(j);
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-12);
        EXPECT_NEAR(Max(fabs(K-K1)),0.0,1e-14);

    }
}

