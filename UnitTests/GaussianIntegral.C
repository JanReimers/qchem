// File: ERIList.C  Test the DFT persistance classes


#include "gtest/gtest.h"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"
#include "Imp/BasisSet/Atom/l/Gaussian_IE.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/l/Gaussian_IBS.H"
#include "Imp/BasisSet/Atom/l/Yl.H"
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
        SMatrix<double> Knum = mintegrator->Grad(*oi); //This give the wrong answer for l>0

        // We need to add the l*(l+1) term that comes from the angular integrals.
        // Lost of dynamic cast just to get at L!
        const Symmetry& qn=oi->GetSymmetry();
        const Yl_Sym& sqn=dynamic_cast<const Yl_Sym& >(qn);
        int l=sqn.GetL();
        const ::Gaussian::IrrepBasisSet* sg=dynamic_cast<const Gaussian::IrrepBasisSet*>(oi);
        assert(sg);
        for (auto i:Knum.rows())
            for (auto j:Knum.cols(i))
                Knum(i,j)+=((l)*(l+1))*Gaussian::Integral(sg->es(i)+sg->es(j),2*l-2)*sg->ns(i)*sg->ns(j);
        EXPECT_NEAR(Max(fabs(K-Knum)),0.0,1e-12);
        
    }
}

