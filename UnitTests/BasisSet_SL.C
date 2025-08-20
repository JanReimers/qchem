// File: UnitTests/BasisSet_SL.C  unit test the Slater_IBS Evaluator
#include "gtest/gtest.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

import BasisSet.Atom.Slater_IBS;
import qchem.Mesh.Integrator;
import qchem.Molecule;


class BasisSet_SL: public ::testing::Test
{
   
public:

    BasisSet_SL() :  es{0.5,1.0,2.0}, N(es.size()), LMax(3)
    , cl(new Molecule())
    {
        for (size_t l=0;l<=LMax;l++)
            evals.push_back(new Slater_IBS(es,l,{}));
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(1,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }

    ~BasisSet_SL()
    {
        for (auto ev:evals) delete ev;
        delete cl;
        delete mintegrator;
    }

    std::valarray<double> es;
    size_t N,LMax;
    std::vector<IBS_Evaluator*> evals;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
};

TEST_F(BasisSet_SL,size)
{
    for (auto ev:evals)
        EXPECT_EQ(ev->size(),N);
}

using omls_t=IBS_Evaluator::omls_t;
TEST_F(BasisSet_SL,AnalyticOverlap)
{
    for (auto ev:evals)
    {
        size_t l=ev->Getl();
        omls_t S=ev->Overlap();
        // cout << S << endl;
        for (auto i:S.rows())
            for (auto j:S.cols(i))
            {
                double a=es[i-1],b=es[j-1];
                EXPECT_NEAR(S(i,j),pow(2.0/(sqrt(a/b)+sqrt(b/a)),2*l+3),1e-15);
            }
    }
        
}

TEST_F(BasisSet_SL,NumericalOverlap)
{
    for (auto ev:evals)
    {
        omls_t S=ev->Overlap();
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        omls_t Snum = mintegrator->Overlap(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,1e-15);
    }
        
}

TEST_F(BasisSet_SL,Grad2)
{
    for (auto ev:evals)
    {
        omls_t S=ev->Grad2();
        omls_t Snum = mintegrator->Grad2(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,2e-15);
    }
        
}

TEST_F(BasisSet_SL,Inv_r1)
{
    for (auto ev:evals)
    {
        omls_t S=ev->Inv_r1();
        omls_t Snum = mintegrator->Inv_r1(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,3e-14);
    }
        
}

TEST_F(BasisSet_SL,Inv_r2)
{
    for (auto ev:evals)
    {
        omls_t S=ev->Inv_r2();
        omls_t Snum = mintegrator->Inv_r2(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,4e-9);
    }
        
}