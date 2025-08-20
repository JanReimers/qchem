// File: UnitTests/BasisSet_SL.C  unit test the Slater_IBS Evaluator
#include "gtest/gtest.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

import BasisSet.Atom.Slater_IBS;
import qchem.Mesh.Integrator;
import qchem.Molecule;
import Common.Constants;


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

double R0(double a, double b, int la, int lb) 
{
    int nab=2*la+2;
    int ncd=2*lb+2;
    double f=FourPi2/(a*b*(a+b));
    if (la==0 && lb==0)
        return 2*f*( 
                    1/(pow(a,1)*pow(b,1)*pow(a+b,0))
                    +1/pow(a+b,2)
                    );
    // if (la==1 && lb==1)
    // return 12*f*( 
    //                  1/(pow(a,1)*pow(b,3)*pow(a+b,0))
    //                 +1/(pow(a,2)*pow(b,2)*pow(a+b,0))
    //                 +2/(pow(a,0)*pow(b,0)*pow(a+b,4))
    //                 -1/(pow(a,0)*pow(b,2)*pow(a+b,2))
    //                 -1/(pow(a,0)*pow(b,3)*pow(a+b,1))
                   
    //                 );
    if (nab==4 && ncd==2)
        return 12*f*( 
                        2/(pow(a,3)*pow(b,1)*pow(a+b,0))
                      + 2/(pow(a,0)*pow(b,0)*pow(a+b,4))
                      + 1/(pow(a,1)*pow(b,0)*pow(a+b,3))
                      - 1/(pow(a,3)*pow(b,0)*pow(a+b,1))
                        );
    if (nab==2 && ncd==4)
        return 12*f*( 
                       2/(pow(a,1)*pow(b,3)*pow(a+b,0))
                     + 2/(pow(a,0)*pow(b,0)*pow(a+b,4))
                     + 1/(pow(a,0)*pow(b,1)*pow(a+b,3))
                     - 1/(pow(a,0)*pow(b,3)*pow(a+b,1))
                    );
    if (nab==4 && ncd==4)
         return 144*f*( 
           1/(pow(a,2)*pow(b,2)*pow(a+b,2))
         + 2/(pow(a,1)*pow(b,1)*pow(a+b,4))
         + 5/(pow(a,0)*pow(b,0)*pow(a+b,6))
         + 1/(pow(a,3)*pow(b,3)*pow(a+b,0))
          );
          
    if (nab==2 && ncd==6)
        return 240*f*( 
                       1/(pow(a,1)*pow(b,5)*pow(a+b,0))
                     + 1/(pow(a,2)*pow(b,4)*pow(a+b,0))
                     + 3/(pow(a,0)*pow(b,0)*pow(a+b,6))
                     - 2/(pow(a,1)*pow(b,0)*pow(a+b,5))
                     - 1/(pow(a,2)*pow(b,0)*pow(a+b,4))
                    );
     if (nab==6 && ncd==2)
        return 240*f*( 
                       1/(pow(a,5)*pow(b,1)*pow(a+b,0))
                     + 1/(pow(a,4)*pow(b,2)*pow(a+b,0))
                     + 3/(pow(a,0)*pow(b,0)*pow(a+b,6))
                     - 2/(pow(a,0)*pow(b,1)*pow(a+b,5))
                     - 1/(pow(a,0)*pow(b,2)*pow(a+b,4))
                    );
    if (nab==4 && ncd==6)
        return 1440*f*( 
                        2/(pow(a,3)*pow(b,5)*pow(a+b,0))
                     +  2/(pow(a,4)*pow(b,4)*pow(a+b,0))
                     + 28/(pow(a,0)*pow(b,0)*pow(a+b,8))
                     -  7/(pow(a,1)*pow(b,0)*pow(a+b,7))
                     - 12/(pow(a,2)*pow(b,0)*pow(a+b,6))
                     -  7/(pow(a,3)*pow(b,0)*pow(a+b,5))
                     -  2/(pow(a,4)*pow(b,0)*pow(a+b,4))
                    );
    if (nab==6 && ncd==4)
        return 1440*f*( 
                        2/(pow(a,5)*pow(b,3)*pow(a+b,0))
                     +  2/(pow(a,4)*pow(b,4)*pow(a+b,0))
                     + 28/(pow(a,0)*pow(b,0)*pow(a+b,8))
                     -  7/(pow(a,0)*pow(b,1)*pow(a+b,7))
                     - 12/(pow(a,0)*pow(b,2)*pow(a+b,6))
                     -  7/(pow(a,0)*pow(b,3)*pow(a+b,5))
                     -  2/(pow(a,0)*pow(b,4)*pow(a+b,4))
                    );
    
    if (nab==6 && ncd==6)
        return 86400*f*( 
                        1/(pow(a,5)*pow(b,5)*pow(a+b,0))
                     +  1/(pow(a,6)*pow(b,4)*pow(a+b,0))
                     + 42/(pow(a,0)*pow(b,0)*pow(a+b,10))
                     - 14/(pow(a,2)*pow(b,0)*pow(a+b,8))
                     - 14/(pow(a,3)*pow(b,0)*pow(a+b,7))
                     -  9/(pow(a,4)*pow(b,0)*pow(a+b,6))
                     -  4/(pow(a,5)*pow(b,0)*pow(a+b,5))
                     -  1/(pow(a,6)*pow(b,0)*pow(a+b,4))
                    );
    
    
                    
    assert(false);
    return 0;

}


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
        int l=ev->Getl();
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

using omlv_t=IBS_Evaluator::omlv_t;
TEST_F(BasisSet_SL,Charge)
{
    for (auto ev:evals)
    {
        omlv_t S=ev->Charge();
        omlv_t Snum = mintegrator->Integrate(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,5e-14);
    }
        
}

using ds_t=IBS_Evaluator::ds_t;
TEST_F(BasisSet_SL,Repulsion)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        if (l>=3) continue;
        omls_t S=ev->Repulsion();
        // cout << "l=" << l << " S=" << S << endl;
        ds_t   ns=ev->Norm();
        for (auto i:S.rows())
            for (auto j:S.cols(i))
            {
                double a=es[i-1],b=es[j-1];
                double rerr=(S(i,j) - R0(a,b,l,l)*ns[i-1]*ns[j-1])/S(i,j);
                EXPECT_NEAR(rerr,0.0,3e-14);
            }
    }
        
}