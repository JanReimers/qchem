// File: UnitTests/BasisSet_SL.C  unit test the Slater_IBS Evaluator
#include "gtest/gtest.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

import BasisSet.Atom.Slater_IBS;
import BasisSet.Atom.Gaussian_IBS;
import qchem.Mesh.Integrator;
import qchem.Molecule;
import Common.Constants;


//----------------------------------------------------------------------------------------
//
//  Testing common to all atom basis set evaluators
//
class BasisSet_Common : public ::testing::Test
{
public:
    BasisSet_Common() : es{0.5,1.0,2.0}, N(es.size()), LMax(3), cl(new Molecule())
    {
        StreamableObject::SetToPretty();
        cl->Insert(new Atom(1,0.0,Vector3D(0,0,0)));
        MeshParams mp({qchem::MHL,200,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    ~BasisSet_Common()
    {
        for (auto ev:evals) delete ev;
        delete cl;
        delete mintegrator;
    }

    void TestOverlap(double eps) const;
    void TestGrad2  (double eps) const;
    void TestInv_r1 (double eps) const;
    void TestInv_r2 (double eps) const;
    void TestCharge (double eps) const;


    std::valarray<double> es;
    size_t N,LMax;
    std::vector<IBS_Evaluator*> evals;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;

};
using omls_t=IBS_Evaluator::omls_t;
using omlv_t=IBS_Evaluator::omlv_t;
void BasisSet_Common::TestOverlap(double eps) const
{
    EXPECT_GT(evals.size(),0);
    for (auto ev:evals)
    {
        omls_t S=ev->Overlap();
        for (auto d:Vector<double>(S.GetDiagonal())) EXPECT_NEAR(d,1.0,1e-15);
        omls_t Snum = mintegrator->Overlap(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,eps);
    }
}
void BasisSet_Common::TestGrad2  (double eps) const
{
    for (auto ev:evals)
    {
        omls_t S=ev->Grad2();
        omls_t Snum = mintegrator->Grad2(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,eps);
    }
        
}
void BasisSet_Common::TestInv_r1 (double eps) const
{
    for (auto ev:evals)
    {
        omls_t S=ev->Inv_r1();
        omls_t Snum = mintegrator->Inv_r1(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,eps);
    }
        
}
void BasisSet_Common::TestInv_r2 (double eps) const
{
    for (auto ev:evals)
    {
        omls_t S=ev->Inv_r2();
        omls_t Snum = mintegrator->Inv_r2(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,eps);
    }
        
}
void BasisSet_Common::TestCharge (double eps) const
{
    for (auto ev:evals)
    {
        omlv_t S=ev->Charge();
        omlv_t Snum = mintegrator->Integrate(*ev);
        EXPECT_NEAR(Max(fabs(S-Snum)),0.0,eps);
    }
        
}




//----------------------------------------------------------------------------------------
//
//  Testing atom Slater basis set evaluators
//

class BasisSet_SL: public BasisSet_Common
{
public:

    BasisSet_SL() : BasisSet_Common()
    {
        for (size_t l=0;l<=LMax;l++)
            evals.push_back(new Slater_IBS(es,l,{}));    
    }
    static double R0(double a, double b, int la, int lb);
};

double BasisSet_SL::R0(double a, double b, int la, int lb) 
{
    double f=FourPi2/(a*b*(a+b));
    if (la==0 && lb==0)
        return 2*f*( 
                    1/(pow(a,1)*pow(b,1)*pow(a+b,0))
                    +1/pow(a+b,2)
                    );
    if (la==1 && lb==1)
         return 144*f*( 
           1/(pow(a,2)*pow(b,2)*pow(a+b,2))
         + 2/(pow(a,1)*pow(b,1)*pow(a+b,4))
         + 5/(pow(a,0)*pow(b,0)*pow(a+b,6))
         + 1/(pow(a,3)*pow(b,3)*pow(a+b,0))
          );
          
    if (la==2 && lb==2)
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

TEST_F(BasisSet_SL,Overlap) {TestOverlap(1e-15);}
TEST_F(BasisSet_SL,Grad2  ) {TestGrad2  (2e-15);}
TEST_F(BasisSet_SL,Inv_r1 ) {TestInv_r1 (3e-14);}
TEST_F(BasisSet_SL,Inv_r2 ) {TestInv_r2 (4e-9);}
TEST_F(BasisSet_SL,Charge ) {TestCharge (5e-14);}

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
TEST_F(BasisSet_SL,AnalyticRepulsion)
{
    using ds_t=IBS_Evaluator::ds_t;
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

//----------------------------------------------------------------------------------------
//
//  Testing atom Gaussian basis set evaluators
//

class BasisSet_SG: public BasisSet_Common
{
public:

    BasisSet_SG() : BasisSet_Common()
    {
        for (size_t l=0;l<=LMax;l++)
            evals.push_back(new Gaussian_IBS(es,l,{}));    
    }
    static double R0(double a, double b, int la, int lb);
};

double BasisSet_SG::R0(double a, double b, int la, int lb) 
{
    double f=FourPi2*Pi12/8*(1/(a*b*sqrt(a+b)));
    if (la==0 && lb==0)
        return f;
    if (la==1 && lb==1)
         return f*3/4*(2/(a*b)+1/((a+b)*(a+b)));
          
    if (la==2 && lb==2)
        return f*15/(4*pow(a+b,4))*(2*(a+b)*(b/(a*a)+a/(b*b)) + 7*(b/a+a/b)+63./4.);
    
    
                    
    assert(false);
    return 0;

}


TEST_F(BasisSet_SG,Overlap) {TestOverlap(1e-15);}
TEST_F(BasisSet_SG,Grad2  ) {TestGrad2  (4e-15);}
TEST_F(BasisSet_SG,Inv_r1 ) {TestInv_r1 (3e-14);}
TEST_F(BasisSet_SG,Inv_r2 ) {TestInv_r2 (4e-9);}
TEST_F(BasisSet_SG,Charge ) {TestCharge (5e-14);}

TEST_F(BasisSet_SG,AnalyticOverlap)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        omls_t S=ev->Overlap();
        for (auto i:S.rows())
            for (auto j:S.cols(i))
            {
                double a=es[i-1],b=es[j-1];
                EXPECT_NEAR(S(i,j),pow(2.0/(sqrt(a/b)+sqrt(b/a)),(2.0*l+3.0)/2.0),1e-15);
            }
    }
        
}

TEST_F(BasisSet_SG,AnalyticRepulsion)
{
    using ds_t=IBS_Evaluator::ds_t;
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

