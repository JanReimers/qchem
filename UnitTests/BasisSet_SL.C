// File: UnitTests/BasisSet_SL.C  unit test the Slater_IBS Evaluator
#include "gtest/gtest.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

import BasisSet.Atom.Slater.NR.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.BS_Evaluator;
import BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import BasisSet.Atom.Gaussian_BS;
import BasisSet.Atom.BSpline_IBS;
import BasisSet.Atom.BSpline_BS;
import qchem.Mesh.Integrator;
import qchem.Molecule;
import Common.Constants;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet;
import qchem.BasisSet.Atom.Slater.NR.BS;
import qchem.BasisSet.Atom.Gaussian.NR.BS;
import qchem.Orbital_HF_IBS;


bool operator==(const ERI4& a, const ERI4& b)
{
    static double eps=1e-16;
    if (a.GetLimits()!=b.GetLimits()) return false;
    for (size_t i=1;i<a.Nab();i++)
        for (size_t j=1;j<a.Nab();j++)
            if (FrobeniusNorm(a(i,j)-b(i,j))>eps) return false;
    return true;
}
//----------------------------------------------------------------------------------------
//
//  Testing common to all atom basis set evaluators
//
class BasisSet_Common : public ::testing::Test
{
public:
    BasisSet_Common(BS_Evaluator* _bseval) 
        : es{0.5,1.0,2.0}
        , N(es.size())
        , LMax(3)
        , cl(new Molecule())
        , bs_eval(_bseval)
        , bs(0)
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
        delete bs_eval;
        delete bs;
    }
    void Insert(IBS_Evaluator* eval)
    {
        evals.push_back(eval);
        bs_eval->Register(eval);
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
    BS_Evaluator* bs_eval;
    BasisSet* bs;
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
        cout.precision(2);
        // cout << S-Snum << endl;
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

    BasisSet_SL() : BasisSet_Common(new Slater_BS)
    {
        for (size_t l=0;l<=LMax;l++)
            Insert(new Slater_IBS(es,l)); 
        bs=new Atoml::Slater::BasisSet(convert(es),LMax);
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

TEST_F(BasisSet_SL,HF_ERIs)
{
    auto a=evals.begin();
    for (auto aibs:bs->Iterate<Orbital_HF_IBS<double>>())
    {
        auto c=evals.begin();
        for (auto cibs:bs->Iterate<Orbital_HF_IBS<double>>())
        {
            if (aibs->GetID()<cibs->GetID())
            {
                ERI4 J1=bs_eval->Direct(*a,*c);
                ERI4 J2=aibs->Direct(*cibs);
                EXPECT_TRUE(J1==J2);
                ERI4 K1=bs_eval->Exchange(*a,*c);
                ERI4 K2=aibs->Exchange(*cibs);
                EXPECT_TRUE(K1==K2);
            }
            ++c;
        }
        ++a;
    }
}
//----------------------------------------------------------------------------------------
//
//  Testing atom Gaussian basis set evaluators
//

class BasisSet_SG: public BasisSet_Common
{
public:

    BasisSet_SG() : BasisSet_Common(new Gaussian_BS)
    {
        for (size_t l=0;l<=LMax;l++)
            Insert(new Gaussian_IBS(es,l));    
        bs=new Atoml::Gaussian::BasisSet(convert(es),LMax);
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
TEST_F(BasisSet_SG,HF_ERIs)
{
    auto a=evals.begin();
    for (auto aibs:bs->Iterate<Orbital_HF_IBS<double>>())
    {
        auto c=evals.begin();
        for (auto cibs:bs->Iterate<Orbital_HF_IBS<double>>())
        {
            if (aibs->GetID()<cibs->GetID())
            {
                ERI4 J1=bs_eval->Direct(*a,*c);
                ERI4 J2=aibs->Direct(*cibs);
                EXPECT_TRUE(J1==J2);
                ERI4 K1=bs_eval->Exchange(*a,*c);
                ERI4 K2=aibs->Exchange(*cibs);
                EXPECT_TRUE(K1==K2);
            }
            ++c;
        }
        ++a;
    }
}

//----------------------------------------------------------------------------------------
//
//  Testing atom Gaussian basis set evaluators
//

class BasisSet_BS: public BasisSet_Common
{
public:

    BasisSet_BS() : BasisSet_Common(new BSpline_BS<6>)
    {
        for (size_t l=0;l<=0;l++)
            Insert(new BSpline_IBS<6>(9+2*l,0.01,20.0,l));    
    }
   
};

TEST_F(BasisSet_BS,Overlap) {TestOverlap(4e-8);}
TEST_F(BasisSet_BS,Grad2  ) {TestGrad2  (7e-3);}
TEST_F(BasisSet_BS,Inv_r1 ) {TestInv_r1 (4e-6);}
// TEST_F(BasisSet_BS,Inv_r2 ) {TestInv_r2 (4e-8);}
TEST_F(BasisSet_BS,Charge ) {TestCharge (2e-5);}
