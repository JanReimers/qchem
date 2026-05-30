// File: UnitTests/BasisSet_Atom.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <blaze/Math.h>
using std::cout;
using std::endl;

import qchem.Unittests.BasisSetPool;
import qchem.BasisSet.Atom.Evaluators.Slater.IBS;
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.Concepts;

import qchem.Mesh.Integrator;
import qchem.Symmetry.Yl;
import qchem.Hamiltonian.Types;

using qchem::Hamiltonian::ohfbs_t;
using qchem::Hamiltonian::obs_t;
using namespace BasisSet::Atom::Evaluators;
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;
//----------------------------------------------------------------------------------------
//
//  Testing common to all atom basis set evaluators
//
template <isHF_NR_Evaluator E> class BasisSet_Common : public ::testing::Test
{
public:
    BasisSet_Common() 
        : es{0.5,1.0,2.0}
        , N(es.size())
        , LMax(3)
        , cl(new Atom(1,0.0,Vector3D(0,0,0)))
        , bs(0)
    {
        
        MeshParams mp({qchem::MHL,500,3,2.0,qchem::Gauss,1,0,0,3});
        mintegrator=new MeshIntegrator<double>(cl->CreateMesh(mp));
    }
    ~BasisSet_Common()
    {
        for (auto ev:evals) delete ev;
        delete cl;
        delete mintegrator;
        delete bs;
    }
    void Insert(E* eval)
    {
        evals.push_back(eval);
    }

    void TestOverlap(double eps) const;
    void TestGrad2  (double eps) const;
    void TestInv_r1 (double eps) const;
    void TestInv_r2 (double eps) const;
    // void TestCharge (double eps) const;


    rvec_t es;
    size_t N,LMax;
    std::vector<E*> evals;
    Cluster* cl;
    MeshIntegrator<double>* mintegrator;
    Real_BS* bs;
};
template <isHF_NR_Evaluator E> void BasisSet_Common<E>::TestOverlap(double eps) const
{
    EXPECT_GT(evals.size(),0);
    size_t index=0;
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Overlap(i,j);
        for (auto d:diagonal(S)) EXPECT_NEAR(d,1.0,1e-15);
        rsmat_t Snum = mintegrator->Overlap(*ev);
        cout.precision(2);
        // cout << S-Snum << endl;
        cout << "l=" << index++ << std::endl;
        EXPECT_NEAR(max(abs(S-Snum)),0.0,eps);
    }
}
template <isHF_NR_Evaluator E> void BasisSet_Common<E>::TestGrad2  (double eps) const
{
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Grad2(i,j);
        rsmat_t Snum = mintegrator->Grad2(*ev);
        EXPECT_NEAR(max(abs(S-Snum)),0.0,eps);
    }
        
}
template <isHF_NR_Evaluator E> void BasisSet_Common<E>::TestInv_r1 (double eps) const
{
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Inv_r1(i,j);
        rsmat_t Snum = mintegrator->Inv_r1(*ev);
        EXPECT_NEAR(max(abs(S-Snum)),0.0,eps);
    }
        
}
template <isHF_NR_Evaluator E> void BasisSet_Common<E>::TestInv_r2 (double eps) const
{
    for (auto ev:evals)
    {
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Inv_r2(i,j);
        rsmat_t Snum = mintegrator->Inv_r2(*ev);
        EXPECT_NEAR(max(abs(S-Snum)),0.0,eps);
    }
        
}
// template <isFit_Evaluator E> void BasisSet_Common<E>::TestCharge (double eps) const
// {
//     for (auto ev:evals)
//     {
//         size_t N=ev->size();
//         rvec_t S(N);
//         for (auto i:iv_t(0,N))
//                 S[i]= ev->Charge(i);
//         rvec_t Snum = mintegrator->Integrate(*ev);
//         EXPECT_NEAR(max(abs(S-Snum)),0.0,eps);
//     }
        
// }

//----------------------------------------------------------------------------------------
//
//  Testing atom Slater basis set evaluators
//

class BasisSet_SL: public BasisSet_Common<Slater::Slater_IBS_Evaluator>
{
public:

    BasisSet_SL() : BasisSet_Common<Slater::Slater_IBS_Evaluator>()
    {
        Atom_EC ec(86); //Radon has f orbtials with no magnetic splitting.
        for (auto ir:ec.GetIrreps())
            Insert(new Slater::Slater_IBS_Evaluator(es,ir)); 
        bs=PoolFactory(BasisSetAccuracy::N3,BasisSet::Atom::Type::Slater,86);
        cout << es << endl << *bs << endl;
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
TEST_F(BasisSet_SL,Grad2  ) {TestGrad2  (3e-15);}
TEST_F(BasisSet_SL,Inv_r1 ) {TestInv_r1 (3e-14);}
TEST_F(BasisSet_SL,Inv_r2 ) {TestInv_r2 (4e-9);}
// TEST_F(BasisSet_SL,Charge ) {TestCharge (9e-14);}

TEST_F(BasisSet_SL,AnalyticOverlap)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Overlap(i,j);
        // cout << S << endl;
        for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                EXPECT_NEAR(S(i,j),pow(2.0/(sqrt(a/b)+sqrt(b/a)),2*l+3),1e-15);
            }
    }
        
}
TEST_F(BasisSet_SL,AnalyticRepulsion)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        if (l>=3) continue;
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Repulsion(i,j);
        // cout << "l=" << l << " S=" << S << endl;
        rvec_t   ns=ev->Norm();
        for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                double rerr=(S(i,j) - R0(a,b,l,l)*ns[i]*ns[j])/S(i,j);
                EXPECT_NEAR(rerr,0.0,3e-14);
            }
    }
        
}

TEST_F(BasisSet_SL,HF_ERIs)
{
    auto a=evals.begin();
    for (auto aibs:bs->Iterate<ohfbs_t>())
    {
        auto c=evals.begin();
        for (auto cibs:bs->Iterate<ohfbs_t>())
        {
            ERI4 J2=aibs->Direct(*cibs);
            ERI4 K2=aibs->Exchange(*cibs);
            ++c;
        }
        ++a;
    }
}

std::string angularIDs[]={"0 {}","1 {}","2 {}","3 {}"};
TEST_F(BasisSet_SL,IDs)
{
    size_t index=0;
    for (auto ibs:bs->Iterate<obs_t>())
    {
        EXPECT_EQ(ibs->Name(),"SL ");
        EXPECT_EQ(ibs->RadialID(),"SL N=3 {0.5 ... 2}");
        EXPECT_EQ(ibs->AngularID(),angularIDs[index++]);
    }
}
//----------------------------------------------------------------------------------------
//
//  Testing atom Gaussian basis set evaluators
//

class BasisSet_SG: public BasisSet_Common<Gaussian::Gaussian_IBS_Evaluator>
{
public:

    BasisSet_SG() : BasisSet_Common()
    {
        // for (size_t l=0;l<=LMax;l++)
        //     Insert(new Gaussian_IBS_Evaluator(es,l));    
        // bs=new AtomBS::Gaussian::BasisSet(es,LMax);

        Atom_EC ec(86); //Radon has f orbtials with no magnetic splitting.
        for (auto ir:ec.GetIrreps())
            Insert(new Gaussian::Gaussian_IBS_Evaluator(es,ir)); 
        bs=PoolFactory(BasisSetAccuracy::N3,BasisSet::Atom::Type::Gaussian,86);
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
TEST_F(BasisSet_SG,Grad2  ) {TestGrad2  (6e-15);}
TEST_F(BasisSet_SG,Inv_r1 ) {TestInv_r1 (3e-14);}
TEST_F(BasisSet_SG,Inv_r2 ) {TestInv_r2 (4e-9);}
// TEST_F(BasisSet_SG,Charge ) {TestCharge (5e-14);}

TEST_F(BasisSet_SG,AnalyticOverlap)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Overlap(i,j);

         for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                EXPECT_NEAR(S(i,j),pow(2.0/(sqrt(a/b)+sqrt(b/a)),(2.0*l+3.0)/2.0),1e-15);
            }
    }
        
}
TEST_F(BasisSet_SG,AnalyticRepulsion)
{
    for (auto ev:evals)
    {
        int l=ev->Getl();
        if (l>=3) continue;
        size_t N=ev->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= ev->Repulsion(i,j);
        // cout << "l=" << l << " S=" << S << endl;
        rvec_t   ns=ev->Norm();
         for (auto i:iv_t(0,S.rows()))
            for (auto j:iv_t(i,S.rows()))
            {
                double a=es[i],b=es[j];
                double rerr=(S(i,j) - R0(a,b,l,l)*ns[i]*ns[j])/S(i,j);
                EXPECT_NEAR(rerr,0.0,3e-14);
            }
    }
        
}
TEST_F(BasisSet_SG,HF_ERIs)
{
    auto a=evals.begin();
    for (auto aibs:bs->Iterate<ohfbs_t>())
    {
        auto c=evals.begin();
        for (auto cibs:bs->Iterate<ohfbs_t>())
        {
            ERI4 J2=aibs->Direct(*cibs);
            ERI4 K2=aibs->Exchange(*cibs);
            ++c;
        }
        ++a;
    }
}

TEST_F(BasisSet_SG,IDs)
{
    size_t index=0;
    for (auto ibs:bs->Iterate<obs_t>())
    {
        EXPECT_EQ(ibs->Name(),"SG ");
        EXPECT_EQ(ibs->RadialID(),"SG N=3 {0.5 ... 2}");
        EXPECT_EQ(ibs->AngularID(),angularIDs[index++]);
    }
}

//----------------------------------------------------------------------------------------
//
//  Testing atom BSpline basis set evaluators
//

class BasisSet_BS: public BasisSet_Common<BSpline::Evaluator<6>>
{
public:

    BasisSet_BS() : BasisSet_Common()
    {
        for (size_t l=0;l<=3;l++)
            Insert(new BSpline::Evaluator<6>(5,0.01,20.0,Irrep_QNs::sym_t(new Yl_Sym(l))));

        bs=PoolFactory(BasisSetAccuracy::N5,BasisSet::Atom::Type::BSpline6,86);
    }
   
};

TEST_F(BasisSet_BS,Overlap) {TestOverlap(3e-11);}
TEST_F(BasisSet_BS,Grad2  ) {TestGrad2  (1e-4);}
TEST_F(BasisSet_BS,Inv_r1 ) {TestInv_r1 (5e-9);}
TEST_F(BasisSet_BS,Inv_r2 ) {TestInv_r2 (5e-3);}
// TEST_F(BasisSet_BS,Charge ) {TestCharge (8e-7);}

TEST_F(BasisSet_BS,IDs)
{
    size_t index=0;
    for (auto ibs:bs->Iterate<obs_t>())
    {
        // cout << "index=" << index << " ibs=" << *ibs << endl;
        EXPECT_EQ(ibs->Name(),"BSpline<6>");
        EXPECT_EQ(ibs->RadialID(),"BSpline<6> grid: N=6 {0,0.25,0.5 ... 4}");
        EXPECT_EQ(ibs->AngularID(),angularIDs[index++]);
    }
}
