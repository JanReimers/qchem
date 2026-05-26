#include <iostream>
#include <iomanip>
#include <blaze/Math.h>
#include "gtest/gtest.h"

import qchem.Unittests.QchemTester;
import qchem.LASolver;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;

class BasisSetPoolTests : public ::testing::Test
{
public:
    void Orthogonality(BasisSet::Atom::Type type)
    {
        const double trunc_tol=0;
        for (size_t Z:{1,10,20,40,80,100})
        {
            cout << "---------------- Z=" << Z << " ---------------"<< endl;
            for (auto acc:{Low,Medium,High})
            {
                
                BasisSet::Real_BS* bs=PoolFactory(acc,type,Z);
                for (auto ibs:bs->Iterate<Real_OIBS>())
                {
                    // cout << BasisSetAccuracyStrs[static_cast<size_t>(acc)] << " " << *ibs;
                    LASolver<double>* las=LASolver<double>::Factory(qchem::Cholsky,trunc_tol);
                    las->SetBasisOverlap(ibs->Overlap());
                    double smin=blaze::min(las->Get_BS_Diagonal());
                    cout << BasisSetAccuracyStrs[static_cast<size_t>(acc)] << " " 
                        << std::setprecision(2) << std::setw(6) << smin << "  " << *ibs;
                    EXPECT_GE(smin,1e-13);
                    delete las;
                }
                delete bs;
            }
        }
    }
};


// qchem::Ortho orthos[] = {qchem::SVD,qchem::Eigen,qchem::Cholsky};
// std::string OrthStrs[]={"Cholsky","Eigen  ","SVD    "};

TEST_F(BasisSetPoolTests,SlaterOrthogonality)
{
    Orthogonality(BasisSet::Atom::Type::Slater);
}
TEST_F(BasisSetPoolTests,GaussianOrthogonality)
{
    Orthogonality(BasisSet::Atom::Type::Gaussian);
}
TEST_F(BasisSetPoolTests,BSplineOrthogonality)
{
    Orthogonality(BasisSet::Atom::Type::BSpline6);
}



using namespace qchem::Hamiltonian;
class A_HF_U : public virtual QchemTester, public ::testing::TestWithParam<size_t>, TestAtom
{
public:
    A_HF_U() : TestAtom(GetParam()) {};
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::HF,Pol::UnPolarized,cluster);
    }
};

class BS_U_High : public A_HF_U {};
TEST_P(BS_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::BSpline6);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-7    ,1e-7 , 2.5e-13      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,1e-9);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_High,::testing::Values(70,2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 

class BS_U_Medium : public A_HF_U {};
TEST_P(BS_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::BSpline6);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   30     ,Z*1e-7    ,1e-7 , 2.5e-7      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,1e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Medium,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 
class BS_U_Low : public A_HF_U {};
TEST_P(BS_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Low,BasisSet::Atom::Type::BSpline6);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   30     ,Z*1e-7    ,1e-7 , 5e-5      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Low,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 


class SG_U_High : public A_HF_U {};
TEST_P(SG_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Gaussian);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-5    ,1e-7 , 1e-5      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),2e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,SG_U_High,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 
class SG_U_Medium : public A_HF_U {};
TEST_P(SG_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Gaussian);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-5    ,1e-7 , 5e-2      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),2e-4); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,SG_U_Medium,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 

class SL_U_High : public A_HF_U {};
TEST_P(SL_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Slater);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   32     ,Z*1e-5    ,1e-7 , 1e-6      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),1e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_High,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 

class SL_U_Medium : public A_HF_U {};
TEST_P(SL_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Slater);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   22     ,Z*1e-4    ,1e-5 , 5e-4      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),20e-6); 
    EXPECT_TRUE(Converged()); 
        
}

INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Medium,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 

class SL_U_Low : public A_HF_U {};
TEST_P(SL_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Low,BasisSet::Atom::Type::Slater);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   30     ,Z*1e-4    ,1e-4 , 5e-1      ,Z*2e-5 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),0.01); //1% 
    EXPECT_TRUE(Converged()); 
        
}

INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Low,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 


// INSTANTIATE_TEST_SUITE_P(SlaterHFGroundStates,A_HF_U,::testing::Values(2,70));//)); 
