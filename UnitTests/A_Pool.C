#include <iostream>
#include <iomanip>
#include "BasisSetPool.H"
#include "gtest/gtest.h"
#include "QchemTester.H"

import qchem.Hamiltonian.Factory;
import qchem.Factory;
import qchem.Cluster;

import qchem.LASolver;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;

class BasisSetPoolTests : public ::testing::Test
{
public:
};

TEST_F(BasisSetPoolTests,SlaterDisplay)
{
    for (size_t Z:{1,2,10,20,40,80,100})
    {
        cout << "---------------- Z=" << Z << " ---------------"<< endl;
        for (auto acc:{Low,Medium,High,Extreme})
        {
            rvec_t es=SlaterExponents(acc,Z);
            cout << BasisSetAccuracyStrs[static_cast<size_t>(acc)] << " " 
                 << "N=" << es.size() << " {" << es[0] << "," << es[1] << " ... " << es[es.size()-1] << "}" << endl;

        }
    }
}

// qchem::Ortho orthos[] = {qchem::SVD,qchem::Eigen,qchem::Cholsky};
// std::string OrthStrs[]={"Cholsky","Eigen  ","SVD    "};

TEST_F(BasisSetPoolTests,SlaterOrthonoality)
{
    const double trunc_tol=0;
    using enum BasisSet::Atom::Type;
    for (size_t Z:{1,10,20,40,80,100})
    {
        cout << "---------------- Z=" << Z << " ---------------"<< endl;
        for (auto acc:{Low,Medium,High,Extreme})
        {
            
            BasisSet::Real_BS* bs=Factory(acc,Slater,Z);
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

TEST_P(A_HF_U,Slater_HFGroundStates_Extreme)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    BasisSet::Real_BS* bs=Factory(Extreme,BasisSet::Atom::Type::Slater,Z); 
    QchemTester::Init(1e-3,bs);
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   32     ,Z*1e-5    ,1e-7 , 1e-6      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),1e-6); 
    EXPECT_TRUE(Converged()); 
        
}

INSTANTIATE_TEST_SUITE_P(SlaterHFGroundStates,A_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 
// INSTANTIATE_TEST_SUITE_P(SlaterHFGroundStates,A_HF_U,::testing::Values(2,70));//)); 
