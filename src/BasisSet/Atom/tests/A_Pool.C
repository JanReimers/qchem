#include <iostream>
#include <iomanip>
#include "gtest/gtest.h"

import qchem.BasisSet.Atom.Factory;
import qchem.LASolver;
import qchem.Blaze;
import qchem.BasisSet.Internal.DB_Cache_RAM;

using std::cout;
using std::endl;

using namespace BasisSet::Atom;
using enum BasisSetAccuracy;
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;

class BasisSetPoolTests : public ::testing::Test
{
public:
    ~BasisSetPoolTests()
    {
        delete BasisSet::theGlobalCache;  
        BasisSet::theGlobalCache=new BasisSet::IntegralsCache_RAM<double>(true); 
    }
    void Orthogonality(BasisSet::Atom::Type type)
    {
        const double trunc_tol=0;
        for (size_t Z:{1,10,20,40,80,100})
        {
            cout << "---------------- Z=" << Z << " ---------------"<< endl;
            for (auto acc:{Low,Medium,High})
            {
                
                BasisSet::Real_BS* bs=Factory(acc,type,Z);
                for (auto ibs:bs->Iterate<Real_OIBS>())
                {
                    // cout << BasisSetAccuracyStrs[static_cast<size_t>(acc)] << " " << *ibs;
                    LASolver<double>* las=LASolver<double>::Factory(qchem::Cholesky,trunc_tol);
                    las->SetBasisOverlap(ibs->Overlap());
                    double smin=blazem::min(las->Get_BS_Diagonal());
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


// qchem::Ortho orthos[] = {qchem::SVD,qchem::Eigen,qchem::Cholesky};
// std::string OrthStrs[]={"Cholesky","Eigen  ","SVD    "};

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



