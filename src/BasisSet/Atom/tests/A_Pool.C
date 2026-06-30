#include <iostream>
#include <iomanip>
#include "gtest/gtest.h"

import qchem.BasisSet.Atom.Factory;
import qchem.LASolver;
import qchem.Blaze;
using namespace qchem;

using std::cout;
using std::endl;

using namespace qchem::BasisSet::Atom;
using enum BasisSetAccuracy;
using BasisSet::Real_BS;
using BasisSet::Real_OIBS;

class BasisSetPoolTests : public ::testing::Test
{
public:
    // The integrals cache is a process-wide construct-on-first-use singleton now; there is no longer a
    // per-fixture cache to tear down and recreate (the Cache4 registry self-invalidates on Register).
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
                    LASolver<double>* las=LASolver<double>::Factory(qchem::Cholesky);
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



