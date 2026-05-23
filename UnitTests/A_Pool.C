#include <iostream>
#include <iomanip>
#include "BasisSetPool.H"
// #include <blaze/Math.h>
#include "gtest/gtest.h"
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
    for (size_t Z:{1,10,20,40,80,100})
    {
        cout << "---------------- Z=" << Z << " ---------------"<< endl;
        for (auto acc:{Low,Medium,High,Extreme})
        {
            rvec_t es=SlaterExponents(acc,Z);
            cout << "N=" << es.size() << " {" << es[0] << "," << es[1] << " ... " << es[es.size()-1] << "}" << endl;

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