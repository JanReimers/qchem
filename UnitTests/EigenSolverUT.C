// file: EigenSolverUT.C  Unit test for the eigen solver

#include <iostream>
#include <iomanip> 
#include <cassert>
#include <cmath>
#include "nlohmann/json.hpp"
#include "gtest/gtest.h"
import qchem.LAParams;
import qchem.LASolver;

import qchem.Factory;
import qchem.Irrep_BS;
import qchem.BasisSet;
import oml;
using std::cout;
using std::endl;

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class OrthogonalizeTests : public ::testing::Test
{
public:
    OrthogonalizeTests()
    : Lmax(0) //Higher L gets easier to orthogonalize
    , bs(0)
    {
        StreamableObject::SetToPretty();
    }
    void Set(int N, LAParams lap)
    {
        if (bs) delete bs;
        nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSetAtom::Factory(js,75);
        bs->Set(lap);
    }    
    
    int Lmax;
    BasisSet* bs;
 };

 const double trunc_tol=1e-12;
 const double alg_tol=1e-15;

LAParams laps[] = { 
    {qchem::Lapack,qchem::SVD    ,trunc_tol,alg_tol},
    {qchem::OML   ,qchem::SVD    ,trunc_tol,alg_tol},
    {qchem::Lapack,qchem::Eigen  ,trunc_tol,alg_tol},
    {qchem::OML   ,qchem::Eigen  ,trunc_tol,alg_tol},
    {qchem::Lapack,qchem::Cholsky,trunc_tol,alg_tol},
    {qchem::OML   ,qchem::Cholsky,trunc_tol,alg_tol},
    };


double Norm(const SMatrix<double>& s)
{
    return sqrt(Sum(DirectMultiply(s,s)));
}

TEST_F(OrthogonalizeTests, Types)
{
    typedef SMatrix<double> SMat;
    int NMax=21;
    for (int N=3;N<=NMax;N++)
    {
        for (auto lap:laps)
        {
            Set(N,lap);
            for (auto ibs:bs->Iterate<TOrbital_IBS<double>>())
            {
                LASolver<double>* las=ibs->CreateSolver();
                const SMat& S=ibs->Overlap();
                SMatrix<double> I=las->Transform(S);
                SMatrix<double> I1(ibs->GetNumFunctions());
                Unit(I1);
                double eps=1e-15*pow(N,3);
                cout << N << " " << Max(fabs(I-I1)) << " " << Norm(I-I1) << " " << eps << endl;
                if (N<9)
                {
                    EXPECT_NEAR(Max(fabs(I-I1)),0.0,eps);
                    EXPECT_NEAR(Norm(I-I1),0.0,eps);
                }
                else if (N<12)
                {
                    EXPECT_NEAR(Max(fabs(I-I1)),0.0,N*eps);
                    EXPECT_NEAR(Norm(I-I1),0.0,N*eps);
                }
            }

        }
        cout << "--------------------------------------------------" << endl;
    }
};
