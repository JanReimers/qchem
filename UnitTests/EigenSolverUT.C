// file: EigenSolverUT.C  Unit test for the eigen solver

#include "gtest/gtest.h"
#include "Imp/BasisSet/Atom/l/Slater_BS.H"
#include <LAParams.H>
#include <LASolver.H>
// #include "oml/numeric/LapackSVDSolver.H"
// #include "oml/diagonalmatrix.h"
// #include "oml/numeric.h"
#include <iostream> 
#include <cassert>

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
        bs=new Atoml::Slater::BasisSet(N,0.1,10,Lmax);
        bs->Set(lap);
    }    
    
    int Lmax;
    Atoml::Slater::BasisSet* bs;
 };

 const double trunc_tol=1e-12;
 const double alg_tol=1e-15;

LAParams laps[] = { 
    {qchem::Lapack,qchem::SVD  ,trunc_tol,alg_tol},
    {qchem::Lapack,qchem::Eigen,trunc_tol,alg_tol},
    {qchem::OML   ,qchem::SVD  ,trunc_tol,alg_tol},
    {qchem::OML   ,qchem::Eigen,trunc_tol,alg_tol},
    {qchem::OML   ,qchem::Cholsky,trunc_tol,alg_tol},
    };


double Norm(const SMatrix<double>& s)
{
    return sqrt(Sum(DirectMultiply(s,s)));
}

TEST_F(OrthogonalizeTests, Types)
{
    typedef LASolver<double>::RSMat SMat;
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
                SMat I=las->Transform(S);
                SMat I1(ibs->GetNumFunctions());
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
