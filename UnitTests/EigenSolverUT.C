// file: EigenSolverUT.C  Unit test for the eigen solver

#include <iostream>
#include <iomanip> 
#include <cassert>
#include <cmath>
#include "nlohmann/json.hpp"
#include "gtest/gtest.h"
#include "blaze/Math.h" 

import qchem.LAParams;
import qchem.LASolver;

import qchem.Factory;
import qchem.IrrepBasisSet;
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
    OrthogonalizeTests() : bs(0), Z(75)
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
        bs=BasisSetAtom::Factory(js,Z);
        bs->Set(lap);
    }
    void Set(int N)
    {
        if (bs) delete bs;
        nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSetAtom::Factory(js,Z);
    }    
    BasisSet* bs;
    int Z;
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
            for (auto ibs:bs->Iterate<Real_OIBS>())
            {
                LASolver<double>* las=ibs->CreateSolver();
                const SMat& S=ibs->Overlap();
                SMatrix<double> I=las->Transform(S);
                SMatrix<double> I1(ibs->GetNumFunctions());
                Unit(I1);
                double eps=1.2e-15*pow(N,3);
                // cout << N << " " << Max(fabs(I-I1)) << " " << Norm(I-I1) << " " << eps << endl;
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
        // cout << "--------------------------------------------------" << endl;
    }
};

typedef blaze::SymmetricMatrix < blaze::DynamicMatrix <double ,blaze::columnMajor >> bSMat;
typedef blaze::DynamicMatrix <double ,blaze::columnMajor > bMat;
typedef blaze::IdentityMatrix<double,blaze::columnMajor> bUnit;
typedef blaze::UpperMatrix< bMat > bU;
typedef blaze::LowerMatrix< bMat > bL;


bSMat to_bSMat(const SMatrix<double>& S)
{
    size_t N=S.GetNumRows();
    bSMat bS(N);
    for (auto i:S.rows())
            for (auto j:S.cols(i))
                bS(i-1,j-1)=S(i,j);
    return bS;
}
bMat to_bMat(const SMatrix<double>& S)
{
    size_t N=S.GetNumRows();
    bMat bM(N,N);
    for (auto i:S.rows())
            for (auto j:S.cols(i))
                bM(j-1,i-1)=bM(i-1,j-1)=S(i,j);
    return bM;
}

void zeroLower(bMat& m)
{
    size_t N=m.rows();
    for (size_t i=0;i<N;i++)
            for (size_t j=i+1;j<N;j++)
                m(j,i)=0;
}


TEST_F(OrthogonalizeTests, Blaze)
{
    
    Set(10);
    for (auto ibs:bs->Iterate<Real_OIBS>())
    {
        size_t N=ibs->GetNumFunctions();
        const SMatrix<double>& S=ibs->Overlap();
        bMat  bM=to_bMat(S);
        
        blaze::potrf( bM, 'U' );
        blaze::trtri( bM, 'U', 'N' );
        // blaze::potri( bM, 'U' ); gives the wrong answer!!
        
        zeroLower(bM);
      
        bU U(bM);
        bL L(blaze::trans(U));
        bSMat bS=to_bSMat(S);
        double err=blaze::norm(L*bS*U-bUnit(N));
        cout << "N=" << N << ", error=" << err << endl;
        cout << "------------------------------------" << endl;
        EXPECT_NEAR(err,0.0,N*N*N*1e-14);
    }
}