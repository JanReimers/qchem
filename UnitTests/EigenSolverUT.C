// file: EigenSolverUT.C  Unit test for the eigen solver

#include <iostream>
#include <iomanip> 
#include <cassert>
#include "nlohmann/json.hpp"
#include "gtest/gtest.h"

import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet;
using Real_OIBS=BasisSet::Real_OIBS;

import qchem.LASolver;

import qchem.Factory;
import qchem.Cluster;
import qchem.Hamiltonian.Factory;
import qchem.Symmetry.Spin;

import qchem.WaveFunction.Types;
import qchem.Math;
import qchem.Blaze;

using qchem::WaveFunction::bs_t;

using std::cout;
using std::endl;

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class OrthogonalizeTests : public ::testing::Test
{
public:
    OrthogonalizeTests() : bs(0)
    {
        
    }
    
    double Set(int N, int Z)
    {
        if (bs) delete bs;
        nlohmann::json js = {
        {"type",abs_t::Slater},
        {"N", N}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSet::Atom::Factory(js,Z);
        return pow(10/.1,1./N);
    }    
    Real_BS* bs;
 };

 const double trunc_tol=1e-12;
 



// typedef blazem::SymmetricMatrix < blazem::DynamicMatrix <double ,blaze::columnMajor >> rsmat_t;
// typedef blazem::DynamicMatrix <double ,blazem::columnMajor > rmat_t;
// typedef blazem::DynamicVector <double ,blazem::columnMajor > rvec_t;
typedef blazem::IdentityMatrix<double> bUnit;
typedef blazem::UpperMatrix< rmat_t > bU;
typedef blazem::LowerMatrix< rmat_t > bL;



void zeroLower(rmat_t& m)
{
    size_t N=m.rows();
    for (size_t i=0;i<N;i++)
            for (size_t j=i+1;j<N;j++)
                m(j,i)=0;
}


TEST_F(OrthogonalizeTests, BlazeDemo)
{
    
    Set(10,75);
    for (auto ibs:bs->Iterate<Real_OIBS>())
    {
        size_t N=ibs->GetNumFunctions();
        rmat_t  bM=ibs->Overlap();
        
        blazem::potrf( bM, 'U' );
        blazem::trtri( bM, 'U', 'N' );
        // blaze::potri( bM, 'U' ); gives the wrong answer!!
        
        zeroLower(bM);
      
        bU U(bM);
        bL L(blazem::trans(U));
        rsmat_t bS=ibs->Overlap();
        double err=blazem::norm(L*bS*U-bUnit(N));
        cout << "N=" << N << ", error=" << err << endl;
        cout << "------------------------------------" << endl;
        EXPECT_NEAR(err,0.0,N*N*N*1e-14);
    }
}


qchem::Ortho orthos[] = {qchem::SVD,qchem::Eigen,qchem::Cholsky};
std::string OrthStrs[]={"Cholsky","Eigen  ","SVD    "};

TEST_F(OrthogonalizeTests, Blaze)
{
    int NMax=21;
    for (int N=3;N<=NMax;N++)
    {
        Set(N,75);
        for (auto ibs:bs->Iterate<Real_OIBS>())
        {
            for (auto ortho:orthos)
            {
                rsmat_t S=ibs->Overlap();
                LASolver<double>* las=LASolver<double>::Factory(ortho,trunc_tol);
                las->SetBasisOverlap(S);
                auto I=las->Transform(S);
                bUnit I1(I.rows());
                
                double eps=1.2e-15*pow(N,3);
                cout << OrthStrs[ortho] << " " << ibs->GetSymmetry() << " " << N << " " << blazem::norm(I-I1) << " " << eps << endl;
                if (N<9)
                {
                    EXPECT_NEAR(blazem::norm(I-I1),0.0,eps);
                    // EXPECT_NEAR(Norm(I-I1),0.0,eps);
                }
                else if (N<12)
                {
                    EXPECT_NEAR(blazem::norm(I-I1),0.0,N*eps);
                    // EXPECT_NEAR(Norm(I-I1),0.0,N*eps);
                }
            }

        }
        // cout << "--------------------------------------------------" << endl;
    }
};

TEST_F(OrthogonalizeTests, BlazeHydrogen)
{
    using namespace qchem::Hamiltonian;
    typedef std::shared_ptr<const Cluster> cl_t;
    Cluster* cla=new Atom(1,0);
        Hamiltonian* Ham=qchem::Hamiltonian::Factory(Model::E1,Pol::Polarized,cl_t(cla));
    double beta=Set(21,1);
    cout << "Beta=" << beta << endl;
    for (auto ibs:bs->Iterate<Real_OIBS>())
    {
        for (auto ortho:orthos)
        {
            LASolver<double>* las=LASolver<double>::Factory(ortho,trunc_tol/10);
            las->SetBasisOverlap(ibs->Overlap());
            auto [U,e]=las->Solve(Ham->GetMatrix(ibs,Spin::Down,0));
            cout << "s=" << las->Get_BS_Diagonal() << endl;
            cout << OrthStrs[ortho] << " " << ibs->GetSymmetry() << " " << e[0]+0.5 << " " << e[1]+0.125 << endl;
            EXPECT_NEAR(e[0],-0.5  ,4e-14);
            EXPECT_NEAR(e[1],-0.125,8e-13);
        }
    }
}