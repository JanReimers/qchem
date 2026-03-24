// file: EigenSolverUT.C  Unit test for the eigen solver

#include <iostream>
#include <iomanip> 
#include <cassert>
#include <cmath>
#include "nlohmann/json.hpp"
#include "gtest/gtest.h"
#include "blaze/Math.h" 

import qchem.LAParams;
import qchem.LASolver_blaze;

import qchem.Factory;
import qchem.IrrepBasisSet;
import qchem.BasisSet;
import qchem.Molecule;
import qchem.Hamiltonian.Factory;
import qchem.Symmetry.Spin;
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
    OrthogonalizeTests() : bs(0)
    {
        StreamableObject::SetToPretty();
    }
    void Set(int N, int Z, LAParams lap)
    {
        Set(N,Z);
        bs->Set(lap);
    }
    void Set(int N, int Z)
    {
        if (bs) delete bs;
        nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.1}, {"emax", 10.0},
        };
        bs=BasisSetAtom::Factory(js,Z);
    }    
    BasisSet* bs;
 };

 const double trunc_tol=1e-12;
 


double Norm(const SMatrix<double>& s)
{
    return sqrt(Sum(DirectMultiply(s,s)));
}


typedef blaze::SymmetricMatrix < blaze::DynamicMatrix <double ,blaze::columnMajor >> bSMat;
typedef blaze::DynamicMatrix <double ,blaze::columnMajor > bMat;
typedef blaze::DynamicVector <double ,blaze::columnMajor > bVec;
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


TEST_F(OrthogonalizeTests, BlazeDemo)
{
    
    Set(10,75);
    for (auto ibs:bs->Iterate<Real_OIBS>())
    {
        size_t N=ibs->GetNumFunctions();
        bMat  bM=ibs->Overlap();
        
        blaze::potrf( bM, 'U' );
        blaze::trtri( bM, 'U', 'N' );
        // blaze::potri( bM, 'U' ); gives the wrong answer!!
        
        zeroLower(bM);
      
        bU U(bM);
        bL L(blaze::trans(U));
        bSMat bS=ibs->Overlap();
        double err=blaze::norm(L*bS*U-bUnit(N));
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
                bSMat S=ibs->Overlap();
                LASolver_blaze<double>* las=LASolver_blaze<double>::Factory(ortho,trunc_tol);
                las->SetBasisOverlap(S);
                auto I=las->Transform(S);
                bUnit I1(I.rows());
                
                double eps=1.2e-15*pow(N,3);
                cout << OrthStrs[ortho] << " " << *ibs->GetSymmetry() << " " << N << " " << blaze::norm(I-I1) << " " << eps << endl;
                if (N<9)
                {
                    EXPECT_NEAR(blaze::norm(I-I1),0.0,eps);
                    // EXPECT_NEAR(Norm(I-I1),0.0,eps);
                }
                else if (N<12)
                {
                    EXPECT_NEAR(blaze::norm(I-I1),0.0,N*eps);
                    // EXPECT_NEAR(Norm(I-I1),0.0,N*eps);
                }
            }

        }
        // cout << "--------------------------------------------------" << endl;
    }
};

TEST_F(OrthogonalizeTests, BlazeHydrogen)
{
    typedef std::shared_ptr<const Cluster> cl_t;
    Cluster* cla=new Molecule();
    cla->Insert(new Atom(1,0));
    // cl_t cl(cla);
    Hamiltonian* Ham=HamiltonianF::Factory(HamiltonianF::Model::E1,HamiltonianF::Pol::Polarized,cl_t(cla));
    Set(21,1);
    for (auto ibs:bs->Iterate<Real_OIBS>())
    {
        for (auto ortho:orthos)
        {
            LASolver_blaze<double>* las=LASolver_blaze<double>::Factory(ortho,trunc_tol/10);
            las->SetBasisOverlap(ibs->Overlap());
            auto [U,e]=las->Solve(to_bSMat(Ham->GetMatrix(ibs,Spin::Down,0)));
            cout << OrthStrs[ortho] << " " << *ibs->GetSymmetry() << " " << e[0]+0.5 << " " << e[1]+0.125 << endl;
            EXPECT_NEAR(e[0],-0.5  ,4e-14);
            EXPECT_NEAR(e[1],-0.125,8e-13);
        }
    }
}