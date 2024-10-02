// file: EigenSolverUT.C  Unit test for the eigen solver

#include "gtest/gtest.h"
#include "HartreeFockTester.H"
#include "Cluster/Atom.H"

#include "DFTDataBase/HeapDB/HeapDB.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"

#include "Misc/PeriodicTable.H"
#include "types.H"
#include "oml/numeric/LapackSVDSolver.H"
#include "oml/diagonalmatrix.h"
#include "oml/numeric.h"
#include <iostream> 
#include <cassert>

using std::cout;
using std::endl;

extern PeriodicTable thePeriodicTable;

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class qchem_EigenSolverTests 
    : public ::testing::TestWithParam<LinearAlgebraParams>
{
public:
    
    qchem_EigenSolverTests() 
        : itsDB(0)
        , itsBasisSet(0)
         {}
    void Init(size_t NBasis, double minexp, double maxexp, size_t L,const LinearAlgebraParams& lap)
    {
        if (itsBasisSet) 
        {
            assert(itsDB);
            delete itsDB;   
            delete itsBasisSet;
        }
        itsDB=new HeapDB<double>();
        itsBasisSet=new SphericalGaussianBS(lap,itsDB, NBasis, minexp, maxexp,L);
        
    }

    
    SMatrix<double> GetOverlap() const
    {
        return itsDB->GetOverlap();
    }
//    void OutIn(const PMStreamableObject* pout,PMStreamableObject* pin,StreamableObject::Mode);
    IntegralDataBase<double>* itsDB;
    BasisSet* itsBasisSet;
};

TEST_F(qchem_EigenSolverTests,MinSVEVTests)
{
    StreamableObject::SetToPretty();
    LinearAlgebraParams lap={qchem::Lapack,qchem::SVD,1e-6,1e-12};
    for (size_t L=0;L<=10;L++)
    {
        Init(20,.01,100.0,L,lap);    
        SMatrix<double> S=GetOverlap();
        auto [U,w]=Diagonalize(S);
        
        auto [U1,s,V]=SVD(S);
        cout << "L=" << L << "   w=" << Min(w) << " min(s)=" << Min(s) << endl;;
    }
}

LinearAlgebraParams laps[] = { 
    {qchem::Lapack,qchem::SVD  ,1e-2,1e-12},
    {qchem::Lapack,qchem::Eigen,1e-2,1e-12},
    {qchem::OML   ,qchem::SVD  ,1e-2,1e-12},
    {qchem::OML   ,qchem::Eigen,1e-2,1e-12},
    };

TEST_F(HartreeFockAtomTester, AtomsHFEigenSolvers)
{
    double eps_ro=1e-5; //Converge criterial for delta ro (charge density)
    double eps_e=2e-4; //Since we are truncating the BS we get higher errors.
    int Z=25,NBasis=20; //Mn with some D valance electrons
    int L=thePeriodicTable.GetMaxL(Z);
    double spin=thePeriodicTable.GetNumUnpairedElectrons(Z);
    Atom* atom=new Atom(Z,0,Vector3D<double>(0,0,0));
    Init(atom);
    for (auto lap:laps)
    {
        std::cout << "Testing atom " << thePeriodicTable.GetSymbol(Z)  << " " << lap << std::endl;
        Init(NBasis,L,spin,lap);
        Iterate(1.0,eps_ro,40);
        double E_HF=thePeriodicTable.GetEnergyHF(Z);
        double error=fabs((E_HF-TotalEnergy())/E_HF);
        std::cout.precision(7);
        std::cout << "E_HF relative error=" << error*100.0 << "%" << std::endl;
        EXPECT_LT(error,eps_e);        
    }
};
