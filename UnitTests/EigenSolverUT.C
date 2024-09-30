// file: EigenSolverUT.C  Unit test for the eigen solver

//#include "AtomTester.H"
//#include "HartreeFockTester.H"
#include "gtest/gtest.h"
//#include "oml/imp/stream.h"
//#include "Cluster/Atom.H"
//#include "Cluster/Molecule.H"
//#include "Hamiltonian/Hamiltonian.H"
//#include "BasisSet/BasisGroup.H"
//#include "WaveFunction/WaveFunction.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"

#include "Misc/PeriodicTable.H"
#include "types.H"
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
class qchem_EigenSolverTests : public ::testing::TestWithParam<int>
{
public:
    qchem_EigenSolverTests() 
        : itsDB(0)
        , itsBasisSet(0)
         {}
    void Init(size_t NBasis, double minexp, double maxexp, size_t L)
    {
        if (itsBasisSet) 
        {
            assert(itsDB);
            delete itsDB;   
            delete itsBasisSet;
        }
        itsDB=new HeapDB<double>();
        itsBasisSet=new SphericalGaussianBS(itsDB, NBasis, minexp, maxexp,L);
        
    }
    SMatrix<double> GetOverlap() const
    {
        return itsDB->GetOverlap();
    }
//    void OutIn(const PMStreamableObject* pout,PMStreamableObject* pin,StreamableObject::Mode);
    
    IntegralDataBase<double>* itsDB;
    BasisSet* itsBasisSet;
};

TEST_F(qchem_EigenSolverTests,EigenSTests)
{
    StreamableObject::SetToPretty();
    for (size_t L=0;L<=10;L++)
    {
        Init(5,.01,100.0,L);
        SMatrix<double> S=GetOverlap();
        auto [U,w]=Diagonalize(S);
        
        auto [U1,s,V]=SVD(S);
        cout << "L=" << L << endl << "   w=" <<w << endl << "   s=" << s << endl;;
        //Matrix<double> w1=Transpose(U)*S*U;
    }
}
