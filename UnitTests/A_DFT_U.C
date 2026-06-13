// File A_DFT_U.C  Atom DFT tests for Unpolarized (closed shell) atoms.  Using libxc for exchange functionals
#include "gtest/gtest.h"
import qchem.Unittests.QchemTester;
import qchem.Hamiltonian.Internal.Libxc_LDA_Exchange;
const bool verbose=true;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;
using namespace qchem::Hamiltonian;
class A_DFT_U : public ::testing::TestWithParam<size_t>, public TestAtom
{
public:
    A_DFT_U() : TestAtom(GetParam()) {};
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        ex=new Libxc_LDA_Exchange(7,Spin::None,GetZ());
        cout << *ex << endl;
        Hamiltonian* H=Factory(Pol::UnPolarized,cluster,ex,GetMeshParams(),itsBasisSet);
        // cout << *H << endl;
        return H;
    }
private:
    mutable ExFunctional* ex;
};


// "2": {
//     "Etot": -2.834836,
//     "Ekin": 2.767922,
//     "Ecoul": 1.99612,
//     "Eenuc": -6.625564,
//     "Exc": -0.973314,
//     "1s": -0.570425,
//     "Symbol": "He",
//     "Z": 2
//   },
// "36": {
//     "Etot": -2750.14794,
//     "Ekin": 2747.813142,
//     "Ecoul": 1171.723688,
//     "Eenuc": -6577.865761,
//     "Exc": -91.819009,
//     "1s": -509.982989,
//     "2s": -66.285953,
//     "2p": -60.017328,
//     "3s": -9.315192,
//     "3p": -7.086634,
//     "3d": -3.074109,
//     "4s": -0.820574,
//     "4p": -0.34634,
//     "Symbol": "Kr",
//     "Z": 36
//   },
class SG_DFT_U_High : public A_DFT_U {};
TEST_P(SG_DFT_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Gaussian,verbose);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-5    ,1e-5 , 1e-1      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeDFTError(),2e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_DFT,SG_DFT_U_High,::testing::Values(36));//)); 

class SL_DFT_U_High : public A_DFT_U {};
TEST_P(SL_DFT_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Slater,verbose);
    //       NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-5    ,1e-7 , 2e-2      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeDFTError(),2e-6); 
    EXPECT_TRUE(Converged()); 
        
}
// INSTANTIATE_TEST_SUITE_P(A_DFT,SL_DFT_U_High,::testing::Values(2,4,10,18,36,54));//)); 
INSTANTIATE_TEST_SUITE_P(A_DFT,SL_DFT_U_High,::testing::Values(36));//)); 