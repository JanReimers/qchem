// File A_HF_U.C  Atom Hartree-Fock tests for Unpolarized (closed shell) atoms.
#include "gtest/gtest.h"
import qchem.Unittests.QchemTester;
using namespace qchem;

const bool verbose=true;

using std::cout;
using std::endl;
using namespace qchem::BasisSet::Atom;
using enum BasisSetAccuracy;
using namespace qchem::Hamiltonian;
class A_HF_U : public ::testing::TestWithParam<size_t>, public TestAtom
{
public:
    A_HF_U() : TestAtom(GetParam()) {};
    virtual qchem::Hamiltonian::Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF,Pol::UnPolarized,structure);
    }
};

#ifdef DEBUG
#define LOW
#else
#define MEDIUM
#define LOW
#endif


#ifdef HIGH
class BS_U_High : public A_HF_U {};
TEST_P(BS_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::BSpline6,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-13, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,1e-9);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged());        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_High,::testing::Values(2,88));
// 2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88
class BSr_U_High : public A_HF_U {};
TEST_P(BSr_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::BSpliner6,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-13, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,1e-9);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged());        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BSr_U_High,::testing::Values(2,88));
class SG_U_High : public A_HF_U {};
TEST_P(SG_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Gaussian,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),2e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,SG_U_High,::testing::Values(2,36));//)); 

class SL_U_High : public A_HF_U {};
TEST_P(SL_U_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Slater,verbose);
    Iterate({.NMaxIter = 32, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-6, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),1e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_High,::testing::Values(2,88));//)); 

#endif //HIGH

#ifdef MEDIUM

class BS_U_Medium : public A_HF_U {};
TEST_P(BS_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::BSpline6,verbose);
    Iterate({.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,1e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Medium,::testing::Values(2,4));//)); 
class BSr_U_Medium : public A_HF_U {};
TEST_P(BSr_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::BSpliner6,verbose);
    Iterate({.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,1e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BSr_U_Medium,::testing::Values(2,4));//)); 

class SG_U_Medium : public A_HF_U {};
TEST_P(SG_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Gaussian,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 5e-2, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),2e-4); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,SG_U_Medium,::testing::Values(2,4));//)); 

class SL_U_Medium : public A_HF_U {};
TEST_P(SL_U_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Slater,verbose);
    Iterate({.NMaxIter = 22, .MinΔρ = Z*1e-4, .MinΔFD = 1e-5, .MinVirial = 5e-4, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),20e-6); 
    EXPECT_TRUE(Converged()); 
        
}

INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Medium,::testing::Values(2,88));//)); 

#endif //MEDIUM

#ifdef LOW

class BS_U_Low : public A_HF_U {};
class BSr_U_Low : public A_HF_U {};
#ifdef DEBUG
TEST_P(BS_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Low,BasisSet::Atom::Type::BSpline6,verbose);
    Iterate({.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 5e-5, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Low,::testing::Values(2,4));//));
TEST_P(BSr_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Low,BasisSet::Atom::Type::BSpliner6,verbose);
    Iterate({.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 5e-5, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BSr_U_Low,::testing::Values(2,4));//)); 
#else
TEST_P(BS_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Low,BasisSet::Atom::Type::BSpline6,verbose);
    Iterate({.NMaxIter = 30, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 5e-5, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    double Eerr=RelativeHFError();
    EXPECT_LT(Eerr,40e-6);
    EXPECT_GT(Eerr,-1e-4);
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,BS_U_Low,::testing::Values(2,4));//)); 
#endif

class SL_U_Low : public A_HF_U {};
TEST_P(SL_U_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Low,BasisSet::Atom::Type::Slater,verbose);
    Iterate({.NMaxIter = 30, .MinΔρ = Z*1e-4, .MinΔFD = 1e-4, .MinVirial = 5e-1, .MinFD = Z*2e-5, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    // cout << "RelativeHFError = " << RelativeHFError() << std::endl;
    EXPECT_LT(RelativeHFError(),0.01); //1% 
    EXPECT_TRUE(Converged()); 
        
}
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Low,::testing::Values(2,4,10));//));
#else
INSTANTIATE_TEST_SUITE_P(A_HF,SL_U_Low,::testing::Values(2,88));//)); 
#endif

#endif //LOW



