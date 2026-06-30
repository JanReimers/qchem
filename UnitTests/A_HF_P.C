// File A_HF_P.C  Atom Hartree-Fock tests for Polarized (open shell) atoms.
#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
import qchem.Unittests.QchemTester;

const bool verbose=true;

using std::cout;
using std::endl;
using enum BasisSet::Atom::BasisSetAccuracy;
using namespace qchem::Hamiltonian;
class A_HF_P : public ::testing::TestWithParam<size_t>, public TestAtom
{
public:
    A_HF_P() : TestAtom(GetParam()) {};
    virtual Hamiltonian* GetHamiltonian(st_t& structure) const
    {
        return Factory(Model::HF,Pol::Polarized,structure);
    }
};

#ifdef DEBUG
#define LOW
#else
// #define HIGH
#define MEDIUM
#define LOW
#endif

#ifdef HIGH
class A_SG_HF_P_High : public A_HF_P {};
TEST_P(A_SG_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Gaussian,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-5, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),2.5e-6); 
    EXPECT_TRUE(Converged()); 
        
}
// INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_High,::testing::Values(1,3,5,6,7,8,9,13,14,16,17,21,22,23,25,26,27,28,37,39,40,41,44,45,47,57,58,59,60,61,62,63,65,66,67,68,69,73,91,92)); 
INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_High,::testing::Values(5,21,92)); 

class A_SL_HF_P_High : public A_HF_P {};
TEST_P(A_SL_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::Slater,verbose);
    Iterate({.NMaxIter = 32, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 1e-6, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),2.5e-6); 
    EXPECT_TRUE(Converged()); 
        
}
// INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_High,::testing::Values(1,3,5,6,7,8,9,13,14,16,17,21,22,23,25,26,27,28,37,39,40,41,44,45,47,57,58,59,60,61,62,63,65,66,67,68,69,73,91,92)); 
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_High,::testing::Values(5,21,92)); 


class A_BS_HF_P_High : public A_HF_P {};
TEST_P(A_BS_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::BSpline6,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-12, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),2.5e-6); 
    EXPECT_TRUE(Converged()); 
        
}
// INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_High,::testing::Values(1,3,5,6,7,8,9,13,14,16,17,21,22,23,25,26,27,28,37,39,40,41,44,45,47,57,58,59,60,61,62,63,65,66,67,68,69,73,91,92)); 
INSTANTIATE_TEST_SUITE_P(A_HF,A_BS_HF_P_High,::testing::Values(5,21)); 

class A_BSr_HF_P_High : public A_HF_P {};
TEST_P(A_BSr_HF_P_High,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(High,BasisSet::Atom::Type::BSpliner6,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-7, .MinΔFD = 1e-6, .MinVirial = 5e-12, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),2.5e-6); 
    EXPECT_TRUE(Converged()); 
        
}
// INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_High,::testing::Values(1,3,5,6,7,8,9,13,14,16,17,21,22,23,25,26,27,28,37,39,40,41,44,45,47,57,58,59,60,61,62,63,65,66,67,68,69,73,91,92)); 
INSTANTIATE_TEST_SUITE_P(A_HF,A_BSr_HF_P_High,::testing::Values(5,21)); 
#endif //HIGH

#ifdef MEDIUM
class A_SG_HF_P_Medium : public A_HF_P {};
TEST_P(A_SG_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Gaussian,verbose);
    Iterate({.NMaxIter = 50, .MinΔρ = Z*1e-5, .MinΔFD = 1e-7, .MinVirial = 5e-2, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),60e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_SG_HF_P_Medium,::testing::Values(3,5)); 

class A_SL_HF_P_Medium : public A_HF_P {};
TEST_P(A_SL_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Slater,verbose);
    Iterate({.NMaxIter = 22, .MinΔρ = Z*1e-4, .MinΔFD = 1e-5, .MinVirial = 5e-4, .MinFD = Z*1e-6, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),20e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_Medium,::testing::Values(3,5,21)); 



class A_BS_HF_P_Medium : public A_HF_P {};
TEST_P(A_BS_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::BSpline6,verbose);
    Iterate({.NMaxIter = 40, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),20e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_BS_HF_P_Medium,::testing::Values(3,5)); 
class A_BSr_HF_P_Medium : public A_HF_P {};
TEST_P(A_BSr_HF_P_Medium,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::BSpliner6,verbose);
    Iterate({.NMaxIter = 40, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 2.5e-7, .MinFD = Z*1e-7, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),20e-6); 
    EXPECT_TRUE(Converged()); 
        
}
INSTANTIATE_TEST_SUITE_P(A_HF,A_BSr_HF_P_Medium,::testing::Values(3,5)); 

#endif //MEDIUM

#ifdef LOW
class A_SL_HF_P_Low : public A_HF_P {};
TEST_P(A_SL_HF_P_Low,A)
{
    size_t Z=GetParam();
    cout << "---------------- Z=" << Z << " ---------------"<< endl;
        
    QchemTester::Init(Medium,BasisSet::Atom::Type::Slater,verbose);
    Iterate({.NMaxIter = 14, .MinΔρ = Z*1e-4, .MinΔFD = 1e-4, .MinVirial = 1.1e-4, .MinFD = Z*1e-5, .StartingRelaxRo = Z<40 ? 0.5 : 0.3, .MergeTol = 1e-7, .Verbose = true});
    EXPECT_LT(RelativeHFError(),12e-6); 
    EXPECT_TRUE(Converged()); 
        
}
#ifdef DEBUG
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_Low,::testing::Values(3,5,21)); 
#else
INSTANTIATE_TEST_SUITE_P(A_HF,A_SL_HF_P_Low,::testing::Values(3,5,21,92)); 
#endif

#endif


// inline SCFParams saito_params_BS(int Z) 
// {
//     // Sc (Z=21) needs a lot of iterations to the virial converged.
// #ifdef DEBUG
//     return {.NMaxIter = 55, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 1e-12, .MinFD = Z*5e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true};
// #else
//     return {.NMaxIter = 55, .MinΔρ = Z*1e-7, .MinΔFD = 1e-7, .MinVirial = 1e-12, .MinFD = Z*5e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-7, .Verbose = true};
// #endif
// }

// class A_BS_saito_HF_P : public A_HF_P
// {

// };

// TEST_P(A_BS_saito_HF_P,Saito)
// {
//     int Z=GetParam();
// #ifdef DEBUG
//     nlohmann::json js = {
//         {"type",BasisSet::Atom::Type::BSpline6},
//         {"N", 20}, {"rmin", 0.2}, {"rmax", 20},
//     };
// #else
//    nlohmann::json js = {
//         {"type",BasisSet::Atom::Type::BSpliner6},
//         {"N", 50}, {"rmin", 0.003}, {"rmax", 40},
//     };
// #endif
//     QchemTester::Init(js);
//     Iterate(saito_params_BS(Z));
//     double Eerr=RelativeHFError();
//     EXPECT_LT(Eerr,1e-9);
//     EXPECT_GT(Eerr,-1e-4);
//     EXPECT_TRUE(Converged());
// }

// #ifdef NDEBUG
// INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(1,2,5,21,57,92)); 

// // // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Range(1,93)); 
// // // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(62,63,65,67,68,71,72)); //Known convergence problems. 
// // // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(1,92)); 
// #else
// // // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Range(2,3)); 
// INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(1,2,5,21,57,92)); 
// #endif
