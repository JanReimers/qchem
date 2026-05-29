// File A_HF_P.C  Atom Hartree-Fock tests for Polarized (open shell) atoms.
#include "gtest/gtest.h"
#include "nlohmann/json.hpp"
import qchem.Unittests.QchemTester;

const bool verbose=true;

using std::cout;
using std::endl;
using enum BasisSetAccuracy;
using namespace qchem::Hamiltonian;
class A_HF_P : public virtual QchemTester, public ::testing::TestWithParam<size_t>, TestAtom
{
public:
    A_HF_P() : TestAtom(GetParam()) {};
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::HF,Pol::Polarized,cluster);
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-5    ,1e-7 , 1e-5      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   32     ,Z*1e-5    ,1e-7 , 1e-6      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-7    ,1e-7 , 2.5e-12      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-7    ,1e-6 , 5e-12      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   50     ,Z*1e-5    ,1e-7 , 5e-2      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   22     ,Z*1e-4    ,1e-5 , 5e-4      ,Z*1e-6 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   40     ,Z*1e-7    ,1e-7 , 2.5e-7      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   40     ,Z*1e-7    ,1e-7 , 2.5e-7      ,Z*1e-7 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
    //       NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo    MergeTol verbose
    Iterate({   14     ,Z*1e-4    ,1e-4 , 1.1e-4      ,Z*1e-5 ,Z<40 ? 0.5 : 0.3   ,1e-7  ,true});
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
// //           NMaxIter MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo MergeTol verbose
// #ifdef DEBUG
//     return {   50     ,Z*1e-7    ,1e-9 , 1e-12   ,Z*5e-7        ,0.5     ,1e-7  ,true};
// #else
//     return {   50     ,Z*1e-7    ,1e-9 , 1e-12   ,Z*5e-7        ,0.5     ,1e-7  ,true};
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

// // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Range(1,93)); 
// // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(62,63,65,67,68,71,72)); //Known convergence problems. 
// // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(1,92)); 
// #else
// // INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Range(2,3)); 
// INSTANTIATE_TEST_SUITE_P(Saito,A_BS_saito_HF_P,::testing::Values(2)); 
// #endif
