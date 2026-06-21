// File A_HF.C  Atom Hartree-Fock tests.
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <filesystem>

import qchem.Unittests.QchemTester;

import qchem.Hamiltonian.Factory;
import qchem.Factory;
import qchem.Cluster;


inline SCFParams dft_scf_params(int Z) 
{
//           NMaxIter MinΔρ MinΔFD MinVirial MinFD StartingRelaxRo MergeTol verbose
    return {   20     ,Z*1e-3    ,1e-10  ,1e-13 ,Z*1e-4        ,0.1      ,1e-8  ,true};
}

using namespace qchem::Hamiltonian;
class A_DFT_U : public ::testing::TestWithParam<size_t>, public TestAtom
{
public:
    A_DFT_U() : TestAtom(GetParam()) {};
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPTold.GetSlaterAlpha(GetZ());
        return Factory(Pol::UnPolarized,cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};
class M_DFT_U : public ::testing::TestWithParam<size_t>, public TestMolecule
{
public:
    M_DFT_U() : TestMolecule(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0))) {};
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPTold.GetSlaterAlpha(GetZ());
        return Factory(Pol::UnPolarized,cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};

//---------------------------------------------------------------------------------------------------------------
//
//  Un-polarized tests.
//
class A_SG_DFT_U : public A_DFT_U {};

class A_SL_DFT_U : public A_DFT_U {};

class A_PG_DFT_U : public M_DFT_U
{
public:
    void Init()
    { 
        nlohmann::json js = { {"basis", "dzvp"} };
        QchemTester::Init(js);
    }
};

TEST_P(A_SG_DFT_U,A)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",abs_t::Gaussian},
        {"N", 20}, {"emin", 0.05}, {"emax", 10000*Z},
    };
    QchemTester::Init(js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),2e-3);
}
INSTANTIATE_TEST_SUITE_P(A,A_SG_DFT_U,::testing::Values(2,4,10,18,36,54)); 

TEST_P(A_SL_DFT_U,Multiple)
{
    int Z=GetParam();
    int N=8;
    if (Z>20) N=10;
    if (Z>50) N=11;
    nlohmann::json js = {
        {"type",abs_t::Slater},
        {"N", N}, {"emin", 0.31}, {"emax", 3*Z},
    };
    QchemTester::Init(js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),2e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_SL_DFT_U,::testing::Values(2,4,10,18,36,54));

TEST_P(A_PG_DFT_U,Multiple)
{
    Init();
    Iterate(dft_scf_params(GetParam()));
    EXPECT_LT(RelativeDFTError(),2.2e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_PG_DFT_U,::testing::Values(2,4,10,18,36));










//---------------------------------------------------------------------------------------------------------------
//
//  Polarized tests.
//
class A_DFT_P : public ::testing::TestWithParam<size_t>, public TestAtom
{
public:
    A_DFT_P() : TestAtom(GetParam()) {};
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPTold.GetSlaterAlpha(GetZ());
        return Factory(Pol::Polarized,cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};
class M_DFT_P : public ::testing::TestWithParam<size_t>, public TestMolecule
{
public:
    M_DFT_P() : TestMolecule(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0))) 
    {
        nlohmann::json js = { {"basis", "dzvp"} };
        QchemTester::Init(js);
    };
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPTold.GetSlaterAlpha(GetZ());
        return Factory(Pol::Polarized,cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};
class A_SG_DFT_P : public A_DFT_P {};
class A_SL_DFT_P : public A_DFT_P {};
class A_PG_DFT_P : public M_DFT_P {};

TEST_P(A_SG_DFT_P,Multiple)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",abs_t::Gaussian},
        {"N", 20}, {"emin", 0.01}, {"emax", 4000*Z},
    };
    QchemTester::Init(js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),1e-3);
}

INSTANTIATE_TEST_SUITE_P(Multiple,A_SG_DFT_P,::testing::Values(1,3,7,37,53)); //,3,5,7,37,53

TEST_P(A_SL_DFT_P,Multiple)
{
    int Z=GetParam();
    int N=8;
    if (Z>20) N=10;
    if (Z>50) N=11;
     nlohmann::json js = {
        {"type",abs_t::Slater},
        {"N", N}, {"emin", 0.31}, {"emax", 3*Z},
    };
    QchemTester::Init(js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),2e-3);
}

INSTANTIATE_TEST_SUITE_P(Multiple,A_SL_DFT_P,::testing::Values(1,3,7,37,53)); 



TEST_P(A_PG_DFT_P,Multiple)
{
    Iterate(dft_scf_params(GetParam()));
    EXPECT_LT(RelativeDFTError(),5.1e-3);
}
INSTANTIATE_TEST_SUITE_P(Multiple,A_PG_DFT_P,::testing::Values(3,5,11,37)); //Z=51 is slow.


