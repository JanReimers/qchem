// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"

//
//  Un-polarized tests.
//
class A_SG_SHF_U : public ::testing::TestWithParam<int>
, public TestAtom, public SG_OBasis, SHFHamiltonian, TestUnPolarized
{
public:
    A_SG_SHF_U() : TestAtom(GetParam()), SHFHamiltonian(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_SL_SHF_U : public ::testing::TestWithParam<int>
, public TestAtom, public SL_OBasis, SHFHamiltonian, TestUnPolarized
{
public:
    A_SL_SHF_U() : TestAtom(GetParam()), SHFHamiltonian(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SG_SHF_U,Multiple)
{
    int Z=GetParam();
    Init(20,0.05,10000*Z,GetLMax(Z));
    Iterate({40,Z*1e-4,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_SHF_U,::testing::Values(2,4,10,18,36,54)); 

TEST_P(A_SL_SHF_U,Multiple)
{
    int Z=GetParam();
    Init(10, 0.7,1.5*Z,GetLMax(Z));
    Iterate({40,Z*1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SL_SHF_U,::testing::Values(2,4,10,18,36,54));


//
//  Polarized tests.
//
class A_SG_SHF_P : public ::testing::TestWithParam<int>
, public TestAtom, public SG_OBasis, PolSHFHamiltonian, TestPolarized
{
public:
    A_SG_SHF_P() : TestAtom(GetParam()), PolSHFHamiltonian(GetParam()),TestPolarized(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SG_SHF_P,Multiple)
{
    int Z=GetParam();
    Init(20,0.01,4000*Z,GetLMax(Z));
    Iterate({40,Z*1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_SHF_P,::testing::Values(1,3,5,7,37,53)); //,3,5,7,37,53

class A_SL_SHF_P : public ::testing::TestWithParam<int>
, public TestAtom, public SL_OBasis, PolSHFHamiltonian, TestPolarized
{
public:
    A_SL_SHF_P() : TestAtom(GetParam()), PolSHFHamiltonian(GetParam()),TestPolarized(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SL_SHF_P,Multiple)
{
    int Z=GetParam();
    Init(10,0.7,1.5*Z,GetLMax(Z));
    Iterate({40,Z*1e-2,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_SHF_P,::testing::Values(1,3,5,7,37,53)); 

