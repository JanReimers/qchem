// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"

//
//  Un-polarized tests.
//
class A_SG_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, public SG_OBasis, HFHamiltonian, TestUnPolarized
{
public:
    void Init(int Z,int N, double emin, double emax, int LMax)
    {
        TestAtom::Init(Z);
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init();
    }
};

class A_SL_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, public SL_OBasis, HFHamiltonian, TestUnPolarized
{
public:
    void Init(int Z,int N, double emin, double emax, int LMax)
    {
        TestAtom::Init(Z);
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init();
    }
};

//
//  High precision HE using optimized exponent ranges.
//
TEST_F(A_SG_HF_U,He)
{
    Init(2,20,.1,10000,0);
    Iterate({40,1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

TEST_F(A_SL_HF_U,He)
{
    Init(2,8,.31,10.9,0);
    Iterate({40,1e-4,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

TEST_P(A_SG_HF_U,Multiple)
{
    int Z=GetParam();
    Init(Z,20,0.05,4000*Z,GetLMax(Z));
    Iterate({40,Z*1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_U,::testing::Values(2,4,10,18,36,54)); 

TEST_P(A_SL_HF_U,Multiple)
{
    int Z=GetParam();
    Init(Z,10,1.0,1.5*Z,GetLMax(Z));
    Iterate({40,1e-1,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_U,::testing::Values(2,4,10,18,36,54));


//
//  Polarized tests.
//
class A_SG_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, public SG_OBasis, PolHFHamiltonian, TestPolarized
{
public:
    void Init(int Z,int N, double emin, double emax, int LMax)
    {
        TestAtom::Init(Z);
        SG_OBasis::Init(N,emin,emax,LMax);
        TestPolarized::Init(GetNumUnpairedElectrons(Z));
        QchemTester::Init();
    }
};

TEST_P(A_SG_HF_P,Multiple)
{
    int Z=GetParam();
    Init(Z,20,0.05,4000*Z,GetLMax(Z));
    Iterate({40,Z*1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_P,::testing::Values(1,3,5,7,37,53)); 

class A_SL_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, public SL_OBasis, PolHFHamiltonian, TestPolarized
{
public:
    void Init(int Z,int N, double emin, double emax, int LMax)
    {
        TestAtom::Init(Z);
        SL_OBasis::Init(N,emin,emax,LMax);
        TestPolarized::Init(GetNumUnpairedElectrons(Z));
        QchemTester::Init();
    }
};

TEST_P(A_SL_HF_P,Multiple)
{
    int Z=GetParam();
    Init(Z,10,0.7,1.5*Z,GetLMax(Z));
    Iterate({40,Z*1e-2,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_P,::testing::Values(1,3,5,7,37,53)); 
