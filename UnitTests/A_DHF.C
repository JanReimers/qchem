// File A_DHF.C  Atom Dirac-Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"

class HF_P : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_P(cluster);
    }
};

class DHF : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_DHF(cluster);
    }
};

//
//  Slater functions
//
class A_SLm_HF_ion : public ::testing::TestWithParam<int>
, public TestAtom, SLm_OBasis, HF_P, TestPolarized
{
public:
    A_SLm_HF_ion() : TestAtom(GetParam(),GetParam()-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLm_HF_ion,Multiple)
{
    int Z=GetParam();
    int N=9;
    // if (Z>12) N=14;
    // if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    Init(N,Z/100.,Z*100.,GetLMax(1));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeError(-0.5*Z*Z),1e-14);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_ion,::testing::Values(1,20,60,86,100)); //37,53


class A_SLmj_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SLmj_OBasis, DHF, TestPolarized
{
public:
    A_SLmj_DHF() : TestAtom(GetParam(),GetParam()-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLmj_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLmj_DHF,Multiple)
{
    int Z=GetParam();
    int N=5;
    // if (Z>12) N=14;
    // if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    Init(N,Z/10.,Z*10.,GetLMax(1));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeError(-0.50000666,true),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLmj_DHF,::testing::Values(1,20,60,86,100)); //37,53

//--------------------------------------------------------------------------------------------
//
// Gaussian functions
//
class A_SG_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SG_RKB_OBasis, DHF, TestPolarized
{
public:
    A_SG_DHF() : TestAtom(GetParam(),GetParam()-1) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_RKB_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SG_DHF,Multiple)
{
    int Z=GetParam();
    int N=15;
    // if (Z>12) N=14;
    // if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    double alpha=0.01024,beta=2.5;
    Init(N,alpha,alpha*pow(beta,N-1),GetLMax(1));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeError(-0.50000666,true),1e-8);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DHF,::testing::Values(1,20,60,86,100)); //37,53
