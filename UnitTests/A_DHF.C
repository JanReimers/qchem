// File A_DHF.C  Atom Dirac-Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"


class DHF : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_DHF(cluster);
    }
};

//
//  Polarized tests.
//
class A_SLmj_DHF : public ::testing::TestWithParam<int>
, public TestAtom, SLmj_OBasis, DHF, TestPolarized
{
public:
    A_SLmj_DHF() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLmj_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLmj_DHF,Multiple)
{
    int Z=GetParam();
    int N=3;
    if (Z>12) N=14;
    if (Z>50) N=16;
    //Init(N,1.0,1.0,GetLMax(Z));
    Init(N,0.1,10*Z,GetLMax(Z));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLmj_DHF,::testing::Values(1)); //37,53
