// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include <MeshParams.H>

class DFT_U : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPT.GetSlaterAlpha(GetZ());
        return new Ham_DFT_U(cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};

//
//  Un-polarized tests.
//
class A_SG_DFT_U : public ::testing::TestWithParam<int>
, public TestAtom, SG_OBasis, DFT_U,TestUnPolarized
{
public:
    A_SG_DFT_U() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_SL_DFT_U : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, DFT_U, TestUnPolarized
{
public:
    A_SL_DFT_U() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_PG_DFT_U : public ::testing::TestWithParam<int>
, public TestMolecule, PG_OBasis, DFT_U, TestUnPolarized
{
public:
    A_PG_DFT_U() : TestMolecule() {};
    void Init()
    { 
        Molecule* m=new Molecule;
        m->Insert(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0)));
        TestMolecule::Init(m);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SG_DFT_U,Multiple)
{
    int Z=GetParam();
    Init(20,0.05,10000*Z,GetLMax(Z));
    Iterate({40,Z*1e-4,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DFT_U,::testing::Values(2,4,10,18,36,54)); 

TEST_P(A_SL_DFT_U,Multiple)
{
    int Z=GetParam();
    Init(9, 0.7,5.5*Z,GetLMax(Z));
    Iterate({40,Z*1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SL_DFT_U,::testing::Values(2,4,10,18,36,54));

TEST_P(A_PG_DFT_U,Multiple)
{
    Init();
    Iterate({40,1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_DFT_U,::testing::Values(2,4,10,18,36));

class DFT_P : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPT.GetSlaterAlpha(GetZ());
        return new Ham_DFT_P(cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};
//
//  Polarized tests.
//
class A_SG_DFT_P : public ::testing::TestWithParam<int>
, public TestAtom, SG_OBasis, DFT_P,TestPolarized
{
public:
    A_SG_DFT_P() : TestAtom(GetParam()),TestPolarized(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SG_DFT_P,Multiple)
{
    int Z=GetParam();
    Init(20,0.01,4000*Z,GetLMax(Z));
    Iterate({40,Z*1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DFT_P,::testing::Values(1,3,5,7,37,53)); //,3,5,7,37,53

class A_SL_DFT_P : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, DFT_P, TestPolarized
{
public:
    A_SL_DFT_P() : TestAtom(GetParam()), TestPolarized(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SL_DFT_P,Multiple)
{
    int Z=GetParam();
    Init(9,0.7,5.5*Z,GetLMax(Z));
    Iterate({40,Z*1e-2,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_DFT_P,::testing::Values(1,3,5,7,37,53)); 

class A_PG_DFT_P : public ::testing::TestWithParam<int>
, public TestMolecule, PG_OBasis, DFT_P, TestPolarized
{
public:
    A_PG_DFT_P() : TestMolecule(), TestPolarized(GetParam())  {};
    void Init()
    { 
        Molecule* m=new Molecule;
        m->Insert(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0)));
        TestMolecule::Init(m);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_PG_DFT_P,Multiple)
{
    Init();
    Iterate({40,1e-3,1.0,0.0,false});
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_DFT_P,::testing::Values(3,5,7,37)); //Z=51 is slow.


