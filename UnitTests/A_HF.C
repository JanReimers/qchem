// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"

SCFIterationParams scf_params(int Z) 
{
//           NMaxIter MinDeltaRo MinDelE MinError StartingRelaxRo verbose
    return {   80     ,Z*1e-5    ,1e-10   ,Z*1e-6        ,0.5       ,true};
}

class HF_U : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_U(cluster);
    }
};

//
//  Un-polarized tests.
//
class A_SG_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, SG_OBasis, HF_U
{
public:
    A_SG_HF_U() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_SGm_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, SGm_OBasis, HF_U
{
public:
    A_SGm_HF_U() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SGm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_SL_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, HF_U
{
public:
    A_SL_HF_U() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_SLm_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, SLm_OBasis, HF_U
{
public:
    A_SLm_HF_U() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

class A_PG_HF_U : public ::testing::TestWithParam<int>
, public TestMolecule, PG_OBasis, HF_U
{
public:
    A_PG_HF_U() : TestMolecule() {};
    void Init()
    { 
        Molecule* m=new Molecule;
        m->Insert(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0)));
        TestMolecule::Init(m);
        QchemTester::Init(1e-3);
    }
};

static std::map<int,size_t> expected_itartion_counts={{2,14},{4,12},{10,12},{12,13},{18,10},{20,11},{30,20},{36,18},{38,18},{46,30},{48,29},{54,20},{56,15},{70,28},{80,30},{86,27},{88,20}};
TEST_P(A_SG_HF_U,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>40) N=20;
    if (Z>70) N=25;
    Init(N,0.05,4000*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
    assert(expected_itartion_counts.find(Z)!=expected_itartion_counts.end());
    size_t ic_expected=expected_itartion_counts[Z];
    EXPECT_LE(GetIterationCount(),ic_expected);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88)); 

TEST_P(A_SGm_HF_U,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>40) N=20;
    if (Z>70) N=25;
    Init(N,0.05,4000*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56));//,70,80,86,88)); 

TEST_P(A_SL_HF_U,Multiple)
{
    int Z=GetParam();
    int N=10;
    if (Z>15) N=14;
    if (Z>50) N=18;
    Init(N,0.3,5*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));

TEST_P(A_SLm_HF_U,Multiple)
{
    int Z=GetParam();
    int N=10;
    if (Z>15) N=15;
    if (Z>50) N=18;
    Init(N,0.3,5*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56));//,70,80,86,88));

TEST_P(A_PG_HF_U,Multiple)
{
    int Z=GetParam();
    Init();
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_U,::testing::Values(2,4,10,18,36));

class A_BS_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, BS_OBasis, HF_U
{
public:
    A_BS_HF_U() : TestAtom(GetParam()) {};
    void Init(int N, double rmin, double rmax, int LMax)
    {
        BS_OBasis::Init(N,rmin,rmax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_BS_HF_U,Multiple)
{
    int Z=GetParam();
    int N=50;
    Init(N,0.1,40,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_BS_HF_U,::testing::Values(2,4)); 


class HF_P : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_P(cluster);
    }
};

//
//  Polarized tests.
//
class A_SG_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SG_OBasis, HF_P
{
public:
    A_SG_HF_P() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SG_HF_P,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>40) N=20;
    if (Z>70) N=25;    
    Init(N,0.05,4000*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_P,::testing::Values(1,3,7,21,37,53,56,64)); //Can't do boron, need SGm basis for that.
//INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_P,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88)); 

class A_SL_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, HF_P
{
public:
    A_SL_HF_P() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SL_HF_P,Multiple)
{
    int Z=GetParam();
    int N=10;
    if (Z>15) N=14;
    if (Z>50) N=18;
    Init(N,0.3,5*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_P,::testing::Values(1,3,7,37,53,56,64)); 
//INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_P,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88)); 

class A_BS_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, BS_OBasis, HF_P
{
public:
    A_BS_HF_P() : TestAtom(GetParam()) {};
    void Init(int N, double rmin, double rmax, int LMax)
    {
        BS_OBasis::Init(N,rmin,rmax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_BS_HF_P,Multiple)
{
    int Z=GetParam();
    int N=20;
    Init(N,1.0/Z,30,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

// INSTANTIATE_TEST_CASE_P(Multiple,A_BS_HF_P,::testing::Values(1,2,3,4,)); 
INSTANTIATE_TEST_CASE_P(Multiple,A_BS_HF_P,::testing::Values(1,3,7,37,53,56,64)); 

class A_SLm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SLm_OBasis, HF_P
{
public:
    A_SLm_HF_P() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SLm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SLm_HF_P,Multiple)
{
    int Z=GetParam();
    int N=11;
    if (Z>12) N=16;
    if (Z>50) N=20;
    Init(N,0.125,8*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Values(1,3,5,6,7,8,9,11,13,14,15,16,17,21,37,53,57,64)); 
// INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Range(1,93)); 

class A_SGm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SGm_OBasis, HF_P
{
public:
    A_SGm_HF_P() : TestAtom(GetParam()) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SGm_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_SGm_HF_P,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>40) N=20;
    if (Z>70) N=25;
    Init(N,0.05,6000*Z,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Values(1,3,5,7,21,37,53,57,64)); 
// INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Range(1,93)); //Only goes to 92?!?

class A_BSm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, BSm_OBasis, HF_P
{
public:
    A_BSm_HF_P() : TestAtom(GetParam()) {};
    void Init(int N, double rmin, double rmax, int LMax)
    {
        BSm_OBasis::Init(N,rmin,rmax,LMax);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_BSm_HF_P,Multiple)
{
    int Z=GetParam();
    int N=30;
    Init(N,1.0/Z,30.,GetLMax(Z));
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_BSm_HF_P,::testing::Values(5,7,21,37,53,57,64)); 
// INSTANTIATE_TEST_CASE_P(Multiple,A_BSm_HF_P,::testing::Range(1,93));


class A_PG_HF_P : public ::testing::TestWithParam<int>
, public TestMolecule, PG_OBasis, HF_P
{
public:
    A_PG_HF_P() : TestMolecule() {};
    void Init()
    { 
        Molecule* m=new Molecule;
        m->Insert(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0)));
        TestMolecule::Init(m);
        QchemTester::Init(1e-3);
    }
};

TEST_P(A_PG_HF_P,Multiple)
{
    int Z=GetParam();
    Init();
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_P,::testing::Values(3,5,21,37)); //7 fails Z=51 is slow
//INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_P,::testing::Range(2,25)); //Z=51 is slow


//
//  Uranium atom test for f-orbitals
//
class A_SL_HF_P_92 : public ::testing::Test
, public TestAtom, SL_OBasis, HF_P
{
public:
    A_SL_HF_P_92() : TestAtom(92) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

