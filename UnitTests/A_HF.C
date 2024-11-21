// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"

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
, public TestAtom, SG_OBasis, HF_U, TestUnPolarized
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
, public TestAtom, SGm_OBasis, HF_U, TestUnPolarized
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
, public TestAtom, SL_OBasis, HF_U, TestUnPolarized
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
, public TestAtom, SLm_OBasis, HF_U, TestUnPolarized
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
, public TestMolecule, PG_OBasis, HF_U, TestUnPolarized
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

TEST_P(A_SG_HF_U,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>40) N=20;
    if (Z>70) N=25;
    Init(N,0.05,4000*Z,GetLMax(Z));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88)); 

TEST_P(A_SGm_HF_U,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>40) N=20;
    if (Z>70) N=25;
    Init(N,0.05,4000*Z,GetLMax(Z));
    Iterate({40,Z*1e-4,1.0,0.0,true});
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
    Iterate({40,1e-4,1.0,0.0,true});
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
    Iterate({40,1e-4,1.0,0.00,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56));//,70,80,86,88));

TEST_P(A_PG_HF_U,Multiple)
{
    Init();
    Iterate({40,1e-4,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_U,::testing::Values(2,4,10,18,36));



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
, public TestAtom, SG_OBasis, HF_P, TestPolarized
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
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_P,::testing::Values(1,3,7,21,37,53,56,64)); //Can't do boron, need SGm basis for that.
//INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_P,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88)); 

class A_SL_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, HF_P, TestPolarized
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
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_P,::testing::Values(1,3,7,37,53,56,64)); 
//INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_P,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88)); 

class A_SLm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SLm_OBasis, HF_P, TestPolarized
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
    int N=10;
    if (Z>12) N=14;
    if (Z>50) N=16;
    Init(N,0.3,6*Z,GetLMax(Z));
    Iterate({40,Z*1e-4,1.0,0.0,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

//INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Range(57,86)); //37,53
//INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Values(57)); //37,53
INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Values(1,3,5,7,21,37,53,57)); //37,53
//INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Range(2,56)); //37,53
//INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56));//,70,80,86,88)); 

class A_SGm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SGm_OBasis, HF_P, TestPolarized
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
    Iterate({40,Z*1e-4,1.0,0.0,false});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Values(1,3,5,7,21,37,53)); //,53,57,64
//INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56));//,70,80,86,88)); 
//INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Range(2,56)); //,53,57,64


class A_PG_HF_P : public ::testing::TestWithParam<int>
, public TestMolecule, PG_OBasis, HF_P, TestPolarized
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
    Init();
    Iterate({40,1e-3,1.0,0.0,true});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_P,::testing::Values(3,5,21,37)); //7 fails Z=51 is slow
//INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_P,::testing::Range(2,25)); //Z=51 is slow


//
//  Uranium atom test for f-orbitals
//
class A_SL_HF_P_92 : public ::testing::Test
, public TestAtom, SL_OBasis, HF_P, TestPolarized
{
public:
    A_SL_HF_P_92() : TestAtom(92) {};
    void Init(int N, double emin, double emax, int LMax)
    {
        SL_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
};

//TEST_F(A_SL_HF_P_92,Unranium)
//{
//    Init(10,1.0,1.5*92,3);
//    Iterate({40,1e-1,1.0,0.0,false});
//    EXPECT_LT(RelativeHFError(),MaxRelErrE);
//}

//class A_SG_HF_P_92 : public ::testing::Test
//, public TestAtom, SG_OBasis, HF_P, TestPolarized
//{
//public:
//    A_SG_HF_P_92() : TestAtom(92), TestPolarized(6.0) {};
//    void Init(int N, double emin, double emax, int LMax)
//    {
//        SG_OBasis::Init(N,emin,emax,LMax);
//        QchemTester::Init(1e-3);
//    }
//};
//
//TEST_F(A_SG_HF_P_92,Unranium)
//{
//    Init(20,0.08,7000*92,3);
//    Iterate({40,1e-1,1.0,0.0,true});
//    EXPECT_LT(RelativeHFError(),MaxRelErrE);
//}

//class A_PG_HF_P_92 : public ::testing::Test
//, public TestAtom, SL_OBasis, HF_P, TestPolarized
//{
//public:
//    A_PG_HF_P_92() : TestAtom(21), TestPolarized(3.0) {};
//    void Init(int N, double emin, double emax, int LMax)
//    {
//        SL_OBasis::Init(N,emin,emax,LMax);
//        QchemTester::Init(1e-3);
//    }
//};
//
//TEST_F(A_PG_HF_P_92,Unranium)
//{
//    Init(15,0.1,80,2  );
//    Iterate({40,1e-1,1.0,0.0,true});
//    EXPECT_LT(RelativeHFError(),MaxRelErrE);
//}
//
