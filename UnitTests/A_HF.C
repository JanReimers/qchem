// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include "Cluster/Atom.H"
#include "Cluster/Molecule.H"
#include <Factory.H> //Just to get the types.

const bool verbose=true;
inline SCFParams scf_params(int Z) 
{
//           NMaxIter MinDeltaRo MinDelE MinError StartingRelaxRo MergeTol verbose
    return {   80     ,Z*1e-5    ,1e-10   ,Z*1e-6        ,0.5     ,1e-7  ,verbose};
}
inline SCFParams scf_params_BS(int Z) 
{
//           NMaxIter MinDeltaRo MinDelE MinError StartingRelaxRo MergeTol verbose
    return {   30     ,Z*1e-5    ,1e-10   ,Z*1e-11        ,0.5     ,1e-7  ,verbose};
}

//---------------------------------------------------------------------------------------------------------------
//
//  Un-polarized tests.
//
class HF_U : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_U(cluster);
    }
};

class A_SG_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, SG_OBasis, HF_U
{
public:
    A_SG_HF_U() : TestAtom(GetParam()) {};
};

class A_SL_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, HF_U
{
public:
    A_SL_HF_U() : TestAtom(GetParam()) {};
};

class A_BS_HF_U : public ::testing::TestWithParam<int>
, public TestAtom, BS_OBasis, HF_U
{
public:
    A_BS_HF_U() : TestAtom(GetParam()) {};
 
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
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-3,js);
    }
};



static std::map<int,size_t> expected_itartion_counts={{2,10},{4,11},{10,12},{12,14},{18,13},{20,16},{30,17},{36,12},{38,15},{46,15},{48,16},{54,13},{56,15},{70,22},{80,16},{86,15},{88,17}};
static std::map<int,size_t> NBasis={{2,20},{4,20},{10,22},{12,22},{18,23},{20,25},{30,25},{36,25},{38,27},{46,27},{48,27},{54,30},{56,30},{70,30},{80,30},{86,30},{88,30}};
TEST_P(A_SG_HF_U,Multiple)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Gaussian},
        {"N", NBasis[Z]}, {"emin", 0.01}, {"emax", 10000*Z*sqrt(Z)},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),1e-6);
    assert(expected_itartion_counts.find(Z)!=expected_itartion_counts.end());
    size_t ic_expected=expected_itartion_counts[Z];
    EXPECT_LE(GetIterationCount(),ic_expected);
    EXPECT_EQ(GetIterationCount(),ic_expected);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));//)); 


TEST_P(A_SL_HF_U,Multiple)
{
    int Z=GetParam();
    int N=10;
    if (Z>15) N=14;
    if (Z>50) N=18;
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.3}, {"emax", 5*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_U,::testing::Values(2,4,10,12,18,20,30,36,38,46,48,54,56,70,80,86,88));

TEST_P(A_BS_HF_U,Multiple)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::BSpline},
        {"N", 40}, {"rmin", 0.01}, {"rmax", 50},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params_BS(Z));
    EXPECT_LT(RelativeHFError(),1e-9);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_BS_HF_U,::testing::Values(2,4,10,12,18)); 

TEST_P(A_PG_HF_U,Multiple)
{
    int Z=GetParam();
    Init();
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_HF_U,::testing::Values(2,4,10,18,36));





















//---------------------------------------------------------------------------------------------------------------
//
//  Polarized tests.
//
class HF_P : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_P(cluster);
    }
};

class A_SG_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SG_OBasis, HF_P
{
public:
    A_SG_HF_P() : TestAtom(GetParam()) {};
   
};

class A_SL_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SL_OBasis, HF_P
{
public:
    A_SL_HF_P() : TestAtom(GetParam()) {};
};

class A_BS_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, BS_OBasis, HF_P
{
public:
    A_BS_HF_P() : TestAtom(GetParam()) {};
};

class A_SLm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SLm_OBasis, HF_P
{
public:
    A_SLm_HF_P() : TestAtom(GetParam()) {};
   
};

class A_SGm_HF_P : public ::testing::TestWithParam<int>
, public TestAtom, SGm_OBasis, HF_P
{
public:
    A_SGm_HF_P() : TestAtom(GetParam()) {};
};

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

TEST_P(A_SG_HF_P,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>20) N=25;
    if (Z>40) N=27;    
    if (Z>70) N=30;
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Gaussian},
        {"N", N}, {"emin", 0.01}, {"emax", 10000*Z*sqrt(Z)},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),1e-6);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_HF_P,::testing::Values(1,3,7,25,37,47,63,75)); //Can't do boron or Sc, need SGm basis for that.


TEST_P(A_SL_HF_P,Multiple)
{
    int Z=GetParam();
    int N=10;
    if (Z>15) N=14;
    if (Z>50) N=18;
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.3}, {"emax", 5*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_HF_P,::testing::Values(1,3,7,37,53,56,64)); 


TEST_P(A_BS_HF_P,Multiple)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::BSpline},
        {"N", 20}, {"rmin", 1.0/Z}, {"rmax", 30},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

// INSTANTIATE_TEST_CASE_P(Multiple,A_BS_HF_P,::testing::Values(1,2,3,4,)); 
INSTANTIATE_TEST_CASE_P(Multiple,A_BS_HF_P,::testing::Values(1,3,7,37,53,56,64)); 


TEST_P(A_SLm_HF_P,Multiple)
{
    int Z=GetParam();
    int N=11;
    if (Z>12) N=16;
    if (Z>50) N=20;
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.125}, {"emax", 8*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Values(1,3,5,6,7,8,9,11,13,14,15,16,17,21,37,53,57,64)); 
// INSTANTIATE_TEST_CASE_P(Multiple,A_SLm_HF_P,::testing::Range(1,93)); 

TEST_P(A_SGm_HF_P,Multiple)
{
    int Z=GetParam();
    int N=20;
    if (Z>20) N=25;
    if (Z>70) N=25;
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Gaussian},
        {"N", N}, {"emin", 0.01}, {"emax", 10000*Z*sqrt(Z)},
    };
    QchemTester::Init(1e-3,js);
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Values(1,3,5,7,21,37,53,58,62)); 
// INSTANTIATE_TEST_CASE_P(Multiple,A_SGm_HF_P,::testing::Range(1,93)); //Only goes to 92?!?


TEST_P(A_BSm_HF_P,Multiple)
{
    int Z=GetParam();
     nlohmann::json js = {
        {"type",BasisSetAtom::Type::BSpline},
        {"N", 30}, {"rmin", 1.0/Z}, {"rmax", 30},
    };
    QchemTester::Init(1e-3,js);
   
    Iterate(scf_params(Z));
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_BSm_HF_P,::testing::Values(5,7,21,37,53,57,64)); 


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
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-3,js);
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

