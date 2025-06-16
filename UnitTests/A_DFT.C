// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"
#include "Cluster/Atom.H"
#include "Cluster/Molecule.H"
#include <Mesh/MeshParams.H>
#include <BasisSet/Factory.H> //Just to get the types.
#include <Hamiltonian/Factory.H>

inline SCFParams dft_scf_params(int Z) 
{
//           NMaxIter MinDeltaRo MinDelE MinError StartingRelaxRo MergeTol verbose
    return {   20     ,Z*1e-3    ,1e-10   ,Z*1e-4        ,0.1      ,1e-8  ,true};
}

using namespace HamiltonianF;
class DFT_U : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPT.GetSlaterAlpha(GetZ());
        return Factory(Pol::UnPolarized,cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};

//---------------------------------------------------------------------------------------------------------------
//
//  Un-polarized tests.
//
class A_SG_DFT_U : public ::testing::TestWithParam<int>
, public TestAtom, DFT_U
{
public:
    A_SG_DFT_U() : TestAtom(GetParam()) {};
};

class A_SL_DFT_U : public ::testing::TestWithParam<int>
, public TestAtom, DFT_U
{
public:
    A_SL_DFT_U() : TestAtom(GetParam()) {};
};

class A_PG_DFT_U : public ::testing::TestWithParam<int>
, public TestMolecule,  DFT_U
{
public:
    A_PG_DFT_U() : TestMolecule() {};
    void Init()
    { 
        Molecule* m=new Molecule;
        m->Insert(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0)));
        TestMolecule::Init(m);
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-3,js);
    }
};

TEST_P(A_SG_DFT_U,Multiple)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Gaussian},
        {"N", 20}, {"emin", 0.05}, {"emax", 10000*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DFT_U,::testing::Values(2,4,10,18,36,54)); 

TEST_P(A_SL_DFT_U,Multiple)
{
    int Z=GetParam();
    int N=8;
    if (Z>20) N=10;
    if (Z>50) N=11;
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.31}, {"emax", 3*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_SL_DFT_U,::testing::Values(2,4,10,18,36,54));

TEST_P(A_PG_DFT_U,Multiple)
{
    Init();
    Iterate(dft_scf_params(GetParam()));
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_DFT_U,::testing::Values(2,4,10,18,36));










//---------------------------------------------------------------------------------------------------------------
//
//  Polarized tests.
//
class DFT_P : public virtual QchemTester
{
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        double alpha_ex=QchemTester::itsPT.GetSlaterAlpha(GetZ());
        return Factory(Pol::Polarized,cluster,alpha_ex,GetMeshParams(),itsBasisSet);
    }
};

class A_SG_DFT_P : public ::testing::TestWithParam<int>
, public TestAtom, DFT_P
{
public:
    A_SG_DFT_P() : TestAtom(GetParam()) {};
};
class A_SL_DFT_P : public ::testing::TestWithParam<int>
, public TestAtom,  DFT_P
{
public:
    A_SL_DFT_P() : TestAtom(GetParam()) {};
};

class A_PG_DFT_P : public ::testing::TestWithParam<int>
, public TestMolecule,  DFT_P
{
public:
    A_PG_DFT_P() : TestMolecule() {};
    void Init()
    { 
        Molecule* m=new Molecule;
        m->Insert(new Atom(GetParam(),0.0,Vector3D<double>(0,0,0)));
        TestMolecule::Init(m);
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-3,js);
    }
};
TEST_P(A_SG_DFT_P,Multiple)
{
    int Z=GetParam();
    nlohmann::json js = {
        {"type",BasisSetAtom::Type::Gaussian},
        {"N", 20}, {"emin", 0.01}, {"emax", 4000*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SG_DFT_P,::testing::Values(1,3,7,37,53)); //,3,5,7,37,53

TEST_P(A_SL_DFT_P,Multiple)
{
    int Z=GetParam();
    int N=8;
    if (Z>20) N=10;
    if (Z>50) N=11;
     nlohmann::json js = {
        {"type",BasisSetAtom::Type::Slater},
        {"N", N}, {"emin", 0.31}, {"emax", 3*Z},
    };
    QchemTester::Init(1e-3,js);
    Iterate(dft_scf_params(Z));
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}

INSTANTIATE_TEST_CASE_P(Multiple,A_SL_DFT_P,::testing::Values(1,3,7,37,53)); 



TEST_P(A_PG_DFT_P,Multiple)
{
    Init();
    Iterate(dft_scf_params(GetParam()));
    EXPECT_LT(RelativeDFTError(),MaxRelErrE);
}
INSTANTIATE_TEST_CASE_P(Multiple,A_PG_DFT_P,::testing::Values(3,5,11,37)); //Z=51 is slow.


