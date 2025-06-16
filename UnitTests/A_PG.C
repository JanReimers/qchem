// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"

#include "Cluster/Molecule.H"
#include "Cluster/Atom.H"
#include "Hamiltonians.H"
#include <Mesh/MeshParams.H>

Molecule* MakeN2()
{
    Molecule* N2=new Molecule();
    N2->Insert(new Atom(7,0,Vector3D<double>(-1.03,0,0)));
    N2->Insert(new Atom(7,0,Vector3D<double>( 1.04,0,0)));
    return N2;
}

double E_N2=-109.251;
double Alpha_N2=0.75197;
//
//  Un-polarized tests.
//
class M_PG_HF_U : public ::testing::Test
, public TestMolecule, public PG_OBasis, TestUnPolarized
{
public:
    M_PG_HF_U() {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-3,js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return new Ham_HF_U(cluster);
    }
};

class M_PG_SHF_U : public ::testing::Test
, public TestMolecule, public PG_OBasis, TestUnPolarized
{
public:
    M_PG_SHF_U() {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-2,js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        MeshParams mp({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0});
        return new Ham_SHF_U(cluster,Alpha_N2,GetMeshParams(),itsBasisSet);
    }
};
class M_PG_DFT_U : public ::testing::Test
, public TestMolecule, public PG_OBasis, TestUnPolarized
{
public:
    M_PG_DFT_U() {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(1e-2,js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        //MeshParams mp({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0});
        return new Ham_DFT_U(cluster,Alpha_N2,GetMeshParams(),itsBasisSet);
    }
};



TEST_F(M_PG_HF_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-4,1e-7,1e-5,1.0,1e-4,false});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}

TEST_F(M_PG_SHF_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-4,1e-7,1e-5,1.0,1e-4,false});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}

TEST_F(M_PG_DFT_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-4,1e-7,1e-5,1.0,1e-4,false});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}
