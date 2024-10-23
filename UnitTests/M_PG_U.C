// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"

#include "Imp/Cluster/Molecule.H"
#include "Imp/Cluster/Atom.H"
#include "Imp/Hamiltonian/Hamiltonians.H"
#include <MeshParams.H>

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
        QchemTester::Init(3e-3);
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
        QchemTester::Init(1e-2);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        //MeshParams mp({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0,3});
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
        QchemTester::Init(1e-2);
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
    Iterate({20,1e-3,1.0,0.0,true});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}

TEST_F(M_PG_SHF_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-3,1.0,0.0,true});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}

TEST_F(M_PG_DFT_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-3,1.0,0.0,true});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}
