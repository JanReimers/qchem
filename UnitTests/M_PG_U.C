// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"

#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"
#include "Imp/Cluster/Molecule.H"

Molecule* MakeN2()
{
    RadialMesh*  rm=new MHLRadialMesh(100,2U,2.0);
    Mesh* am=new GaussAngularMesh(12);
    Mesh* atom_mesh=new AtomMesh(*rm,*am);
    Atom* a1=new Atom(7,0,Vector3D<double>(-1.03,0,0));
    Atom* a2=new Atom(7,0,Vector3D<double>( 1.04,0,0));
    a1->SetMesh(atom_mesh);
    a2->SetMesh(atom_mesh);
    Molecule* N2=new Molecule();
    N2->Insert(a1);
    N2->Insert(a2);
    return N2;
}

double E_N2=-109.251;
double Alpha_N2=0.75197;
//
//  Un-polarized tests.
//
class M_PG_HF_U : public ::testing::Test
, public TestMolecule, public PG_OBasis, HFHamiltonian, TestUnPolarized
{
public:
    M_PG_HF_U() {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        QchemTester::Init(1e-2);
    }
};

class M_PG_SHF_U : public ::testing::Test
, public TestMolecule, public PG_OBasis, SHFHamiltonian, TestUnPolarized
{
public:
    M_PG_SHF_U() : SHFHamiltonian (Alpha_N2) {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        QchemTester::Init(1e-2);
    }
};
class M_PG_DFT_U : public ::testing::Test
, public TestMolecule, public PG_OBasis, DFTHamiltonian, TestUnPolarized
{
public:
    M_PG_DFT_U() : DFTHamiltonian (Alpha_N2) {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        QchemTester::Init(1e-2);
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
