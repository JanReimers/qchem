// File A_HF.C  Atom Hartree-Fock tests.

#include "QchemTester.H"

#include "Cluster/Molecule.H"
#include "Cluster/Atom.H"
#include <Hamiltonian/Factory.H>


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
using namespace HamiltonianF;

class M_PG_HF_U : public ::testing::Test
, public TestMolecule
{
public:
    M_PG_HF_U() {};
    void Init(Molecule* m)
    {
        TestMolecule::Init(m);
        nlohmann::json js = { {"filepath","../../../BasisSetData/dzvp.bsd"} };
        QchemTester::Init(3e-3,js);
    }
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const
    {
        return Factory(Model::HF,Pol::UnPolarized,cluster);
    }
};

class M_PG_DFT_U : public ::testing::Test
, public TestMolecule
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
        return Factory(Pol::UnPolarized,cluster,Alpha_N2,GetMeshParams(),itsBasisSet);
    }
};

bool verbose=false;

TEST_F(M_PG_HF_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-4,1e-7,1e-5,1.0,1e-4,verbose});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}

TEST_F(M_PG_DFT_U,N2)
{
    Init(MakeN2());
    Iterate({20,1e-4,1e-7,1e-5,1.0,1e-4,verbose});
    double rerr=fabs((TotalEnergy()-E_N2)/E_N2);
    EXPECT_LT(rerr,MaxRelErrE);
}
