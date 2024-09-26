// File: MoleculeTests.C  Test the SCF calculation for molecules

#include "DFTTester.H"
#include "HartreeFockTester.H"
#include "DFTDataBase/DFTDataBase.H"
#include "Cluster/Atom.H"
#include "Cluster/Molecule.H"
#include "Misc/PeriodicTable.H"
#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"

#include <iostream>
#include <fstream>

using std::cout;
using std::endl;

//DFTDataBase  theDataBase("Atom.db");
extern PeriodicTable thePeriodicTable;
extern double eps_ro; //Converge criterial for delta ro (charge density)
double m_eps_e=1e-2;


TEST_P(DFTMoleculeTester, MoleculesDFTPolarized)
{
    int Z=GetParam();
    std::cout << "Testing atom " << thePeriodicTable.GetSymbol(Z) << std::endl;
    RadialMesh*  rm=new MHLRadialMesh(200,2U,5.0);
    AngularMesh* am=new GaussAngularMesh(12);
    Mesh* atom_mesh=new AtomMesh(*rm,*am);
    Atom* a1=new Atom(Z,0,Vector3D<double>(0,0,0));
    a1->SetMesh(atom_mesh);
    Molecule* Hm=new Molecule();
    Hm->Insert(a1);
    Init(Hm,thePeriodicTable.GetSlaterAlpha(Z),thePeriodicTable.GetNumUnpairedElectrons(Z));
    Iterate(1.0,eps_ro,40);
    double E_DFT=thePeriodicTable.GetEnergyDFT(Z);
    double error=fabs((E_DFT-TotalEnergy())/E_DFT);
    std::cout.precision(5);
    std::cout << "E_DFT relative error=" << error*100.0 << "%" << std::endl;
    EXPECT_LT(error,m_eps_e);
}


TEST_P(HartreeFockMoleculeTester, MoleculesHFPolarized)
{
    int Z=GetParam();
    std::cout << "Testing atom " << thePeriodicTable.GetSymbol(Z) << std::endl;
    Atom* a1=new Atom(Z,0,Vector3D<double>(0,0,0));
    Molecule* Hm=new Molecule();
    Hm->Insert(a1);
    Init(Hm,thePeriodicTable.GetNumUnpairedElectrons(Z));
    Iterate(1.0,eps_ro,40);
    double E_HF=thePeriodicTable.GetEnergyHF(Z);
    double error=fabs((E_HF-TotalEnergy())/E_HF);
    std::cout.precision(5);
    std::cout << "E_HF relative error=" << error*100.0 << "%" << std::endl;
    EXPECT_LT(error,m_eps_e);
}

INSTANTIATE_TEST_CASE_P(MoleculesDFTPolarized,
                        DFTMoleculeTester,
                        ::testing::Values(1,4,5,7));

INSTANTIATE_TEST_CASE_P(MoleculesHFPolarized,
                        HartreeFockMoleculeTester,
                        ::testing::Values(1,4,5,7,10,12,25));

//--------------------------------------------------------------------------------
//
//  N2 tests.
//
double E_N2=-109.251;

TEST_F(DFTMoleculeTester, N2)
{
    RadialMesh*  rm=new MHLRadialMesh(100,2U,2.0);
    AngularMesh* am=new GaussAngularMesh(12);
    Mesh* atom_mesh=new AtomMesh(*rm,*am);
    Atom* a1=new Atom(7,0,Vector3D<double>(-1.03,0,0));
    Atom* a2=new Atom(7,0,Vector3D<double>( 1.04,0,0));
    a1->SetMesh(atom_mesh);
    a2->SetMesh(atom_mesh);
    Molecule* N2=new Molecule();
    N2->Insert(a1);
    N2->Insert(a2);
    Init(N2,0.75197,0.0);
    Iterate(0.5,eps_ro,60);
    EXPECT_LT(fabs((E_N2-TotalEnergy())/E_N2), m_eps_e);
}

TEST_F(SemiHartreeFockMoleculeTester, N2)
{
    RadialMesh*  rm=new MHLRadialMesh(100,2U,2.0);
    AngularMesh* am=new GaussAngularMesh(12);
    Mesh* atom_mesh=new AtomMesh(*rm,*am);
    Atom* a1=new Atom(7,0,Vector3D<double>(-1.03,0,0));
    Atom* a2=new Atom(7,0,Vector3D<double>( 1.04,0,0));
    a1->SetMesh(atom_mesh);
    a2->SetMesh(atom_mesh);
    Molecule* N2=new Molecule();
    N2->Insert(a1);
    N2->Insert(a2);
    Init(N2,0.75197,0.0);
    Iterate(0.5,eps_ro,60);
    EXPECT_LT(fabs((E_N2-TotalEnergy())/E_N2), m_eps_e);
}

TEST_F(HartreeFockMoleculeTester, N2)
{
    Atom* a1=new Atom(7,0,Vector3D<double>(-1.03,0,0));
    Atom* a2=new Atom(7,0,Vector3D<double>( 1.04,0,0));
    Molecule* N2=new Molecule();
    N2->Insert(a1);
    N2->Insert(a2);
    Init(N2,0.0);
    Iterate(0.5,eps_ro,60);
    EXPECT_LT(fabs((E_N2-TotalEnergy())/E_N2), m_eps_e);
}
