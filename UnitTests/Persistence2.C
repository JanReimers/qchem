#include "AtomTester.H"
#include "HartreeFockTester.H"
#include "gtest/gtest.h"
#include "oml/imp/stream.h"
#include "Cluster/Atom.H"
#include "Cluster/Molecule.H"
#include "Hamiltonian/Hamiltonian.H"
#include "BasisSet/BasisGroup.H"
#include "WaveFunction/WaveFunction.H"
#include "Misc/PeriodicTable.H"
#include <iostream> 
#include <fstream>

using std::cout;
using std::endl;

extern PeriodicTable thePeriodicTable;

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class qchem_PersistanceTests : public AtomTester, public HartreeFockTester
{
public:
    qchem_PersistanceTests() 
        : AtomTester(17,.01,20000.)
        , file_name("streamtest2.dat") 
    {
        
    }
    
    void Init(Atom* atom)
    {
        AtomTester::Init(atom);
    }
    void Init(Atom* atom, int Lmax, double spin)
    {
        AtomTester::Init(atom,Lmax,spin);
    }
    void Init(int NBasis, int Lmax, double spin, const LinearAlgebraParams& lap)
    {
        AtomTester::Init(NBasis,Lmax,spin,lap);
    }
//    void OutIn(const PMStreamableObject* pout,PMStreamableObject* pin,StreamableObject::Mode);
    
    std::string file_name;
};


template <class T> T* OutIn(const PMStreamableObject* pout,c_str filename,StreamableObject::Mode m)
{
    StreamableObject::Mode current=StreamableObject::SetOutputMode(m);
    std::ofstream out(filename);
    out << *pout;
    EXPECT_TRUE(out);
    out.close();
    
    std::ifstream in(filename);
    T* pin=T::Factory(in);
    in >> *pin;
    EXPECT_TRUE(in);
    in.close();
    StreamableObject::SetOutputMode(current);
    return pin;
}

    
TEST_F(qchem_PersistanceTests,AtomBasisSets)
{
    
    int Z=51;
    int L=thePeriodicTable.GetMaxL(Z),NBasis=6;
    double spin=thePeriodicTable.GetNumUnpairedElectrons(Z);
    Atom* atom=new Atom(Z,0,Vector3D<double>(0,0,0));
    Init(atom);
    LinearAlgebraParams lap={qchem::Lapack,qchem::SVD,1e-6,1e-12};
    Init(NBasis,L,spin,lap);
    StreamableObject::SetToPretty();
    
    Cluster* cl=GetCluster();
    Cluster* cl1=::OutIn<Cluster>(cl,file_name.c_str(), StreamableObject::ascii);
    cout << *cl << endl << *cl1 << endl;
    Cluster* cl2=::OutIn<Cluster>(cl,file_name.c_str(), StreamableObject::binary);
    cout << *cl2 << endl;
    delete cl1;
    delete cl2;
    
    
    Hamiltonian* h=GetHamiltonian();
    Hamiltonian* h1=::OutIn<Hamiltonian>(h,file_name.c_str(), StreamableObject::ascii);
    Hamiltonian* h2=::OutIn<Hamiltonian>(h,file_name.c_str(), StreamableObject::binary);

    StreamableObject::SetToPretty();
    cout << *h << *h1 << *h2 << endl;
    
    BasisGroup* bg=GetBasisGroup();
    BasisGroup* bg1=::OutIn<BasisGroup>(bg,file_name.c_str(), StreamableObject::ascii);
    BasisGroup* bg2=::OutIn<BasisGroup>(bg,file_name.c_str(), StreamableObject::binary);
    cout << *bg << *bg1 << *bg2 << endl;
    
    WaveFunction* wf=GetWaveFunction();
    WaveFunction* wf1=::OutIn<WaveFunction>(wf,file_name.c_str(), StreamableObject::ascii);
    WaveFunction* wf2=::OutIn<WaveFunction>(wf,file_name.c_str(), StreamableObject::binary);
    cout << *wf << *wf1 << *wf2 << endl;
}
