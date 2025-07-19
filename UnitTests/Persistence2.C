#include "QchemTester.H"

#include <Cluster/Cluster.H>
#include <Hamiltonian/Hamiltonian.H>
#include <iostream> 
#include <fstream>

using std::cout;
using std::endl;

int Z=2;
//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class qchem_PersistanceTests : public ::testing::Test
, public TestAtom, public SG_OBasis, HFHamiltonian, TestUnPolarized
{
public:
    qchem_PersistanceTests() 
        : TestAtom(::Z)
        , file_name("streamtest2.dat") 
    {}
    
    void Init(int N, double emin, double emax, int LMax)
    {
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init(1e-3);
    }
    
    std::string file_name;
};


template <class T> T* OutIn(const Streamable* pout,c_str filename,StreamableObject::Mode m)
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
    int L=itsPT.GetMaxL(::Z);
    Init(6,.01,10000.0,L);
    StreamableObject::SetToPretty();
    
    Cluster* cl=GetCluster();
    Cluster* cl1=::OutIn<Cluster>(cl,file_name.c_str(), StreamableObject::ascii);
    cout << *cl << endl << *cl1 << endl;
    Cluster* cl2=::OutIn<Cluster>(cl,file_name.c_str(), StreamableObject::binary);
    cout << *cl2 << endl;
    delete cl1;
    delete cl2;
    
    
    Hamiltonian* h=itsHamiltonian;
    Hamiltonian* h1=::OutIn<Hamiltonian>(h,file_name.c_str(), StreamableObject::ascii);
    Hamiltonian* h2=::OutIn<Hamiltonian>(h,file_name.c_str(), StreamableObject::binary);

    StreamableObject::SetToPretty();
    cout << *h << *h1 << *h2 << endl;
    delete h1;
    delete h2;
    
//    BasisGroup* bg=GetBasisGroup();
//    BasisGroup* bg1=::OutIn<BasisGroup>(bg,file_name.c_str(), StreamableObject::ascii);
//    BasisGroup* bg2=::OutIn<BasisGroup>(bg,file_name.c_str(), StreamableObject::binary);
//    cout << *bg << *bg1 << *bg2 << endl;
//    
//    WaveFunction* wf=GetWaveFunction();
//    WaveFunction* wf1=::OutIn<WaveFunction>(wf,file_name.c_str(), StreamableObject::ascii);
//    WaveFunction* wf2=::OutIn<WaveFunction>(wf,file_name.c_str(), StreamableObject::binary);
//    cout << *wf << *wf1 << *wf2 << endl;
}
