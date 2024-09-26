// File AtomTester.C Member functions for the atom tester class.

#include "MoleculeTester.H"
#include "Cluster/Molecule.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBS.H"
#include "BasisSet/BasisGroup.H"
#include "Gaussian94RFR/Gaussian94RFR.H"
#include "Mesh/MoleculeMesh.H"


MoleculeTester::MoleculeTester()
{
}

void MoleculeTester::Init(Molecule* m,double spin)
{
    itsCluster.reset(m);
    BaseTester::Init(spin);

}
void MoleculeTester::LoadOrbitalBasisSet()
{
    Gaussian94RFR reader("../BasisSetData/dzvp.bsd");
    BasisSet* bs = new PolarizedGaussianBS(new HeapDB<double>, &reader,itsCluster.get());
    itsBasisGroup->Insert(bs);
}

BasisSet* MoleculeTester::GetCbasisSet() const
{
    Gaussian94RFR reader("../BasisSetData/A2_coul.bsd");
    BasisSet* bs = new PolarizedGaussianBS(new HeapDB<double>, &reader,itsCluster.get());
    return bs;
}

BasisSet* MoleculeTester::GetXbasisSet() const
{
    Mesh* mesh=GetIntegrationMesh();
    assert(mesh);
    Gaussian94RFR reader("../BasisSetData/A1_exch.bsd");
    BasisSet* bs = new PolarizedGaussianBS(new HeapDB<double>, &reader,itsCluster.get(),mesh);
    return bs;
}

Mesh* MoleculeTester::GetIntegrationMesh() const
{
   	Mesh*	me = new MoleculeMesh(*itsCluster,3);
    return me;
}

