// File AtomTester.C Member functions for the atom tester class.

#include "MoleculeTester.H"
#include "Cluster/Molecule.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBS.H"
#include "BasisSet.H"
#include "BasisSetImplementation/Gaussian94RFR.H"
#include "Mesh/MoleculeMesh.H"


MoleculeTester::MoleculeTester() 
: BaseTester() 
, itsSCFIParams({40,1e-3,1.0,0.0})
{};

void MoleculeTester::Init(Molecule* m,double spin)
{
    itsCluster.reset(m);
    BaseTester::Init(spin);
}

void MoleculeTester::Init(Molecule* m,double spin,const LinearAlgebraParams& lap)
{
    itsCluster.reset(m);
    BaseTester::Init(spin,lap);
}

void MoleculeTester::LoadOrbitalBasisSet()
{
    Gaussian94RFR reader("../BasisSetData/dzvp.bsd");
    IrrepBasisSet* bs = new PolarizedGaussianBS(itsLAParams,GetDatabase(), &reader,itsCluster.get());
    itsBasisGroup->Insert(bs);
}

IrrepBasisSet* MoleculeTester::GetCbasisSet() const
{
    Gaussian94RFR reader("../BasisSetData/A2_coul.bsd");
    IrrepBasisSet* bs = new PolarizedGaussianBS(itsLAParams,new HeapDB<double>, &reader,itsCluster.get());
    return bs;
}

IrrepBasisSet* MoleculeTester::GetXbasisSet() const
{
    Gaussian94RFR reader("../BasisSetData/A1_exch.bsd");
    IrrepBasisSet* bs = new PolarizedGaussianBS(itsLAParams,new HeapDB<double>, &reader,itsCluster.get());
    return bs;
}

Mesh* MoleculeTester::GetIntegrationMesh() const
{
   	Mesh*	me = new MoleculeMesh(*itsCluster,3);
    return me;
}

