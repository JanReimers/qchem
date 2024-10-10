// File AtomTester.C Member functions for the atom tester class.

#include "MoleculeTester.H"
#include "Cluster/Molecule.H"
#include "DFTDataBase/HeapDB/HeapDB.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Gaussian94.H"
#include "Mesh/MoleculeMesh.H"


MoleculeTester::MoleculeTester() 
: BaseTester() 
, itsSCFIParams({40,1e-3,1.0,0.0})
{};

void MoleculeTester::Init(Molecule* m,double spin)
{
    itsCluster.reset(m);
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
    auto bg=new BasisGroup();
    IrrepBasisSet* bs = new PolarizedGaussian::BasisSet(itsLAParams,bg->GetDataBase(), &reader,itsCluster.get());
    bg->Insert(bs); 
    BaseTester::Init(bg,spin);
}

void MoleculeTester::Init(Molecule* m,double spin,const LinearAlgebraParams& lap)
{
    itsCluster.reset(m);
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
    auto bg=new BasisGroup();
    IrrepBasisSet* bs = new PolarizedGaussian::BasisSet(itsLAParams,bg->GetDataBase(), &reader,itsCluster.get());
    bg->Insert(bs); 
    BaseTester::Init(bg,spin,lap);
}

void MoleculeTester::LoadOrbitalBasisSet()
{
//    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
//    IrrepBasisSet* bs = new PolarizedGaussian::BasisSet(itsLAParams,GetDatabase(), &reader,itsCluster.get());
//    itsBasisGroup->Insert(bs);
}

IrrepBasisSet* MoleculeTester::GetCbasisSet() const
{
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/A2_coul.bsd");
    IrrepBasisSet* bs = new PolarizedGaussian::BasisSet(itsLAParams,new HeapDB<double>, &reader,itsCluster.get());
    return bs;
}

IrrepBasisSet* MoleculeTester::GetXbasisSet() const
{
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/A1_exch.bsd");
    IrrepBasisSet* bs = new PolarizedGaussian::BasisSet(itsLAParams,new HeapDB<double>, &reader,itsCluster.get());
    return bs;
}

Mesh* MoleculeTester::GetIntegrationMesh() const
{
   	Mesh*	me = new MoleculeMesh(*itsCluster,3);
    return me;
}

