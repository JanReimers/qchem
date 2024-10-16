// File AtomTester.C Member functions for the atom tester class.

#include "DFTTesterSlater.H"
#include "HamiltonianImplementation/CDFittedVee.H"
#include "Mesh/Mesh.H"
#include "Cluster.H"
#include "BasisSet.H"

DFTTesterSlater::DFTTesterSlater()
    : itsCBasisSet()
{
};

DFTTesterSlater::~DFTTesterSlater() {};

void DFTTesterSlater::Init(double Exchange)
{
    SemiHartreeFockTester::Init(Exchange);
}

HamiltonianTerm* DFTTesterSlater::GetVee() const
{
    if(!itsCBasisSet.get()) itsCBasisSet.reset(GetCbasisSet());
    rc_ptr<Mesh>   mesh=GetIntegrationMesh();
    return new CDFittedVee(itsCBasisSet,mesh,itsCluster->GetNumElectrons());
}
