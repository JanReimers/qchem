// File AtomTester.C Member functions for the atom tester class.

#include "DFTTester.H"
#include "HamiltonianImplementation/CDFittedVee.H"
#include "Mesh/Mesh.H"
#include "Cluster.H"
#include "BasisSet.H"

DFTTester::DFTTester()
    : itsCBasisSet()
{
};

DFTTester::~DFTTester() {};

void DFTTester::Init(double Exchange)
{
    SemiHartreeFockTester::Init(Exchange);
}

HamiltonianTerm* DFTTester::GetVee() const
{
    if(!itsCBasisSet.get()) itsCBasisSet.reset(GetCbasisSet());
    rc_ptr<Mesh>   mesh=GetIntegrationMesh();
    return new CDFittedVee(itsCBasisSet,mesh,itsCluster->GetNumElectrons());
}
