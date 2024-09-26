// File AtomTester.C Member functions for the atom tester class.

#include "DFTTester.H"
#include "HamiltonianImplementation/CDFittedVee.H"
#include "Cluster/Cluster.H"
#include "BasisSet/BasisSet.H"

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
    return new CDFittedVee(itsCBasisSet,itsCluster->GetNumElectrons());
}
