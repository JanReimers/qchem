// File AtomTester.C Member functions for the atom tester class.

#include "SemiHFTester.H"
#include "HamiltonianImplementation/FittedVxc.H"
//#include "HamiltonianImplementation/LDAVxc.H"
#include "HamiltonianImplementation/PolarizedFittedVxc.H"
#include "HamiltonianImplementation/SlaterExchange.H"
#include "HamiltonianImplementation/ExactVee.H"
#include "Mesh/Mesh.H"

SemiHartreeFockTester::SemiHartreeFockTester()
    : itsExchange(0.0)
{
};

SemiHartreeFockTester::~SemiHartreeFockTester() {};

void SemiHartreeFockTester::Init(double Exchange)
{
    itsExchange=Exchange;
}

HamiltonianTerm* SemiHartreeFockTester::GetVee() const
{
    return new ExactVee;
}

HamiltonianTerm* SemiHartreeFockTester::GetVxc(double spin) const
{
    if(!itsXBasisSet.get()) itsXBasisSet.reset(GetXbasisSet());
    
    HamiltonianTerm* ret=0;
    if (spin==0.0)
        ret=new FittedVxc(itsXBasisSet, new SlaterExchange(itsExchange),GetIntegrationMesh());
    else
        ret=new PolarizedFittedVxc(itsXBasisSet, new SlaterExchange(itsExchange,Spin(Spin::Up)),GetIntegrationMesh());
    return ret;
}
