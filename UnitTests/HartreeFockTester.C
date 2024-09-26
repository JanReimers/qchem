// File AtomTester.C Member functions for the atom tester class.

#include "HartreeFockTester.H"
#include "HamiltonianImplementation/HartreeFockVxc.H"
#include "HamiltonianImplementation/PolarizedHartreeFockVxc.H"
#include "HamiltonianImplementation/ExactVee.H"

HamiltonianTerm* HartreeFockTester::GetVee() const
{
    return new ExactVee;
}

HamiltonianTerm* HartreeFockTester::GetVxc(double spin) const
{
    HamiltonianTerm* ret=0;
    if (spin==0.0)
        ret=new HartreeFockVxc();
    else
        ret=new PolarizedHartreeFockVxc();
    return ret;
}

void HartreeFockAtomTester::Init(Atom* atom, int Lmax, double spin)
{
    AtomTester::Init(atom,Lmax,spin);
}
