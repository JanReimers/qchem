#ifndef DFTTESTER_H
#define DFTTESTER_H

#include "SemiHFTester.H"

class DFTTester : public virtual BaseTester, private SemiHartreeFockTester
{
public:
    DFTTester();
    ~DFTTester();
    void Init(double Exchange);

    virtual HamiltonianTerm* GetVee() const;

private:
    mutable rc_ptr<BasisSet>       itsCBasisSet; //Charge density fitting basis set.
};

#include "AtomTester.H"
class DFTAtomTester : public AtomTester, public DFTTester
{
public:
    DFTAtomTester() : AtomTester(17,.01,200000.0) {};
    void Init(Atom* atom, double Exchange, int Lmax, double spin)
    {
        DFTTester::Init(Exchange);
        AtomTester::Init(atom,Lmax,spin);
    }

};

#include "MoleculeTester.H"
class DFTMoleculeTester : public MoleculeTester, public DFTTester
{
public:
   void Init(Molecule* m,double Exchange, double spin)
    {
        DFTTester::Init(Exchange);
        MoleculeTester::Init(m,spin);
    }

};

#endif // DFTTESTER_H

