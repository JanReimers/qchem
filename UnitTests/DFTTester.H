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
    DFTAtomTester() : AtomTester(37,.01,200000.0) {};
    void Init(Atom* atom)
    {
        AtomTester::Init(atom);
    }
    void Init(Atom* atom, double Exchange, int Lmax, double spin)
    {
        DFTTester::Init(Exchange);
        AtomTester::Init(atom,Lmax,spin);
    }
    void Init(double Exchange, int NBasis, int Lmax, double spin, const LinearAlgebraParams& lap)
    {
        DFTTester::Init(Exchange);
        AtomTester::Init(NBasis,Lmax,spin,lap);
    }

};

#include "MoleculeTester.H"
class DFTMoleculeTester : public MoleculeTester, public DFTTester
{
public:
    DFTMoleculeTester() : MoleculeTester() {};
    void Init(Molecule* m,double Exchange, double spin)
    {
        DFTTester::Init(Exchange);
        MoleculeTester::Init(m,spin);
    }
    void Init(Molecule* m,double Exchange, double spin, const LinearAlgebraParams& lap)
    {
        DFTTester::Init(Exchange);
        MoleculeTester::Init(m,spin,lap);
    }

};

#endif // DFTTESTER_H

