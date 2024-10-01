#ifndef SEMIHFTESTER_H
#define SEMIHFTESTER_H

#include "BaseTester.H"

class SemiHartreeFockTester : public virtual BaseTester
{
public:
    SemiHartreeFockTester();
    virtual ~SemiHartreeFockTester();

    void Init(double Exchange);

    virtual HamiltonianTerm* GetVee() const;
    virtual HamiltonianTerm* GetVxc(double spin) const;

private:
    mutable rc_ptr<BasisSet>  itsXBasisSet; //Exchange HamiltonianTerm fitting basis set.
    double                    itsExchange;
};

#include "AtomTester.H"
// Use HF for coulomb potential and DFT for exchange potential
class SemiHartreeFockAtomTester : public AtomTester, public SemiHartreeFockTester
{
public:
    SemiHartreeFockAtomTester() : AtomTester(17,.01,200000.0) {};
    void Init(Atom* atom, double Exchange, int Lmax, double spin)
    {
        SemiHartreeFockTester::Init(Exchange);
        AtomTester::Init(atom,Lmax,spin);
    }
};

#include "MoleculeTester.H"
class SemiHartreeFockMoleculeTester : public MoleculeTester, public SemiHartreeFockTester
{
public:
   void Init(Molecule* m,double Exchange, double spin)
    {
        SemiHartreeFockTester::Init(Exchange);
        MoleculeTester::Init(m,spin);
    }

};

#endif // SEMIHFTESTER_H

