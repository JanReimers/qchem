#ifndef HartreeFockTESTER_H
#define HartreeFockTESTER_H

#include "BaseTester.H"

// Use HF for coulomb potential and DFT for exchange potential
class HartreeFockTester : public virtual BaseTester
{
public:
    virtual HamiltonianTerm* GetVee() const;
    virtual HamiltonianTerm* GetVxc(double spin) const;
};

#include "AtomTester.H"
class HartreeFockAtomTester
    : public HartreeFockTester
    , public AtomTester

{
public:
  void Init(Atom* atom, int Lmax, double spin);
};

#include "MoleculeTester.H"
class HartreeFockMoleculeTester : public MoleculeTester, public HartreeFockTester
{
public:
   void Init(Molecule* m, double spin)
    {
        MoleculeTester::Init(m,spin);
    }

};

#endif // HartreeFockTESTER_H

