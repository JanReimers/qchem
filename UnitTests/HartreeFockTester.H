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

//#include "AtomTester.H"
//class HartreeFockAtomTester
//    : public HartreeFockTester
//    , public AtomTester
//
//{
//public:
//  HartreeFockAtomTester() : AtomTester(27,.01,20000.0) {};
//  void Init(Atom* atom);
//  void Init(Atom* atom, int Lmax, double spin);
//  void Init(int NBasis, int Lmax, double spin, const LAParams&);
//};
//
//#include "AtomTesterSlater.H"
//class HartreeFockAtomTesterSlater
//    : public HartreeFockTester
//    , public AtomTesterSlater
//
//{
//public:
//  HartreeFockAtomTesterSlater() : AtomTesterSlater(10,.5,30.0) {};
//  void Init(Atom* atom);
//  void Init(Atom* atom, int Lmax, double spin);
//  void Init(int NBasis, int Lmax, double spin, const LAParams&);
//};


#include "MoleculeTester.H"
class HartreeFockMoleculeTester : public MoleculeTester, public HartreeFockTester
{
public:
    HartreeFockMoleculeTester() : MoleculeTester() {}
    void Init(Molecule* m,double spin)
    {
        MoleculeTester::Init(m,spin);
    }
    void Init(Molecule* m,double spin, const LAParams& lap)
    {
        MoleculeTester::Init(m,spin,lap);
    }

};

#endif // HartreeFockTESTER_H

