#ifndef DFTTESTERSlater_H
#define DFTTESTERSlater_H

#include "SemiHFTester.H"

class DFTTesterSlater : public virtual BaseTester, private SemiHartreeFockTester
{
public:
    DFTTesterSlater();
    ~DFTTesterSlater();
    void Init(double Exchange);

    virtual HamiltonianTerm* GetVee() const;

private:
    mutable rc_ptr<IrrepBasisSet> itsCBasisSet; //Charge density fitting basis set.
};

#include "AtomTesterSlater.H"
class DFTAtomTesterSlater : public AtomTesterSlater, public DFTTesterSlater
{
public:
    DFTAtomTesterSlater() : AtomTesterSlater(10,.5,30.0) {};
    void Init(Atom* atom)
    {
        AtomTesterSlater::Init(atom);
    }
    void Init(Atom* atom, double Exchange, int Lmax, double spin)
    {
        DFTTesterSlater::Init(Exchange);
        AtomTesterSlater::Init(atom,Lmax,spin);
    }
    void Init(double Exchange, int NBasis, int Lmax, double spin, const LAParams& lap)
    {
        DFTTesterSlater::Init(Exchange);
        AtomTesterSlater::Init(NBasis,Lmax,spin,lap);
    }

};

#endif // DFTTESTERSlater_H

