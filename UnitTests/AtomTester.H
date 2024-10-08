#ifndef ATOMTESTER_H
#define ATOMTESTER_H

#include "BaseTester.H"

class Atom;

class AtomTester : public virtual BaseTester
{
public:
    AtomTester(size_t n, double emin, double emax);
    void Init(Atom* atom); //Only call this once
    void Init(Atom* atom, int Lmax,double spin);
    void Init(int NBasis, int Lmax,double spin, const LinearAlgebraParams&);
    virtual IrrepBasisSet*  GetXbasisSet() const;
    virtual IrrepBasisSet*  GetCbasisSet() const;
    virtual Mesh*      GetIntegrationMesh() const;
    virtual void LoadOrbitalBasisSet();

    int    itsNbasis;
    double itsEmin,itsEmax;
    int    itsLmax;
    SCFIterationParams itsSCFIParams;
};

#endif // ATOMTESTER_H

