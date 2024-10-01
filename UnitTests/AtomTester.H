#ifndef ATOMTESTER_H
#define ATOMTESTER_H

#include "BaseTester.H"

class Atom;

class AtomTester : public virtual BaseTester
{
public:
    AtomTester(size_t n, double emin, double emax);
    void Init(Atom* atom, int Lmax, double spin);
    void Init(Atom* atom, int NBasis, int Lmax, double spin);
    virtual BasisSet*  GetXbasisSet() const;
    virtual BasisSet*  GetCbasisSet() const;
    virtual Mesh*      GetIntegrationMesh() const;
    virtual void LoadOrbitalBasisSet();

    int    itsNbasis;
    double itsEmin,itsEmax;
    int    itsLmax;
};

#endif // ATOMTESTER_H

