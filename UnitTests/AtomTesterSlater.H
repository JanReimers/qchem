#ifndef ATOMTESTERSlater_H
#define ATOMTESTERSlater_H

#include "BaseTester.H"

class Atom;

class AtomTesterSlater : public virtual BaseTester
{
public:
    AtomTesterSlater(size_t n, double emin, double emax);
    void Init(Atom* atom); //Only call this once
    void Init(Atom* atom, int Lmax,double spin);
    void Init(int NBasis, int Lmax,double spin, const LAParams&);
    virtual IrrepBasisSet*  GetXbasisSet() const;
    virtual IrrepBasisSet*  GetCbasisSet() const;
    virtual rc_ptr<Mesh>    GetIntegrationMesh() const;
    virtual void LoadOrbitalBasisSet() {};

    int    itsNbasis;
    double itsEmin,itsEmax;
    int    itsLmax;
    SCFIterationParams itsSCFIParams;
};

#endif // ATOMTESTERSlater_H

