#ifndef MOLECULETESTER_H
#define MOLECULETESTER_H

#include "BaseTester.H"

class Molecule;

class MoleculeTester : public virtual BaseTester
{
public:
    MoleculeTester();
    void Init(Molecule*, double spin);
    virtual BasisSet*  GetXbasisSet() const;
    virtual BasisSet*  GetCbasisSet() const;
    virtual Mesh*      GetIntegrationMesh() const;
    virtual void LoadOrbitalBasisSet();
};

#endif // MOLECULETESTER_H

