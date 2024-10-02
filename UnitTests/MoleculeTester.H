#ifndef MOLECULETESTER_H
#define MOLECULETESTER_H

#include "BaseTester.H"

class Molecule;

class MoleculeTester : public virtual BaseTester
{
public:
    MoleculeTester();
    void Init(Molecule*, double spin);
    void Init(Molecule*, double spin, const LinearAlgebraParams& lap);
    
    virtual BasisSet*  GetXbasisSet() const;
    virtual BasisSet*  GetCbasisSet() const;
    virtual Mesh*      GetIntegrationMesh() const;
    virtual void LoadOrbitalBasisSet();
    
    SCFIterationParams itsSCFIParams;

};

#endif // MOLECULETESTER_H

