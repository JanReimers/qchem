#ifndef BASETESTER_H
#define BASETESTER_H

#include "gtest/gtest.h"
#include "Misc/rc_ptr.h"

class HamiltonianTerm;
class BasisSet;
class Mesh;
class Hamiltonian;
class WaveFunction;
class SCFIterator;
class Cluster;
class BasisGroup;

//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class BaseTester
//: public virtual ::testing::Test
    : public ::testing::TestWithParam<int>
{
public:
    BaseTester();
    virtual ~BaseTester();
    void Init(double spin);

    virtual HamiltonianTerm* GetVee() const=0;
    virtual HamiltonianTerm* GetVxc(double spin) const=0;
    virtual BasisSet*  GetXbasisSet() const=0;
    virtual BasisSet*  GetCbasisSet() const=0;
    virtual Mesh*      GetIntegrationMesh() const=0;
    virtual void       LoadOrbitalBasisSet()=0;
    
    Cluster*      GetCluster     () const {return &*itsCluster;}
    BasisGroup*   GetBasisGroup  () const {return itsBasisGroup;}
    Hamiltonian*  GetHamiltonian () const {return itsHamiltonian;}
    WaveFunction* GetWaveFunction() const {return itsWaveFunction;}

    void   Iterate(double relax, double epsilon, int Niter);
    double TotalEnergy() const;

protected:
    rc_ptr<Cluster>        itsCluster;
    BasisGroup*            itsBasisGroup;
private:
    Hamiltonian*           itsHamiltonian;
    WaveFunction*          itsWaveFunction;
    SCFIterator*           itsSCFIterator;
};


#endif // BASETESTER_H

