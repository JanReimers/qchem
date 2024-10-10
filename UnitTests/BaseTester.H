#ifndef BASETESTER_H
#define BASETESTER_H

#include "gtest/gtest.h"
#include "SCFIterator/IterationParams.H"
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
    void Init(BasisGroup*, double spin);
    void Init(BasisGroup*, double spin,const LinearAlgebraParams&);

    virtual HamiltonianTerm* GetVee() const=0;
    virtual HamiltonianTerm* GetVxc(double spin) const=0;
    virtual IrrepBasisSet*  GetXbasisSet() const=0;
    virtual IrrepBasisSet*  GetCbasisSet() const=0;
    virtual Mesh*      GetIntegrationMesh() const=0;
    virtual void       LoadOrbitalBasisSet()=0;
    
    Cluster*      GetCluster     () const {return &*itsCluster;}
    BasisGroup*   GetBasisGroup  () const {return itsBasisGroup;}
    Hamiltonian*  GetHamiltonian () const {return itsHamiltonian;}
    WaveFunction* GetWaveFunction() const {return itsWaveFunction;}
    IntegralDataBase<double>* GetDatabase() const;

    void   Iterate(const SCFIterationParams&);
    double TotalEnergy() const;

protected:
    LinearAlgebraParams    itsLAParams;
    rc_ptr<Cluster>        itsCluster;
    BasisGroup*            itsBasisGroup;
private:
    Hamiltonian*           itsHamiltonian;
    WaveFunction*          itsWaveFunction;
    SCFIterator*           itsSCFIterator;
};


#endif // BASETESTER_H

