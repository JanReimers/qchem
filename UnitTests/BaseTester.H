#ifndef BASETESTER_H
#define BASETESTER_H

#include "gtest/gtest.h"
#include "SCFIterator/IterationParams.H"
#include "Imp/LAParams.H"
#include "Misc/rc_ptr.H"

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
    void Init(BasisSet*, double spin);
    void Init(BasisSet*, double spin,const LAParams&);

    virtual HamiltonianTerm* GetVee() const=0;
    virtual HamiltonianTerm* GetVxc(double spin) const=0;
    virtual IrrepBasisSet*  GetXbasisSet() const=0;
    virtual IrrepBasisSet*  GetCbasisSet() const=0;
    virtual rc_ptr<Mesh>  GetIntegrationMesh() const=0;
    virtual void       LoadOrbitalBasisSet()=0;
    
    Cluster*      GetCluster     () const {return &*itsCluster;}
    BasisSet*   GetBasisSet  () const {return itsBasisSet;}
    Hamiltonian*  GetHamiltonian () const {return itsHamiltonian;}
    WaveFunction* GetWaveFunction() const {return itsWaveFunction;}
    IntegralDataBase<double>* GetDatabase() const;

    void   Iterate(const SCFIterationParams&);
    double TotalEnergy() const;

protected:
    LAParams               itsLAParams;
    rc_ptr<Cluster>        itsCluster;
    BasisSet*              itsBasisSet;
private:
    Hamiltonian*           itsHamiltonian;
    WaveFunction*          itsWaveFunction;
    SCFIterator*           itsSCFIterator;
};


#endif // BASETESTER_H

