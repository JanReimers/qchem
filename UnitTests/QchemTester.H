#ifndef QchemTester_H
#define QchemTester_H

#include "gtest/gtest.h"
#include "types.H"
#include "Misc/rc_ptr.H"
#include "SCFIterator/IterationParams.H"
#include "Imp/LAParams.H"
#include "Misc/PeriodicTable.H"


class QchemTester
{
public:
    QchemTester();
    virtual ~QchemTester() {};
    void   Init();
    void   Iterate(const SCFIterationParams&);
    double TotalEnergy() const;
    double RelativeHFError(bool quiet=false) const;

    
protected:
    // Atom of Molecule functions
    virtual Cluster*     GetCluster    () const=0;
    // Orbital Basis Set functions SG, PG, Slater
    virtual BasisSet*    GetBasisSet   (const Cluster*) const=0;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(const rc_ptr<Cluster>& cl) const=0;
    // Polarized or un-polarized
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const=0;
private:
    rc_ptr<Cluster> itsCluster;
    BasisSet*     itsBasisSet;
    Hamiltonian*  itsHamiltonian;
    WaveFunction* itsWaveFunction;
    SCFIterator*  itsSCFIterator;
public:
    PeriodicTable itsPT;
    double MaxRelErrE;
};

class TestHamiltonian : public virtual QchemTester
{
public:
    virtual Hamiltonian* GetHamiltonian(const rc_ptr<Cluster>& cl) const;
    virtual HamiltonianTerm* GetVee() const=0;
    virtual HamiltonianTerm* GetVxc() const=0;
};

class HFHamiltonian : public TestHamiltonian
{
public:
    virtual HamiltonianTerm* GetVee() const;
    virtual HamiltonianTerm* GetVxc() const;
};

class PolHFHamiltonian : public HFHamiltonian
{
public:
    virtual HamiltonianTerm* GetVxc() const;
};

class SG_OBasis : public virtual QchemTester
{
public:
    void Init(int _N, double _emin, double _emax, int _Lmax, const LAParams& _lap)
    {
        N=_N;Lmax=_Lmax;emin=_emin;emax=_emax;
        lap=_lap;
    }
    void Init(int _N, double _emin, double _emax, int _Lmax) 
    {
        Init(_N,_emin,_emax,_Lmax,LAParams({qchem::Lapack,qchem::SVD,1e-6,1e-12}));
    }
    
private:
    virtual BasisSet* GetBasisSet (const Cluster*) const;
    int N,Lmax;
    double emin,emax;
    LAParams lap;
};

class TestAtom : public virtual QchemTester
{
public:
    void Init(int _Z, int _q=0) {Z=_Z; q=_q;}
    virtual Cluster* GetCluster() const;
    
private:
    int Z,q; //Atomic # and charge.

};

class TestUnPolarized  : public virtual QchemTester
{
public:
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
};

class TestPolarized  : public virtual QchemTester
{
public:
    void Init(double _spin) {spin=_spin;}
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
private:
    double spin;
};

class A_SG_HF_U 
: public ::testing::Test, public TestAtom, public SG_OBasis, HFHamiltonian, TestUnPolarized
{
public:
    void Init(int Z,int N, double emin, double emax, int LMax)
    {
        TestAtom::Init(Z);
        SG_OBasis::Init(N,emin,emax,LMax);
        QchemTester::Init();
    }
};

TEST_F(A_SG_HF_U,He)
{
    Init(2,25,.1,10000,0);
    Iterate({40,1e-3,1.0,0.0});
    EXPECT_LT(RelativeHFError(),MaxRelErrE);
}

#endif //QchemTester_H
