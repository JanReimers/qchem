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
    virtual ~QchemTester();
    void   Init();
    void   Iterate(const SCFIterationParams&);
    double TotalEnergy() const;
    double RelativeHFError(bool quiet=false) const;
    double RelativeDFTError(bool quiet=false) const;
    int    GetLMax                (int Z) const {return itsPT.GetMaxL(Z);}
    double GetNumUnpairedElectrons(int Z) const {return itsPT.GetNumUnpairedElectrons(Z);} 
    double GetExParam             (int Z) const {return itsPT.GetSlaterAlpha(Z);}
    
protected:
    // Atom of Molecule functions
    virtual Cluster*     GetCluster        () const=0;
    virtual Mesh*        GetIntegrationMesh() const=0;
    // Orbital Basis Set functions SG, PG, Slater
    virtual BasisSet*      GetBasisSet   (const Cluster*) const=0;
    virtual IrrepBasisSet* GetXBasisSet  () const=0;
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

//Why??
#include "Hamiltonian/ExchangeFunctional.H" 
#include "BasisSet.H" 

class SHFHamiltonian : public HFHamiltonian
{
public:
    virtual ~SHFHamiltonian();
    void Init(double exparam);
    virtual HamiltonianTerm* GetVxc() const;
protected:
    rc_ptr<ExchangeFunctional> XcFunct;
};

class PolSHFHamiltonian : public SHFHamiltonian
{
public:
    void Init(double exparam);
    virtual HamiltonianTerm* GetVxc() const;
private:
    rc_ptr<ExchangeFunctional> XcFunct;
};

//
//  Basis sets
//
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
    virtual BasisSet*      GetBasisSet (const Cluster*) const;
    virtual IrrepBasisSet* GetXBasisSet() const;
    int N,Lmax;
    double emin,emax;
    LAParams lap;
    mutable IntegralDataBase<double>* idb;
};

class SL_OBasis : public virtual QchemTester
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
    virtual BasisSet*      GetBasisSet (const Cluster*) const;
    virtual IrrepBasisSet* GetXBasisSet() const;
    int N,Lmax;
    double emin,emax;
    LAParams lap;
    mutable IntegralDataBase<double>* idb;
};

class TestAtom : public virtual QchemTester
{
public:
    void Init(int _Z, int _q=0) {Z=_Z; q=_q;}
    virtual Cluster* GetCluster() const;
    virtual Mesh*    GetIntegrationMesh() const;

    
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
    TestPolarized() : spin(-1.0) {};
    void Init(double _spin) {spin=_spin;}
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
private:
    double spin;
};



#endif //QchemTester_H
