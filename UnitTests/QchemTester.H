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
    void   Init(double eps);
    void   Iterate(const SCFIterationParams&);
    double TotalEnergy() const;
    double RelativeHFError(bool quiet=false) const;
    double RelativeDFTError(bool quiet=false) const;
    int    GetLMax(int Z) const {return itsPT.GetMaxL(Z);}
    
protected:
    // Atom of Molecule functions
    virtual Cluster*     GetCluster        () const=0;
    virtual Mesh*        GetIntegrationMesh() const=0;
    virtual int          GetZ              () const;
    // Orbital Basis Set functions SG, PG, Slater
    virtual BasisSet*      GetBasisSet   (const Cluster*) const=0;
    virtual IrrepBasisSet* GetCBasisSet  () const=0;
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
    static PeriodicTable itsPT;
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
    SHFHamiltonian(int Z);
    virtual ~SHFHamiltonian();
    virtual HamiltonianTerm* GetVxc() const;
protected:
    rc_ptr<ExchangeFunctional> XcFunct;
};

class PolSHFHamiltonian : public HFHamiltonian
{
public:
    PolSHFHamiltonian(int Z);
    virtual HamiltonianTerm* GetVxc() const;
private:
    rc_ptr<ExchangeFunctional> XcFunct;
};

class DFTHamiltonian : public SHFHamiltonian
{
public:
    DFTHamiltonian(int Z) : SHFHamiltonian(Z) {};
    virtual HamiltonianTerm* GetVee() const;
};

class PolDFTHamiltonian : public PolSHFHamiltonian
{
public:
    PolDFTHamiltonian(int Z) : PolSHFHamiltonian(Z) {};
    virtual HamiltonianTerm* GetVee() const;
};

//------------------------------------------------------------------------------------
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
        Init(_N,_emin,_emax,_Lmax,LAParams({qchem::Lapack,qchem::SVD,1e-4,1e-12}));
    }
     
private:
    virtual BasisSet*      GetBasisSet (const Cluster*) const;
    virtual IrrepBasisSet* GetCBasisSet() const;
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
        Init(_N,_emin,_emax,_Lmax,LAParams({qchem::Lapack,qchem::SVD,1e-4,1e-12}));
    }
    
private:
    virtual BasisSet*      GetBasisSet (const Cluster*) const;
    virtual IrrepBasisSet* GetCBasisSet() const;
    virtual IrrepBasisSet* GetXBasisSet() const;
    int N,Lmax;
    double emin,emax;
    LAParams lap;
    mutable IntegralDataBase<double>* idb;
};

class PG_OBasis : public virtual QchemTester
{
public:
    PG_OBasis() : lap({qchem::Lapack,qchem::SVD,1e-4,1e-12}) {};     
private:
    virtual BasisSet*      GetBasisSet (const Cluster*) const;
    virtual IrrepBasisSet* GetCBasisSet() const;
    virtual IrrepBasisSet* GetXBasisSet() const;
    LAParams lap;
    mutable IntegralDataBase<double>* idb;
};

//----------------------------------------------------------------------------------
//
//  Atoms and Molecules
//
class TestAtom : public virtual QchemTester
{
public:
    TestAtom(int _Z, int _q=0) : Z(_Z), q(_q) {};
    virtual Cluster* GetCluster() const;
    virtual Mesh*    GetIntegrationMesh() const;
    virtual int      GetZ              () const {return Z;}
private:
    int Z,q; //Atomic # and charge.

};

class TestMolecule : public virtual QchemTester
{
public:
    TestMolecule() : itsCluster(0) {};
    void Init(Molecule*);
    
    virtual Cluster* GetCluster() const;
    virtual Mesh*    GetIntegrationMesh() const;
private:
    Molecule* itsCluster;
};


//--------------------------------------------------------------------------------
//
//  Wave functions
//
class TestUnPolarized  : public virtual QchemTester
{
public:
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
};

class TestPolarized  : public virtual QchemTester
{
public:
    TestPolarized(int Z);
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
private:
    double spin;
};



#endif //QchemTester_H
