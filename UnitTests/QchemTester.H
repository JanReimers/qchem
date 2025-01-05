#ifndef QchemTester_H
#define QchemTester_H

#include "gtest/gtest.h"
#include <types.H>
#include <IterationParams.H>
#include <LAParams.H>
#include "Imp/Misc/PeriodicTable.H"
#include <memory>

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
    typedef std::shared_ptr<const Cluster> cl_t;

    // Atom of Molecule functions
    virtual const Cluster*  GetCluster   () const {return itsCluster.get();}
    virtual MeshParams      GetMeshParams() const=0;
    virtual int             GetZ         () const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const=0;

    // Orbital Basis Set functions SG, PG, Slater
    virtual BasisSet*      GetBasisSet   () const=0;
    virtual IrrepBasisSet* GetCBasisSet  () const;
    virtual IrrepBasisSet* GetXBasisSet  () const;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const=0;
    // Polarized or un-polarized
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const=0;
protected:
    cl_t          itsCluster;
    BasisSet*     itsBasisSet;
    Hamiltonian*  itsHamiltonian;
    WaveFunction* itsWaveFunction;
    SCFIterator*  itsSCFIterator;
public:
    static PeriodicTable itsPT;
    double MaxRelErrE;
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
     
protected:
    virtual BasisSet*      GetBasisSet () const;

    int N,Lmax;
    double emin,emax;
    LAParams lap;
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
    
protected:
    virtual BasisSet*      GetBasisSet () const;
    int N,Lmax;
    double emin,emax;
    LAParams lap;
};

class SLm_OBasis : public SL_OBasis
{
 protected:
    virtual BasisSet*      GetBasisSet () const;   
};


class SGm_OBasis : public SG_OBasis
{
 protected:
    virtual BasisSet*      GetBasisSet () const;   
};

class PG_OBasis : public virtual QchemTester
{
public:
    PG_OBasis() : lap({qchem::Lapack,qchem::SVD,1e-4,1e-12}), N(0) {};     
    void Init(int _N, double _emin, double _emax, int _LMax)
    {
      N=_N;LMax=_LMax;emin=_emin;emax=_emax;
    }
private:
    virtual BasisSet*      GetBasisSet () const;
    LAParams lap;
    int N,LMax;
    double emin,emax;
};

// For Dirac wave functions
class SLmj_OBasis : public SL_OBasis
{
 protected:
    virtual BasisSet*      GetBasisSet () const;   
};


#include "Imp/WaveFunction/ElectronConfiguration.H"
//----------------------------------------------------------------------------------
//
//  Atoms and Molecules
//
class TestAtom : public virtual QchemTester
{
public:
    TestAtom(int _Z, int _q=0);
    virtual MeshParams GetMeshParams() const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return &ec;}
private:
    AtomElectronConfiguration ec;
};

class TestMolecule : public virtual QchemTester
{
public:
    TestMolecule() {};
    void Init(Molecule*);
    virtual MeshParams  GetMeshParams() const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return &ec;}
private:
    MoleculeElectronConfiguration ec;
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
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
};

class TestDirac : public virtual QchemTester
{
    public:
    virtual WaveFunction* GetWaveFunction(const BasisSet*) const;
};  


#endif //QchemTester_H
