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
     
private:
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
    
private:
    virtual BasisSet*      GetBasisSet () const;
    int N,Lmax;
    double emin,emax;
    LAParams lap;
};

class PG_OBasis : public virtual QchemTester
{
public:
    PG_OBasis() : lap({qchem::Lapack,qchem::SVD,1e-4,1e-12}) {};     
private:
    virtual BasisSet*      GetBasisSet () const;
    LAParams lap;
};

//----------------------------------------------------------------------------------
//
//  Atoms and Molecules
//
class TestAtom : public virtual QchemTester
{
public:
    TestAtom(int _Z, int _q=0);
    virtual MeshParams GetMeshParams() const;
};

class TestMolecule : public virtual QchemTester
{
public:
    TestMolecule() {};
    void Init(Molecule*);
    virtual MeshParams  GetMeshParams() const;
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
