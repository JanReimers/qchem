#ifndef QchemTester_H
#define QchemTester_H

#include "gtest/gtest.h"
#include <types.H>
#include <BasisSet.H>
#include <IterationParams.H>
#include <LAParams.H>
#include "Imp/Misc/PeriodicTable.H"
#include <TotalEnergy.H>

#include <memory>

class QchemTester
{
public:
    typedef BasisSet::symv_t symv_t;
    QchemTester();
    virtual ~QchemTester();
    void   Init(double eps);
    void   Iterate(const SCFIterationParams&);
    double TotalEnergy() const;
    EnergyBreakdown GetEnergyBreakdown() const;
    double TotalCharge() const;
    const Orbitals* GetOrbitals(const Irrep_QNs& qns) const;
    const Orbital* GetOrbital(size_t index, const Irrep_QNs& qns) const;
    double RelativeError(double expected,bool quiet=false) const;
    double RelativeHFError(bool quiet=false) const;
    double RelativeDFTError(bool quiet=false) const;
    int    GetLMax(int Z) const {return itsPT.GetMaxL(Z);}
    symv_t GetSymmetries() const;
    size_t GetIterationCount() const;
    
protected:
    typedef std::shared_ptr<const Cluster> cl_t;

    // Atom of Molecule functions
    virtual const Cluster*  GetCluster   () const {return itsCluster.get();}
    virtual MeshParams      GetMeshParams() const=0;
    virtual int             GetZ         () const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const=0;

    // Orbital Basis Set functions SG, PG, Slater
    virtual BasisSet* GetBasisSet   () const=0;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const=0;
protected:
    cl_t          itsCluster;
    BasisSet*     itsBasisSet;
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
        Init(_N,_emin,_emax,_Lmax,LAParams({qchem::Lapack,qchem::SVD,1e-9,1e-12}));
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
        Init(_N,_emin,_emax,_Lmax,LAParams({qchem::Lapack,qchem::SVD,1e-9,1e-12}));
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

class SG_RKB_OBasis : public SG_OBasis
{
 protected:
    virtual BasisSet*      GetBasisSet () const;   
};

class BS_OBasis : public virtual QchemTester
{
public:
    void Init(int _N, double _rmin, double _rmax, int _LMax, const LAParams& _lap)
    {
        N=_N;LMax=_LMax;rmin=_rmin;rmax=_rmax;
        lap=_lap;
    }
    void Init(int _N, double _rmin, double _rmax, int _LMax) 
    {
        Init(_N,_rmin,_rmax,_LMax,LAParams({qchem::Lapack,qchem::Cholsky,1e-9,1e-12}));
    }
 protected:
    virtual BasisSet*      GetBasisSet () const;   
    LAParams lap;
    int N,LMax;
    double rmin,rmax;
};

class BSm_OBasis : public virtual BS_OBasis
{
 protected:
    virtual BasisSet*      GetBasisSet () const;  

};



#include "Imp/BasisSet/Atom/EC.H"
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
    Atom_EC ec;
};

#include "Imp/BasisSet/Molecule/EC.H"

class TestMolecule : public virtual QchemTester
{
public:
    TestMolecule() {};
    void Init(Molecule*);
    virtual MeshParams  GetMeshParams() const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return &ec;}
private:
    Molecule_EC ec;
};



#endif //QchemTester_H
