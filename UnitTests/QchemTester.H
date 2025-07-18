#ifndef QchemTester_H
#define QchemTester_H

#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
#include <memory>

#include <Hamiltonian/fwd.H>
#include <Orbitals/fwd.H>
#include <SCFIteratorfwd.H>

#include <BasisSet/BasisSet.H>
#include <SCFParams.H>
#include <LASolver/LAParams.H>
#include <Hamiltonian/TotalEnergy.H>

import Common.PeriodicTable;
import Mesh; //To get Meshparams
import qchem.Molecule;
import qchem.Symmetry.Irrep;
import qchem.Symmetry.AtomEC;

class QchemTester
{
public:
    typedef BasisSet::symv_t symv_t;
    QchemTester();
    virtual ~QchemTester();
    void   Init(double eps, const nlohmann::json&, bool verbose=false, LAParams lap={qchem::Lapack,qchem::Cholsky,1e-9,1e-12});
    void   Iterate(const SCFParams&);
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
    virtual BasisSet* GetBasisSet   (const nlohmann::json&) const=0;
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
    virtual BasisSet* GetBasisSet (const nlohmann::json&) const;
    Atom_EC ec;
    int itsZ;
};

import qchem.Symmetry.MoleculeEC;
class TestMolecule : public virtual QchemTester
{
public:
    TestMolecule() {};
    void Init(Molecule*);
    virtual MeshParams  GetMeshParams() const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return &ec;}
private:
    virtual BasisSet* GetBasisSet (const nlohmann::json&) const;
    Molecule_EC ec;
};



#endif //QchemTester_H
