#ifndef QchemTester1_H
#define QchemTester1_H

#include "gtest/gtest.h"
#include <memory>
#include <nlohmann/json.hpp>

import qchem.SCFIterator;
import qchem.LAParams;
import qchem.Hamiltonian;
import qchem.BasisSet1;
import Common.PeriodicTable;
import qchem.Cluster;
import qchem.Symmetry.AtomEC;
import qchem.Orbitals;
import qchem.Mesh; //To get Meshparams

class QchemTester1
{
public:
    typedef BasisSet1::irrepv_t irrepv_t;
    QchemTester1();
    virtual ~QchemTester1();
    void   Init(double eps, const nlohmann::json&, bool verbose=false, LAParams lap={qchem::Cholsky,1e-12});
    void   Iterate(const SCFParams&);
    double TotalEnergy() const;
    EnergyBreakdown GetEnergyBreakdown() const;
    double TotalCharge() const;
    const Orbitals* GetOrbitals(const Irrep_QNs& qns) const;
    const Orbital* GetOrbital(size_t index, const Irrep_QNs& qns) const;
    double RelativeError(double expected,bool quiet=false) const;
    double RelativeHFError(bool quiet=false) const;
    double RelativeDFTError(bool quiet=false) const;
    bool   Converged() const;
    int    GetLMax(int Z) const {return itsPT.GetMaxL(Z);}
    size_t GetIterationCount() const;
    irrepv_t GetIrreps(const Spin& ms) const;
protected:
    typedef std::shared_ptr<const Cluster> cl_t;

    // Atom of Molecule functions
    virtual const Cluster*  GetCluster   () const {return itsCluster.get();}
    virtual MeshParams      GetMeshParams() const=0;
    virtual int             GetZ         () const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const=0;

    // Orbital Basis Set functions SG, PG, Slater
    virtual BasisSet1* GetBasisSet   (const nlohmann::json&) const=0;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const=0;
protected:
    cl_t          itsCluster;
    BasisSet1*     itsBasisSet;
    SCFIterator*  itsSCFIterator;
public:
    static PeriodicTable itsPT;
    double MaxRelErrE;
};


//----------------------------------------------------------------------------------
//
//  Atoms and Molecules
//
class TestAtom : public virtual QchemTester1
{
public:
    TestAtom(int _Z, int _q=0);
    virtual MeshParams GetMeshParams() const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return &ec;}
private:
    virtual BasisSet1* GetBasisSet (const nlohmann::json&) const;
    Atom_EC ec;
    int itsZ;
};

// import qchem.Symmetry.MoleculeEC;
// class TestMolecule : public virtual QchemTester
// {
// public:
//     TestMolecule() {};
//     void Init(Cluster*);
//     virtual MeshParams  GetMeshParams() const;
//     virtual const ElectronConfiguration* GetElectronConfiguration() const {return &ec;}
// private:
//     virtual BasisSet1* GetBasisSet (const nlohmann::json&) const;
//     Molecule_EC ec;
// };



#endif //QchemTester1_H
