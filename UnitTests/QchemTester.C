// File: UnitTests/QchemTester.C  Heper class for doing SCF calculations in unit tests.
module;
#include "gtest/gtest.h"
#include <memory>
#include <nlohmann/json.hpp>

export module qchem.Unittests.QchemTester;
export import qchem.Unittests.BasisSetPool;
export import qchem.Hamiltonian;
export import qchem.Hamiltonian.Factory;
export import qchem.SCFParams;
import qchem.SCFIterator;
import qchem.LAParams;
import qchem.Cluster;
import qchem.Orbitals;
import qchem.Mesh; //To get Meshparams
import qchem.Factory;
import qchem.Symmetry.AtomEC;
import qchem.Symmetry.MoleculeEC;
import Common.PeriodicTable;

export typedef BasisSet::irrepv_t irrepv_t;

export using qchem::Orbitals::Orbital;
export using qchem::Orbitals::Orbitals;
using qchem::Hamiltonian::Hamiltonian;
using qchem::SCFIterator::SCFIterator;

export class QchemTester
{
public:
    QchemTester(ElectronConfiguration* ec);
    virtual ~QchemTester();
    void   Init(const nlohmann::json&, bool verbose=false, LAParams lap={qchem::Cholsky,1e-12});
    void   Init(Real_BS*, bool verbose=false, LAParams lap={qchem::Cholsky,1e-12});
    void   Init(BasisSetAccuracy acc, BasisSet::Atom::Type type,bool verbose=false,LAParams lap={qchem::Cholsky,1e-12});
    void   Iterate(const SCFParams&);
    
    double          TotalEnergy() const;
    EnergyBreakdown GetEnergyBreakdown() const;
    double          TotalCharge() const;
    const Real_BS*  GetBasisSet() const {return itsBasisSet;}
    Hamiltonian* GetHamiltonian() const {return itsHamiltonian;}
    const Orbitals* GetOrbitals(const Irrep& qns) const;
    const Orbital*  GetOrbital(size_t index, const Irrep& qns) const;
    double          RelativeError(double expected,bool quiet=false) const;
    double          RelativeHFError(bool quiet=false) const;
    double          RelativeDFTError(bool quiet=false) const;
    bool            Converged() const;
    
    int             GetLMax(int Z) const {return itsPT.GetMaxL(Z);}
    size_t          GetIterationCount() const;
    irrepv_t        GetIrreps(const Spin& ms) const;
protected:
    typedef std::shared_ptr<const Cluster> cl_t;

    // Atom of Molecule functions
    virtual const Cluster*  GetCluster   () const {return itsCluster.get();}
    virtual MeshParams      GetMeshParams() const=0;
    virtual int             GetZ         () const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return itsEC;}

    // Orbital Basis Set functions SG, PG, Slater
    virtual Real_BS* GetBasisSet   (const nlohmann::json&) const=0;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(cl_t& cluster) const=0;
protected:
    cl_t                   itsCluster;
    ElectronConfiguration* itsEC;
    Real_BS*               itsBasisSet;
    Hamiltonian*           itsHamiltonian;
    SCFIterator*           itsSCFIterator;
public:
    static PeriodicTableSaito itsPT;
    static PeriodicTable  itsPTold;
};


//----------------------------------------------------------------------------------
//
//  Atoms and Molecules
//
export class TestAtom : public QchemTester
{
public:
    TestAtom(int _Z, int _q=0);
    virtual MeshParams GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
    int itsZ;
};

export class TestMolecule : public QchemTester
{
public:
    TestMolecule(Cluster*);
    virtual MeshParams  GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
};



