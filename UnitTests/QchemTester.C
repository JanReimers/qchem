// File: UnitTests/QchemTester.C  Heper class for doing SCF calculations in unit tests.
module;
#include "gtest/gtest.h"
#include <memory>
#include <nlohmann/json.hpp>

export module qchem.Unittests.QchemTester;
export import qchem.BasisSet.Atom.Factory;
export import qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
export import qchem.SCFParams;
import qchem.SCFIterator;
import qchem.Structure;
import qchem.Orbitals;
import qchem.Mesh; //To get Meshparams
import qchem.Factory;
import qchem.ElectronConfiguration;
import Common.PeriodicTable;

export typedef BasisSet::irrepv_t irrepv_t;

export using qchem::Orbitals::Orbital;
export using qchem::Orbitals::Orbitals;
using qchem::Hamiltonian::Hamiltonian;
using qchem::SCFIterator::SCFIterator;
using namespace BasisSet::Atom;

export class QchemTester
{
public:
    QchemTester(ElectronConfiguration* ec);
    virtual ~QchemTester();
    void   Init(const nlohmann::json&, bool verbose=false);
    void   Init(Real_BS*, bool verbose=false);
    void   Init(BasisSetAccuracy acc, BasisSet::Atom::Type type, bool verbose=false);
    void   Iterate(const SCFParams&);
    // Choose the SCF accelerator via JSON, e.g. {"type":"Ladder","floor":1e-4,"stall":5}.
    // "type" is "DIIS" (default), "GDM", or "Ladder"; other keys override the defaults.
    void   SetAcceleratorConfig(const nlohmann::json& j) {itsAccConfig=j;}

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
    double          RelativeDHFError(bool quiet=false) const;
    bool            Converged() const;
    
    int             GetLMax(int Z) const {return itsPT.GetMaxL(Z);}
    size_t          GetIterationCount() const;
    irrepv_t        GetIrreps(const Spin& ms) const;
protected:
    typedef std::shared_ptr<const Structure> cl_t;

    // Atom of Molecule functions
    virtual const Structure*  GetStructure   () const {return itsStructure.get();}
    virtual MeshParams      GetMeshParams() const=0;
    virtual int             GetZ         () const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return itsEC;}

    // Orbital Basis Set functions SG, PG, Slater
    virtual Real_BS* GetBasisSet   (const nlohmann::json&) const=0;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(cl_t& structure) const=0;
protected:
    cl_t                   itsStructure;
    ElectronConfiguration* itsEC;
    Real_BS*               itsBasisSet;
    Hamiltonian*           itsHamiltonian;
    SCFIterator*           itsSCFIterator;
    nlohmann::json         itsAccConfig; //SCF accelerator config (empty => DIIS defaults)
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

export class TestDiracAtom : public QchemTester
{
public:
    TestDiracAtom(int _Z, int _q=0);
    virtual MeshParams GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
    int itsq;
};

export class TestMolecule : public QchemTester
{
public:
    TestMolecule(Structure*);
    virtual MeshParams  GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
};



