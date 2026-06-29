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
export import qchem.ChargeDensity;       // DM_CD (sampling rho(r) for the SAD atomic-density generator)
export import qchem.ChargeDensity.Seed;  // SeedStrategy (the SCF seed knob)
import qchem.SCFIterator;
import qchem.Structure;
import qchem.Orbitals;
import qchem.Mesh; //qcMesh::MeshParams
import qchem.Factory;
import qchem.ElectronConfiguration;
import qchem.PeriodicTable;

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
    // Choose the SCF seed strategy (call BEFORE Init; default resolves to CoreGuess for molecules).
    // SAD (superposition of atomic densities) is a DFT-only seed -- do not use it with HF.
    void   SetSeedStrategy(qchem::ChargeDensity::SeedStrategy s) {itsSeed=s;}

    double          TotalEnergy() const;
    EnergyBreakdown GetEnergyBreakdown() const;
    double          TotalCharge() const;
    //! Converged charge density (caller OWNS the returned heap density).  A ScalarFunction<double> you can
    //! sample rho(r) from -- used by the SAD atomic-density generator (scfrun --out) to dump a radial density.
    qchem::ChargeDensity::DM_CD* GetChargeDensity() const;
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
    typedef std::shared_ptr<const Structure> st_t;

    // Atom of Molecule functions
    virtual const Structure*  GetStructure   () const {return itsStructure.get();}
    virtual qcMesh::MeshParams      GetMeshParams() const=0;
    virtual int             GetZ         () const;
    virtual const ElectronConfiguration* GetElectronConfiguration() const {return itsEC;}

    // Orbital Basis Set functions SG, PG, Slater
    virtual Real_BS* GetBasisSet   (const nlohmann::json&) const=0;
    // Hamiltonian functions HF,semi HF, DFT all Pol or un-polarized.
    virtual Hamiltonian* GetHamiltonian(st_t& structure) const=0;
protected:
    st_t                   itsStructure;
    ElectronConfiguration* itsEC;
    Real_BS*               itsBasisSet;
    Hamiltonian*           itsHamiltonian;
    SCFIterator*           itsSCFIterator;
    nlohmann::json         itsAccConfig; //SCF accelerator config (empty => DIIS defaults)
    qchem::ChargeDensity::SeedStrategy itsSeed = qchem::ChargeDensity::SeedStrategy::Default;
public:
    static PeriodicTableSaito itsPT;
};


//----------------------------------------------------------------------------------
//
//  Atoms and Molecules
//
export class TestAtom : public QchemTester
{
public:
    TestAtom(int _Z, int _q=0);
    virtual qcMesh::MeshParams GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
    int itsZ;
};

export class TestDiracAtom : public QchemTester
{
public:
    TestDiracAtom(int _Z, int _q=0);
    virtual qcMesh::MeshParams GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
    int itsq;
};

export class TestMolecule : public QchemTester
{
public:
    TestMolecule(Structure*);
    virtual qcMesh::MeshParams  GetMeshParams() const;
private:
    virtual Real_BS* GetBasisSet (const nlohmann::json&) const;
};



