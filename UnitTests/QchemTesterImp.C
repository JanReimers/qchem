// File: UnitTests/QchemTester.C  Heper class for doing SCF calculations in unit tests.
module;
#include <memory>
#include <nlohmann/json.hpp>
#include <iostream>
module qchem.Unittests.QchemTester;
import qchem.ElectronConfiguration.AtomNR;
import qchem.ElectronConfiguration.AtomDirac;
import qchem.ElectronConfiguration.Molecule;
import qchem.SCFAccelerator.Factory;
import qchem.Math;


using qchem::SCFAccelerators::SCFAccelerator;

PeriodicTableSaito QchemTester::itsPT;

QchemTester::QchemTester(ElectronConfiguration* ec)
: itsStructure(0)
, itsEC(ec)
, itsBasisSet(0)
, itsSCFIterator(0)
{
    //Cannot call virtual functions from here.
}

QchemTester::~QchemTester()
{
    delete itsEC;
    delete itsBasisSet;
    delete itsSCFIterator;   // owns and deletes itsHamiltonian (-> its terms -> the DFT fit bases)
}

void QchemTester::Init(const nlohmann::json& js, bool verbose)
{
    Init(GetBasisSet(js), verbose);
}

void QchemTester::Init(BasisSetAccuracy acc, BasisSet::Atom::Type type, bool verbose)
{
    Init(Factory(acc, type, GetZ()), verbose);
}


void QchemTester::Init(Real_BS* bs, bool verbose)
{
    itsBasisSet=bs;
    
    assert(itsStructure);
    assert(&*itsStructure);
    assert(itsBasisSet);
    // if (verbose)
    {
        std::cout << " " << *itsBasisSet << std::endl;
    }
    int Z=GetZ();
    itsHamiltonian=GetHamiltonian(itsStructure);
    nlohmann::json jsacc={{"NProj",4},{"EMax",Z*Z*0.1/32},{"EMin",1e-7},{"SVTol",5e-9}};
    using qchem::SCFAccelerators::Type;
    std::string ts="DIIS";
    if (itsAccConfig.is_object())   //caller supplied an accelerator config (else default DIIS)
    {
        for (auto& [k,v]:itsAccConfig.items()) jsacc[k]=v;  //overrides incl. heuristic params
        ts=itsAccConfig.value("type","DIIS");
    }
    bool directmin = (ts=="directmin" || ts=="DirectMin");
    if (directmin) jsacc["EMax"]=1e10; //GDM always steps; the direct-min loop seeds via diagonalize
    Type type = directmin ? Type::GDM
              : ts=="Ladder" ? Type::Ladder : ts=="GDM" ? Type::GDM : Type::DIIS;
    SCFAccelerator* acc=qchem::SCFAccelerators::Factory(type,jsacc);
    delete itsSCFIterator;
    itsSCFIterator=new SCFIterator(itsBasisSet,GetElectronConfiguration(),itsHamiltonian,acc);
    if (directmin) itsSCFIterator->SetDirectMin(true);
    assert(itsSCFIterator);
}

void QchemTester::Iterate(const SCFParams& ipar)
{
    assert(itsSCFIterator);
    itsSCFIterator->Iterate(ipar);
}


double QchemTester::TotalEnergy() const
{
    return GetEnergyBreakdown().GetTotalEnergy();
}

EnergyBreakdown QchemTester::GetEnergyBreakdown() const
{
    return itsSCFIterator->GetEnergy();
}

double QchemTester::TotalCharge() const
{
    return itsSCFIterator->GetWaveFunction()->GetChargeDensity()->GetTotalCharge();
}

const Orbitals* QchemTester::GetOrbitals(const Irrep& qns) const
{
    return itsSCFIterator->GetWaveFunction()->GetOrbitals(qns);
}

const Orbital* QchemTester::GetOrbital(size_t index, const Irrep& qns) const
{
    const Orbitals* orbs=GetOrbitals(qns);
    assert(index<(size_t)orbs->GetNumOrbitals());
    const Orbital* o=0;
    for (auto oi:orbs->Iterate<Orbital>())
    {
        if (index==0) 
        {
            o=oi;
            break;
        }
        index--;
    }
    return o;
}

size_t QchemTester::GetIterationCount() const 
{
    return itsSCFIterator->GetIterationCount();
}

bool   QchemTester::Converged() const
{
    return itsSCFIterator->Converged();
}

double QchemTester::RelativeError(double E,bool quiet) const
{
    double error=(E-TotalEnergy())/E;
    if (!quiet)
    {
        std::cout.precision(9);
        std::cout << "E relative error=" << error*100.0 << "%, ";
        std::cout.precision(2);
        if (fabs(error)>1e-7)
            std::cout << error*1e6 << "(ppm)" << std::endl;
        else if (fabs(error)>1e-10)
            std::cout << error*1e9 << "(ppb)" << std::endl;
        else
            std::cout << error*1e12 << "(ppt)" << std::endl;
    }
    return error;
}



double QchemTester::RelativeHFError(bool quiet) const
{
    double E_HF=itsPT.GetEnergyHF(itsStructure->GetNuclearCharge());
    return RelativeError(E_HF,quiet);
}

double QchemTester::RelativeDFTError(bool quiet) const
{
    double E_DFT=itsPT.GetEnergyDFT(itsStructure->GetNuclearCharge());
    return RelativeError(E_DFT,quiet);
}

double QchemTester::RelativeDHFError(bool quiet) const
{
    double E_DHF=itsPT.GetEnergyDHF(itsStructure->GetNuclearCharge());
    return RelativeError(E_DHF,quiet);
}

int QchemTester::GetZ() const
{
    return GetStructure()->GetNuclearCharge();
}

irrepv_t QchemTester::GetIrreps(const Spin& ms) const
{
    return itsBasisSet->GetIrreps(ms);
}



TestAtom::TestAtom(int Z, int q) : QchemTester(new Atom_EC(Z-q)) //Pass in the electr4on config.
{
    itsZ=Z-q;
    itsStructure=st_t(new Atom(Z,q,Vector3D<double>(0,0,0)));
};

qcMesh::MeshParams TestAtom::GetMeshParams() const
{
    return {.radial=qcMesh::RadialKind::MHL, .nRadial=50, .mhl_m=3, .mhl_alpha=2.0,
            .angular=qcMesh::AngularKind::Gauss, .nAngular=1, .beckeOrder=2};
}

Real_BS* TestAtom::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSet::Atom::Factory(js,itsZ);
}

TestDiracAtom::TestDiracAtom(int Z, int q) : QchemTester(new AtomDirac_EC(Z-q)), itsq(q)
{
    itsStructure=st_t(new Atom(Z,q,Vector3D<double>(0,0,0)));
};

qcMesh::MeshParams TestDiracAtom::GetMeshParams() const
{
    return {.radial=qcMesh::RadialKind::MHL, .nRadial=50, .mhl_m=3, .mhl_alpha=2.0,
            .angular=qcMesh::AngularKind::Gauss, .nAngular=1, .beckeOrder=2};
}

Real_BS* TestDiracAtom::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSet::Atom::Factory(js,GetZ()-itsq);
}

TestMolecule::TestMolecule(Structure* m) : QchemTester(new Molecule_EC(m->GetNumElectrons())) 
{
    itsStructure=st_t(m);
};

qcMesh::MeshParams TestMolecule::GetMeshParams() const
{
    return {.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=3, .mhl_alpha=2.0,
            .angular=qcMesh::AngularKind::Gauss, .nAngular=12, .beckeOrder=2};
}
Real_BS* TestMolecule::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSet::Molecule::Factory(js,GetStructure());
}


