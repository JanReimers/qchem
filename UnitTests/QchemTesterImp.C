// File: UnitTests/QchemTester.C  Heper class for doing SCF calculations in unit tests.
module;
#include <memory>
#include <cmath>
#include <nlohmann/json.hpp>
#include <iostream>
module qchem.Unittests.QchemTester;
import qchem.SCFAccelerator.Factory;


using qchem::SCFAccelerators::SCFAccelerator;

PeriodicTable QchemTester::itsPT;

QchemTester::QchemTester()
: itsCluster(0)
, itsBasisSet(0)
, itsSCFIterator(0)
{
    //Cannot call virtual functions from here.
}

QchemTester::~QchemTester()
{
    delete itsBasisSet;
    delete itsSCFIterator;
}

void QchemTester::Init(const nlohmann::json& js, bool verbose,LAParams lap)
{
    Init(GetBasisSet(js),verbose,lap);
}

void QchemTester::Init(BasisSetAccuracy acc, BasisSet::Atom::Type type,bool verbose,LAParams lap)
{
    Init(PoolFactory(acc,type,GetZ()),verbose,lap); 
}
    

void QchemTester::Init(Real_BS* bs, bool verbose,LAParams lap)
{
    itsBasisSet=bs;
    assert(eps>0.0);
    
    assert(itsCluster);
    assert(&*itsCluster);
    assert(itsBasisSet);
    // if (verbose)
    {
        std::cout << " " << *itsBasisSet << std::endl;
    }
    int Z=GetZ();
    nlohmann::json jsacc={{"NProj",4},{"EMax",Z*Z*0.1/32},{"EMin",1e-7},{"SVTol",5e-9}};
    SCFAccelerator* acc=qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS,jsacc);
    // SCFAccelerator* acc=new SCFAcceleratorDIIS({8,Z*Z*0.1/16,1e-7,1e-9});
    delete itsSCFIterator;
    itsSCFIterator=new SCFIterator(itsBasisSet,GetElectronConfiguration(),GetHamiltonian(itsCluster),acc);
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

const Orbitals* QchemTester::GetOrbitals(const Irrep_QNs& qns) const
{
    return itsSCFIterator->GetWaveFunction()->GetOrbitals(qns);
}

const Orbital* QchemTester::GetOrbital(size_t index, const Irrep_QNs& qns) const
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
    double E_HF=itsPT.GetEnergyHF(itsCluster->GetNuclearCharge());
    return RelativeError(E_HF,quiet);
}

double QchemTester::RelativeDFTError(bool quiet) const
{
    double E_DFT=itsPT.GetEnergyDFT(itsCluster->GetNuclearCharge());
    return RelativeError(E_DFT,quiet);
}

int QchemTester::GetZ() const
{
    return GetCluster()->GetNuclearCharge();
}

irrepv_t QchemTester::GetIrreps(const Spin& ms) const
{
    return itsBasisSet->GetIrreps(ms);
}



TestAtom::TestAtom(int Z, int q) : ec(Z-q) //Pass in # of electrons.
{
    itsZ=Z-q;
    itsCluster=cl_t(new Atom(Z,q,Vector3D<double>(0,0,0)));
};

MeshParams TestAtom::GetMeshParams() const
{
    return MeshParams({qchem::MHL,50,3,2.0,qchem::Gauss,1,0,0,2});
}

Real_BS* TestAtom::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSet::Atom::Factory(js,itsZ);
}

void TestMolecule::Init(Cluster* m)
{
    assert(m);
    itsCluster=cl_t(m);
    ec=Molecule_EC(m->GetNumElectrons());
}

MeshParams TestMolecule::GetMeshParams() const
{
    return MeshParams({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0,2});
}
Real_BS* TestMolecule::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSet::Molecule::Factory(js,GetCluster());
}


