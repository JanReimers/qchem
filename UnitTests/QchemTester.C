
#include <memory>
#include <cmath>
#include <memory>

#include <BasisSet/Factory.H>
#include "QchemTester.H"
#include <SCFAccelerator/SCFAccelerator.H>
#include <SCFAccelerator/Factory.H>
#include <WaveFunction/WaveFunction.H>
#include <Hamiltonian/Hamiltonian.H>
#include <BasisSet/BasisSet.H>
#include <BasisSet/Irrep_BS.H>
#include <Hamiltonian/TotalEnergy.H>
#include <Orbitals/Orbitals.H>
#include <ChargeDensity/ChargeDensity.H>
#include <SCFIterator.H>
#include "Cluster/Molecule.H"

PeriodicTable QchemTester::itsPT;

import qchem.Cluster;
import qchem.Atom;


QchemTester::QchemTester()
: itsCluster(0)
, itsBasisSet(0)
, itsSCFIterator(0)
, MaxRelErrE(0)
{
    //Cannot call virtual functions from here.
}

QchemTester::~QchemTester()
{
    delete itsBasisSet;
    delete itsSCFIterator;
}

void QchemTester::Init(double eps,const nlohmann::json& js, bool verbose,LAParams lap)
{
    assert(eps>0.0);
    MaxRelErrE=eps;
    
    assert(itsCluster);
    assert(&*itsCluster);
    itsBasisSet=GetBasisSet(js); //SG, PG, Slater
    assert(itsBasisSet);
    itsBasisSet->Set(lap);
    if (verbose)
    {
        std::cout << " " << *itsBasisSet << std::endl;
    }
    int Z=GetZ();
    nlohmann::json jsacc={{"NProj",8},{"EMax",Z*Z*0.1/16},{"EMin",1e-7},{"SVTol",1e-9}};
    SCFAccelerator* acc=SCFAcceleratorF::Factory(SCFAcceleratorF::Type::DIIS,jsacc);
    // SCFAccelerator* acc=new SCFAccelerator_DIIS({8,Z*Z*0.1/16,1e-7,1e-9});
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
        else
            std::cout << error*1e9 << "(ppb)" << std::endl;
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

QchemTester::symv_t QchemTester::GetSymmetries() const
{
    return itsBasisSet->GetSymmetries();
}



TestAtom::TestAtom(int Z, int q) : ec(Z-q) //Pass in # of electrons.
{
    itsZ=Z-q;
    Cluster* cl=new Molecule;
    cl->Insert(new Atom(Z,q,Vector3D<double>(0,0,0)));
    itsCluster=cl_t(cl);
};

MeshParams TestAtom::GetMeshParams() const
{
    return MeshParams({qchem::MHL,50,3,2.0,qchem::Gauss,1,0,0,2});
}

BasisSet* TestAtom::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSetAtom::Factory(js,itsZ);
}

void TestMolecule::Init(Molecule* m)
{
    assert(m);
    itsCluster=cl_t(m);
    ec=Molecule_EC(m->GetNumElectrons());
}

MeshParams TestMolecule::GetMeshParams() const
{
    return MeshParams({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0,2});
}
BasisSet* TestMolecule::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSetMolecule::Factory(js,GetCluster());
}


