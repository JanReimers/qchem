
#include <memory>
#include <cmath>
#include <memory>

#include "QchemTester1.H"

import qchem.SCFIterator;
import qchem.SCFAccelerator.Factory;
import qchem.WaveFunction;
import qchem.ChargeDensity;
import qchem.Factory;
// import qchem.IrrepBasisSet;
import qchem.Cluster;
import qchem.Streamable;

PeriodicTable QchemTester1::itsPT;

QchemTester1::QchemTester1()
: itsCluster(0)
, itsBasisSet(0)
, itsSCFIterator(0)
, MaxRelErrE(0)
{
    //Cannot call virtual functions from here.
}

QchemTester1::~QchemTester1()
{
    delete itsBasisSet;
    delete itsSCFIterator;
}

void QchemTester1::Init(double eps,const nlohmann::json& js, bool verbose,LAParams lap)
{
    assert(eps>0.0);
    MaxRelErrE=eps;
    
    assert(itsCluster);
    assert(&*itsCluster);
    itsBasisSet=GetBasisSet(js); //SG, PG, Slater
    assert(itsBasisSet);
    // if (verbose)
    {
        std::cout << " " << *itsBasisSet << std::endl;
    }
    int Z=GetZ();
    nlohmann::json jsacc={{"NProj",4},{"EMax",Z*Z*0.1/32},{"EMin",1e-7},{"SVTol",1e-9}};
    SCFAccelerator* acc=SCFAcceleratorF::Factory(SCFAcceleratorF::Type::DIIS,jsacc);
    // SCFAccelerator* acc=new SCFAcceleratorDIIS({8,Z*Z*0.1/16,1e-7,1e-9});
    itsSCFIterator=new SCFIterator(itsBasisSet,GetElectronConfiguration(),GetHamiltonian(itsCluster),acc);
    assert(itsSCFIterator);
}

void QchemTester1::Iterate(const SCFParams& ipar)
{
    assert(itsSCFIterator);
    itsSCFIterator->Iterate(ipar);
}


double QchemTester1::TotalEnergy() const
{
    return GetEnergyBreakdown().GetTotalEnergy();
}

EnergyBreakdown QchemTester1::GetEnergyBreakdown() const
{
    return itsSCFIterator->GetEnergy();
}

double QchemTester1::TotalCharge() const
{
    return itsSCFIterator->GetWaveFunction()->GetChargeDensity()->GetTotalCharge();
}

const Orbitals* QchemTester1::GetOrbitals(const Irrep_QNs& qns) const
{
    return itsSCFIterator->GetWaveFunction()->GetOrbitals(qns);
}

const Orbital* QchemTester1::GetOrbital(size_t index, const Irrep_QNs& qns) const
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

size_t QchemTester1::GetIterationCount() const 
{
    return itsSCFIterator->GetIterationCount();
}

bool   QchemTester1::Converged() const
{
    return itsSCFIterator->Converged();
}

double QchemTester1::RelativeError(double E,bool quiet) const
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



double QchemTester1::RelativeHFError(bool quiet) const
{
    double E_HF=itsPT.GetEnergyHF(itsCluster->GetNuclearCharge());
    return RelativeError(E_HF,quiet);
}

double QchemTester1::RelativeDFTError(bool quiet) const
{
    double E_DFT=itsPT.GetEnergyDFT(itsCluster->GetNuclearCharge());
    return RelativeError(E_DFT,quiet);
}

int QchemTester1::GetZ() const
{
    return GetCluster()->GetNuclearCharge();
}

QchemTester1::irrepv_t QchemTester1::GetIrreps(const Spin& ms) const
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

BasisSet1* TestAtom::GetBasisSet (const nlohmann::json& js) const
{
    return BasisSetAtom1::Factory(js,itsZ);
}

// void TestMolecule::Init(Cluster* m)
// {
//     assert(m);
//     itsCluster=cl_t(m);
//     ec=Molecule_EC(m->GetNumElectrons());
// }

// MeshParams TestMolecule::GetMeshParams() const
// {
//     return MeshParams({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0,2});
// }
// BasisSet* TestMolecule::GetBasisSet (const nlohmann::json& js) const
// {
//     return BasisSetMolecule::Factory(js,GetCluster());
// }


