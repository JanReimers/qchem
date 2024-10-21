
#include "QchemTester.H"
#include <MeshParams.H>
#include <SCFIterator.H>
#include <WaveFunction.H>
#include <Hamiltonian.H>
#include <Cluster.H>
#include <BasisSet.H>
#include <TotalEnergy.H>
#include <memory>

PeriodicTable QchemTester::itsPT;

QchemTester::QchemTester()
: itsCluster(0)
, itsBasisSet(0)
, itsHamiltonian(0)
, itsWaveFunction(0)
, itsSCFIterator(0)
, MaxRelErrE(0)
{
    //Cannot call virtual functions from here.
}

QchemTester::~QchemTester()
{
    delete itsBasisSet;
    delete itsWaveFunction;
    delete itsSCFIterator;
}

void QchemTester::Init(double eps)
{
    assert(eps>0.0);
    MaxRelErrE=eps;
    
    assert(itsCluster);
    assert(&*itsCluster);
    itsBasisSet=GetBasisSet(); //SG, PG, Slater
    assert(itsBasisSet);
    itsHamiltonian=GetHamiltonian(itsCluster); //HF,semi HF, DFT all Pol or un-polarized.
    assert(itsHamiltonian);
    itsWaveFunction=GetWaveFunction(itsBasisSet); //Polarized or un-polarized
    assert(itsWaveFunction);
    ChargeDensity* GuessCD=itsWaveFunction->GetChargeDensity();
    itsSCFIterator=itsWaveFunction->
                   MakeIterator(itsHamiltonian,GuessCD,itsCluster->GetNumElectrons(),0.0,false); //show plot flag.

}

void QchemTester::Iterate(const SCFIterationParams& ipar)
{
    assert(itsSCFIterator);
    itsSCFIterator->Iterate(ipar);
}


double QchemTester::TotalEnergy() const
{
    assert(itsHamiltonian);
    return itsHamiltonian->GetTotalEnergy().GetTotalEnergy();
}

#include <cmath> //fabs
double QchemTester::RelativeHFError(bool quiet) const
{
    double E_HF=itsPT.GetEnergyHF(itsCluster->GetNuclearCharge());
    double error=fabs((E_HF-TotalEnergy())/E_HF);
    if (!quiet)
    {
        std::cout.precision(6);
        std::cout << "E_HF relative error=" << error*100.0 << "%, ";
        std::cout.precision(2);
        std::cout << error*1e6 << "(ppm)" << std::endl;            
    }
    return error;
}

double QchemTester::RelativeDFTError(bool quiet) const
{
    double E_DFT=itsPT.GetEnergyDFT(itsCluster->GetNuclearCharge());
    double error=fabs((E_DFT-TotalEnergy())/E_DFT);
    if (!quiet)
    {
        std::cout.precision(6);
        std::cout << "E_DFT relative error=" << error*100.0 << "%, ";
        std::cout.precision(2);
        std::cout << error*1e6 << "(ppm)" << std::endl;            
    }
    return error;
}

int QchemTester::GetZ() const
{
    return GetCluster()->GetNuclearCharge();
}

IrrepBasisSet* QchemTester::GetCBasisSet() const
{
    return itsBasisSet->CreateCDFitBasisSet(itsCluster.get());
}

IrrepBasisSet* QchemTester::GetXBasisSet() const
{
    return itsBasisSet->CreateVxcFitBasisSet(itsCluster.get());    
}
    

#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
BasisSet* SG_OBasis::GetBasisSet () const
{
    return  new SphericalGaussian::BasisSet(lap,N,emin,emax,Lmax);
}

#include "Imp/BasisSet/Slater/BasisSet.H"
BasisSet* SL_OBasis::GetBasisSet () const
{
    return new Slater::BasisSet(lap,N,emin,emax,Lmax);
}

#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Gaussian94.H"

BasisSet* PG_OBasis::GetBasisSet () const
{
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
    return new PolarizedGaussian::BasisSet(lap, &reader,GetCluster());
}


#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
TestAtom::TestAtom(int Z, int q)  
{
    Cluster* cl=new Molecule;
    cl->Insert(new Atom(Z,q,Vector3D<double>(0,0,0)));
    itsCluster=cl_t(cl);
};

MeshParams TestAtom::GetMeshParams() const
{
    return MeshParams({qchem::MHL,50,3,2.0,qchem::Gauss,1,0,0});
}

void TestMolecule::Init(Molecule* p)
{
    assert(p);
    itsCluster=cl_t(p);
}

MeshParams TestMolecule::GetMeshParams() const
{
    return MeshParams({qchem::MHL,30,3,2.0,qchem::Gauss,12,0,0});
}

    

#include "Imp/WaveFunction/MasterUnPolarizedWF.H"
WaveFunction* TestUnPolarized::GetWaveFunction(const BasisSet* bs) const
{
    return new MasterUnPolarizedWF(bs);
}

#include "Imp/WaveFunction/MasterPolarizedWF.H"

TestPolarized::TestPolarized(int Z) : spin(QchemTester::itsPT.GetNumUnpairedElectrons(Z)) {};

WaveFunction* TestPolarized::GetWaveFunction(const BasisSet* bs) const
{
    assert(spin>=0);
    return new MasterPolarizedWF(bs,spin);
}

