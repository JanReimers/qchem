
#include "QchemTester.H"
#include "SCFIterator/SCFIterator.H"
#include <WaveFunction.H>
#include <Hamiltonian.H>
#include <Cluster.H>
#include <BasisSet.H>
#include <TotalEnergy.H>

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
    
    itsCluster=GetCluster(); //Atom or Molecule
    assert(&*itsCluster);
    itsBasisSet=GetBasisSet(&*itsCluster); //SG, PG, Slater
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
    



#include "HamiltonianImplementation/HamiltonianImplementation.H"
#include "HamiltonianImplementation/Kinetic.H"
#include "HamiltonianImplementation/HartreeFockVxc.H"
#include "HamiltonianImplementation/ExactVnn.H"
#include "HamiltonianImplementation/ExactVen.H"
Hamiltonian* TestHamiltonian::GetHamiltonian(const rc_ptr<Cluster>& cl) const
{
    Hamiltonian* H=new HamiltonianImplementation();
    
    H->Add(new Kinetic);
    H->Add(new ExactVnn(cl));
    H->Add(new ExactVen(cl));
    H->Add(GetVee());
    H->Add(GetVxc()); //pol or un pol
    return H;
}

#include "HamiltonianImplementation/ExactVee.H"
#include "HamiltonianImplementation/HartreeFockVxc.H"

HamiltonianTerm* HFHamiltonian:: GetVee() const
{
    return new ExactVee;
}

HamiltonianTerm* HFHamiltonian:: GetVxc() const
{
    return new HartreeFockVxc;
}

#include "HamiltonianImplementation/SlaterExchange.H"
#include "HamiltonianImplementation/FittedVxc.H"
#include "HamiltonianImplementation/PolarizedFittedVxc.H"
#include "Hamiltonian/ExchangeFunctional.H" 

SHFHamiltonian::SHFHamiltonian(int Z)
{
    XcFunct=new SlaterExchange(QchemTester::itsPT.GetSlaterAlpha(Z));
}
SHFHamiltonian::SHFHamiltonian(double ex)
{
    XcFunct=new SlaterExchange(ex);
}

SHFHamiltonian::~SHFHamiltonian()
{}


HamiltonianTerm* SHFHamiltonian:: GetVxc() const
{
    rc_ptr<IrrepBasisSet> XFitBasis=GetXBasisSet();    
    return new FittedVxc(XFitBasis, XcFunct,GetIntegrationMesh());
}

PolSHFHamiltonian::PolSHFHamiltonian(int Z)
{
    XcFunct=new SlaterExchange(QchemTester::itsPT.GetSlaterAlpha(Z),Spin(Spin::Up));
}
PolSHFHamiltonian::PolSHFHamiltonian(double ex)
{
    XcFunct=new SlaterExchange(ex);
}

HamiltonianTerm* PolSHFHamiltonian:: GetVxc() const
{
    rc_ptr<IrrepBasisSet> XFitBasis=GetXBasisSet();    
    return new PolarizedFittedVxc(XFitBasis,XcFunct ,GetIntegrationMesh());
}


#include "HamiltonianImplementation/PolarizedHartreeFockVxc.H"
HamiltonianTerm* PolHFHamiltonian:: GetVxc() const
{
    return new PolarizedHartreeFockVxc;
}

#include "HamiltonianImplementation/CDFittedVee.H"

HamiltonianTerm* DFTHamiltonian::GetVee() const
{
    rc_ptr<IrrepBasisSet> CFitBasis=GetCBasisSet(); 
    return new CDFittedVee(CFitBasis,GetIntegrationMesh(),GetZ());
}
HamiltonianTerm* PolDFTHamiltonian::GetVee() const
{
    rc_ptr<IrrepBasisSet> CFitBasis=GetCBasisSet(); 
    return new CDFittedVee(CFitBasis,GetIntegrationMesh(),GetZ());
}

#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
BasisSet* SG_OBasis::GetBasisSet (const Cluster*) const
{
    BasisSet* bs=new SphericalGaussian::BasisSet(lap,N,emin,emax,Lmax);
    idb=bs->GetDataBase();
    return  bs;
}

IrrepBasisSet* SG_OBasis::GetCBasisSet () const
{
    return new SphericalGaussian::IrrepBasisSet(lap,idb,N,emin*2.0,emax*2.0,0);
}

IrrepBasisSet* SG_OBasis::GetXBasisSet () const
{
    return new SphericalGaussian::IrrepBasisSet(lap,idb,N,emin*2.0/3.0,emax*2.0/3.0,0);
}

#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
BasisSet* SL_OBasis::GetBasisSet (const Cluster*) const
{
    BasisSet* bs=new Slater::BasisSet(lap,N,emin,emax,Lmax);
    idb=bs->GetDataBase();
    return  bs;
}

IrrepBasisSet* SL_OBasis::GetCBasisSet () const
{
    return new Slater::IrrepBasisSet(lap,idb,N,emin*2.0,emax*2.0,0);
}

IrrepBasisSet* SL_OBasis::GetXBasisSet () const
{
    return new Slater::IrrepBasisSet(lap,idb,N,emin*2.0/3.0,emax*2.0/3.0,0);
}


#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Gaussian94.H"

BasisSet* PG_OBasis::GetBasisSet (const Cluster* cl) const
{
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/dzvp.bsd");
    auto bs=new PolarizedGaussian::BasisSet(lap, &reader,cl);
    idb=bs->GetDataBase();
    return  bs;
}

IrrepBasisSet* PG_OBasis::GetCBasisSet () const
{
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/A2_coul.bsd");
    return new PolarizedGaussian::IrrepBasisSet(lap,idb, &reader,GetCluster());
}

IrrepBasisSet* PG_OBasis::GetXBasisSet () const
{
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/A1_exch.bsd");
    return new PolarizedGaussian::IrrepBasisSet(lap,idb, &reader,GetCluster());
}


#include "Cluster/Molecule.H"

Cluster* TestAtom::GetCluster() const
{
    Cluster* cl=new Molecule;
    cl->Insert(new Atom(Z,q,Vector3D<double>(0,0,0)));
    return cl;
}

#include "Mesh/RadialMesh/MHLRadialMesh.H"
#include "Mesh/AngularMesh/GaussAngularMesh.H"
#include "Mesh/AtomMesh.H"

Mesh* TestAtom::GetIntegrationMesh() const
{
    RadialMesh*            rm=new MHLRadialMesh(50,2U,2.0); //mem leak
    AngularMesh*           am=new GaussAngularMesh(1);      //mem leak
    return new AtomMesh(*rm,*am); //why not own?
}

void TestMolecule::Init(Molecule* p)
{
    assert(p);
    itsCluster=p;
}

Cluster* TestMolecule::GetCluster() const {return itsCluster;}

#include "Mesh/MoleculeMesh.H"
Mesh*    TestMolecule::GetIntegrationMesh() const
{
    return  new MoleculeMesh(*itsCluster,3);
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

