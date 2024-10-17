
#include "QchemTester.H"
#include "SCFIterator/SCFIterator.H"
#include <WaveFunction.H>
#include <Hamiltonian.H>
#include <Cluster.H>
#include <TotalEnergy.H>


QchemTester::QchemTester()
: itsCluster(0)
, itsBasisSet(0)
, itsHamiltonian(0)
, itsWaveFunction(0)
, itsSCFIterator(0)
, MaxRelErrE(1e-3)
{
    //Cannot call virtual functions from here.
}

void QchemTester::Init()
{
    itsCluster=GetCluster(); //Atom or Molecule
    assert(itsCluster);
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
    double E_HF=itsPT.GetEnergyHF(2);
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

#include "HamiltonianImplementation/PolarizedHartreeFockVxc.H"
HamiltonianTerm* PolHFHamiltonian:: GetVxc() const
{
    return new PolarizedHartreeFockVxc;
}

#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"

BasisSet* SG_OBasis::GetBasisSet (const Cluster*) const
{
    return new SphericalGaussian::BasisSet(lap,N,emin,emax,Lmax); 
}

#include "Cluster/Molecule.H"

Cluster* TestAtom::GetCluster() const
{
    Cluster* cl=new Molecule;
    cl->Insert(new Atom(Z,q,Vector3D<double>(0,0,0)));
    return cl;
}

#include "Imp/WaveFunction/MasterUnPolarizedWF.H"
WaveFunction* TestUnPolarized::GetWaveFunction(const BasisSet* bs) const
{
    return new MasterUnPolarizedWF(bs);
}

#include "Imp/WaveFunction/MasterPolarizedWF.H"
WaveFunction* TestPolarized::GetWaveFunction(const BasisSet* bs) const
{
    return new MasterPolarizedWF(bs,spin);
}

