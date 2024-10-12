// File AtomTester.C Member functions for the atom tester class.

#include "BaseTester.H"
#include "BasisSet.H"
#include "TotalEnergy.H"
#include "HamiltonianImplementation/HamiltonianImplementation.H"
#include "HamiltonianImplementation/ExactVen.H"
#include "HamiltonianImplementation/ExactVnn.H"
#include "HamiltonianImplementation/Kinetic.H"
#include "Imp/WaveFunction/MasterUnPolarizedWF.H"
#include "Imp/WaveFunction/MasterPolarizedWF.H"
#include "SCFIterator/SCFIterator.H"
#include "Cluster/Molecule.H"

BaseTester::BaseTester()
    : itsLAParams({qchem::Lapack,qchem::SVD,1e-6,1e-12})
    , itsCluster(new Molecule())
    , itsBasisSet(0)
    , itsHamiltonian(0)
    , itsWaveFunction(0)
    , itsSCFIterator(0)
{
};

BaseTester::~BaseTester()
{
    //delete itsHamiltonian; //SCFIterator owns the Hamiltonian
    delete itsWaveFunction;
    delete itsBasisSet;
    delete itsSCFIterator;
    
}

void BaseTester::Init(BasisSet* bs,double spin,const LAParams& lap)
{
    itsLAParams=lap;
    Init(bs,spin);
}

void BaseTester::Init(BasisSet* bs,double spin)
{
    assert(bs); //Derived should already have created this.
    assert(itsCluster->GetNumAtoms()>0);
    //if(itsHamiltonian) delete itsHamiltonian; SCFIterator owns the Hamiltonian
    if (itsBasisSet) delete itsBasisSet;
    if (itsWaveFunction) delete itsWaveFunction;
    if (itsSCFIterator) delete itsSCFIterator;
//    itsFittedChargeDensity=new FittedCDImplementation<double>(itsCBasisSet,itsCluster->GetNumElectrons ());
    itsBasisSet=bs;
    LoadOrbitalBasisSet();

    itsHamiltonian=new HamiltonianImplementation();
    
    itsHamiltonian->Add(new Kinetic);
    itsHamiltonian->Add(new ExactVnn(itsCluster));
    itsHamiltonian->Add(new ExactVen(itsCluster));
    itsHamiltonian->Add(GetVee());
    itsHamiltonian->Add(GetVxc(spin));
    if (spin==0.0)
    {
        itsWaveFunction=new MasterUnPolarizedWF(itsBasisSet);
    }
    else
    {
        itsWaveFunction=new MasterPolarizedWF(itsBasisSet,spin);
    }

    ChargeDensity* GuessCD=itsWaveFunction->GetChargeDensity();
    itsSCFIterator=itsWaveFunction->
                   MakeIterator(itsHamiltonian,GuessCD,itsCluster->GetNumElectrons(),0.0,false); //show plot flag.
}
IntegralDataBase<double>*   BaseTester::GetDatabase() const 
{
    return itsBasisSet->GetDataBase();
}

void BaseTester::Iterate(const SCFIterationParams& ipar)
{
    itsSCFIterator->Iterate(ipar);
}

double BaseTester::TotalEnergy() const
{
    return itsHamiltonian->GetTotalEnergy().GetTotalEnergy();
}
