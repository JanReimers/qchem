// File AtomTester.C Member functions for the atom tester class.

#include "BaseTester.H"
#include "BasisSet/BasisGroup.H"
#include "HamiltonianImplementation/HamiltonianImplementation.H"
#include "HamiltonianImplementation/ExactVen.H"
#include "HamiltonianImplementation/ExactVnn.H"
#include "HamiltonianImplementation/Kinetic.H"
#include "WaveFunctionImp/MasterWF/MasterUnPolarizedWF.H"
#include "WaveFunctionImp/MasterWF/MasterPolarizedWF.H"
#include "WaveFunction/SCFIterator.H"
#include "Cluster/Molecule.H"

BaseTester::BaseTester()
    : itsCluster(new Molecule())
    , itsBasisGroup(new BasisGroup)
    , itsHamiltonian(new HamiltonianImplementation)
{
};

BaseTester::~BaseTester()
{
    delete itsHamiltonian;
    delete itsBasisGroup;
}

void BaseTester::Init(double spin)
{
    assert(itsCluster->GetNumAtoms()>0);
//    itsFittedChargeDensity=new FittedCDImplementation<double>(itsCBasisSet,itsCluster->GetNumElectrons ());
    LoadOrbitalBasisSet();

    itsHamiltonian->Add(new Kinetic);
    itsHamiltonian->Add(new ExactVnn(itsCluster));
    itsHamiltonian->Add(new ExactVen(itsCluster));
    itsHamiltonian->Add(GetVee());
    itsHamiltonian->Add(GetVxc(spin));
    if (spin==0.0)
    {
        itsWaveFunction=new MasterUnPolarizedWF(itsBasisGroup);
    }
    else
    {
        itsWaveFunction=new MasterPolarizedWF(itsBasisGroup,spin);
    }

    ChargeDensity* GuessCD=itsWaveFunction->GetChargeDensity();
    itsSCFIterator=itsWaveFunction->
                   MakeIterator(itsHamiltonian,GuessCD,itsCluster->GetNumElectrons(),0.0,false); //show plot flag.
}

void BaseTester::Iterate(double relax, double epsilon, int Niter)
{
    itsSCFIterator->Iterate(relax,epsilon,Niter,0.0);
}

double BaseTester::TotalEnergy() const
{
    return itsHamiltonian->GetTotalEnergy().GetTotalEnergy();
}
