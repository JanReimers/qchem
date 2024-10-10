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
    , itsBasisGroup(0)
    , itsHamiltonian(0)
    , itsWaveFunction(0)
    , itsSCFIterator(0)
{
};

BaseTester::~BaseTester()
{
    //delete itsHamiltonian;
    delete itsBasisGroup;
}

void BaseTester::Init(BasisGroup* bg,double spin,const LinearAlgebraParams& lap)
{
    itsLAParams=lap;
    Init(bg,spin);
}

void BaseTester::Init(BasisGroup* bg,double spin)
{
    assert(bg); //Derived should already have created this.
    assert(itsCluster->GetNumAtoms()>0);
    //if(itsHamiltonian) delete itsHamiltonian; SCFIterator owns the Hamiltonian
    if(itsBasisGroup) delete itsBasisGroup;
    if (itsWaveFunction) delete itsWaveFunction;
    if (itsSCFIterator) delete itsSCFIterator;
//    itsFittedChargeDensity=new FittedCDImplementation<double>(itsCBasisSet,itsCluster->GetNumElectrons ());
    itsBasisGroup=bg;
    LoadOrbitalBasisSet();

    itsHamiltonian=new HamiltonianImplementation();
    
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
IntegralDataBase<double>*   BaseTester::GetDatabase() const 
{
    return itsBasisGroup->GetDataBase();
}

void BaseTester::Iterate(const SCFIterationParams& ipar)
{
    itsSCFIterator->Iterate(ipar);
}

double BaseTester::TotalEnergy() const
{
    return itsHamiltonian->GetTotalEnergy().GetTotalEnergy();
}
