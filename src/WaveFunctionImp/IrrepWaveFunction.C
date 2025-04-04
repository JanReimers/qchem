// File: IrrepWaveFunction.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <Hamiltonian.H>
#include <Irrep_BS.H>
#include <Symmetry.H>
#include "oml/imp/binio.h"
#include <cassert>
#include <iomanip>



IrrepWaveFunction::IrrepWaveFunction(const TOrbital_IBS<double>* bs, const Spin& ms)
    : itsOrbitals(new  TOrbitalsImp<double>(bs,ms))
    , itsSpin    (ms )
    , itsQN      (&bs->GetQuantumNumber())
    , itsQNs     (ms,&bs->GetQuantumNumber())
{
    assert(itsOrbitals);
};

IrrepWaveFunction::~IrrepWaveFunction()
{
    delete itsOrbitals;
}

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  
//
void IrrepWaveFunction::DoSCFIteration(Hamiltonian& ham)
{
    assert(itsOrbitals);
    itsOrbitals->UpdateOrbitals(ham,itsSpin);
}

Exact_CD* IrrepWaveFunction::GetChargeDensity() const
{
    assert(itsOrbitals);
    return itsOrbitals->GetChargeDensity();
}

Orbitals* IrrepWaveFunction::GetOrbitals(const Irrep_QNs& qns) const
{
    assert(itsOrbitals);
    assert(qns==itsQNs);

    return itsOrbitals;
}
//
//  There are three steps here:
//
const EnergyLevels& IrrepWaveFunction::FillOrbitals(const ElectronConfiguration* ec)
{
    // Step one: How many electron for this Irrep(qn,spin) ?
    double ne=ec->GetN(*itsQN,itsSpin);
    //std::cout << "ne=" << ne << " QN=" << *itsQN << std::endl;
    //  Loop over orbitals and consume the electrons quota.
    for (auto& o:*itsOrbitals)
    {
        ne=o->TakeElectrons(ne);
        if (ne<=0.0) break;
    }
    //  Now update the list of energy levels.
    itsELevels.clear();
    for (auto o:*itsOrbitals)
        itsELevels.insert(o->MakeEnergyLevel(itsSpin));
    
    //  Display the occupied orbitals with eigen vectors.
    // for (auto o:*itsOrbitals)
    //     if (o->GetOccupation()>0.0)
    //         std::cout << *o << std::endl;

    return itsELevels;
}

void  IrrepWaveFunction::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
}

