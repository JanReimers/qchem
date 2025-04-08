// File: Irrep_WF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <Irrep_BS.H>
#include <cassert>
#include <iostream>



Irrep_WF::Irrep_WF(const TOrbital_IBS<double>* bs, const Spin& ms)
    : itsOrbitals(new  TOrbitalsImp<double>(bs,ms))
    , itsQNs     (ms,&bs->GetSymmetry())
{
    assert(itsOrbitals);
};

Irrep_WF::~Irrep_WF()
{
    delete itsOrbitals;
}

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  
//
void Irrep_WF::DoSCFIteration(Hamiltonian& ham,const DM_CD* cd)
{
    assert(itsOrbitals);
    itsOrbitals->UpdateOrbitals(ham,itsQNs.ms,cd);
}

DM_CD* Irrep_WF::GetChargeDensity() const
{
    assert(itsOrbitals);
    return itsOrbitals->GetChargeDensity();
}

const Orbitals* Irrep_WF::GetOrbitals() const
{
    assert(itsOrbitals);
    return itsOrbitals;
}
//
//  There are three steps here:
//
const EnergyLevels& Irrep_WF::FillOrbitals(const ElectronConfiguration* ec)
{
    // Step one: How many electron for this Irrep(qn,spin) ?
    double ne=ec->GetN(*itsQNs.sym,itsQNs.ms);
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
        itsELevels.insert(o->MakeEnergyLevel(itsQNs.ms));
    
    //  Display the occupied orbitals with eigen vectors.
    // for (auto o:*itsOrbitals)
    //     if (o->GetOccupation()>0.0)
    //         std::cout << *o << std::endl;

    return itsELevels;
}

void  Irrep_WF::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
}

