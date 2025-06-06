// File: Irrep_WF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/Irrep_WF.H"
#include <ElectronConfiguration.H>
#include "Imp/Orbitals/TOrbitals.H"
#include <Irrep_BS.H>
#include <cassert>
#include <iostream>



Irrep_WF::Irrep_WF(const TOrbital_IBS<double>* bs, const Spin& ms,SCFIrrepAccelerator* acc)
    : itsOrbitals(new  TOrbitalsImp<double>(bs,ms,acc))
    , itsIrrep     (ms,&bs->GetSymmetry())
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
    itsOrbitals->UpdateOrbitals(ham,cd);
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
Orbitals* Irrep_WF::GetOrbitals() 
{
    assert(itsOrbitals);
    return itsOrbitals;
}
//
//  There are three steps here:
//
const EnergyLevels& Irrep_WF::FillOrbitals(const ElectronConfiguration* ec)
{
    
    double ne=ec->GetN(itsIrrep); // Step one: How many electron for this Irrep(qn,spin) ?
    ne=itsOrbitals->TakeElectrons(ne); // Step two: Dump electrons into the orbitals
    assert(ne==0.0); 
    
    // Step three: Make a list of energy levels.  Degenerate levels should get merged.
    itsELevels.clear();
    for (auto o:itsOrbitals->Iterate<Orbital>())
        itsELevels.insert(EnergyLevel(o));
    
    //  Display the occupied orbitals with eigen vectors.
    // for (auto o:itsOrbitals->Iterate<Orbital>())
    //     if (o->GetOccupation()>0.0)
    //         std::cout << *o << std::endl;
    // for (auto el:itsELevels)
    //     if (el.second.occ>0)
    //     {
    //         el.second.Report(std::cout);
    //         std::cout << std::endl;
    //     }

    return itsELevels;
}

void  Irrep_WF::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
}

