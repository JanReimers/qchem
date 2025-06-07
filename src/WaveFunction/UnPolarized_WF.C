// File: UnPolarized_WF.C  Wave function for an unpolarized atom.

#include "Imp/WaveFunction/UnPolarized_WF.H"
#include "Imp/WaveFunction/Irrep_WF.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/Orbitals/TOrbitals.H"
#include "Imp/SCFAccelerator.H"
#include <BasisSet.H>
#include <Irrep_BS.H>
#include <cassert>
#include <Spin.H>
#include <Orbital_QNs.H>
#include <iostream>

UnPolarized_WF::UnPolarized_WF()
   
{};

UnPolarized_WF::UnPolarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator& acc)
    : Composite_WF(bs,ec)
{
    MakeIrrep_WFs(acc,Spin::None);
};


DM_CD* UnPolarized_WF::GetChargeDensity() const
{
    return GetChargeDensity(Spin::None);
}

void UnPolarized_WF::DisplayEigen() const
{
    StreamableObject::SetToPretty();

    std::cout << "Alpha+Beta spin :" << std::endl;
    GetEnergyLevels().Report(std::cout);
}
