// File: UnPolarized_WF.C  Wave function for an unpolarized atom.

#include "Imp/WaveFunction/UnPolarized_WF.H"
#include <iostream>

UnPolarized_WF::UnPolarized_WF(const BasisSet* bs,const ElectronConfiguration* ec,SCFAccelerator* acc)
    : Composite_WF(bs,ec,acc)
{
    MakeIrrep_WFs(Spin::None);
};


void UnPolarized_WF::DisplayEigen() const
{
    StreamableObject::SetToPretty();

    std::cout << "Alpha+Beta spin :" << std::endl;
    GetEnergyLevels().Report(std::cout);
}
