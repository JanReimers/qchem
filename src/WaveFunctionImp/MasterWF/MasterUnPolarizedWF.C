// File: MasterUnPolarizedWF.C  Wave function for an unpolarized atom.



#include "WaveFunctionImp/MasterWF/MasterUnPolarizedWF.H"
#include "WaveFunctionImp/MasterWF/UnPolarizedSCFIterator.H"
#include "WaveFunctionImp/WaveFunctionGroup/WaveFunctionGroup.H"
#include "Misc/Spin.H"
#include "Misc/ptr_vector.h"
#include "Misc/ptrvector_io.h"
#include <cassert>

MasterUnPolarizedWF::MasterUnPolarizedWF()
    : itsGroup(0)
{};

MasterUnPolarizedWF::MasterUnPolarizedWF(const BasisGroup* bg)
    :itsGroup(new WaveFunctionGroup(bg,Spin(Spin::None)))
{
    assert(itsGroup);
};

MasterUnPolarizedWF::~MasterUnPolarizedWF()
{
    delete itsGroup;
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the ElectronDumper
//  to fill up the orbitals with electrons.
//
void MasterUnPolarizedWF::DoSCFIteration(Hamiltonian& ham)
{
    assert(itsGroup);
    itsGroup->DoSCFIteration(ham);
}

ChargeDensity* MasterUnPolarizedWF::GetChargeDensity(Spin s) const
{
    assert(itsGroup);
    assert(s==Spin::None);
    return itsGroup->GetChargeDensity(s);
}

SCFIterator* MasterUnPolarizedWF::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons, double kT, bool showplot)
{
    return new UnPolarizedSCFIterator(this, H, cd,nElectrons, kT,showplot);
}

void MasterUnPolarizedWF::UpdateElectronDumper(ElectronDumper& ed)
{
    assert(itsGroup);
    itsGroup->UpdateElectronDumper(ed);
}

std::ostream& MasterUnPolarizedWF::Write(std::ostream& os) const
{
    assert(itsGroup);
    os << itsGroup;
    return os;
}

std::istream& MasterUnPolarizedWF::Read (std::istream& is)
{
    delete itsGroup;
    itsGroup=WaveFunction::Factory(is);
    assert(itsGroup);
    is >> itsGroup;
    return is;
}

