// File: MasterUnPolarizedWF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/MasterUnPolarizedWF.H"
#include "Imp/WaveFunction/WaveFunctionGroup.H"
#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include <Spin.H>
#include <cassert>
#include <iostream>

MasterUnPolarizedWF::MasterUnPolarizedWF()
    : itsGroup(0)
    , itsEC(0)
{};

MasterUnPolarizedWF::MasterUnPolarizedWF(const BasisSet* bg,const ElectronConfiguration* ec)
    :itsGroup(new WaveFunctionGroup(bg,Spin(Spin::None)))
    ,itsEC(ec)
{
    assert(itsGroup);
    assert(itsEC);
};

MasterUnPolarizedWF::~MasterUnPolarizedWF()
{
    delete itsGroup;
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
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

Orbitals* MasterUnPolarizedWF::GetOrbitals(const QuantumNumber& qn,Spin s) const
{
    assert(itsGroup);
    assert(s==Spin::None);
    return itsGroup->GetOrbitals(qn,s);
}

SCFIterator* MasterUnPolarizedWF::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons)
{
    return new SCFIteratorUnPol(this, H, cd,nElectrons);
}

const EnergyLevels& MasterUnPolarizedWF::FillOrbitals(const ElectronConfiguration*)
{
    assert(itsGroup);
    return itsGroup->FillOrbitals(itsEC);
}

void MasterUnPolarizedWF::DisplayEigen() const
{
    std::cout << "Alpha+Beta spin :" << std::endl;
    itsGroup->DisplayEigen();
}


std::ostream& MasterUnPolarizedWF::Write(std::ostream& os) const
{
    assert(itsGroup);
    os << *itsGroup;
    return os;
}

std::istream& MasterUnPolarizedWF::Read (std::istream& is)
{
    delete itsGroup;
    itsGroup=WaveFunction::Factory(is);
    assert(itsGroup);
    is >> *itsGroup;
    return is;
}

