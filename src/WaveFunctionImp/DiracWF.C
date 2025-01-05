// File: DiracWF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/DiracWaveFunction.H"
#include "Imp/WaveFunction/MasterPolarizedWF.H"
#include "Imp/SCFIterator/SCFIteratorPol.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include <ChargeDensity.H>
#include <Spin.H>
#include <QuantumNumber.H>
#include "oml/imp/binio.h"
#include <cassert>

DiracWF::DiracWF()
    : itsLargeComponent(0)
    , itsSmallComponent(0)
    , itsEC      (0)
{};

DiracWF::DiracWF(const BasisSet* bg,const ElectronConfiguration* ec)
    : itsLargeComponent(new MasterPolarizedWF(bg,ec))
    , itsSmallComponent(new MasterPolarizedWF(bg,ec))
    , itsEC           (ec) //Electron cofiguration
{
    assert(itsLargeComponent);
    assert(itsSmallComponent);
    assert(itsEC);
};

DiracWF::~DiracWF()
{
    delete itsLargeComponent;
    delete itsSmallComponent;
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void DiracWF::DoSCFIteration(Hamiltonian& ham)
{
    assert(itsLargeComponent);
    assert(itsSmallComponent);
    itsLargeComponent->DoSCFIteration(ham);
    itsSmallComponent->DoSCFIteration(ham);
}

ChargeDensity* DiracWF::GetChargeDensity(Spin s) const
{
    assert(itsLargeComponent);
    assert(itsSmallComponent);
    ChargeDensity* large=itsLargeComponent->GetChargeDensity(s);
    ChargeDensity* small=itsSmallComponent->GetChargeDensity(s);
    large->MixIn(*small,0.5);
    large->ReScale(2.0);
    return large;
}


const EnergyLevels& DiracWF::FillOrbitals(const ElectronConfiguration*)
{
    assert(itsLargeComponent);
    assert(itsSmallComponent);
    itsLargeLevels=itsLargeComponent->FillOrbitals(itsEC);
    itsSmallLevels=itsSmallComponent->FillOrbitals(itsEC);
    return itsLargeLevels;
}


SCFIterator* DiracWF::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons)
{
    return new SCFIteratorPol(this,H,cd,nElectrons,0.0);
}

WaveFunction* DiracWF::GetWaveFunction(const Spin& S)
{
    return itsLargeComponent->GetWaveFunction(S);
}

void DiracWF::DisplayEigen() const
{
    itsLargeComponent->DisplayEigen();
}

std::ostream& DiracWF::Write(std::ostream& os) const
{
    assert(itsLargeComponent);
    assert(itsSmallComponent);
    if (Pretty())
        os << "Dirac wave function :" << std::endl;
    os << *itsLargeComponent << *itsSmallComponent;
    
    return os;
}

std::istream& DiracWF::Read (std::istream& is)
{
    return is;
}


