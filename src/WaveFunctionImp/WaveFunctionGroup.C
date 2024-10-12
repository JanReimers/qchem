// File: WaveFunctionGroup.C  Wave function for an unpolarized atom.



#include "ChargeDensity.H"
#include "BasisSet.H"
#include "Imp/WaveFunction/WaveFunctionGroup.H"
#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "ChargeDensityImplementation/CompositeCD/CompositeCD.H"
#include "Imp/SCFIterator/UnPolarizedSCFIterator.H"
#include "Imp/Containers/ptr_vector_io.h"
#include <cassert>


WaveFunctionGroup::WaveFunctionGroup()
{};

WaveFunctionGroup::WaveFunctionGroup(const BasisSet* bg, const Spin& S)
: itsBasisGroup(bg)
{
    assert(itsBasisGroup);
    for (auto b:*itsBasisGroup)
    {
        const TIrrepBasisSet<double>* tbs=dynamic_cast<const TIrrepBasisSet<double>*>(b); //TODO avoid casting here?
        assert(tbs);
        itsIrrepWFs.push_back(new IrrepWaveFunction(tbs,S));
    }
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the ElectronDumper
//  to fill up the orbitals with electrons.
//
void WaveFunctionGroup::DoSCFIteration(Hamiltonian& ham)
{
    for (auto w:itsIrrepWFs) w->DoSCFIteration(ham);
}

ChargeDensity* WaveFunctionGroup::GetChargeDensity(Spin s) const
{
    CompositeCD* cd = new CompositeCD();
    for (auto w:itsIrrepWFs) cd->Insert(w->GetChargeDensity(s));
    return cd;
}

void WaveFunctionGroup::UpdateElectronDumper(ElectronDumper& ed)
{
    for (auto w:itsIrrepWFs) w->UpdateElectronDumper(ed);
}

SCFIterator* WaveFunctionGroup::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double NElectrons, double kT, bool showplot)
{
    return new UnPolarizedSCFIterator(this, H, cd,NElectrons,kT,showplot);
}

std::ostream& WaveFunctionGroup::Write(std::ostream& os) const
{
    os << itsIrrepWFs;
    return os;
}

std::istream& WaveFunctionGroup::Read (std::istream& is)
{
    is >> itsIrrepWFs;
    return is;
}


