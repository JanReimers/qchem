// File: WaveFunctionGroup.C  Wave function for an unpolarized atom.



#include "ChargeDensity.H"
#include <BasisSet.H>
#include <Irrep_BS.H>
#include "Imp/WaveFunction/WaveFunctionGroup.H"
#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/ChargeDensity/CompositeCD.H"
#include "Imp/Orbitals/TOrbitals.H"
#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include "Imp/Containers/ptr_vector_io.h"
#include <QuantumNumber.H>
#include <cassert>


WaveFunctionGroup::WaveFunctionGroup()
{};

WaveFunctionGroup::WaveFunctionGroup(const BasisSet* bg, const Spin& S)
: itsBasisSet(bg)
{
    assert(itsBasisSet);
    assert(itsBasisSet->GetNumFunctions()>0);
    for (auto b:*itsBasisSet)
    {
        auto tbs=dynamic_cast<const TOrbital_IBS<double>*>(b); //TODO avoid casting here?
        assert(tbs);
        itsIrrepWFs.push_back(new IrrepWaveFunction(tbs,S));
    }
};

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
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

Orbitals* WaveFunctionGroup::GetOrbitals(const QuantumNumber& qn, Spin s) const
{

    auto w=itsIrrepWFs.begin();
    for (auto b:*itsBasisSet)
    {
        if (qn==b->GetQuantumNumber()) 
            return (*w)->GetOrbitals(qn,s); 
        w++;
    }
    assert(false);
    return NULL;
}

const EnergyLevels& WaveFunctionGroup::FillOrbitals(const ElectronConfiguration* ec)
{
    itsELevels.clear();
    for (auto w:itsIrrepWFs) 
        itsELevels.merge(w->FillOrbitals(ec),0.0001);
    return itsELevels;
}


SCFIterator* WaveFunctionGroup::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double NElectrons)
{
    return new SCFIteratorUnPol(this, H, cd,NElectrons);
}

void WaveFunctionGroup::DisplayEigen() const
{
    itsELevels.Report(std::cout);
//    for (auto w:itsIrrepWFs) w->DisplayEigen();
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


