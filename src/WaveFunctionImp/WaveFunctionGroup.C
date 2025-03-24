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
    for (auto b:itsBasisSet->Iterate<TOrbital_IBS<double> >())
        itsIrrepWFs.push_back(new IrrepWaveFunction(b,S));
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

Exact_CD* WaveFunctionGroup::GetChargeDensity(Spin s) const
{
    Composite_Exact_CD* cd = new Composite_Exact_CD();
    for (auto w:itsIrrepWFs) cd->Insert(w->GetChargeDensity(s));
    return cd;
}

Orbitals* WaveFunctionGroup::GetOrbitals(const QuantumNumber& qn, Spin s) const
{

    auto w=itsIrrepWFs.begin();
    for (auto b:itsBasisSet->Iterate<Orbital_IBS>())
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


SCFIterator* WaveFunctionGroup::MakeIterator(Hamiltonian* H, Exact_CD* cd, double NElectrons)
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


