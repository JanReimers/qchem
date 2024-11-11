// File: IrrepWaveFunction.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/WaveFunction/ElectronDumper.H"
#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <Hamiltonian.H>
#include "oml/imp/binio.h"
#include <cassert>

IrrepWaveFunction::IrrepWaveFunction()
    : itsOrbitals(0)
    , itsSpin    ( )
{};

IrrepWaveFunction::IrrepWaveFunction(const IrrepBasisSet* bs, const Spin& S)
    : itsOrbitals(new  TOrbitalsImp<double>(dynamic_cast<const TIrrepBasisSet<double>*>(bs)))
    , itsSpin    (S )
{
    assert(itsOrbitals);
};

IrrepWaveFunction::~IrrepWaveFunction()
{
    delete itsOrbitals;
}

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the ElectronDumper
//  to fill up the orbitals with electrons.
//
void IrrepWaveFunction::DoSCFIteration(Hamiltonian& ham)
{
    assert(itsOrbitals);
    itsOrbitals->UpdateOrbitals(ham,itsSpin);
}

ChargeDensity* IrrepWaveFunction::GetChargeDensity(Spin s) const
{
    assert(itsOrbitals);
    return itsOrbitals->GetChargeDensity(s);
}

void IrrepWaveFunction::FillOrbitals(const ElectronConfiguration* ec, const Spin& s)
{
    ElectronDumper ed(0.0001,0.0);
    ed.Add(itsOrbitals);
    ed.MakeEnergyLevels();
    Orbital* o1=*itsOrbitals->begin();
    ed.DumpInElectrons(ec->GetN(o1->GetQuantumNumber(),s));
}

void  IrrepWaveFunction::DisplayEigen() const
{
    itsOrbitals->DisplayEigen();
}

SCFIterator* IrrepWaveFunction::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons)
{
    return new SCFIteratorUnPol(this, H, cd,nElectrons);
}

std::string spin_strs[]={"Down","None","Up"};
std::ostream& IrrepWaveFunction::Write(std::ostream& os) const
{
    assert(itsOrbitals);
    if (Pretty()) 
        os << "    Irreducible rep. Wave finction, spin=" << spin_strs[itsSpin.itsState] << std::endl;
    else
        os << itsSpin;
        
    os << *itsOrbitals;
    if (Pretty()) os << "        ";

    return os;
}

std::istream& IrrepWaveFunction::Read (std::istream& is)
{
    is >> itsSpin;

    delete itsOrbitals;
    itsOrbitals = Orbitals::Factory(is);
    assert(itsOrbitals);
    is >> *itsOrbitals;

    return is;
}

