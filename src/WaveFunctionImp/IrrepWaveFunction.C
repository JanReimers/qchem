// File: IrrepWaveFunction.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/IrrepWaveFunction.H"
#include "Imp/WaveFunction/ElectronConfiguration.H"
#include "Imp/SCFIterator/SCFIteratorUnPol.H"
#include "Imp/Orbitals/TOrbitals.H"
#include <Hamiltonian.H>
#include <Irrep_BS.H>
#include <QuantumNumber.H>
#include "oml/imp/binio.h"
#include <cassert>
#include <iomanip>

IrrepWaveFunction::IrrepWaveFunction()
    : itsOrbitals(0)
    , itsSpin    ( )
{};

IrrepWaveFunction::IrrepWaveFunction(const TOrbital_IBS<double>* bs, const Spin& S)
    : itsOrbitals(new  TOrbitalsImp<double>(bs))
    , itsSpin    (S )
    , itsQN      (&bs->GetQuantumNumber())
{
    assert(itsOrbitals);
};

IrrepWaveFunction::~IrrepWaveFunction()
{
    delete itsOrbitals;
}

//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  
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

Orbitals* IrrepWaveFunction::GetOrbitals(const QuantumNumber& qn, Spin s) const
{
    assert(itsOrbitals);
    assert(qn==*itsQN);

    return itsOrbitals;
}
//
//  There are three steps here:
//
const EnergyLevels& IrrepWaveFunction::FillOrbitals(const ElectronConfiguration* ec)
{
    // Step one: How many electron for this Irrep(qn,spin) ?
    double ne=ec->GetN(*itsQN,itsSpin);
    //std::cout << "ne=" << ne << " QN=" << *itsQN << std::endl;
    //  Loop over orbitals and consume the electrons quota.
    for (auto& o:*itsOrbitals)
    {
        ne=o->TakeElectrons(ne);
        if (ne<=0.0) break;
    }
    //  Now update the list of energy levels.
    itsELevels.clear();
    for (auto o:*itsOrbitals)
        itsELevels.insert(o->MakeEnergyLevel(itsSpin));
    
    //  Display the occupied orbitals with eigen vectors.
    // for (auto o:*itsOrbitals)
    //     if (o->GetOccupation()>0.0)
    //         std::cout << *o << std::endl;

    return itsELevels;
}

void  IrrepWaveFunction::DisplayEigen() const
{
    itsELevels.Report(std::cout);
   
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

