// File: IrrepWaveFunction.C  Wave function for an unpolarized atom.



#include "WaveFunctionImp/IrrepWaveFunction/IrrepWaveFunction.H"
#include "WaveFunctionImp/MasterWF/UnPolarizedSCFIterator.H"
#include "Hamiltonian.H"
#include "Orbital.H"
#include "Orbital/ElectronDumper.H"
#include "ChargeDensityImplementation/ExactIrrepCD/ExactIrrepCD.H"
#include "oml/imp/binio.h"
#include "oml/smatrix.h"
#include <cassert>

IrrepWaveFunction::IrrepWaveFunction()
    : itsOrbitals(0)
    , itsBasisSet(0)
    , itsSpin    ( )
{};

IrrepWaveFunction::IrrepWaveFunction(const BasisSet* bs, const Spin& S)
    : itsOrbitals(0 )
    , itsBasisSet(bs)
    , itsSpin    (S )
{};

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
    if (itsOrbitals) delete itsOrbitals;
    itsOrbitals = itsBasisSet->CreateOrbitals(itsBasisSet,&ham,itsSpin);
}

ChargeDensity* IrrepWaveFunction::GetChargeDensity(Spin s) const
{
    ChargeDensity* cd=0;
    if (itsOrbitals)
    {
        cd=itsOrbitals->GetChargeDensity(s);
    }
    else
    {
        int n=itsBasisSet->GetNumFunctions();
        SMatrix<double> D(n,n);
        Fill(D,0.0);
        cd=new ExactIrrepCD<double>(D,itsBasisSet,s);
    }
    assert(cd);
    return cd;
}

void IrrepWaveFunction::UpdateElectronDumper(ElectronDumper& ed)
{
    assert(itsOrbitals);
    ed.Add(itsOrbitals);
}

SCFIterator* IrrepWaveFunction::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons, double kT, bool showplot)
{
    return new UnPolarizedSCFIterator(this, H, cd,nElectrons, kT, showplot);
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
    os << *itsBasisSet; 

    return os;
}

std::istream& IrrepWaveFunction::Read (std::istream& is)
{
    is >> itsSpin;

    delete itsOrbitals;
    itsOrbitals = OrbitalGroup::Factory(is);
    assert(itsOrbitals);
    is >> *itsOrbitals;

    BasisSet* temp=BasisSet::Factory(is);
    is >> *temp;
    itsBasisSet.reset(temp);
    FixUpPointer(itsOrbitals,itsBasisSet);


    return is;
}

