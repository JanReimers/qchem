// File: WaveFunctionGroup.C  Wave function for an unpolarized atom.



#include "WaveFunctionImp/WaveFunctionGroup/WaveFunctionGroup.H"
#include "ChargeDensity/ChargeDensity.H"
#include "BasisSet/TBasisSet.H"
#include "BasisSet/BasisGroup.H"
#include "WaveFunctionImp/IrrepWaveFunction/IrrepWaveFunction.H"
#include "ChargeDensityImplementation/CompositeCD/CompositeCD.H"
#include "Hamiltonian/TotalEnergy.H"
#include "WaveFunctionImp/MasterWF/UnPolarizedSCFIterator.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "Misc/ptrvector_io.h"
#include "Misc/ptr_vector.h"
#include <cassert>


WaveFunctionGroup::WaveFunctionGroup()
{};

WaveFunctionGroup::WaveFunctionGroup(const BasisGroup* bg, const Spin& S)
: itsBasisGroup(bg)
{
    assert(itsBasisGroup);
    for (auto b:*itsBasisGroup)
    {
        const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(b); //TODO avoid casting here?
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
    for (optr_vector<WaveFunction*>::iterator i(itsIrrepWFs.begin()); i!=itsIrrepWFs.end(); i++)
        i->DoSCFIteration(ham);
}

ChargeDensity* WaveFunctionGroup::GetChargeDensity(Spin s) const
{
    CompositeCD* cd = new CompositeCD();
    optr_vector<WaveFunction*>::const_iterator sb(itsIrrepWFs.begin());
    for (; sb!=itsIrrepWFs.end(); sb++)
        cd->Insert(sb->GetChargeDensity(s));
    return cd;
}

void WaveFunctionGroup::UpdateElectronDumper(ElectronDumper& ed)
{
    for (optr_vector<WaveFunction*>::iterator i(itsIrrepWFs.begin()); i!=itsIrrepWFs.end(); i++)
        i->UpdateElectronDumper(ed);
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


