// File: MasterPolarizedWF.C  Wave function for an unpolarized atom.



#include "Imp/WaveFunction/MasterPolarizedWF.H"
#include "Imp/WaveFunction/WaveFunctionGroup.H"
#include "Imp/SCFIterator/SCFIteratorPol.H"
#include "Imp/ChargeDensity/PolarizedCD.H"
#include <ChargeDensity.H>
#include <Spin.H>
#include "oml/imp/binio.h"
#include <cassert>

MasterPolarizedWF::MasterPolarizedWF()
    : itsSpinUpGroup  (0)
    , itsSpinDnGroup(0)
    , itsEC      (0)
{};

MasterPolarizedWF::MasterPolarizedWF(const BasisSet* bg,const ElectronConfiguration* ec)
    : itsSpinUpGroup  (new WaveFunctionGroup(bg,Spin(Spin::Up  )))
    , itsSpinDnGroup(new WaveFunctionGroup(bg,Spin(Spin::Down)))
    , itsEC           (ec) //Electron cofiguration
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDnGroup);
    assert(itsEC);
};

MasterPolarizedWF::~MasterPolarizedWF()
{
    delete itsSpinUpGroup;
    delete itsSpinDnGroup;
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
//  to fill up the orbitals with electrons.
//
void MasterPolarizedWF::DoSCFIteration(Hamiltonian& ham)
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDnGroup);
    itsSpinUpGroup  ->DoSCFIteration(ham);
    itsSpinDnGroup->DoSCFIteration(ham);
}

ChargeDensity* MasterPolarizedWF::GetChargeDensity(Spin s) const
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDnGroup);
    assert(s==Spin::None);
    ChargeDensity* up=itsSpinUpGroup->GetChargeDensity(Spin::Up);
    ChargeDensity* dn=itsSpinDnGroup->GetChargeDensity(Spin::Down);
    return new PolarizedCDImp(up,dn);
}


const EnergyLevels& MasterPolarizedWF::FillOrbitals(const ElectronConfiguration*)
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDnGroup);
    itsUpELevels=itsSpinUpGroup->FillOrbitals(itsEC);
    itsDnELevels=itsSpinDnGroup->FillOrbitals(itsEC);
    return itsUpELevels;
}


SCFIterator* MasterPolarizedWF::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons)
{
    return new SCFIteratorPol(this,H,cd,nElectrons,0.0);
}

WaveFunction* MasterPolarizedWF::GetWaveFunction(const Spin& S)
{
    assert(S.itsState!=Spin::None);
    WaveFunction* ret=0;
    if (S.itsState==Spin::Up  ) ret=itsSpinUpGroup  ;
    if (S.itsState==Spin::Down) ret=itsSpinDnGroup;
    return ret;
}

void MasterPolarizedWF::DisplayEigen() const
{
    std::cout << "  Spin:     up                 down" << std::endl;
    auto iup=itsUpELevels.begin();
    auto idn=itsDnELevels.begin();
    while (iup!=itsUpELevels.end() && idn!=itsDnELevels.end())
    {
        if (iup!=itsUpELevels.end() && iup->first<=0.0)
            iup->second.Report(std::cout);
        else
            std::cout << "                                     ";
        
        if (idn!=itsDnELevels.end() && idn->first<=0.0)
            idn->second.Report(std::cout);
        
        
            
        std::cout << std::endl;
        if (iup!=itsUpELevels.end()) iup++;
        if (idn!=itsDnELevels.end()) idn++;
        if (iup->first>0.0 && idn->first>0.0) break;
    }
   
}

std::ostream& MasterPolarizedWF::Write(std::ostream& os) const
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDnGroup);
    if (Pretty())
        os << "Polarized wave function :" << std::endl;
    os << *itsSpinUpGroup << *itsSpinDnGroup;
    
    return os;
}

std::istream& MasterPolarizedWF::Read (std::istream& is)
{
    delete itsSpinUpGroup;
    itsSpinUpGroup=WaveFunction::Factory(is);
    assert(itsSpinUpGroup  );
    is >> *itsSpinUpGroup;

    delete itsSpinDnGroup;
    itsSpinDnGroup=WaveFunction::Factory(is);
    assert(itsSpinDnGroup  );
    is >> *itsSpinDnGroup;

    return is;
}


