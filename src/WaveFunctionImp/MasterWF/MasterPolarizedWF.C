// File: MasterPolarizedWF.C  Wave function for an unpolarized atom.



#include "WaveFunctionImp/MasterWF/MasterPolarizedWF.H"
#include "WaveFunctionImp/MasterWF/PolarizedSCFIterator.H"
#include "WaveFunctionImp/WaveFunctionGroup/WaveFunctionGroup.H"
#include "ChargeDensity/PolarizedCD.H"
#include "ChargeDensityImplementation/PolarizedCDImplementation.H"
#include "Misc/Spin.H"
#include "oml/imp/binio.h"
#include <cassert>

MasterPolarizedWF::MasterPolarizedWF()
    : itsSpinUpGroup  (0)
    , itsSpinDownGroup(0)
    , itsNetSpin      (0)
{};

MasterPolarizedWF::MasterPolarizedWF(const BasisGroup* bg,double spin)
    : itsSpinUpGroup  (new WaveFunctionGroup(bg,Spin(Spin::Up  )))
    , itsSpinDownGroup(new WaveFunctionGroup(bg,Spin(Spin::Down)))
    , itsNetSpin      (spin)
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
};

MasterPolarizedWF::~MasterPolarizedWF()
{
    delete itsSpinUpGroup;
    delete itsSpinDownGroup;
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the ElectronDumper
//  to fill up the orbitals with electrons.
//
void MasterPolarizedWF::DoSCFIteration(Hamiltonian& ham)
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    itsSpinUpGroup  ->DoSCFIteration(ham);
    itsSpinDownGroup->DoSCFIteration(ham);
}

ChargeDensity* MasterPolarizedWF::GetChargeDensity(Spin s) const
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    assert(s==Spin::None);
    return new PolarizedCDImplementation(itsSpinUpGroup->GetChargeDensity(Spin::Up),itsSpinDownGroup->GetChargeDensity(Spin::Down));
}

void MasterPolarizedWF::UpdateElectronDumper(ElectronDumper& ed)
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    itsSpinUpGroup  ->UpdateElectronDumper(ed);
    itsSpinDownGroup->UpdateElectronDumper(ed);
}


SCFIterator* MasterPolarizedWF::MakeIterator(Hamiltonian* H, ChargeDensity* cd, double nElectrons, double kT, bool showplot)
{
    return new PolarizedSCFIterator(this,H,cd,nElectrons,itsNetSpin,kT,showplot);
}

WaveFunction* MasterPolarizedWF::GetWaveFunction(const Spin& S)
{
    assert(S.itsState!=Spin::None);
    WaveFunction* ret=0;
    if (S.itsState==Spin::Up  ) ret=itsSpinUpGroup  ;
    if (S.itsState==Spin::Down) ret=itsSpinDownGroup;
    return ret;
}

std::ostream& MasterPolarizedWF::Write(std::ostream& os) const
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    if (Pretty())
        os << "Polarized wave function with net spin=" << itsNetSpin << ":" << std::endl;
    os << *itsSpinUpGroup << *itsSpinDownGroup;
    
    if (Binary()) BinaryWrite(itsNetSpin,os);
    if (Ascii ()) os << itsNetSpin << " ";
    return os;
}

std::istream& MasterPolarizedWF::Read (std::istream& is)
{
    delete itsSpinUpGroup;
    itsSpinUpGroup=WaveFunction::Factory(is);
    assert(itsSpinUpGroup  );
    is >> *itsSpinUpGroup;

    delete itsSpinDownGroup;
    itsSpinDownGroup=WaveFunction::Factory(is);
    assert(itsSpinDownGroup  );
    is >> *itsSpinDownGroup;

    if (Binary())
        BinaryRead(itsNetSpin,is);
    else
        is >> itsNetSpin;

    return is;
}


