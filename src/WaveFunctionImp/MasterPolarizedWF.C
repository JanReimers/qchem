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
    , itsSpinDownGroup(0)
    , itsEC      (0)
{};

MasterPolarizedWF::MasterPolarizedWF(const BasisSet* bg,const ElectronConfiguration* ec)
    : itsSpinUpGroup  (new WaveFunctionGroup(bg,Spin(Spin::Up  )))
    , itsSpinDownGroup(new WaveFunctionGroup(bg,Spin(Spin::Down)))
    , itsEC           (ec) //Electron cofiguration
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    assert(itsEC);
};

MasterPolarizedWF::~MasterPolarizedWF()
{
    delete itsSpinUpGroup;
    delete itsSpinDownGroup;
}
//----------------------------------------------------------------------------
//
//  This function will creat EMPTY orbtials.  One must use the FillOrbitals member function
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
    ChargeDensity* up=itsSpinUpGroup->GetChargeDensity(Spin::Up);
    ChargeDensity* dn=itsSpinDownGroup->GetChargeDensity(Spin::Down);
    return new PolarizedCDImp(up,dn);
}


void MasterPolarizedWF::FillOrbitals(const ElectronConfiguration*)
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    itsSpinUpGroup  ->FillOrbitals(itsEC);
    itsSpinDownGroup->FillOrbitals(itsEC);
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
    if (S.itsState==Spin::Down) ret=itsSpinDownGroup;
    return ret;
}

void MasterPolarizedWF::DisplayEigen() const
{
    std::cout << "Alpha spin :" << std::endl;
    itsSpinUpGroup->DisplayEigen();
    std::cout << "Beta spin :" << std::endl;
    itsSpinDownGroup->DisplayEigen();
}

std::ostream& MasterPolarizedWF::Write(std::ostream& os) const
{
    assert(itsSpinUpGroup  );
    assert(itsSpinDownGroup);
    if (Pretty())
        os << "Polarized wave function :" << std::endl;
    os << *itsSpinUpGroup << *itsSpinDownGroup;
    
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

    return is;
}


