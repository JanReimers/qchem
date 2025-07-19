// File: PolarizedCD.C  Interface for the charge density.
// module;
#include <cstddef>
#include <cassert>
#include "PolarizedCD.H"
import qchem.ChargeDensity;
import qchem.Symmetry.Spin;

//---------------------------------------------------------------------------------
//
//  Construction zone.
//
Polarized_CDImp::Polarized_CDImp()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{}; // No UT coverage

Polarized_CDImp::Polarized_CDImp(DM_CD* up, DM_CD* down)
    : itsSpinUpCD  (up  )
    , itsSpinDownCD(down)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

Polarized_CDImp::~Polarized_CDImp()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}

//-------------------------------------------------------------------
//
//  Access to individual components.
//
DM_CD* Polarized_CDImp::GetChargeDensity(const Spin& s)
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    DM_CD* ret=0;
    if(s==Spin::Up  ) ret=itsSpinUpCD  ;
    if(s==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const DM_CD* Polarized_CDImp::GetChargeDensity(const Spin& s) const
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const DM_CD* ret=0;
    if(s==Spin::Up  ) ret=itsSpinUpCD  ;
    if(s==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

