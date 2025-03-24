// File: PolarizedCD.C  Interface for the charge density.

#include "Imp/ChargeDensity/PolarizedCD.H"
#include <cassert>

//---------------------------------------------------------------------------------
//
//  Construction zone.
//
Polarized_CDImp::Polarized_CDImp()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{}; // No UT coverage

Polarized_CDImp::Polarized_CDImp(Exact_CD* up, Exact_CD* down)
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
Exact_CD* Polarized_CDImp::GetChargeDensity(const Spin& S)
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    Exact_CD* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const Exact_CD* Polarized_CDImp::GetChargeDensity(const Spin& S) const
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const Exact_CD* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}




