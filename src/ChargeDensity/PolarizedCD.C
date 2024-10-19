// File: PolarizedCD.C  Interface for the charge density.

#include "Imp/ChargeDensity/PolarizedCD.H"
#include <cassert>

//---------------------------------------------------------------------------------
//
//  Construction zone.
//
PolarizedCDImp::PolarizedCDImp()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{};

PolarizedCDImp::PolarizedCDImp(ChargeDensity* up, ChargeDensity* down)
    : itsSpinUpCD  (up  )
    , itsSpinDownCD(down)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

PolarizedCDImp::~PolarizedCDImp()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}

//-------------------------------------------------------------------
//
//  Access to individual components.
//
ChargeDensity* PolarizedCDImp::GetChargeDensity(const Spin& S)
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    ChargeDensity* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const ChargeDensity* PolarizedCDImp::GetChargeDensity(const Spin& S) const
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const ChargeDensity* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}



//--------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& PolarizedCDImp::Write(std::ostream& os) const
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    return os << *itsSpinUpCD << *itsSpinDownCD;
}

std::istream& PolarizedCDImp::Read (std::istream& is)
{
    delete itsSpinUpCD;
    itsSpinUpCD=ChargeDensity::Factory(is);
    assert(itsSpinUpCD);
    is >> *itsSpinUpCD;

    delete itsSpinDownCD;
    itsSpinDownCD=ChargeDensity::Factory(is);
    assert(itsSpinDownCD);
    is >> *itsSpinDownCD;

    return is;
}




