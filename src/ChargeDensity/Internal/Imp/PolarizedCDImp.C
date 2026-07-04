// File: PolarizedCD.C  Interface for the charge density.
module;
#include <cassert>
module qchem.ChargeDensity.Imp.PolarizedCD;
import qchem.ChargeDensity;
import qchem.Symmetry.Spin;

namespace qchem::ChargeDensity
{

//---------------------------------------------------------------------------------
//
//  Construction zone.
//
Polarized_CDImp::Polarized_CDImp()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{}; // No UT coverage

Polarized_CDImp::Polarized_CDImp(rDM_CD* up, rDM_CD* down)
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
rDM_CD* Polarized_CDImp::GetChargeDensity(const Spin& s)
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    rDM_CD* ret=0;
    if(s==Spin::Up  ) ret=itsSpinUpCD  ;
    if(s==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const rDM_CD* Polarized_CDImp::GetChargeDensity(const Spin& s) const
{
    assert(s!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const rDM_CD* ret=0;
    if(s==Spin::Up  ) ret=itsSpinUpCD  ;
    if(s==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

} //namespace
