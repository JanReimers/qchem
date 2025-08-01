// File: RadialCommon.C  Partial implementation for the radial part of a basis function.
module;
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <iomanip>
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.Common;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite1;

namespace PolarizedGaussian
{
//#######################################################################
//
//   gaussian implementation
//

RadialCommon::RadialCommon()
    : itsCenter(0,0,0)
    , itsL     (0)
    , itsH1    (0)
{};

RadialCommon::RadialCommon(const RVec3& theCenter, int theL)
    : itsCenter(theCenter)
    , itsL     (theL     )
    , itsH1    (0)
{};

RadialCommon::RadialCommon(const RadialCommon& rfi)
    : itsCenter(rfi.itsCenter)
    , itsL     (rfi.itsL)
    , itsH1    (0)
{};

RadialCommon::~RadialCommon()
{
    delete itsH1;
}

bool RadialCommon::operator==(const RadialFunction& other) const
{
    bool ret=false;
    const RadialCommon* rfi = dynamic_cast<const RadialCommon*>(&other);
    assert(rfi);
    if (rfi)
    {
        ret = norm(itsCenter-rfi->itsCenter) < 0.01; //0.01 Atom units.
//    ret= distance && (itsL==rfi->itsL);
    }
    return ret;
}

const Hermite1& RadialCommon::GetH1() const
{
//    Hermite1* temp=itsH1; //This convoluted due to gcc-2.7.2 optimization BUG for mutables.
    if(itsH1==0) itsH1=MakeH1();
    return *itsH1;
}

std::ostream& RadialCommon::Write(std::ostream& os) const
{
    UniqueIDImp::Write(os);
    os << itsCenter << " " << itsL << " ";

    return os;
}


} //namespace PolarizedGaussian
