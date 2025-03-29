// File: RadialCommon.C  Partial implementation for the radial part of a basis function.



#include "Imp/BasisSet/Molecule/PolarizedGaussian/Radial/Common.H"
#include "Imp/BasisSet/Molecule/PolarizedGaussian/MnD/Hermite1.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include <cmath>
#include <iostream>
#include <cassert>
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
    if ( StreamableObject::Binary())
    {
        os << itsCenter;
        BinaryWrite(itsL      ,os);
    }
    if (StreamableObject::Ascii ())
        os << itsCenter              << " "
        << itsL                   << " ";

    return os;
}

std::istream& RadialCommon::Read(std::istream& is)
{
    UniqueIDImp::Read(is);
    if (StreamableObject::Binary())
    {
        is >> itsCenter;
        BinaryRead(itsL    ,is);
    }
    else
    {
        is >> itsCenter >> itsL;
    }
    itsH1=0; //old one is not longer valid.
    return is;
}

} //namespace PolarizedGaussian
