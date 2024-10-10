// File: RadialFunctionImplementation.C  Partial implementation for the radial part of a basis function.



#include "Imp/BasisSet/PolarizedGaussian/Radial/Common.H"
#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite1.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include <cmath>
#include <iostream>
#include <cassert>

//#######################################################################
//
//   gaussian implementation
//

RadialFunctionImplementation::RadialFunctionImplementation()
    : itsCenter(0,0,0)
    , itsL     (0)
    , itsH1    (0)
{};

RadialFunctionImplementation::RadialFunctionImplementation(const RVec3& theCenter, int theL)
    : itsCenter(theCenter)
    , itsL     (theL     )
    , itsH1    (0)
{};

RadialFunctionImplementation::RadialFunctionImplementation(const RadialFunctionImplementation& rfi)
    : itsCenter(rfi.itsCenter)
    , itsL     (rfi.itsL)
    , itsH1    (0)
{};

RadialFunctionImplementation::~RadialFunctionImplementation()
{
    delete itsH1;
}

bool RadialFunctionImplementation::operator==(const RadialFunction& other) const
{
    bool ret=false;
    const RadialFunctionImplementation* rfi = dynamic_cast<const RadialFunctionImplementation*>(&other);
    assert(rfi);
    if (rfi)
    {
        ret = norm(itsCenter-rfi->itsCenter) < 0.01; //0.01 Atom units.
//    ret= distance && (itsL==rfi->itsL);
    }
    return ret;
}

const Hermite1& RadialFunctionImplementation::GetH1() const
{
//    Hermite1* temp=itsH1; //This convoluted due to gcc-2.7.2 optimization BUG for mutables.
    if(itsH1==0) itsH1=MakeH1();
    return *itsH1;
}

std::ostream& RadialFunctionImplementation::Write(std::ostream& os) const
{
    UniqueID::Write(os);
//  ScalarFunctionBuffer::Write(os);
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

std::istream& RadialFunctionImplementation::Read(std::istream& is)
{
    UniqueID::Read(is);
//  ScalarFunctionBuffer::Read(is);
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

