// File: BasisFunction.C  Polarized Gaussian in 3D space.


#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

#include "PolarizedGaussian/BasisFunction.H"

namespace PolarizedGaussian
{

//#######################################################################
//
//  Polarized gaussian implementation
//

BasisFunction::BasisFunction()
    : itsRadial       (0)
    , itsNormalization(0)
{};

BasisFunction::BasisFunction(const RadialFunction* theRF,const Polarization& thePol, double norm)
    : itsRadial(theRF )
    , itsPol   (thePol)
    , itsNormalization(norm)
{
    assert(itsRadial);
};

bool BasisFunction::operator==(const ::Real_BF& bf) const
{
    assert(itsRadial);
    const BasisFunction& pgbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&pgbf);
    return (*itsRadial)==(*pgbf.itsRadial) && (itsPol==pgbf.itsPol);
}

//------------------------------------------------------------------------
//
//  These basis functions are really just proxys, they dont own the radial.
//  So no information needs be pickled.
//
std::ostream& BasisFunction::Write(std::ostream& os) const
{
    assert(itsRadial);
    
    os.setf(std::ios::fixed,std::ios::floatfield);
    os <<  std::setw(5) << std::setprecision(2) << itsRadial->GetCenter() << " "
    << itsPol << " " << *itsRadial;
    return os;
}

std::istream& BasisFunction::Read(std::istream& is)
{
    return is;
}

void BasisFunction::Insert(const RadialFunction* theRF,const Polarization& thePol)
{
    itsRadial       =theRF;
    itsPol          =thePol;
    assert(itsRadial);
}

double BasisFunction::operator()(const RVec3& r) const
{
    assert(itsRadial);
    const RadialFunction& rf=*itsRadial;
    RVec3 dr=r-rf.GetCenter();
    return itsNormalization*itsPol(dr) * rf(r);
}

RVec3 BasisFunction::Gradient(const RVec3& r) const
{
    assert(itsRadial);
    const RadialFunction& rf=*itsRadial;
    RVec3 dr=r-rf.GetCenter();
    return itsNormalization*(itsPol.Gradient(dr) * rf(r) + itsPol(dr) * rf.Gradient(r));
}

BasisFunction* BasisFunction::Clone() const
{
    assert(itsRadial);
    return new  BasisFunction(*this);
}

} //namespace PolarizedGaussian
