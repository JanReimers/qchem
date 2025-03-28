// File: Atom/ml/Slater_BF.C  r^l exp(-ar)*Y_lm type basis function 

#include "Imp/BasisSet/Atom/ml/Slater_BF.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
#include <iostream>
#include <cassert>

namespace Atom_ml
{
namespace Slater
{


BasisFunction::BasisFunction(double e, int n, int l, int _ml, double norm)
: Base(e,n,l,norm)
, ml(_ml)
{};

bool BasisFunction::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction& sgbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&sgbf);
    return  Base::operator==(bf) && (ml==sgbf.ml);
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

}} //namespace
