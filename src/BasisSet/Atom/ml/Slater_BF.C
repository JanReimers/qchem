// File: Atom/ml/Slater_BF.C  r^l exp(-ar)*Y_lm type basis function 

#include <iostream>
#include <cassert>
#include "ml/Slater_BF.H"

import Common.IntPow;
import oml;

namespace Atom_ml
{
namespace Slater
{


BasisFunction::BasisFunction(double e, int n, int l, int _ml, double norm)
: Base(e,n,l,norm)
, ml(_ml)
{};


BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

}} //namespace
