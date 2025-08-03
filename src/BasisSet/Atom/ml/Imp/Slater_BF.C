// File: Atom/ml/Slater_BF.C  r^l exp(-ar)*Y_lm type basis function 
module;
#include <iostream>
#include <cassert>
module qchem.BasisSet.Atom.Internal.ml.SlaterBS;
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


}} //namespace
