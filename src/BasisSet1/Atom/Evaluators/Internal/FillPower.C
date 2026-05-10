// File: BasisSet/Atom/Internal/FillPower.C  power series to create tempered basis set exponents.
module;
#include <cmath>
#include <cassert>
#include <blaze/math/DynamicVector.h>

export module  qchem.BasisSet.Atom.Internal.FillPower;
import qchem.Types;

export template <class T> void FillPower(vec_t<T>& arr,T start, T stop)
{
  size_t N=arr.size();
  assert(N>0);
  double beta=N>1 ? std::pow(stop/start,1.0/(N-1)) : 1.0;
  for (auto& a:arr) 
  {
    a=start;
    start*=beta;
  }
}
