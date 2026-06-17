// File: BasisSet/Atom/Evaluators/Internal/FillPower.C  power series to create tempered basis set exponents.
module;
#include <cassert>

export module  qchem.BasisSet.Atom.Internal.FillPower;
import qchem.Math;
import qchem.Blaze;

export template <class T> void FillPower(vec_t<T>& arr,T start, T stop)
{
  size_t N=arr.size();
  assert(N>0);
  double beta=N>1 ? pow(stop/start,1.0/(N-1)) : 1.0;
  for (auto& a:arr) 
  {
    a=start;
    start*=beta;
  }
}
