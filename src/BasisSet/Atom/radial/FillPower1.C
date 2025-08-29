// File: BasisSet/Atom/radial/FillPower1.C  power series to create tempered basis set exponents.
module;
#include <cmath>
#include <valarray>
#include <cassert>
export module qchem.BasisSet.Atom.Internal.radial.FillPower1;

export template <class T> void FillPower(std::valarray<T>& arr,T start, T stop)
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

