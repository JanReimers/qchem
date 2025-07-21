// File: BasisSet/Atom/radial/FillPower.C  power series to create tempered basis set exponents.
module;
#include <cmath>
export module qchem.BasisSet.Atom.radial.FillPower;
export import oml.Vector;

export template <class T> void FillPower(Vector<T>& arr,T start, T stop)
{
  double del=0.5*(start+stop); //n=1 case
  if (arr.size()>1)
    del=(std::log(stop/start))/(double)(arr.size()-1);
  typename Vector<T>::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i=T(start*std::exp(n*del));
}

