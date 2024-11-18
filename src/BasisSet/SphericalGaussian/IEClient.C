
#include "Imp/BasisSet/SphericalGaussian/IEClient.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"

template <class T> inline void FillPower(Vector<T>& arr,T start, T stop)
{
  double del=(std::log(stop/start))/(double)(arr.size()-1);
  typename Vector<T>::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i=T(start*std::exp(n*del));
}

namespace SphericalGaussian
{

double IrrepIEClient::Norm(double e, size_t l) const
{
    return GaussianNorm(e,l);
}


} //namespace
