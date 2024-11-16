
#include "Imp/BasisSet/SphericalGaussian/IEClient.H"
#include "Imp/Integrals/GaussianIntegrals.H"

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

const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}


} //namespace
