
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

const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}

Vector<double>  IEClient::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
    const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double>  IEClient::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
    const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
    return cd->ExchangeRk(la,lc);
}

} //namespace
