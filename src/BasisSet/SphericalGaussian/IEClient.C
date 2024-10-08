
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
    
void IrrepIEClient::Init(double minexp,double maxexp,size_t L)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ls,L);
      for (auto i:es.indices())  ns(i)=GaussianNorm(es(i),L);
}

} //namespace
