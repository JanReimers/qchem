
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

void IEClient::Append(const IrrepIEClient* ic)
{
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ls.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ls(j)=ic->Ls(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
    }

}

} //namespace
