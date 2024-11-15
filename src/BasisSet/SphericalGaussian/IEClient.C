
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

double IrrepIEClient::Norm(double e, double l)
{
    return GaussianNorm(e,l);
}
    
//void IrrepIEClient::Init(double minexp,double maxexp,size_t L)
//{
//    
//      FillPower(es,minexp,maxexp);
//      Fill(Ls,L);
//      for (auto i:es.indices())  ns(i)=GaussianNorm(es(i),L);
//}

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
        BFGrouper::Append(es(j),Ls(j),j);
    }

}

const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}


} //namespace
