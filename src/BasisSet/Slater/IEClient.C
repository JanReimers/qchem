
#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/Integrals/SlaterIntegrals.H"

template <class T> void FillPower(Vector<T>& arr,T start, T stop)
{
  double del=0.5*(start+stop); //n=1 case
  if (arr.size()>1)
    del=(std::log(stop/start))/(double)(arr.size()-1);
  typename Vector<T>::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i=T(start*std::exp(n*del));
}

template void FillPower(Vector<double>& arr,double start, double stop);

namespace Slater
{
    
void IrrepIEClient::Init(double minexp,double maxexp,size_t L)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      for (auto i:es.indices())  ns(i)=SlaterNorm(es(i),Ns(i));
}

void IEClient::Append(const IrrepIEClient* ic)
{
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ns.SetLimits(N,true);
    Ls.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ns(j)=ic->Ns(i);
        Ls(j)=ic->Ls(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
        BFGrouper::Append(es(j),Ls(j),j);
    }

}

const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}


} //namespace
