
#include "Imp/BasisSet/Slater_m/IEClient.H"
#include "Imp/Integrals/SlaterIntegrals.H"

template <class T> void FillPower(Vector<T>& arr,T start, T stop);
//template void FillPower(Vector<double>& arr,double start, double stop);

namespace Slater_m
{
    
void IrrepIEClient::Init(double minexp,double maxexp,size_t L, int m)
{
    
      FillPower(es,minexp,maxexp);
      Fill(Ns,L+1);
      Fill(Ls,L);
      Fill(Ms,m);
      for (auto i:es.indices())  ns(i)=SlaterNorm(es(i),Ns(i));
}

void IEClient::Append(const IrrepIEClient* ic)
{
    size_t j=size()+1;
    size_t N=size()+ic->size();
    Ns.SetLimits(N,true);
    Ls.SetLimits(N,true);
    Ms.SetLimits(N,true);
    es.SetLimits(N,true);
    ns.SetLimits(N,true);
    for (size_t i=1;i<=ic->size();i++,j++)
    {
        Ns(j)=ic->Ns(i);
        Ls(j)=ic->Ls(i);
        Ms(j)=ic->Ms(i);
        es(j)=ic->es(i);
        ns(j)=ic->ns(i);
    }

}

} //namespace
