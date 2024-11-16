
#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/SlaterCD.H"

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
    
double IrrepIEClient::Norm(double e, size_t l) const
{
     return SlaterNorm(e,l+1);  
}



const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}


Vector<double> IEClient::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double> IEClient::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->ExchangeRk(la,lc);
}


} //namespace
