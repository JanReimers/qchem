// File: ExponentScaler.C  Rescale Slater exponents based in angular momentum L.

#include "radial/Slater/ExponentScaler.H"


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



ExponentScaler::ExponentScaler(size_t N, double emin, double emax, size_t LMax)
    : itsN(N)
    , itsLMax(LMax)
    , itsemin(emin)
    , itsemax(emax)
    , es(N)
{
    FillPower(es,itsemin,itsemax);
};

ExponentScaler::RVec   ExponentScaler::Get_es(size_t L) const
{
    if (L==0) return es;
    int N=itsN-1*L;
    if (N<1) N=1;
    if (N+L>itsN) L=itsN-N;
    RVec esL(N);
    for (auto i:esL.indices()) esL(i)=es(i+L);
    return esL;
}

} //namespace