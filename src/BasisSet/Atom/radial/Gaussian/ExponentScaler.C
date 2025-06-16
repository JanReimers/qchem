// File: Gaussian/ExponentScaler.C  Rescale Gaussian exponents based in angular momentum L.

#include "radial/Gaussian/ExponentScaler.H"

template <class T> void FillPower(Vector<T>& arr,T start, T stop);

namespace Gaussian
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
    int N=itsN-4*L;
    if (N<1) N=1;
    RVec esL(N);
    for (auto i:esL.indices()) esL(i)=es(i);
    return esL;
}

} //namespace