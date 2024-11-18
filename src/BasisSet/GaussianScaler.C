// File: GaussianScaler.C  Rescale Gaussian exponents based in angular momentum L.

#include "Imp/BasisSet/GaussianScaler.H"
template <class T> void FillPower(Vector<T>& arr,T start, T stop);

 GaussianScaler::GaussianScaler(size_t N, double emin, double emax, size_t LMax)
    : itsN(N)
    , itsLMax(LMax)
    , itsemin(emin)
    , itsemax(emax)
    , es(N)
{
    FillPower(es,itsemin,itsemax);
};

GaussianScaler::RVec   GaussianScaler::Get_es(size_t L) const
{
    if (L==0) return es;
    int N=itsN-2*L;
    if (N<1) N=1;
    RVec esL(N);
    for (auto i:esL.indices()) esL(i)=es(i+L);
    return esL;
}
