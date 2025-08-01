// File: ExponentScaler.C  Rescale Slater exponents based in angular momentum L.
module;
#include <cmath>
module qchem.BasisSet.Atom.Internal.radial.Slater.ExponentScaler;
import qchem.BasisSet.Atom.Internal.radial.FillPower;

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

RVec   ExponentScaler::Get_es(size_t L) const
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