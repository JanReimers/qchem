// File: Gaussian/ExponentScaler.C  Rescale Gaussian exponents based in angular momentum L.
module;
#include <valarray>
module qchem.BasisSet.Atom.Internal.radial.Gaussian.ExponentScaler; 
import qchem.BasisSet.Atom.Internal.radial.FillPower;

namespace Gaussian
{


 ExponentScaler::ExponentScaler(size_t N, double emin, double emax, size_t LMax)
    : itsN(N)
    , itsLMax(LMax)
    , itsemin(emin)
    , itsemax(emax)
    , es(N)
{
    ::FillPower(es,itsemin,itsemax);
};

ExponentScaler::ds_t ExponentScaler::Get_es(size_t L) const
{
    if (L==0) return es;
    int N=itsN-4*L;
    if (N<1) N=1;
    return es[std::slice(0,N,1)];
}

} //namespace