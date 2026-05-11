// File: src/BasisSet1/Atom/Evaluators/Slater/Internal/Imp/ExponentScaler.C  Rescale Slater exponents based in angular momentum L.
module;
#include <cmath>
#include <blaze/math/Subvector.h>
module qchem.BasisSet1.Atom.Evaluators.Slater.Internal.ExponentScaler;
import qchem.BasisSet1.Atom.Internal.FillPower;
import qchem.Symmetry.Yl;

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

rvec_t ExponentScaler::Get_es(size_t L) const
{
    if (L==0) return es;
    int N=itsN-1*L;
    if (N<1) N=1;
    if (N+L>itsN) L=itsN-N;
    return blaze::subvector(es,L,N);
}

rvec_t ExponentScaler::Get_es (const Irrep_QNs::sym_t& ir) const
{
    const Yl_Sym* yl=dynamic_cast<const Yl_Sym*>(ir.get());
    assert(yl);
    return Get_es(yl->GetL());
}
} //namespace