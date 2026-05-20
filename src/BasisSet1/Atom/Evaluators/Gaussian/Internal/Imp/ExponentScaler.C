// File: Gaussian/ExponentScaler.C  Rescale Gaussian exponents based in angular momentum L.
module;
#include <blaze/math/views/Subvector.h>
module qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.ExponentScaler; 
import qchem.BasisSet.Atom.Internal.FillPower;
import qchem.Symmetry.Yl;

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

rvec_t ExponentScaler::Get_es(size_t L) const
{
    if (L==0) return es;
    int N=itsN-4*L;
    if (N<1) N=1;
    return blaze::subvector(es,0,N);
}
rvec_t ExponentScaler::Get_es (const Irrep_QNs::sym_t& ir) const
{
    const Yl_Sym* yl=dynamic_cast<const Yl_Sym*>(ir.get());
    assert(yl);
    return Get_es(yl->GetL());
}

} //namespace