// File: src/BasisSet/Atom/Evaluators/Slater/Internal/Imp/ExponentScaler.C  Rescale Slater exponents based in angular momentum L.
module;
module qchem.BasisSet.Atom.Evaluators.Slater.Internal.ExponentScaler;
import qchem.BasisSet.Atom.Internal.FillPower;
import qchem.Symmetry.Spherical;
import qchem.Math;
import qchem.Blaze;

namespace qchem::Slater
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
    return blazem::subvector(es,L,N);
}

rvec_t ExponentScaler::Get_es (const sym_t& ir) const
{
    return Get_es(Getl(ir));
}
} //namespace