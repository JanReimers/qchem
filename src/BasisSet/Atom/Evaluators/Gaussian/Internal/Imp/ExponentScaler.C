// File: Gaussian/ExponentScaler.C  Rescale Gaussian exponents based in angular momentum L.
module;
module qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.ExponentScaler; 
import qchem.BasisSet.Atom.Internal.FillPower;
import qchem.Symmetry.Atom.Spherical;
import qchem.Blaze;

namespace qchem::Gaussian
{


 ExponentScaler::ExponentScaler(size_t N, double emin, double emax, size_t LMax)
    : itsN(N)
    , itsLMax(LMax)
    , itsemin(emin)
    , itsemax(emax)
    , es(N)
{
    qchem::FillPower(es,itsemin,itsemax);
};

rvec_t ExponentScaler::Get_es(size_t L) const
{
    if (L==0) return es;
    int N=itsN-4*L;
    if (N<1) N=1;
    return blazem::subvector(es,0,N);
}
rvec_t ExponentScaler::Get_es (const sym_t& ir) const
{
    return Get_es(Symmetry::Atom::Getl(ir));
}

} //namespace