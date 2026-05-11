// File: BasisSet1/Atom/Evaluators/Gaussian/Internal/ExponentScaler.C  Rescale Gaussian exponents based in angular momentum L.
module;
export module qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.ExponentScaler; 
import qchem.Types;
import qchem.Symmetry.Irrep;

export namespace Gaussian
{

class ExponentScaler
{
public:
    ExponentScaler(size_t N, double emin, double emax, size_t LMax);
    rvec_t Get_es (size_t L) const;
    rvec_t Get_es (const Irrep_QNs::sym_t& ir) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    rvec_t es;
};

} //namespace

