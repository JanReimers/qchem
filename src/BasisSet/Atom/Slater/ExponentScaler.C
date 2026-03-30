// File: ExponentScaler.H  Rescale Slater exponents based in angular momentum L.
module;
export module qchem.BasisSet.Atom.Slater.ExponentScaler;
import qchem.Types;

export namespace Slater
{

class ExponentScaler
{
public:
    ExponentScaler(size_t N, double emin, double emax, size_t LMax);
    rvec_t Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    rvec_t es;
};

} //namespace

