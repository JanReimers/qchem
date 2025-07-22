// File: GaussianScaler.H  Rescale Gaussian exponents based in angular momentum L.
module;
export module qchem.BasisSet.Atom.radial.Gaussian.ExponentScaler; 
import oml.Vector;
export import qchem.Types;

export namespace Gaussian
{

class ExponentScaler
{
public:
    ExponentScaler(size_t N, double emin, double emax, size_t LMax);
    RVec   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    RVec es;
};

} //namespace

