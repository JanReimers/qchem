// File: GaussianScaler.H  Rescale Gaussian exponents based in angular momentum L.
module;
#include <valarray>
export module qchem.BasisSet.Atom.Internal.radial.Gaussian.ExponentScaler; 

export namespace Gaussian
{

class ExponentScaler
{
public:
    using ds_t=std::valarray<double>;
    ExponentScaler(size_t N, double emin, double emax, size_t LMax);
    ds_t   Get_es (size_t L) const;
private:        
    size_t itsN,itsLMax;
    double itsemin,itsemax;
    ds_t es;
};

} //namespace

