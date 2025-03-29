// File: Atom/ml/Gaussian_BS.H  r^l exp(-ar^2)*Y_lm type basis set.
#ifndef _Atom_ml_Gaussian_BS_H_
#define _Atom_ml_Gaussian_BS_H_

#include "../../BS_Common.H"
#include "../IE.H"
#include "Imp/BasisSet/Atom/l/Gaussian_BS.H"

namespace PolarizedGaussian
{
    class Reader;    
}

namespace Atom_ml
{
namespace Gaussian
{

class BasisSet 
: public ::Gaussian::BS_Common
, public ::AtomIE_BS_2E_lm //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
    BasisSet(PolarizedGaussian::Reader*, const Atom*);
};

} //namespace 
} //namespace 

#endif // _Atom_ml_Gaussian_BS_H_
