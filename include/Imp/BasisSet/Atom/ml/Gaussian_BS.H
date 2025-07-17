// File: Atom/ml/Gaussian_BS.H  r^l exp(-ar^2)*Y_lm type basis set.
#ifndef _Atom_ml_Gaussian_BS_H_
#define _Atom_ml_Gaussian_BS_H_

#include <vector>
#include "IE_HF_Angular.H"
#include "../../BS_Common.H"
#include "../IE.H"
#include "../radial/Gaussian/BS_Common.H"

namespace Atom_ml
{
namespace Gaussian
{

class BasisSet 
: public ::Gaussian::BS_Common
, public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
};

} //namespace 
} //namespace 

#endif // _Atom_ml_Gaussian_BS_H_
