// File: Atom/l/Gaussian_BS.H Gaussian Basis Set for atoms.
#ifndef _Atom_l_Gaussian_BS_H_
#define _Atom_l_Gaussian_BS_H_

#include "IE_HF_Angular.H"
#include "../radial/Gaussian/BS_Common.H"
#include "../../BS_Common.H"

namespace Atoml
{
namespace Gaussian
{


class BasisSet 
: public ::Gaussian::BS_Common
, public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
};

} //namespace
} //namespace

#endif // _Atom_l_Gaussian_BS_H_
