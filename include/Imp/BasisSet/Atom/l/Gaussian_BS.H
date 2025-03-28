// File: Atom/l/Gaussian_BS.H Gaussian Basis Set for atoms.
#ifndef _Atom_l_Gaussian_BS_H_
#define _Atom_l_Gaussian_BS_H_

#include "../radial/Gaussian/BS_Common.H"
#include "../../TCommon.H"

namespace Atoml
{
namespace Gaussian
{


class BasisSet 
: public ::Gaussian::BS_Common
, public ::AtomIE_BS_2E_l //Pick angular integrals.
{
protected:
    typedef ::AtomIE_BS_2E_Angular::RVec RVec;
    typedef ::BasisSet Base;
    typedef Base::bs_t bs_t;
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
};

} //namespace
} //namespace

#endif // _Atom_l_Gaussian_BS_H_
