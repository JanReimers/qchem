// File: Atom/l/BSpline_BS.H BSpline Basis Set for atoms.
#ifndef _Atom_l_BSpline_BS_H_
#define _Atom_l_BSpline_BS_H_

#include "../radial/BSpline/BS_Common.H"
#include "IE_HF_Angular.H"
#include "../../BS_Common.H"

namespace Atoml
{
namespace BSpline
{

template <size_t K> class BasisSet 
    : public ::BSpline::BS_Common<K>
    , public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet(size_t N, double rmin, double rmax, size_t Lmax); 
};

}} //namespace Atoml::BSpline

#endif // _Atom_l_BSpline_BS_H_
