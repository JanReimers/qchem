// File: Atom/l/Slater_BS.H Slater Basis Set for atoms.
#ifndef _Atom_l_Slater_BS_H_
#define _Atom_l_Slater_BS_H_

#include "../radial/Slater/BS_Common.H"
#include "../../TCommon.H"

namespace Atoml
{
namespace Slater
{

class BasisSet 
    : public ::Slater::BS_Common
    , public ::AtomIE_BS_2E_l //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax); 
};

}} //namespace Slater::l

#endif // _Atom_l_Slater_BS_H_
