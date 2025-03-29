// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Slater_BS_H_
#define _Atom_kappa_Slater_BS_H_

#include "../../BS_Common.H"
#include "Imp/DataBase/HeapDB.H"

namespace Atom_kappa
{
namespace Slater
{


class BasisSet 
    : public ::BS_Common
    , public DB_cache<double>
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lMax);
    
};

}} //namespace Atom_kappa::Slater

#endif // _Atom_kappa_Slater_BS_H_
