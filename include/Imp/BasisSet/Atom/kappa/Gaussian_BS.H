// File: Atom/kappa/Gaussian_BS.H  Restricted Kinetic Balance (RKB) Basis Set (BS).
#ifndef _Atom_kappa_Gaussian_BS_H_
#define _Atom_kappa_Gaussian_BS_H_

#include "../../TCommon.H"
#include "Imp/DataBase/HeapDB.H"

namespace Atom_kappa
{
namespace Gaussian
{

class BasisSet 
    : public ::BasisSetImp
    , public DB_cache<double>
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lmax);
    
};

}} //namespace

#endif // _Atom_kappa_Gaussian_BS_H_
