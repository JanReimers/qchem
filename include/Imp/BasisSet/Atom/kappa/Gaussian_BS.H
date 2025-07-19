// File: Atom/kappa/Gaussian_BS.H  Restricted Kinetic Balance (RKB) Basis Set (BS).
#ifndef _Atom_kappa_Gaussian_BS_H_
#define _Atom_kappa_Gaussian_BS_H_

import qchem.BasisSet.Imp.HeapDB;
import qchem.BasisSet.Common;

namespace Atom_kappa
{
namespace Gaussian
{

class BasisSet 
    : public ::BS_Common
    , public DB_cache<double>
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lmax);
    
};

}} //namespace

#endif // _Atom_kappa_Gaussian_BS_H_
