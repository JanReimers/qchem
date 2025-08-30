// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
module;
#include <iosfwd>
#include <valarray>
class DiracIntegralTests;

export module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import BasisSet.Atom.Slater_IBS;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.HeapDB;

namespace Atom_kappa
{
namespace Slater
{


//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r) as:
//
//                                  -e*r^l*exp(-e*r), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)g = {
//                                  ((2l+1)/r-e)*r^lexp(-e*r), kappa>0
//
export class Orbital_RKB_IBS
    : public IrrepBasisSet_Common<double>
    , public Orbital_RKB_IBS_Common<double> 
{
    using ds_t=std::valarray<double>;
public:
    Orbital_RKB_IBS(const DB_cache<double>* db,const ds_t& exponents, int kappa);
    virtual size_t size() const {return Orbital_RKB_IBS_Common<double>::size();}
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    friend class ::DiracIntegralTests;

};

export template <class T> class Orbital_RKBL_IBS
    : Slater_IBS
    , public Atom::IrrepBasisSet //Use NR slater basis
    , public Atom::Orbital_RKBL_IBS<T>
{
    using ds_t=std::valarray<T>;
public:
    Orbital_RKBL_IBS(const DB_cache<T>*,const ds_t& exponents, int kappa);
};

export template <class T> class Orbital_RKBS_IBS
    : public Slater_RKBS_IBS
    , public Atom::IrrepBasisSet
    , public Atom::Orbital_RKBS_IBS<T> 
{
    using Orbital_RKBS_IBS_Common<T>::large;
    using ds_t=std::valarray<T>;
public:
    Orbital_RKBS_IBS(const DB_cache<double>*,const ds_t& exponents, int kappa);
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    using Orbital_RKBS_IBS_Common<T>::kappa;
};

// Full basis set
export class BasisSet 
    : public ::BS_Common
    , public DB_cache<double>
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lMax);
    
};

}} //namespace Atom_kappa::Slater

