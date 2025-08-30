// File: Atom/kappa/Gaussian_BS.H  Restricted Kinetic Balance (RKB) Basis Set (BS).
module;
#include <iosfwd>
#include <valarray>
class DiracIntegralTests;
export module qchem.BasisSet.Atom.Internal.kappa.GaussianBS;
import BasisSet.Atom.Gaussian_IBS;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Internal.HeapDB;

namespace Atom_kappa
{
namespace Gaussian
{ 


//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r^2) as:
//
//                                   { -2e*r^(l+1)*exp(-e*r^2), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)P(r) = {
//                                   { ((2l+1)/r-2er)*r^l*exp(-e*r^2), kappa>0
//
export class Orbital_RKB_IBS
    : public IrrepBasisSet_Common<double>
    , public Orbital_RKB_IBS_Common<double>
{
    using ds_t=std::valarray<double>;
public:
    Orbital_RKB_IBS(const DB_cache<double>*, const ds_t& exponents, int kappa);
    virtual size_t size() const {return Orbital_RKB_IBS_Common<double>::size();}
    virtual std::ostream&  Write(std::ostream&    ) const;

private:
    friend class ::DiracIntegralTests;

};

export template <class T> class Orbital_RKBL_IBS
    : public Gaussian_IBS
    , public Atom::IrrepBasisSet //Use NR Gaussian basis
    , public Atom::Orbital_RKBL_IBS<T>
{
    using ds_t=std::valarray<T>;
public:
    Orbital_RKBL_IBS(const DB_cache<T>*,const ds_t& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return Gaussian_IBS::size();}
};

export template <class T> class Orbital_RKBS_IBS
    : public Gaussian_RKBS_IBS
    , public Atom::IrrepBasisSet //Use NR Gaussian basis
    , public Atom::Orbital_RKBS_IBS<T> 
{
    using Orbital_RKBS_IBS_Common<T>::large;
    using ds_t=std::valarray<T>;
public:
    Orbital_RKBS_IBS(const DB_cache<T>*,const ds_t& exponents,int kappa);
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
    BasisSet(size_t N, double minexp, double maxexp, size_t lmax);
    
};

}} //namespace

