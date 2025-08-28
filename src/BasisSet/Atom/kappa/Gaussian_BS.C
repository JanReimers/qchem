// File: Atom/kappa/Gaussian_BS.H  Restricted Kinetic Balance (RKB) Basis Set (BS).
module;
#include <iosfwd>
class DiracIntegralTests;
export module qchem.BasisSet.Atom.Internal.kappa.GaussianBS;
import BasisSet.Atom.Gaussian_IBS;
import qchem.BasisSet.Atom.IBS;

import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;

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

// Irrep Basis set
// All integrals are handled at the Orbital_RKB_IBS_Common.  i.e. they are not Gaussian
// specific.
export class Orbital_RKB_IBS
    : public IrrepBasisSet_Common<double>
    , public Orbital_RKB_IBS_Common<double>
{
public:
    Orbital_RKB_IBS(const DB_cache<double>*, const Vector<double>& exponents, int kappa);
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
public:
    Orbital_RKBL_IBS(const DB_cache<T>*,const Vector<T>& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return Gaussian_IBS::size();}
};

export template <class T> class Orbital_RKBS_IBS
    : public IrrepBasisSet_Common<T> 
    , public Gaussian_RKBS_IBS
    , public Atom::Orbital_RKBS_IBS<T> 
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    using Orbital_RKBS_IBS_Common<T>::large;
public:
    Orbital_RKBS_IBS(const DB_cache<T>*,const Vector<T>& exponents,int kappa);

    virtual size_t  size() const {return Gaussian_RKBS_IBS::size();}
    virtual size_t  GetNumFunctions() const {return Gaussian_RKBS_IBS::size();}
    virtual Vec     operator() (const RVec3& r) const {return Gaussian_RKBS_IBS::operator()(r);}
    virtual Vec3Vec Gradient   (const RVec3& r) const {return Gaussian_RKBS_IBS::Gradient(r);}
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

