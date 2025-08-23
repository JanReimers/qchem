// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
module;
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import BasisSet.Atom.Slater_IBS;

import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;

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


// Irrep basis set

// All integrals are handled at the Orbital_RKB_IBS_Common.  i.e. they are not Slater function
// specific.
export class Orbital_RKB_IBS
    : public IrrepBasisSet_Common<double>
    , public Orbital_RKB_IBS_Common<double> 
{
public:
    Orbital_RKB_IBS(const DB_cache<double>* db,const Vector<double>& exponents, int kappa);
    virtual size_t size() const {return Orbital_RKB_IBS_Common<double>::size();}
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    friend class ::DiracIntegralTests;

};

export template <class T> class Orbital_RKBL_IBS
    : public ::Slater::IrrepBasisSet //Use NR slater basis
    , Slater_IBS
    , public Atom::Orbital_RKBL_IBS<T>
{
public:
    Orbital_RKBL_IBS(const DB_cache<T>*,const Vector<T>& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return Slater_IBS::size();}
};

export template <class T> class Orbital_RKBS_IBS
    : public IrrepBasisSet_Common<T> 
    , public Slater_RKBS_IBS
    , public Atom::Orbital_RKBS_IBS<T> 
    , public AtomIrrepIEClient
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    using Orbital_RKBS_IBS_Common<T>::large;
public:
    Orbital_RKBS_IBS(const DB_cache<double>*,const Vector<T>& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return Slater_RKBS_IBS::size();}

    // using Slater_RKBS_IBS::operator();
    // using Slater_RKBS_IBS::Gradient;
    virtual Vec     operator() (const RVec3& r) const {return Slater_RKBS_IBS::operator()(r);}
    virtual Vec3Vec Gradient   (const RVec3& r) const {return Slater_RKBS_IBS::Gradient(r);}
    virtual std::ostream&  Write(std::ostream&    ) const;

    

private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
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

