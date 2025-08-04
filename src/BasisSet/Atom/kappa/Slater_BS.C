// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
module;
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.Common;
import qchem.Irrep_BS;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.Internal.radial.Slater.IE_Primatives;
import qchem.BasisSet.Atom.IE;

export namespace Atom_kappa
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


// Integral engine
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public virtual ::Slater::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
    , public virtual ::Slater::IE_Primatives
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db) : AtomIE_RKBS<T>(db) {};
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
};


// Irrep basis set

// All integrals are handled at the Orbital_RKB_IBS_Common.  i.e. they are not Slater function
// specific.
class Orbital_IBS
    : public virtual Real_OIBS
    , public         Orbital_RKB_IBS_Common1<double> 
{
public:
    Orbital_IBS(const DB_cache<double>* db,const Vector<double>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;

private:
    friend class ::DiracIntegralTests;

};

template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public     Orbital_RKBL_IBS_Common1<T> 
    , public     Orbital_RKBL_IE<T>
    , public     AtomIrrepIEClient
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
public:
    using AtomIrrepIEClient::es;
    using AtomIrrepIEClient::ns;
    
    Large_Orbital_IBS(const DB_cache<T>*, const Vector<T>& exponents, int kappa);
    virtual size_t  GetNumFunctions() const {return size();}

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;

//private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBL_IBS_Common1<T>::kappa;
};

template <class T> class Small_Orbital_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     Orbital_RKBS_IBS_Common1<T> 
    , public     Orbital_RKBS_IE<T>
    , public     AtomIrrepIEClient
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    using Orbital_RKBS_IBS_Common1<T>::large;
public:
    Small_Orbital_IBS(const DB_cache<double>*, const Vector<T>& exponents, int kappa);
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l);
    virtual size_t  GetNumFunctions() const {return size();}

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;

private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBS_IBS_Common1<T>::kappa;
};

// Full basis set
class BasisSet 
    : public ::BS_Common
    , public DB_cache<double>
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lMax);
    
};

}} //namespace Atom_kappa::Slater

