// File: Atom/kappa/Gaussian_BS.H  Restricted Kinetic Balance (RKB) Basis Set (BS).
module;
#include <iosfwd>
class DiracIntegralTests;
export module qchem.BasisSet.Atom.Internal.kappa.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.IE_Primatives;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;

export namespace Atom_kappa
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

// Integral engine
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public virtual ::Gaussian::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db,const ::IE_Primatives* pie) : AtomIE_RKBL<T>(db,pie) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
    , public virtual ::Gaussian::IE_Primatives
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db,const ::IE_Primatives* pie) : AtomIE_RKBS<T>(db,pie) {};
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
};

// Irrep Basis set
// All integrals are handled at the Orbital_RKB_IBS_Common.  i.e. they are not Gaussian
// specific.
class Orbital_RKB_IBS
    : public virtual Real_OIBS
    , public Orbital_RKB_IBS_Common<double>
{
public:
    Orbital_RKB_IBS(const DB_cache<double>*, const ::IE_Primatives*, const Vector<double>& exponents, int kappa);
   
    virtual std::ostream&  Write(std::ostream&    ) const;

private:
    friend class ::DiracIntegralTests;

};

template <class T> class Orbital_RKBL_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public Orbital_RKBL_IBS_Common<T> 
    , public Orbital_RKBL_IE<T>
    , public AtomIrrepIEClient
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
public:
    Orbital_RKBL_IBS(const DB_cache<T>*,const IE_Primatives* pie, const Vector<T>& exponents, int kappa);

    virtual size_t  GetNumFunctions() const {return size();}
    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;
//private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBL_IBS_Common<T>::kappa;
};

template <class T> class Orbital_RKBS_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     Orbital_RKBS_IBS_Common<T> 
    , public     Orbital_RKBS_IE<T>
    , public     AtomIrrepIEClient
{
    typedef typename VectorFunction<T>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<T>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    using Orbital_RKBS_IBS_Common<T>::large;
public:
    Orbital_RKBS_IBS(const DB_cache<T>*,const IE_Primatives* pie,const Vector<T>& exponents,int kappa);

    virtual size_t  GetNumFunctions() const {return size();}
    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual const SMatrix<T>& Overlap() const {return DB_Kinetic<T>::Kinetic();}
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBS_IBS_Common<T>::kappa;
};



// Full basis set
class BasisSet 
    : public ::BS_Common
    , public DB_cache<double>
    , public ::Gaussian::IE_Primatives
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t lmax);
    
};

}} //namespace

