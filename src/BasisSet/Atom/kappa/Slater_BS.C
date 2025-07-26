// File: Atom/kappa/Slater_BS.H  Slater Basis Set (BS) with Restricted Kinetic Balance (RKB).
module;
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.Common;
import qchem.Irrep_BS;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisFunction;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.Internal.radial.Slater.IE_Primatives;
import qchem.BasisSet.Atom.IE;

export namespace Atom_kappa
{
namespace Slater
{

// Basis function
class Small_BasisFunction;
    
class Large_BasisFunction
    : public virtual ::Real_BF
{
public:
    Large_BasisFunction() {};
    Large_BasisFunction(double ex, int kappa, double mj, double norm);
    
    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::Real_BF* Clone(        ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;
private:
    friend class Small_BasisFunction;
    double itsExponent;
    int kappa;
    double mj;
    int l;
    double itsNormalization;
};

//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r) as:
//
//                                  -e*r^l*exp(-e*r), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)g = {
//                                  ((2l+1)/r-e)*r^lexp(-e*r), kappa>0
//
class Small_BasisFunction
    : public virtual ::Real_BF
{
public:
    Small_BasisFunction();
    Small_BasisFunction(const Large_BasisFunction*,double norm);
    
    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::Real_BF* Clone(        ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3   Gradient  (const RVec3&) const;
private:
    const Large_BasisFunction* Pr;
    double itsNormalization;
};


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
    : public virtual TOrbital_IBS<double>
    , public         Orbital_RKB_IBS_Common<double> 
{
public:
    Orbital_IBS(const DB_cache<double>* db,const Vector<double>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};

template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public     Orbital_RKBL_IBS_Common<T> 
    , public     Orbital_RKBL_IE<T>
    , public     AtomIrrepIEClient
{
    public:
    Large_Orbital_IBS(const DB_cache<T>*, const Vector<T>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBL_IBS_Common<T>::kappa;
};

template <class T> class Small_Orbital_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     Orbital_RKBS_IBS_Common<T> 
    , public     Orbital_RKBS_IE<T>
    , public     AtomIrrepIEClient
{
public:
    Small_Orbital_IBS(const DB_cache<double>*, const Vector<T>& exponents, int kappa);
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBS_IBS_Common<T>::kappa;
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

