// File: Atom/kappa/Gaussian_BS.H  Restricted Kinetic Balance (RKB) Basis Set (BS).
module;
#include <iosfwd>
class DiracIntegralTests;
export module qchem.BasisSet.Atom.kappa.GaussianBS;
import qchem.BasisSet.Atom.radial.IE_Primatives;
import qchem.BasisSet.Imp.HeapDB;
import qchem.BasisSet.Common;
import qchem.BasisSet.IBS_Common;
import qchem.BasisFunction;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;

export namespace Atom_kappa
{
namespace Gaussian
{

    // Basis functions
    class Small_BasisFunction;

class Large_BasisFunction
    : public virtual ::Real_BF
{
public:
    Large_BasisFunction(                             );
    Large_BasisFunction(double theExponent, int kappa, double norm);
    
    virtual std::ostream&  Write(std::ostream&) const;
    virtual Real_BF* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    friend class Small_BasisFunction;

    double itsExponent;
    int    kappa,l; // l is redundant but convenient for r^l
    double itsNormalization;
};

//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r^2) as:
//
//                                   { -2e*r^(l+1)*exp(-e*r^2), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)P(r) = {
//                                   { ((2l+1)/r-2er)*r^l*exp(-e*r^2), kappa>0
//
class Small_BasisFunction
    : public virtual ::Real_BF
{
public:
    Small_BasisFunction();
    Small_BasisFunction(const Large_BasisFunction*,double norm);

    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::Real_BF* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;
private:
    const Large_BasisFunction* Pr;
    double itsNormalization;
};

// Integral engine
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public virtual ::Gaussian::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
    , public virtual ::Gaussian::IE_Primatives
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db) : AtomIE_RKBS<T>(db) {};
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
};

// Irrep Basis set
// All integrals are handled at the Orbital_RKB_IBS_Common.  i.e. they are not Gaussian
// specific.
class Orbital_IBS
    : public virtual TOrbital_IBS<double>
    , public Orbital_RKB_IBS_Common<double>
{
public:
    using RVec3=::IrrepBasisSet::RVec3;
    Orbital_IBS(const DB_cache<double>*, const Vector<double>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};

template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public Orbital_RKBL_IBS_Common<T> 
    , public Orbital_RKBL_IE<T>
    , public AtomIrrepIEClient
{
    using RVec3=::IrrepBasisSet::RVec3;
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
    using RVec3=::IrrepBasisSet::RVec3;
public:
    Small_Orbital_IBS(const DB_cache<T>*,const Vector<T>& exponents,int kappa);
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
    BasisSet(size_t N, double minexp, double maxexp, size_t lmax);
    
};

}} //namespace

