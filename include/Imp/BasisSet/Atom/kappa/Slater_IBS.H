// File: Atom/kappa/Slater_IBS.H  Slater Irrep Basis Set (IBS) with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Slater_IBS_H_
#define _Atom_kappa_Slater_IBS_H_



#include "Slater_IE.H"
#include "../IEC.H"
#include "../../DHF_IBS_Common.H"

class DiracIntegralTests;

namespace Atom_kappa
{
namespace Slater
{
   
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
    using ::IrrepBasisSet::RVec3;
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
    using ::IrrepBasisSet::RVec3;
public:
    Small_Orbital_IBS(const DB_cache<double>*, const Vector<T>& exponents, int kappa);
    virtual void InsertBasisFunctions(const Orbital_RKBL_IBS<T>* l);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
    using Orbital_RKBS_IBS_Common<T>::kappa;
};



}} //namespace Slater_mj

#endif //_Atom_kappa_Slater_IBS_H_
