// File: Atom/kappa/Slater_IBS.H  Slater Irrep Basis Set (IBS) with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Slater_IBS_H_
#define _Atom_kappa_Slater_IBS_H_



#include "Slater_IEC.H"
#include "Slater_IE.H"
#include "../../DHF_IBS_Common.H"

class DiracIntegralTests;

namespace Atom_kappa
{
namespace Slater
{
   
template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public     Orbital_RKBL_IBS_Common<T> 
    , public     Orbital_RKBL_IE<T>
    , public     IrrepIEClient
{
    public:
    Large_Orbital_IBS() {};
    Large_Orbital_IBS(const DB_cache<T>*, const Vector<T>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;

};

template <class T> class Small_Orbital_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     Orbital_RKBS_IBS_Common<T> 
    , public     Orbital_RKBS_IE<T>
    , public     Small_IrrepIEClient
{
public:
    Small_Orbital_IBS() {};
    Small_Orbital_IBS(const DB_cache<double>*,const Large_Orbital_IBS<T>*);
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;

};

class Dirac_IrrepBasisSet
    : public virtual TOrbital_IBS<double>
    , public         Orbital_RKB_IBS_Common<double> 
    , public         Dirac_IrrepIEClient
{
    typedef Integrals_Base<double>::SMat SMat;
public:
    Dirac_IrrepBasisSet(const DB_cache<double>* db,const Vector<double>& exponents, int kappa);

    virtual size_t size() const {return Orbital_RKB_IBS_Common<double>::size();}
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};


}} //namespace Slater_mj

#endif //_Atom_kappa_Slater_IBS_H_
