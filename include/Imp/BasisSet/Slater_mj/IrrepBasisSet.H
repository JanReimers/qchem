// File: Slater_mj/IrrepBasisSet.H  Spherical Slater basis set with orbital angular momentum kappa,mj.
#ifndef _Slater_mj_IrrepBasisSet_H_
#define _Slater_mj_IrrepBasisSet_H_



#include "IEClient.H"
#include "Imp/BasisSet/Dirac_IBS.H"

#include "IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"

class DiracIntegralTests;

namespace Slater_mj
{
   
template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public     IrrepBasisSetCommon
    , protected  TIrrepBasisSetCommon<T> 
    , public     Orbital_RKBL_IE<T>
    , public     IrrepIEClient
{
    public:
    Large_Orbital_IBS() {};
    Large_Orbital_IBS(const DB_cache<T>*, const Vector<T>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

template <class T> class Small_Orbital_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     IrrepBasisSetCommon
    , protected  TIrrepBasisSetCommon<T> 
    , public     Orbital_RKBS_IE<T>
    , public     Small_IrrepIEClient
{
public:
    Small_Orbital_IBS() {};
    Small_Orbital_IBS(const DB_cache<double>*,const Large_Orbital_IBS<T>*);
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

class Dirac_IrrepBasisSet
    : public virtual TOrbital_IBS<double>
    , public Dirac::IrrepBasisSet<double>
    , public         Dirac_IrrepIEClient
{
    typedef Integrals_Base<double>::SMat SMat;
public:
    //Dirac_IrrepBasisSet() {};
    Dirac_IrrepBasisSet(const DB_cache<double>* db,const Vector<double>& exponents, int kappa);

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};


} //namespace Slater_mj

#endif //_Slater_mj_IrrepBasisSet_H_
