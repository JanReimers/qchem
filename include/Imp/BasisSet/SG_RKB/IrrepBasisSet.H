// File: SphericalGaussian_RKB_IrrepBasisSet.H  .
#ifndef _SphericalGaussian_RKB_BS_H_
#define _SphericalGaussian_RKB_BS_H_



#include "IEClient.H"
#include "Imp/BasisSet/Dirac_IBS.H"
#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../SphericalGaussian/IntegralEngine.H"
#include "oml/vector.h"
class DiracIntegralTests;


namespace SphericalGaussian_RKB
{

template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public     IrrepBasisSetCommon
    , protected  TIrrepBasisSetCommon<T> 
    , public     DB_RKBL<T>
    , public     IrrepIEClient
{
    public:
    Large_Orbital_IBS() {};
    Large_Orbital_IBS(const LAParams&,
        const Vector<T>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

template <class T> class Small_Orbital_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     IrrepBasisSetCommon
    , protected  TIrrepBasisSetCommon<T> 
    , public     DB_RKBS<T>
    , public     Small_IrrepIEClient
{
public:
    Small_Orbital_IBS() {};
    Small_Orbital_IBS(const LAParams&,const Large_Orbital_IBS<T>*);
    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const;
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};


class Dirac_IrrepBasisSet
    : public virtual TOrbital_IBS<double>
    , public Dirac::IrrepBasisSet<double>
    , public         Dirac_IrrepIEClient
{
    typedef DB_1E<double>::SMat SMat;
public:
    //Dirac_IrrepBasisSet() {};
    Dirac_IrrepBasisSet(const LAParams&,const DB_cache<double>*, const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
   
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};


} //namespace

#endif //_SphericalGaussian_RKB_BS_H_
