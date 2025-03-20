// File: SphericalGaussian_RKB_IrrepBasisSet.H  .
#ifndef _SphericalGaussian_RKB_BS_H_
#define _SphericalGaussian_RKB_BS_H_



#include "IEClient.H"
#include "IntegralEngine.H"
#include "../Dirac_IBS.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
class DiracIntegralTests;


namespace SphericalGaussian_RKB
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
    Large_Orbital_IBS(const LAParams&,const DB_cache<T>*, const Vector<T>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
   
    

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
    Small_Orbital_IBS(const LAParams&,const DB_cache<T>*,const Large_Orbital_IBS<T>*);
    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
  

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};


class Dirac_IrrepBasisSet
    : public virtual TOrbital_IBS<double>
    , public Dirac::IrrepBasisSet<double>
    , public         Dirac_IrrepIEClient
{
public:
    Dirac_IrrepBasisSet(const LAParams&,const DB_cache<double>*, const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
   
    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};


} //namespace

#endif //_SphericalGaussian_RKB_BS_H_
