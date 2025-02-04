// File: SphericalGaussian_RKB_IrrepBasisSet.H  .
#ifndef _SphericalGaussianBS_H_
#define _SphericalGaussianBS_H_



#include "IEClient.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "oml/vector.h"
class DiracIntegralTests;

namespace SphericalGaussian_RKB
{
 
class Large_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         TIrrepBasisSetCommon<double>
    , public         IrrepIEClient
{
public:
    Large_IrrepBasisSet();
    Large_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,
        const Vector<double>& exponents, int kappa);

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

class Small_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         TIrrepBasisSetCommon<double>
    , public         Small_IrrepIEClient
{
public:
    Small_IrrepBasisSet() {};
    Small_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Large_IrrepBasisSet*);

    virtual ::IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

class Dirac_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         TIrrepBasisSetCommon<double>
    , public         Dirac_IrrepIEClient
{
public:
    Dirac_IrrepBasisSet() {};
    Dirac_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Vector<double>& exponents, int kappa);

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    friend class ::DiracIntegralTests;
    
    Large_IrrepBasisSet* itsLargeBS;
    Small_IrrepBasisSet* itsSmallBS;
};


} //namespace

#endif //_SphericalGaussianBS_H_
