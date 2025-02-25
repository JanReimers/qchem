// File: SphericalGaussian_RKB_IrrepBasisSet.H  .
#ifndef _SphericalGaussian_RKB_BS_H_
#define _SphericalGaussian_RKB_BS_H_



#include "IEClient.H"
#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../SphericalGaussian/IntegralEngine.H"
#include "oml/vector.h"
class DiracIntegralTests;


namespace SphericalGaussian_RKB
{

class Large_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         Orbital_IBS_Common<double>
    , public         SphericalGaussian::IntegralEngine1
    , public         IrrepIEClient
{
public:
    Large_IrrepBasisSet();
    Large_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,
        const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}

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
    , public         Orbital_IBS_Common<double>
    , public         Small_IrrepIEClient
    , public         SphericalGaussian_RKB::Small_IntegralEngine1

{
public:
    Small_IrrepBasisSet() {};
    Small_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Large_IrrepBasisSet*);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
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
    , public         Orbital_IBS_Common<double>
    , public         DB_1E<double>
    , public         Dirac_IrrepIEClient
{
    typedef DB_1E<double>::SMat SMat;
public:
    Dirac_IrrepBasisSet() {};
    Dirac_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
    virtual SMat MakeIntegrals(qchem::IType,const Cluster*) const;

private:
    friend class ::DiracIntegralTests;
    
    Large_IrrepBasisSet* itsLargeBS;
    Small_IrrepBasisSet* itsSmallBS;
};


} //namespace

#endif //_SphericalGaussian_RKB_BS_H_
