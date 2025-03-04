// File: Slater_mj/IrrepBasisSet.H  Spherical Slater basis set with orbital angular momentum kappa,mj.
#ifndef _Slater_mj_IrrepBasisSet_H_
#define _Slater_mj_IrrepBasisSet_H_



#include "IEClient.H"
#include "IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../Slater/IntegralEngine.H"
#include "oml/vector.h"

class DiracIntegralTests;

namespace Slater_mj
{
   
class Large_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         Orbital_IBS_Common<double>
    , public         Slater::IntegralEngine1
    , public         IrrepIEClient
{
public:
    Large_IrrepBasisSet() {};
    Large_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual ::IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

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
    , public         Slater_mj::Small_IntegralEngine1
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
    //, public         DB_2E<double>    
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

    virtual SMat MakeOverlap() const;
    virtual SMat MakeKinetic() const;
    virtual SMat MakeNuclear(const Cluster*) const;
    typedef typename DB_2E<double>::bs_t bs_t;
    // virtual ERI4 MakeDirect  (const bs_t& c) const;
    // virtual ERI4 MakeExchange(const bs_t& b) const;
private:
    friend class ::DiracIntegralTests;
    
    Large_IrrepBasisSet* itsLargeBS;
    Small_IrrepBasisSet* itsSmallBS;
};


} //namespace Slater_mj

#endif //_Slater_mj_IrrepBasisSet_H_
