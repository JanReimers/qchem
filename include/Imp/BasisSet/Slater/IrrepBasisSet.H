// File: SlaterBS.H  Spherical Slater basis set with orbital angular momentum l.
#ifndef _SlaterBS_H_
#define _SlaterBS_H_



#include "IEClient.H"
#include "IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "oml/vector.h"


namespace Slater
{
   
// Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         IrrepBasisSetCommon
    , public         IrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t L);
    virtual std::ostream&  Write(std::ostream&    ) const;
};

class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Orbital_IE
{
public:
    Orbital_IBS(const LAParams& lap,const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , Orbital_IBS_Common<double>(lap)
    , Orbital_IE(db)
    {};


    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};

class Fit_IBS 
: public virtual ::Fit_IBS
, public virtual FitIntegrals_1E<double>
, public IrrepBasisSet
, public TIrrepBasisSetCommon<double>
, public Fit_IBS_Common
, public Fit_IE

{
public:
    Fit_IBS(const LAParams& lap,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , TIrrepBasisSetCommon<double>(lap)
    {};
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
};

} //namespace Slater

#endif //_SlaterBS_H_
