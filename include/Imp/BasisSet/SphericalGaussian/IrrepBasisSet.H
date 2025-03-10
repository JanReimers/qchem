// File: SphericalGaussianBS.H  Sphericalic gaussian basis set with orbital angular momentum.
#ifndef _SphericalGaussianBS_H_
#define _SphericalGaussianBS_H_



#include "IEClient.H"
#include "IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "oml/vector.h"

namespace SphericalGaussian
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
    , public virtual TOrbital_DFT_IBS<double>
    , public         IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         IntegralEngine1

{
public:
    Orbital_IBS(const LAParams& lap,IntegralDataBase<double>* theDB,const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , Orbital_IBS_Common<double>(lap,theDB)
    , IntegralEngine1(db)
    {};


    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};

class Fit_IBS 
: public virtual ::Fit_IBS
, public         IrrepBasisSet
, public         TIrrepBasisSetCommon<double>
, public         IntegralEngine1
{
public:
    Fit_IBS(const LAParams& lap,IntegralDataBase<double>* theDB,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , TIrrepBasisSetCommon<double>(lap,theDB)
    , IntegralEngine1(0)
    {};
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
};


} //namespace

#endif //_SphericalGaussianBS_H_
