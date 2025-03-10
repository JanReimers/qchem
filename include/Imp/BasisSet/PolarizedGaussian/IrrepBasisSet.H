// File: PolarizedGaussian/IrrepBasisSet.H  Polarized Gaussian basis set, for MO calculations.
#ifndef _PolarizedGaussian_IrrepBasisSet_H_
#define _PolarizedGaussian_IrrepBasisSet_H_

#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "Polarization.H"
#include "IntegralEngine.H"
#include "IEClient.H"
#include "Block.H"
#include "Imp/Containers/ptr_vector.h"


namespace PolarizedGaussian
{
class Reader;

// Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         IrrepBasisSetCommon
    , public         IrrepIEClient
{
public:
    typedef typename ::IrrepBasisSet::RVec RVec;
    //IrrepBasisSet();
    IrrepBasisSet(Reader*, const Cluster*);
    IrrepBasisSet(const Vector<double>& exponents, size_t L, const Cluster*);
    IrrepBasisSet(const Vector<double>& exponents, size_t L);

 
    virtual std::ostream&  Write(std::ostream&    ) const;

private:
    IrrepBasisSet(const IrrepBasisSet* bs,const optr_vector1<Block*>&);
    void MakeBasisFunctions(const RVec& norms);
   
    optr_vector1<Block*> itsBlocks;
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
    Orbital_IBS(const LAParams&, IntegralDataBase<double>*,const DB_BS_2E<double>* db, Reader*, const Cluster*);
    Orbital_IBS(const LAParams&, IntegralDataBase<double>*,const DB_BS_2E<double>* db, const Vector<double>& exponents, size_t L, const Cluster*);
    Orbital_IBS(const LAParams&, IntegralDataBase<double>*,const DB_BS_2E<double>* db, const Vector<double>& exponents, size_t L);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* Clone(const RVec3&) const;

};

class Fit_IBS 
    : public virtual ::Fit_IBS
    , public         IrrepBasisSet
    , public         TIrrepBasisSetCommon<double>
    , public         IntegralEngine1
{
    typedef DB_Fit::SMat SMat;
public:
    Fit_IBS(const LAParams&, IntegralDataBase<double>*, Reader*, const Cluster*);
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
    virtual UniqueID::IDtype GetID() const {return UniqueID::GetID();}
};

} //namespace PolarizedGaussian

#endif //_PolarizedGaussian_IrrepBasisSet_H_
