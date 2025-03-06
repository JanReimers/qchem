// File: Slater_m/IrrepBasisSet.H  Spherical Slater basis set with orbital angular momentum l,m.
#ifndef _Slater_m_IrrepBasisSet_H_
#define _Slater_m_IrrepBasisSet_H_



#include "IEClient.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../Slater/IntegralEngine.H"
#include "oml/vector.h"

namespace Slater_m
{
   
class IrrepBasisSet
    : public virtual TOrbital_HF_IBS<double>
    , public         IrrepBasisSetCommon
    , public         Orbital_IBS_Common<double>
    , public         Slater::IntegralEngine1
    , public         IrrepIEClient
{
public:
    //IrrepBasisSet();
    IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const DB_BS_2E<double>* db, const Vector<double>& exponents, size_t L,int m);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

} //namespace Slater_m

#endif //_Slater_m_IrrepBasisSet_H_
