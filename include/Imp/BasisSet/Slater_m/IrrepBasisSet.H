// File: Slater_m/IrrepBasisSet.H  Spherical Slater basis set with orbital angular momentum l,m.
#ifndef _Slater_m_IrrepBasisSet_H_
#define _Slater_m_IrrepBasisSet_H_



#include "IEClient.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "oml/vector.h"

namespace Slater_m
{
   
class IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         TIrrepBasisSetCommon<double>
    , public         IrrepIEClient
{
public:
    IrrepBasisSet();
    IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,size_t size, double minexp, double maxexp, size_t L,int m);

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

} //namespace Slater_m

#endif //_Slater_m_IrrepBasisSet_H_
