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
   
class IrrepBasisSet
    : public virtual Orbital_IBS<double>
    , protected         IrrepBasisSetCommon
    , protected         Orbital_IBS_Common<double>
    , protected         IntegralEngine1
    , protected         IrrepIEClient
{
public:
    IrrepBasisSet();
    IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Vector<double>& exponents, size_t L);

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

} //namespace Slater

#endif //_SlaterBS_H_
