// File: Slater_mj/BasisSet.H
#ifndef _Slater_mj_BasisSet_H_
#define _Slater_mj_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"
#include "Imp/Containers/ERI4.H"

namespace Slater_mj
{


class DiracBasisSet 
    : public ::BasisSetImp
{
public:
    DiracBasisSet() {};
    DiracBasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t lMax);
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const
    {
        return itsIE->MakeDirect(a,c);
    }
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const
    {
        return itsIE->MakeExchange(a,b);
    }
};

} //namespace Slater_mj

#endif // _Slater_mj_BasisSet_H_
