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
    
};

} //namespace Slater_mj

#endif // _Slater_mj_BasisSet_H_
