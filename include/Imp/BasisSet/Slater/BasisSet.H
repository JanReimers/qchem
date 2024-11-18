// File Slater/BasisSet.H
#ifndef _Slater_BasisSet_H_
#define _Slater_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"

namespace Slater
{

class BasisSet 
    : public ::BasisSetImp
{
public:
    BasisSet() {};
    BasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    
};

} //namespace Slater

#endif // _Slater_BasisSet_H_
