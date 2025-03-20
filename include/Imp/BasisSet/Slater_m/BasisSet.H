// File: Slater_m/BasisSet.H
#ifndef _Slater_m_BasisSet_H_
#define _Slater_m_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"
#include "../Slater/BasisSet.H"

namespace Slater_m
{

class BasisSet 
    : public Slater::BS_Common
    , public AtomIE_BS_2E_lm
{
public:
    BasisSet() {};
    BasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    
};

} //namespace Slater_m

#endif // _Slater_m_BasisSet_H_
