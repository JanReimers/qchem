// File: Slater_m/BasisSet.H
#ifndef _Slater_m_BasisSet_H_
#define _Slater_m_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"
#include "../Slater/BasisSet.H"

namespace Slater_m
{

class BasisSet 
    : public Slater::BasisSet
{
    typedef Slater::BasisSet::RVec RVec;
public:
    BasisSet() {};
    BasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;

};

} //namespace Slater_m

#endif // _Slater_m_BasisSet_H_
