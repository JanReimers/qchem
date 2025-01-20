// File: Slater_m/IntegralEngine.H  Integral Engine for Slater function with l,m QNs.
#ifndef _Slater_m_IntegralEngine_H_
#define _Slater_m_IntegralEngine_H_


#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/IEClient.H"


namespace Slater_m
{
    
class IntegralEngine
    : public virtual AnalyticIE<double>
    , public Slater::IntegralEngine  
{

public:


protected:
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;

};

} //namespace

#endif //_Slater_m_IntegralEngine_H_
