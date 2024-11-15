// File: Slater_m/IntegralEngine.H  Integral Engine for Slater function with l,m QNs.
#ifndef _Slater_m_IntegralEngine_H_
#define _Slater_m_IntegralEngine_H_


#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/Slater_m/IEClient.H"


namespace Slater_m
{
    
class IEClient;
    
class IntegralEngine
    : public virtual AnalyticIE<double>
    , public Slater::IntegralEngine  
{

public:
    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock


private:
    static double FourPi2;

};

} //namespace

#endif //_Slater_m_IntegralEngine_H_
