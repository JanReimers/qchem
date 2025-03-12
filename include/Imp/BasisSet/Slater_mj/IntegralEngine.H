// File: Slater_mj/IntegralEngine.H  Integral Engine for Slater function with kappa,mj QNs.
#ifndef _Slater_mj_IntegralEngine_H_
#define _Slater_mj_IntegralEngine_H_

#include "Imp/BasisSet/Dirac_IBS.H"

class DiracIntegralTests;

namespace Slater_mj
{
    

class Small_IntegralEngine1
    : public virtual Integrals_RKBS<double>
    , public DB_RKBS<double>
{
public:
    virtual  double Integral(qchem::IType,double ea , double eb,size_t l) const;
};

} //namespace

#endif //_Slater_m_IntegralEngine_H_
