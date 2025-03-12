// File: SphericalGaussian_RKB_IE.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussian_RKB_IE1_H_
#define _SphericalGaussian_RKB_IE1_H_


#include <AnalyticIE.H>
#include "IEClient.H"
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"

class DiracIntegralTests;

namespace SphericalGaussian_RKB
{
    
class IrrepIEClient;
    
class Small_IntegralEngine1
    : public virtual Integrals_1E<double>
    , public AtomIE_1E<double>
    , public DB_RKBS<double>
{
public:
    virtual  double Integral(qchem::IType,double ea , double eb,size_t l) const;
};



} //namespace SphericalGaussian_RKB


#endif //_SphericalGaussian_RKB_IE_H_
