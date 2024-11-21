// File: SphericalGaussian_mIE.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussian_mIE_H_
#define _SphericalGaussian_mIE_H_

#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "Imp/BasisSet/SphericalGaussian_m/IEClient.H" 

namespace SphericalGaussian_m
{

    
class IntegralEngine
    : public virtual AnalyticIE<double>
    , public SphericalGaussian::IntegralEngine
{

public:
   

private:
    typedef AtomIrrepIEClient iriec;
   
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    
};

} //namespace

#endif //_SphericalGaussian_mIE_H_
