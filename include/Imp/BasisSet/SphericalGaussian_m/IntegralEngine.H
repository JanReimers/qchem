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
   
    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock

private:
    typedef AtomIrrepIEClient iriec;
    virtual ERIJ MakeDirect  (const iriec* a, const iriec* c,const AtomIEClient*) const;
    virtual ERIK MakeExchange(const iriec* a, const iriec* b,const AtomIEClient*) const;
   
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    
    static double FourPi2;
};

} //namespace

#endif //_SphericalGaussian_mIE_H_
