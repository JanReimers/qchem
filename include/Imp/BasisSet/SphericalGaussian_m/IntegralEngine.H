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
    virtual void MakeDirect  (erij_t& Jab, const ::IEClient*) const; //Used for Hartree-Fock
    virtual void MakeExchange(erik_t& Kab, const ::IEClient*) const; //Used for Hartree-Fock

private:
    typedef AtomIrrepIEClient iriec;
    ERIJ MakeDirect  (const iriec* a, const iriec* c,const IEClient*) const;
    ERIK MakeExchange(const iriec* a, const iriec* b,const IEClient*) const;
   
    static double FourPi2;
};

} //namespace

#endif //_SphericalGaussian_mIE_H_
