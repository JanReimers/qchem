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
    typedef AnalyticIE<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;

            
    IntegralEngine() {};
    
    virtual size_t size() const {return -1;}

    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock
    virtual void MakeDirect  (erij_t& Jab, const ::IEClient*) const; //Used for Hartree-Fock
    virtual void MakeExchange(erik_t& Kab, const ::IEClient*) const; //Used for Hartree-Fock

    virtual void Report(std::ostream&) const;
    virtual IE* Clone() const;

private:
    static const IrrepIEClient* dcast(iec_t*);
    ERIJ MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c,const IEClient*) const;
    ERIK MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b,const IEClient*) const;
   
    static double FourPi2;
};

} //namespace

#endif //_SphericalGaussian_mIE_H_
