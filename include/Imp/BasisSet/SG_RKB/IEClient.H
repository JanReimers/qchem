// File: SphericalGaussian_RKB/IEClient.H
#ifndef _SphericalGaussian_RKB_IEClient_H_
#define _SphericalGaussian_RKB_IEClient_H_

#include "Imp/BasisSet/AtomIEClient.H"

namespace SphericalGaussian_RKB
{
    
struct IrrepIEClient //Use this for the large component.
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N, int kappa);
    
    int kappa;
    double j;
    
    void Init(const Vector<double>& exponents);
    virtual bool Large() const {return true;}
private:
    virtual double Norm(double e, size_t l) const;
};

struct Small_IrrepIEClient
    : public virtual ::IrrepIEClient
    , public IrrepIEClient
{
public:
    Small_IrrepIEClient() {}; 
    Small_IrrepIEClient(size_t N,int kappa) : IrrepIEClient(N,kappa) {};
    virtual bool Large() const {return false;}

private:
    virtual double Norm(double e, size_t l) const;
    static double Kinetic(double ea, double eb,size_t l);     
};

struct Dirac_IrrepIEClient
    : public virtual ::IrrepIEClient
{
    Dirac_IrrepIEClient() : itsLargeIEC(0), itsSmallIEC(0) {};
    virtual size_t size() const;

    void Init(const SphericalGaussian_RKB::IrrepIEClient* large,const SphericalGaussian_RKB::IrrepIEClient* small);
    
    const SphericalGaussian_RKB::IrrepIEClient* itsLargeIEC;
    const SphericalGaussian_RKB::IrrepIEClient* itsSmallIEC;
        
};


} //namespace

#endif // _SphericalGaussian_RKB_IEClient_H_
