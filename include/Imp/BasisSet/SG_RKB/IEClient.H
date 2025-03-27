// File: SphericalGaussian_RKB/IEClient.H
#ifndef _SphericalGaussian_RKB_IEClient_H_
#define _SphericalGaussian_RKB_IEClient_H_

#include "Imp/BasisSet/Atom/radial/Gaussian/IE_Primatives.H"
#include "Imp/BasisSet/AtomIEClient.H"

namespace SphericalGaussian_RKB
{
    
struct IrrepIEClient //Use this for the large component.
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
    , public virtual Gaussian::IE_Primatives
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N, int kappa);
    
    int kappa;
    double j;
    
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
