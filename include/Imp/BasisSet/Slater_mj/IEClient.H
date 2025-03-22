// File: Slater_mj/IEClient.H
#ifndef _Slater_mj_IEClient_H_
#define _Slater_mj_IEClient_H_

#include "../Slater/IE_Primatives.H"
#include "Imp/BasisSet/AtomIEClient.H"
    
namespace Slater_mj
{
 
struct IrrepIEClient //Use this for the large component.
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient() {}; 
    IrrepIEClient(size_t N,int _kappa);

    int kappa;
    double j;
    
    void Init(const Vector  <double>& exponents);
    virtual bool Large() const {return true;}
    
private:
    virtual double Norm(double e, size_t l) const;
};

struct Small_IrrepIEClient
    : public virtual ::IrrepIEClient
    , public IrrepIEClient
    , public virtual Slater::IE_Primatives
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

    void Init(const Slater_mj::IrrepIEClient* large,const Slater_mj::IrrepIEClient* small);
    
    const Slater_mj::IrrepIEClient* itsLargeIEC;
    const Slater_mj::IrrepIEClient* itsSmallIEC;
        
};


} //namespace Slater_mj

#endif // _Slater_mj_IEClient_H_
