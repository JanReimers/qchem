// File: Atom/kappa/Slater_IEC.H  Inegral Engine Client (IEC) for Slater basis set with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Slater_IEC_H_
#define _Atom_kappa_Slater_IEC_H_

#include "../radial/Slater/IE_Primatives.H"
#include "../IEC.H"
    
namespace Atom_kappa
{
namespace Slater
{
 
struct IrrepIEClient //Use this for the large component.
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient(size_t N,int _kappa);

    int kappa;
    double j;
    
    // void Init(const Vector  <double>& exponents);
    virtual bool Large() const {return true;}
    
private:
    virtual double Norm(double e, size_t l) const;
};

struct Small_IrrepIEClient
    : public virtual ::IrrepIEClient
    , public IrrepIEClient
    , public virtual ::Slater::IE_Primatives
{
public:
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

    void Init(const IrrepIEClient* large,const IrrepIEClient* small);
    
    const IrrepIEClient* itsLargeIEC;
    const IrrepIEClient* itsSmallIEC;
        
};


}} //namespace Atom_kappa::Slater

#endif // _Atom_kappa_Slater_IEC_H_
