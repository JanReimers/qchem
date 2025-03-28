// File: Atom/kappa/Gaussian_IEC.H  Integral Engine Client (IEC) for RKB Gaussians.
#ifndef _Atom_kappa_Gaussian_IEC_H_
#define _Atom_kappa_Gaussian_IEC_H_

#include "Imp/BasisSet/Atom/radial/Gaussian/IE_Primatives.H"
#include "Imp/BasisSet/AtomIEClient.H"

namespace Atom_kappa
{
namespace Gaussian
{
    
struct IrrepIEClient //Use this for the large component.
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
    , public virtual ::Gaussian::IE_Primatives
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



}} //namespace

#endif // _Atom_kappa_Gaussian_IEC_H_
