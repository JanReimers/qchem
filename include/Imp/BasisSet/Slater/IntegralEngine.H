// File: SlaterIE.H  Integral Engine for Spherical gaussians.
#ifndef _SlaterIE_H_
#define _SlaterIE_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace Slater
{
    
class IrrepIEClient;
    
class IntegralEngine
    : public virtual AnalyticIE<double>
    , public AtomIE
{
public:
    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock

private:
    virtual double Overlap(double ea, double eb,size_t l) const;
    virtual double Kinetic(double ea ,double eb,size_t l) const;
    virtual double Nuclear(double ea, double eb,size_t l) const;
    virtual double Charge (double ea,           size_t l) const;
    virtual double Repulsion(double eab, double ec,size_t la,size_t lc) const;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    
    static double FourPi2;

};

} //namespace

#endif //_SlaterIE_H_
