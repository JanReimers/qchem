// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace SphericalGaussian
{
    
class IrrepIEClient;
    

class IntegralEngine1
: public virtual AtomIE_FitDFT_Common
, public virtual AtomIE_1EFit_Common
, public AtomIE_1E<double>
{
public:
virtual double Integral(qchem::IType,double ea , double eb,size_t l) const;
virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;

};

class Orbital_IE
: public IntegralEngine1
, public DB_2E<double>
, public AtomIE_DFT<double>
{
public:
virtual double Overlap (double ea , double eb,size_t l) const;
protected:
Orbital_IE(const DB_BS_2E<double>* db) : DB_2E<double>(db) {};

};

class Fit_IE
: public IntegralEngine1
, public AtomIE_Fit

{
public:
virtual double Charge  (double ea ,           size_t l) const;
};
} //namespace

#endif //_SphericalGaussianIE1_H_
