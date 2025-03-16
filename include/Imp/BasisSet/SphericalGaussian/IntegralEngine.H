// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace SphericalGaussian
{
    
class IrrepIEClient;
    

class IE_Common
: public virtual AtomIE_FitDFT_Common
, public virtual AtomIE_1EFit_Common
, public AtomIE_1E<double>
{
public:
    virtual double Overlap(double ea , double eb,size_t l_total) const;
    virtual double Kinetic(double ea , double eb,size_t la, size_t lb) const;
    virtual double Nuclear(double ea , double eb,size_t l_total) const;
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;
protected:
    IE_Common(const DB_cache<double>* db) : AtomIE_1E<double>(db) {};

};

class Orbital_IE
: public IE_Common
, public DB_2E<double>
, public AtomIE_DFT<double>
{
public:
protected:
Orbital_IE(const DB_BS_2E<double>* db) : IE_Common(db), DB_2E<double>(db), AtomIE_DFT<double>(db) {};

};

class Fit_IE
: public IE_Common
, public AtomIE_Fit

{
public:
    virtual double Charge  (double ea ,           size_t l) const;
protected:
    Fit_IE(const DB_cache<double>* db) : IE_Common(db), AtomIE_Fit(db) {};
};
} //namespace

#endif //_SphericalGaussianIE1_H_
