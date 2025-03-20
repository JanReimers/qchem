// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace SphericalGaussian
{
    
class IrrepIEClient;
    
class IE_Primatives 
    : public virtual  Primative_Overlap<double>
    , public virtual  Primative_Kinetic<double>
    , public virtual  Primative_Nuclear<double>
    , public virtual Primative_Repulsion<double>
{
    virtual double Overlap(double ea , double eb,size_t l_total) const;
    virtual double Kinetic(double ea , double eb,size_t la, size_t lb) const;
    virtual double Nuclear(double ea , double eb,size_t l_total) const;
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;
};

class IE_Common
: public IE_Primatives
, public AtomIE_Overlap<double>
, public AtomIE_Kinetic<double>
, public AtomIE_Nuclear<double>
{
public:
    
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;
protected:
    IE_Common(const DB_cache<double>* db) : AtomIE_Overlap<double>(db), AtomIE_Kinetic<double>(db), AtomIE_Nuclear<double>(db) {};

};

class Orbital_IE
: public IE_Common
, public DB_2E<double>
, public AtomIE_DFT<double>
{
protected:
    Orbital_IE(const DB_BS_2E<double>* db) : IE_Common(db), DB_2E<double>(db), AtomIE_DFT<double>(db) {};

};

class Fit_IE
: public AtomIE_Fit
, public AtomIE_Overlap<double>
, public IE_Primatives
{
public:
    virtual double Charge  (double ea ,           size_t l) const;
protected:
    Fit_IE(const DB_cache<double>* db) : AtomIE_Fit(db), AtomIE_Overlap<double>(db) {};
};
} //namespace

#endif //_SphericalGaussianIE1_H_
