// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


// #include <AnalyticIE.H>
#include "IE_Primatives.H"
#include "../AtomIE.H"
#include "oml/vector.h"


namespace SphericalGaussian
{
    
class IE_Common
: public IE_Primatives
, public AtomIE_Overlap<double>
, public AtomIE_Grad2<double>
, public AtomIE_Nuclear<double>
{
public:
    
    virtual double Repulsion(double ea, double ec,size_t la,size_t lc) const;
protected:
    IE_Common(const DB_cache<double>* db) : AtomIE_Overlap<double>(db), AtomIE_Grad2<double>(db), AtomIE_Nuclear<double>(db) {};

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
