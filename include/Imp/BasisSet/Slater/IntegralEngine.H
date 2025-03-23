// File: SlaterIE.H  Integral Engine for Spherical gaussians.
#ifndef _SlaterIE_H_
#define _SlaterIE_H_

#include "../AtomIE.H"
#include "IE_Primatives.H"

namespace Slater
{

class IE_Common
    : public virtual IE_Primatives
    , public AtomIE_Overlap<double>
    , public AtomIE_Grad2<double>
    , public AtomIE_Nuclear<double>
{
public:
protected:
    IE_Common(const DB_cache<double>* db) : AtomIE_Overlap<double>(db), AtomIE_Grad2<double>(db), AtomIE_Nuclear<double>(db) {};
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
: public AtomIE_Fit
, public AtomIE_Overlap<double>
, public virtual IE_Primatives

{
public:
    virtual double Charge  (double ea ,           size_t l) const;
protected:
    Fit_IE(const DB_cache<double>* db) : AtomIE_Fit(db), AtomIE_Overlap<double>(db) {};
 

};

} //namespace

#endif //_SlaterIE_H_
