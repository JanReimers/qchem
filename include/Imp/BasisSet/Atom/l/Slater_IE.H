// File: Atom/l/Slater_IE.H Integral Engine for Slater atom basis functions.
#ifndef _Atom_l_Slater_IE_H_
#define _Atom_l_Slater_IE_H_

#include "../../Atom_IE_HF.H"
#include "../../Atom_IE_DHF.H"
#include "../../Atom_IE_DFT.H"
#include "../../Atom_IE_Fit.H"
#include "../radial/Slater/IE_Primatives.H"

namespace Atoml
{
namespace Slater
{

class IE_Common
    : public virtual ::Slater::IE_Primatives
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
, public virtual ::Slater::IE_Primatives

{
protected:
    Fit_IE(const DB_cache<double>* db) : AtomIE_Fit(db), AtomIE_Overlap<double>(db) {};
};

}} //namespace

#endif //_Atom_l_Slater_IE_H_
