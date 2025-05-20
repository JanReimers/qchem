// File: Atom/l/BSpline_IE.H Integral Engine for BSpline basis functions.
#ifndef _Atom_l_BSpline_IE_H_
#define _Atom_l_BSpline_IE_H_

#include "../radial/BSpline/IE_HF.H"
#include "../radial/BSpline/IE_DHF.H"
#include "../radial/BSpline/IE_DFT.H"
#include "../radial/BSpline/IE_Fit.H"
#include "../radial/BSpline/IE_Primatives.H"

namespace Atoml
{
namespace BSpline
{

template <size_t K> class IE_Common
    : public virtual ::BSpline::IE_Primatives<K>
    , public ::BSpline::IE_Overlap<double,K>
    , public ::BSpline::IE_Kinetic  <double,K>
    , public ::BSpline::IE_Inv_r1<double,K>
{
public:
protected:
    IE_Common(const DB_cache<double>* db) : ::BSpline::IE_Overlap<double,K>(db), ::BSpline::IE_Kinetic<double,K>(db), ::BSpline::IE_Inv_r1<double,K>(db) {};
};

template <size_t K> class Orbital_IE
: public IE_Common<K>
, public DB_2E<double>
, public ::BSpline::IE_DFT<double,K>
{
public:
protected:
    Orbital_IE(const DB_BS_2E<double>* db) : IE_Common<K>(db), DB_2E<double>(db), ::BSpline::IE_DFT<double,K>(db) {};

};

template <size_t K> class Fit_IE
: public ::BSpline::IE_Fit<K>
, public ::BSpline::IE_Overlap<double,K>
, public virtual ::BSpline::IE_Primatives<K>

{
protected:
    Fit_IE(const DB_cache<double>* db) : ::BSpline::IE_Fit<K>(db), ::BSpline::IE_Overlap<double,K>(db) {};
};

}} //namespace

#endif //_Atom_l_BSpline_IE_H_
